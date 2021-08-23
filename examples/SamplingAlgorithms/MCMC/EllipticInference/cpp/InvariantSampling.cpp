#include "FlowEquation.h"

#include "MUQ/Approximation/GaussianProcesses/CovarianceKernels.h"
#include "MUQ/Approximation/GaussianProcesses/GaussianProcess.h"

#include "MUQ/Modeling/ModPiece.h"
#include "MUQ/Modeling/WorkGraph.h"
#include "MUQ/Modeling/ModGraphPiece.h"
#include "MUQ/Modeling/LinearAlgebra/SliceOperator.h"
#include "MUQ/Modeling/Distributions/DensityProduct.h"
#include "MUQ/Modeling/CwiseOperators/CwiseUnaryOperator.h"

#include "MUQ/SamplingAlgorithms/SamplingProblem.h"
#include "MUQ/SamplingAlgorithms/DILIKernel.h"
#include "MUQ/SamplingAlgorithms/SingleChainMCMC.h"
#include "MUQ/SamplingAlgorithms/SampleCollection.h"

#include "MUQ/Utilities/RandomGenerator.h"

#include <boost/property_tree/ptree.hpp>

using namespace muq::Modeling;
using namespace muq::SamplingAlgorithms;
using namespace muq::Approximation;
using namespace muq::Utilities;

struct Discretization
{   
    Discretization(unsigned int numCellsIn)
    {     
        numCells = numCellsIn;
        numNodes = numCells+1;

        nodeLocs = Eigen::VectorXd::LinSpaced(numCells+1,0,1);

        cellLocs.resize(1,numCells);
        cellLocs.row(0) = 0.5*(nodeLocs.tail(numCells) + nodeLocs.head(numCells));
    };

    unsigned int numCells;
    unsigned int numNodes;

    // Node locations 
    Eigen::VectorXd nodeLocs;

    // Cell locations in a row vector
    Eigen::MatrixXd cellLocs;
};


Eigen::VectorXd GetTrueLogConductivity(Discretization const& mesh)
{
    return (10.0*mesh.cellLocs.row(0)).array().cos();
}

Eigen::VectorXd GenerateData(Discretization const& mesh, unsigned int obsThin, double obsVar)
{   
    // Generate the data
    unsigned numRefine = 1;
    Discretization fineMesh(numRefine*mesh.numCells);
    Eigen::VectorXd trueCond = GetTrueLogConductivity(fineMesh).array().exp();

    // Create the model with twice the mesh resolution
    Eigen::VectorXd recharge = Eigen::VectorXd::Ones(fineMesh.numCells);
    auto mod = std::make_shared<FlowEquation>(recharge);

    // Solve the forward problem with the true conductivity
    Eigen::VectorXd trueSol = mod->Evaluate( trueCond ).at(0);

    // Take every N node as an "observation"
    auto slicer = std::make_shared<SliceOperator>(fineMesh.numNodes,0,fineMesh.numCells,numRefine*obsThin);

    return slicer->Evaluate(trueSol).at(0) + std::sqrt(obsVar)*RandomGenerator::GetNormal(slicer->outputSizes(0));
}

std::shared_ptr<Gaussian> CreatePrior(Discretization const& mesh)
{
    // Define the prior distribution
    double priorVar = 1.0;
    double priorLength = 0.05;
    double priorNu = 3.0/2.0;
    
    auto covKernel = std::make_shared<MaternKernel>(1, priorVar, priorLength, priorNu); // The first argument "1" specifies we are working in 1d
    
    auto meanFunc = std::make_shared<ZeroMean>(1,1); // dimension of x, components in k(x) if it was vector-valued

    auto priorGP = std::make_shared<GaussianProcess>(meanFunc,covKernel);

    return priorGP->Discretize(mesh.cellLocs);
}

std::shared_ptr<ModPiece> DefinePosterior(Discretization const& mesh, Eigen::VectorXd const& data, unsigned int obsThin, double obsVar)
{   
    // Define the forward model
    Eigen::VectorXd recharge = Eigen::VectorXd::Ones(mesh.numCells);
    auto forwardMod = std::make_shared<FlowEquation>(recharge);
    
    auto graph = std::make_shared<WorkGraph>();
    graph->AddNode(std::make_shared<IdentityOperator>(mesh.numCells), "Log Conductivity");
    graph->AddNode(std::make_shared<ExpOperator>(mesh.numCells), "Conductivity");
    graph->AddEdge("Log Conductivity", 0, "Conductivity", 0);

    graph->AddNode(forwardMod, "Forward Model");
    graph->AddEdge("Conductivity", 0, "Forward Model", 0);
    graph->AddNode(std::make_shared<SliceOperator>(mesh.numNodes,0,mesh.numCells,obsThin), "Observables");
    graph->AddEdge("Forward Model", 0, "Observables", 0);
    
    
    auto priorDist = CreatePrior(mesh);
    auto likelihood = std::make_shared<Gaussian>(data, obsVar*Eigen::VectorXd::Ones(data.size()));

    graph->AddNode(likelihood->AsDensity(), "Likelihood");
    graph->AddEdge("Observables", 0, "Likelihood", 0);
    
    graph->AddNode(priorDist->AsDensity(), "Prior");
    graph->AddEdge("Log Conductivity", 0, "Prior", 0);
    
    graph->AddNode(std::make_shared<DensityProduct>(2), "Posterior");
    graph->AddEdge("Prior",0,"Posterior",0);
    graph->AddEdge("Likelihood",0,"Posterior",1);

    graph->Visualize("WorkGraph.pdf");

    return graph->CreateModPiece("Posterior");
}


std::shared_ptr<SampleCollection> SampleDILI(std::shared_ptr<ModPiece> const& posterior, Eigen::VectorXd const& startPt, unsigned int numSamps)
{
    boost::property_tree::ptree pt;
    pt.put("NumSamples",numSamps);
    pt.put("BurnIn", 5000);
    pt.put("PrintLevel",3);
    pt.put("HessianType","Exact");
    pt.put("Adapt Interval", 10000);
    pt.put("Initial Weight", 1);
    pt.put("Prior Node", "Prior");
    pt.put("Likelihood Node", "Likelihood");

    pt.put("Eigensolver Block", "EigOpts");
    pt.put("EigOpts.NumEigs",100); // Maximum number of generalized eigenvalues to compute (e.g., maximum LIS dimension)
    pt.put("EigOpts.ExpectedRank", 5);
    pt.put("EigOpts.OversamplingFactor", 2);
    //pt.put("EigOpts.AbsoluteTolerance", );
    pt.put("EigOpts.Verbosity",3);

    pt.put("LIS Block", "LIS");
    pt.put("LIS.Method", "MHKernel");
    pt.put("LIS.Proposal","MyProposal");
    pt.put("LIS.MyProposal.Method","MALAProposal");
    pt.put("LIS.MyProposal.ProposalVariance", 5.0);
    pt.put("LIS.MyProposal.StepSize", 0.05);

    pt.put("CS Block", "CS");
    pt.put("CS.Method", "MHKernel");
    pt.put("CS.Proposal","MyProposal");
    pt.put("CS.MyProposal.Method", "CrankNicolsonProposal");
    pt.put("CS.MyProposal.Beta",0.5);
    pt.put("CS.MyProposal.PriorNode","Prior");

    // create a sampling problem
    auto problem = std::make_shared<SamplingProblem>(posterior);

    std::vector<std::shared_ptr<TransitionKernel>> kernels(1);
    kernels.at(0) = std::make_shared<DILIKernel>(pt, problem);

    auto sampler = std::make_shared<SingleChainMCMC>(pt, kernels);

    return sampler->Run(startPt); // Use a true posterior sample to avoid burnin
}


std::shared_ptr<SampleCollection> SamplePCN(std::shared_ptr<ModPiece> const& posterior, Eigen::VectorXd const& startPt, unsigned int numSamps)
{   
    boost::property_tree::ptree pt;
    pt.put("NumSamples", numSamps); // number of Monte Carlo samples
    pt.put("BurnIn",5000);
    pt.put("PrintLevel",3);
    pt.put("KernelList", "Kernel1"); // the transition kernel
    pt.put("Kernel1.Method","MHKernel");
    pt.put("Kernel1.Proposal", "MyProposal"); // the proposal
    pt.put("Kernel1.MyProposal.Method", "CrankNicolsonProposal");
    pt.put("Kernel1.MyProposal.Beta", 0.3);
    pt.put("Kernel1.MyProposal.PriorNode", "Prior"); // The node in the WorkGraph containing the prior density

    // create a sampling problem
    auto problem = std::make_shared<SamplingProblem>(posterior);

    auto sampler = std::make_shared<SingleChainMCMC>(pt,problem);

    return sampler->Run(startPt); // Use a true posterior sample to avoid burnin
}

int main(){

    
    // Define the mesh
    unsigned int numCells = 100;
    Discretization mesh(numCells);

    // Generate synthetic "truth" data
    unsigned int obsThin = 10;
    double obsVar = 0.01*0.01;
    auto data = GenerateData(mesh, obsThin, obsVar);

    auto posterior = DefinePosterior(mesh, data, obsThin, obsVar);

    // Use a random draw from the prior as a starting point for MCMC
    Eigen::VectorXd startPt = GetTrueLogConductivity(mesh);//CreatePrior(mesh)->Sample();

    unsigned int numSamps = 100000;
    auto diliSamps = SampleDILI(posterior, startPt, numSamps);
    std::cout << "DILI Min ESS: " << diliSamps->ESS().minCoeff() << std::endl;

    auto pcnSamps = SamplePCN(posterior, startPt, numSamps);
    std::cout << "PCN Min ESS: " << pcnSamps->ESS().minCoeff() << std::endl;
    return 0;
}