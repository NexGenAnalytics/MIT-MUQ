#include "MUQ/Modeling/WorkGraph.h"
#include "MUQ/Modeling/LinearAlgebra/IdentityOperator.h"
#include "MUQ/Modeling/WorkGraphPiece.h"
#include "MUQ/Modeling/ModGraphPiece.h"

#include "MUQ/Modeling/Distributions/Gaussian.h"
#include "MUQ/Modeling/Distributions/DensityProduct.h"

#include "MUQ/SamplingAlgorithms/SingleChainMCMC.h"
#include "MUQ/SamplingAlgorithms/MHKernel.h"
#include "MUQ/SamplingAlgorithms/InfMALAProposal.h"

namespace pt = boost::property_tree;
using namespace muq::Modeling;
using namespace muq::SamplingAlgorithms;
using namespace muq::Utilities;

class Model : public ModPiece {
public:
  static const unsigned int dim = 2;

  inline Model() : ModPiece(Eigen::VectorXi::Constant(1, dim), Eigen::VectorXi::Constant(1, dim)) {}

  virtual ~Model() = default;
private:

  inline virtual void EvaluateImpl(ref_vector<Eigen::VectorXd> const& inputs) override {
    const Eigen::VectorXd& in = inputs[0].get();

    outputs.resize(1);
    outputs[0] = Eigen::VectorXd::Zero(dim);
    for( unsigned int i=0; i<dim; ++i ) {
      outputs[0][i] = 0.5*in[i]*in[i] + std::sin(M_PI*in[(i+1)%dim]);
    }
  }

  inline virtual void JacobianImpl(unsigned int const inwrt, unsigned int const outwrt, ref_vector<Eigen::VectorXd> const& inputs) override {
    const Eigen::VectorXd& in = inputs[0].get();

    // since we only have one in/output
    assert(inwrt==0); assert(outwrt==0);

    jacobian = Eigen::MatrixXd::Zero(dim, dim);
    for( unsigned int i=0; i<dim; ++i ) {
      jacobian(i, i) = in[i];
      jacobian(i, (i+1)%dim) = M_PI*std::cos(M_PI*in[(i+1)%dim]);
    }
  }
};

int main() {
  // create the forward model
  auto model = std::make_shared<Model>();

  // check the forward model implementation
  {
    // example input
    //const Eigen::VectorXd input = Eigen::VectorXd::Random(Model::dim);
    const Eigen::VectorXd input = Eigen::VectorXd::Zero(Model::dim);
    std::cout << "input vector: " << input.transpose() << std::endl;

    // evaluate the model
    const Eigen::VectorXd output = model->Evaluate(input) [0];
    std::cout << "output vector: " << output.transpose() << std::endl;

    const Eigen::MatrixXd jac = model->Jacobian(0, 0, input);
    std::cout << "Jacobian: " << std::endl << jac << std::endl;
    const Eigen::MatrixXd jacfd = model->JacobianByFD(0, 0, input);
    std::cout << "Jacobian (using finite differences): " << std::endl << jacfd << std::endl;

    const Eigen::VectorXd sens = Eigen::VectorXd::Random(Model::dim);
    const Eigen::VectorXd grad = model->Gradient(0, 0, input, sens);
    std::cout << "Gradient: " << grad.transpose() << std::endl;
    const Eigen::VectorXd gradfd = model->GradientByFD(0, 0, std::vector<Eigen::VectorXd>({input}), sens);
    std::cout << "Gradient (using finite differences): " << gradfd.transpose() << std::endl;

    const Eigen::VectorXd vec = Eigen::VectorXd::Random(Model::dim);
    const Eigen::VectorXd jacapp = model->ApplyJacobian(0, 0, input, vec);
    std::cout << "Jacobian action: " << jacapp.transpose() << std::endl;
    const Eigen::VectorXd jacappfd = model->ApplyJacobianByFD(0, 0, std::vector<Eigen::VectorXd>({input}), vec);
    std::cout << "Jacobian action: (using finite differences): " << jacappfd.transpose() << std::endl;

    std::cout << std::endl << "-----------" << std::endl << std::endl;
  }

  // number of Monte Carlo samples
  const unsigned int N = 1e6;

  // parameters for the sampler
  pt::ptree pt;
  pt.put("NumSamples", N);
  pt.put("BurnIn", (unsigned int)(0.1*N));
  pt.put("PrintLevel", 2);
  pt.put("KernelList", "Kernel1");
  pt.put("Kernel1.Method","MHKernel");
  pt.put("Kernel1.Proposal", "MyProposal");
  pt.put("Kernel1.MyProposal.Method", "InfMALAProposal");
  pt.put("Kernel1.MyProposal.Beta", 0.5);
  pt.put("Kernel1.MyProposal.StepSize", 0.5);
  pt.put("Kernel1.MyProposal.PriorNode", "Prior");

  // create a Gaussian distribution---the sampling problem is built around characterizing this distribution
  const Eigen::VectorXd mu = Eigen::VectorXd::Zero(Model::dim);
  const Eigen::MatrixXd priorCov = Eigen::MatrixXd::Identity(Model::dim, Model::dim);
  auto prior = std::make_shared<Gaussian>(mu, priorCov);

  Eigen::VectorXd data = Eigen::VectorXd::Zero(Model::dim) + 0.01*Eigen::VectorXd::Random(Model::dim);
  Eigen::MatrixXd likeCov = 0.1*Eigen::MatrixXd::Identity(Model::dim, Model::dim);
  auto likelihood = std::make_shared<Gaussian>(data, likeCov);

  auto graph = std::make_shared<WorkGraph>();
  graph->AddNode(std::make_shared<IdentityOperator>(Model::dim), "Parameters");
  graph->AddNode(prior->AsDensity(), "Prior");
  graph->AddNode(likelihood->AsDensity(), "Likelihood");
  graph->AddNode(model, "Model");

  graph->AddEdge("Parameters", 0, "Prior", 0);
  graph->AddEdge("Parameters", 0, "Model", 0);
  graph->AddEdge("Model", 0, "Likelihood", 0);

  graph->AddNode(std::make_shared<DensityProduct>(2), "Posterior");
  graph->AddEdge("Prior", 0, "Posterior", 0);
  graph->AddEdge("Likelihood", 0, "Posterior", 1);

  auto posterior = graph->CreateModPiece("Posterior");

  // create a sampling problem
  auto problem = std::make_shared<SamplingProblem>(posterior);

  // starting point for the MCMC chain
  std::vector<Eigen::VectorXd> start(1);
  start.at(0) = mu;

  // create an instance of MCMC and do the sampling
  auto mcmc = std::make_shared<SingleChainMCMC>(pt, problem);
  std::shared_ptr<SampleCollection> samps = mcmc->Run(start);

  const Eigen::VectorXd sampMean = samps->Mean();
  const Eigen::MatrixXd sampCov = samps->Covariance();

  std::cout << std::endl;
  std::cout << "Mean: " << sampMean.transpose() << std::endl;
  std::cout << "Cov: " << std::endl << sampCov << std::endl;
  std::cout << std::endl;
  std::cout << std::endl;
}
