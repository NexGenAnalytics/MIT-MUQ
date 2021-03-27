#include <iostream>
#include <string>

#include <boost/property_tree/ptree.hpp>

#include <MUQ/Modeling/OneStepCachePiece.h>

#include <MUQ/SamplingAlgorithms/SamplingProblem.h>
#include <MUQ/SamplingAlgorithms/SingleChainMCMC.h>
#include <MUQ/SamplingAlgorithms/MCMCFactory.h>

#include "HTTPModPiece.h"

namespace pt = boost::property_tree;

using namespace muq::Utilities;
using namespace muq::Modeling;
using namespace muq::SamplingAlgorithms;


int main(int argc, char* argv[])
{

  auto modelModPiece = std::make_shared<HTTPModPiece>("localhost:4242");

  auto samplingProblem = std::make_shared<SamplingProblem>(modelModPiece);

  pt::ptree pt;
  pt.put("NumSamples", 1e3); // number of MCMC steps
  pt.put("BurnIn", 5e1);
  pt.put("PrintLevel",3);
  pt.put("KernelList", "Kernel1"); // Name of block that defines the transition kernel
  pt.put("Kernel1.Method","MHKernel");  // Name of the transition kernel class
  pt.put("Kernel1.Proposal", "MyProposal"); // Name of block defining the proposal distribution
  pt.put("Kernel1.MyProposal.Method", "AMProposal"); // Name of proposal class
  pt.put("Kernel1.MyProposal.InitialVariance", 1);
  pt.put("Kernel1.MyProposal.AdaptSteps", 250);
  pt.put("Kernel1.MyProposal.AdaptStart", 250);

  Eigen::VectorXd startPt = Eigen::VectorXd::Zero(1);
  auto mcmc = MCMCFactory::CreateSingleChain(pt, samplingProblem);

  mcmc->Run(startPt);

  std::cout << "Sample mu: " << std::endl << mcmc->GetSamples()->Mean() << std::endl;
  std::cout << "Sample sigma: " << std::endl << mcmc->GetSamples()->Variance() << std::endl;

  return 0;
}