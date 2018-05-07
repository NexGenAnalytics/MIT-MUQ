#include <gtest/gtest.h>

#include <boost/property_tree/ptree.hpp>

#include "MUQ/Modeling/Distributions/Gaussian.h"

#include "MUQ/SamplingAlgorithms/SingleChainMCMC.h"
#include "MUQ/SamplingAlgorithms/MHKernel.h"
#include "MUQ/SamplingAlgorithms/MHProposal.h"
#include "MUQ/SamplingAlgorithms/AMProposal.h"

namespace pt = boost::property_tree;
using namespace muq::Modeling;
using namespace muq::SamplingAlgorithms;


TEST(MCMC, MHKernel_MHProposal) {

  const unsigned int N = 1e4;

  // parameters for the sampler
  pt::ptree pt;
  pt.put("MyMCMC.NumSamples", N); // number of Monte Carlo samples
  pt.put("MyMCMC.KernelList", "Kernel1"); // the transition kernel
  pt.put("MyMCMC.Kernel1.Method","MHKernel");
  pt.put("MyMCMC.Kernel1.Proposal", "MyProposal"); // the proposal
  pt.put("MyMCMC.Kernel1.MyProposal.Method", "MHProposal");
  pt.put("MyMCMC.Kernel1.MyProposal.ProposalVariance", 1e-3); // the variance of the isotropic MH proposal

  // create a Gaussian distribution---the sampling problem is built around characterizing this distribution
  const Eigen::VectorXd mu = Eigen::VectorXd::Ones(2);
  auto dist = std::make_shared<Gaussian>(mu); // standard normal Gaussian

  // create a sampling problem
  std::vector<int> blockSizes(1,mu.size());
  auto problem = std::make_shared<SamplingProblem>(dist, blockSizes);

  // starting point
  const Eigen::VectorXd start = Eigen::VectorXd::Random(2);

  // evaluate
  // create an instance of MCMC
  auto mcmc = std::make_shared<SingleChainMCMC>(pt.get_child("MyMCMC"),problem);

  // Make sure the kernel and proposal are correct
  std::shared_ptr<TransitionKernel> kernelBase = mcmc->Kernels().at(0);
  std::shared_ptr<MHKernel> kernelMH = std::dynamic_pointer_cast<MHKernel>(kernelBase);
  EXPECT_TRUE(kernelMH);

  std::shared_ptr<MCMCProposal> proposalBase = kernelMH->Proposal();
  std::shared_ptr<MHProposal> proposalMH = std::dynamic_pointer_cast<MHProposal>(proposalBase);
  EXPECT_TRUE(proposalMH);

  mcmc->Run(start);
}

TEST(MCMC, MHKernel_AMProposal) {


  const unsigned int N = 2.0e5;

  // parameters for the sampler
  pt::ptree pt;
  pt.put("MyMCMC.NumSamples", N); // number of Monte Carlo samples
  pt.put("MyMCMC.KernelList", "Kernel1"); // the transition kernel
  pt.put("MyMCMC.Kernel1.Method","MHKernel");
  pt.put("MyMCMC.Kernel1.Proposal", "MyProposal"); // the proposal
  pt.put("MyMCMC.Kernel1.MyProposal.Method", "AMProposal");
  pt.put("MyMCMC.Kernel1.MyProposal.ProposalVariance", 1e-3); // the variance of the isotropic MH proposal
  pt.put("MyMCMC.Kernel1.MyProposal.AdaptSteps",200);
  pt.put("MyMCMC.Kernel1.MyProposal.AdaptStart",2000);
  pt.put("MyMCMC.Kernel1.MyProposal.AdaptScale",1e-1);

  // create a Gaussian distribution---the sampling problem is built around characterizing this distribution
  const Eigen::VectorXd mu = Eigen::VectorXd::Ones(2);
  auto dist = std::make_shared<Gaussian>(mu); // standard normal Gaussian

  // create a sampling problem
  std::vector<int> blockSizes(1,mu.size());
  auto problem = std::make_shared<SamplingProblem>(dist, blockSizes);

  // starting point
  const Eigen::VectorXd start = Eigen::VectorXd::Random(2);

  // evaluate
  // create an instance of MCMC
  auto mcmc = std::make_shared<SingleChainMCMC>(pt.get_child("MyMCMC"), problem);
  mcmc->Run(start);

  std::shared_ptr<TransitionKernel> kernelBase = mcmc->Kernels().at(0);
  std::shared_ptr<MHKernel> kernelMH = std::dynamic_pointer_cast<MHKernel>(kernelBase);
  EXPECT_TRUE(kernelMH);

  std::shared_ptr<MCMCProposal> proposalBase = kernelMH->Proposal();
  std::shared_ptr<AMProposal> proposalAM = std::dynamic_pointer_cast<AMProposal>(proposalBase);
  EXPECT_TRUE(proposalAM);
}
