#include "MUQ/SamplingAlgorithms/SLMCMC.h"
#include "MUQ/SamplingAlgorithms/GreedyMLMCMC.h"
#include "MUQ/SamplingAlgorithms/MIMCMC.h"

#include "MUQ/Modeling/Distributions/Gaussian.h"
#include "MUQ/Modeling/Distributions/Density.h"

#include "MUQ/SamplingAlgorithms/MHKernel.h"
#include "MUQ/SamplingAlgorithms/MHProposal.h"
#include "MUQ/SamplingAlgorithms/CrankNicolsonProposal.h"
#include "MUQ/SamplingAlgorithms/SamplingProblem.h"
#include "MUQ/SamplingAlgorithms/SubsamplingMIProposal.h"

#include "MUQ/SamplingAlgorithms/MIComponentFactory.h"
#include "MUQ/SamplingAlgorithms/Diagnostics.h"

#include <boost/property_tree/ptree.hpp>

namespace pt = boost::property_tree;
using namespace muq::Modeling;
using namespace muq::SamplingAlgorithms;
using namespace muq::Utilities;


#include "Problem.h"


int main(){

  pt::ptree pt;

  pt.put("NumSamples", 1e4); // number of samples for single level MCMC
  pt.put("NumSamples_0_0", 1e5);
  pt.put("NumSamples_0_1", 1e5);
  pt.put("NumSamples_0_2", 1e4);
  pt.put("NumSamples_1_0", 1e5);
  pt.put("NumSamples_1_1", 1e4);
  pt.put("NumSamples_1_2", 1e3);
  pt.put("NumSamples_2_0", 1e4);
  pt.put("NumSamples_2_1", 1e3);
  pt.put("NumSamples_2_2", 1e3);
  pt.put("MLMCMC.Subsampling_0_0", 5);
  pt.put("MLMCMC.Subsampling_0_1", 5);
  pt.put("MLMCMC.Subsampling_0_2", 5);
  pt.put("MLMCMC.Subsampling_1_0", 5);
  pt.put("MLMCMC.Subsampling_1_1", 5);
  pt.put("MLMCMC.Subsampling_1_2", 5);
  pt.put("MLMCMC.Subsampling_2_0", 5);
  pt.put("MLMCMC.Subsampling_2_1", 5);
  pt.put("MLMCMC.Subsampling_2_2", 5);
  auto componentFactory = std::make_shared<MyMIComponentFactory>(pt);


  unsigned int numChains = 5;
  std::vector<std::shared_ptr<MultiIndexEstimator>> estimators(numChains);

  for(int chainInd=0; chainInd<numChains; ++chainInd){
    Eigen::VectorXd x0 = RandomGenerator::GetNormal(2);
    auto componentFactory = std::make_shared<MyMIComponentFactory>(x0, pt);

    std::cout << "\n=============================\n";
    std::cout << "Running MIMCMC Chain " << chainInd << ": \n";
    std::cout << "-----------------------------\n";

    MIMCMC mimcmc (pt, componentFactory);
    estimators.at(chainInd) = mimcmc.Run();

    std::cout << "mean QOI: " << estimators.at(chainInd)->Mean().transpose() << std::endl;

    std::stringstream filename;
    filename << "MultilevelGaussianSampling_Chain" << chainInd << ".h5";
    greedymlmcmc.WriteToFile(filename.str());
  }
  
  std::cout << "\n=============================\n";
  std::cout << "Multilevel Summary: \n";
  std::cout << "-----------------------------\n";
  std::cout << "  Rhat:               " << Diagnostics::Rhat(estimators).transpose() << std::endl;
  std::cout << "  Mean (chain 0):     " << estimators.at(0)->Mean().transpose() << std::endl;
  std::cout << "  MCSE (chain 0):     " << estimators.at(0)->StandardError().transpose() << std::endl;
  std::cout << "  ESS (chain 0):      " << estimators.at(0)->ESS().transpose() << std::endl;
  std::cout << "  Variance (chain 0): " << estimators.at(0)->Variance().transpose() << std::endl;
  std::cout << std::endl;



  std::cout << "\n=============================\n";
  std::cout << "Running Single Level Chain" << ": \n";
  std::cout << "-----------------------------\n";

  SLMCMC slmcmc (pt, componentFactory);
  auto samps = slmcmc.Run();

  std::cout << "\n=============================\n";
  std::cout << "Single Level Summary: \n";
  std::cout << "-----------------------------\n";
  std::cout << "  Mean:               " << samps->Mean().transpose() << std::endl;
  std::cout << "  MCSE:               " << samps->StandardError().transpose() << std::endl;
  std::cout << "  ESS:                " << samps->ESS().transpose() << std::endl;
  std::cout << "  Variance:           " << samps->Variance().transpose() << std::endl;
  std::cout << std::endl;

  return 0;
}
