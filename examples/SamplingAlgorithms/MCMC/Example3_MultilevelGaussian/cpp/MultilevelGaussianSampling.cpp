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

#include <boost/property_tree/ptree.hpp>

namespace pt = boost::property_tree;
using namespace muq::Modeling;
using namespace muq::SamplingAlgorithms;
using namespace muq::Utilities;

#include "Problem.h"


int main(){

  pt::ptree pt;

  pt.put("NumSamples", 1e2); // number of samples for single level
  pt.put("NumInitialSamples", 1e3); // number of initial samples for greedy MLMCMC
  pt.put("GreedyTargetVariance", 0.1); // estimator variance to be achieved by greedy algorithm
  pt.put("verbosity", 1); // show some output
  pt.put("MLMCMC.Subsampling_0", 8);
  pt.put("MLMCMC.Subsampling_1", 4);
  pt.put("MLMCMC.Subsampling_2", 2);
  pt.put("MLMCMC.Subsampling_3", 0);


  auto componentFactory = std::make_shared<MyMIComponentFactory>(pt);


  std::cout << std::endl << "*************** greedy multilevel chain" << std::endl << std::endl;

  GreedyMLMCMC greedymlmcmc (pt, componentFactory);
  greedymlmcmc.Run();
  std::cout << "mean QOI: " << greedymlmcmc.MeanQOI().transpose() << std::endl;

  greedymlmcmc.WriteToFile("MultilevelGaussianSampling.h5");


  std::cout << std::endl << "*************** single chain reference" << std::endl << std::endl;

  SLMCMC slmcmc (pt, componentFactory);
  slmcmc.Run();
  std::cout << "mean QOI: " << slmcmc.MeanQOI().transpose() << std::endl;

  return 0;
}
