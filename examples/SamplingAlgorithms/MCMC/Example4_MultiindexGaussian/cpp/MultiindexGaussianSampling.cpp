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


  std::cout << std::endl << "*************** multiindex chain" << std::endl << std::endl;

  MIMCMC mimcmc (pt, componentFactory);
  mimcmc.Run();
  mimcmc.Draw(false);
  std::cout << "mean QOI: " << mimcmc.MeanQOI().transpose() << std::endl;

  std::cout << std::endl << "*************** single chain reference" << std::endl << std::endl;

  SLMCMC slmcmc (pt, componentFactory);
  slmcmc.Run();
  std::cout << "mean QOI: " << slmcmc.MeanQOI().transpose() << std::endl;

  return 0;
}
