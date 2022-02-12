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

#include "MultilevelProblem.h"


int main(int argc, char* argv[]){
  if (argc <= 1) {
    std::cerr << "Call with: ./binary HOSTNAME:PORT" << std::endl;
    exit(-1);
  }
  std::string host(argv[1]);

  pt::ptree pt;

  pt.put("NumSamples", 1e2); // number of samples for single level
  pt.put("verbosity", 1); // show some output
  pt.put("MLMCMC.Subsampling_0", 25);
  pt.put("MLMCMC.Subsampling_1", 0);
  //pt.put("MLMCMC.Subsampling_2", 2);
  //pt.put("MLMCMC.Subsampling_3", 0);

  pt.put("NumSamples_0", 1e4);
  pt.put("NumSamples_1", 1e2);
  //pt.put("NumSamples_2", 1e2);
  //pt.put("NumSamples_3", 1e1);
  pt.put("RegressionOptions.NumNeighbors", 10);
  pt.put("NumNeighbors", 10);

  //pt.put("GammaScale", 1e-4);
  //pt.put("TailCorrection", 1);
  //pt.put("LyapunovScale", 1e-2);
  pt.put("GammaScale", .1);
  pt.put("TailCorrection", 0);
  pt.put("LyapunovScale", 1.0);

  auto componentFactory = std::make_shared<MyMIComponentFactory>(pt, host);


  std::cout << std::endl << "*************** multilevel chain" << std::endl << std::endl;

  MIMCMC mimcmc (pt, componentFactory);
  mimcmc.Run();
  for (int boxindex = 0; boxindex < mimcmc.GetIndices()->Size(); boxindex++) {
    auto index = mimcmc.GetIndices()->at(boxindex);
    std::cout << mimcmc.GetBox(index)->FinestChain()->GetSamples()->size() << " samples on level " << *index << std::endl;
    /*auto samplingProblem = std::dynamic_pointer_cast<ExpensiveSamplingProblem>(mimcmc.GetBox(index)->GetFinestProblem());
    std::cout << *index << ": " << samplingProblem->cumbeta + samplingProblem->cumgamma + samplingProblem->cumkappa << std::endl;*/
  }
  //std::cout << "mean QOI: " << mimcmc.GetQOIs()->Mean() << std::endl;

  mimcmc.WriteToFile("MultilevelGaussianSampling.h5");


  /*std::cout << std::endl << "*************** single chain reference" << std::endl << std::endl;

  SLMCMC slmcmc (pt, componentFactory);
  slmcmc.Run();
  std::cout << "mean QOI: " << slmcmc.MeanQOI().transpose() << std::endl;*/

  return 0;
}
