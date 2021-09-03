/***
## Overview
This example shows how to use the high-level (basic) API to MUQ's Multilevel MCMC algorithms.  The actual sampling 
problem is quite simple: we want to draw samples from a multivariate Gaussian distribution with 
mean 
\f[
  \mu = \left[ \begin{array}{c} 1\\ 2\end{array}\right]
\f]
and covariance 
\f[
\Sigma = \left[\begin{array}{cc} 0.7& 0.6\\ 0.6 & 1.0\end{array}\right].
\f]
It's of course possible to sample this distribution directly, but we will use Multilevel MCMC methods in 
this example to illustrate their use without introducing unnecessary complexity to the problem definition.   

*/

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

#include "MUQ/SamplingAlgorithms/Diagnostics.h"

#include "MUQ/SamplingAlgorithms/MIComponentFactory.h"

#include "MUQ/Utilities/RandomGenerator.h"

#include <boost/property_tree/ptree.hpp>
#include <sstream>

namespace pt = boost::property_tree;
using namespace muq::Modeling;
using namespace muq::SamplingAlgorithms;
using namespace muq::Utilities;


/***
  ## Define the Target Distributions
  To apply MLMCMC, we need to define 
*/
std::vector<std::shared_ptr<ModPiece>> ConstructDensities()
{
  unsigned int numLevels = 4;
  std::vector<std::shared_ptr<ModPiece>> logDensities(numLevels);

  Eigen::VectorXd tgtMu(2);
  Eigen::MatrixXd tgtCov(2,2);
  tgtMu  << 1.0, 2.0;
  tgtCov << 0.7, 0.6,
            0.6, 1.0;

  // Define the problem at the coarsest level
  Eigen::VectorXd levelMu  = 0.8*tgtMu;
  Eigen::MatrixXd levelCov = 2.0*tgtCov;
  logDensities.at(0) = std::make_shared<Gaussian>(levelMu, levelCov)->AsDensity();

  // Define the second coarsest level
  levelMu  = 0.9*tgtMu;
  levelCov = 1.5*tgtCov;
  logDensities.at(1) = std::make_shared<Gaussian>(levelMu, levelCov)->AsDensity();

  // Define the second finest level
  levelMu = 0.99*tgtMu;
  levelCov = 1.1*tgtCov;
  logDensities.at(2) = std::make_shared<Gaussian>(levelMu, levelCov)->AsDensity();

  // Deifne the finest level.  This should be the target distribution. 
  logDensities.at(3) = std::make_shared<Gaussian>(tgtMu, tgtCov)->AsDensity();

  return logDensities;
}

int main(){

  std::vector<std::shared_ptr<ModPiece>> logDensities = ConstructDensities();


  pt::ptree options;
  
  options.put("NumSamples", 1e2); // number of samples for single level
  options.put("NumInitialSamples", 1e3); // number of initial samples for greedy MLMCMC
  options.put("GreedyTargetVariance", 0.1); // estimator variance to be achieved by greedy algorithm
  options.put("verbosity", 1); // show some output
  options.put("MLMCMC.Subsampling_0", 8);
  options.put("MLMCMC.Subsampling_1", 4);
  options.put("MLMCMC.Subsampling_2", 2);
  options.put("MLMCMC.Subsampling_3", 0);

  options.put("Proposal.Method", "MHProposal");
  options.put("Proposal.ProposalVariance", 0.5);


  unsigned int numChains = 5;
  std::vector<std::shared_ptr<MultiIndexEstimator>> estimators(numChains);

  for(int chainInd=0; chainInd<numChains; ++chainInd){
    Eigen::VectorXd x0 = RandomGenerator::GetNormal(2);
    
    std::cout << std::endl << "*************** greedy multilevel chain" << std::endl << std::endl;
   
    GreedyMLMCMC sampler(options, x0, logDensities);
    estimators.at(chainInd) = sampler.Run();

    std::cout << "mean QOI: " << estimators.at(chainInd)->Mean().transpose() << std::endl;

    std::stringstream filename;
    filename << "MultilevelGaussianSampling_Chain" << chainInd << ".h5";
    sampler.WriteToFile(filename.str());
  }

  std::cout << "Rhat Convergence Diagnostic: " << Diagnostics::Rhat(estimators) << std::endl;

  std::cout << std::endl << "*************** single chain reference" << std::endl << std::endl;

  // Eigen::VectorXd x0 = RandomGenerator::GetNormal(numChains);
  
  // SLMCMC slmcmc (pt, componentFactory);
  // slmcmc.Run();
  // std::cout << "mean QOI: " << slmcmc.GetQOIs()->Mean().transpose() << std::endl;

  return 0;
}
