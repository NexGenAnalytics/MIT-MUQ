#include "MUQ/SamplingAlgorithms/Diagnostics.h"

#include <stan/math/fwd/scal.hpp>
#include <algorithm>

using namespace muq::SamplingAlgorithms;

Eigen::VectorXd muq::SamplingAlgorithms::Diagnostics::Rhat(std::vector<std::shared_ptr<SampleCollection>> const& origChains,
                                                           boost::property_tree::ptree options)
{

  std::vector<std::shared_ptr<SampleCollection>> chains;

  unsigned int chainLength = origChains.at(0)->size();
  unsigned int numChains = origChains.size();
  const unsigned int dim = origChains.at(0)->at(0)->TotalDim();

  // Split the original chains
  if(options.get("Split",true)){

    // Figure out how long the split chains will be
    chainLength = std::floor(0.5*chainLength);
    numChains = 2*origChains.size();

    chains.resize(2*origChains.size());

    // Extract the split chains
    for(int i=0; i<origChains.size();++i){
      chains.at(2*i) = origChains.at(i)->head(chainLength);
      chains.at(2*i+1) = origChains.at(i)->tail(chainLength);
    }

  }else{
    chains.insert(chains.begin(), origChains.begin(), origChains.end());
  }

  // The total number of samples available using all chains
  const unsigned int totalSamps = numChains*chainLength;

  // A matrix of the chain means.
  Eigen::MatrixXd chainMeans(dim,numChains);
  Eigen::MatrixXd chainVars(dim,numChains);

  // Apply the rank normalization if desired
  if(options.get("Normalize",true)){

    for(unsigned int i=0; i<dim; ++i){
      // Compute the ranks
      std::vector<Eigen::VectorXd> ranks = ComputeRanks(chains,i);

      // Apply a normal transformation to the ranks and compute chain means and variances.  See eqn. (14) in https://arxiv.org/pdf/1903.08008.pdf
      for(unsigned int chainInd=0; chainInd<ranks.size(); ++chainInd){
        ranks.at(chainInd) = ( (ranks.at(chainInd).array()+0.625)/(totalSamps + 0.25) ).unaryExpr([](double v){return stan::math::inv_Phi(v);});
        chainMeans(i, chainInd) = ranks.at(chainInd).mean();
        chainVars(i, chainInd) = (ranks.at(chainInd).array() - chainMeans(i, chainInd)).pow(2).sum() / (chainLength -1);
      }
    }

  }else{

    for(int i=0; i<numChains; ++i){
      chainMeans.col(i) = chains.at(i)->Mean();
      chainVars.col(i) = chains.at(i)->Variance();
    }
  }


  // /////////////////////////////////////////////////////////////////////
  // Now we're good to compute the standard Rhat estimator of Gelman 2013

  Eigen::VectorXd globalMean = chainMeans.rowwise().mean();
  Eigen::VectorXd W = chainVars.rowwise().mean(); // average within-chain variance

  // Compute the between-chain variance
  Eigen::VectorXd B = (chainLength/(numChains-1)) * (chainMeans.colwise()-globalMean).eval().array().pow(2).rowwise().sum();

  // Get an estimate of the marginal posterior variance
  Eigen::VectorXd rhat = ((chainLength-1.0)/chainLength + B.array()/(W.array()*chainLength) ).sqrt();
  return rhat;

}



std::vector<Eigen::VectorXd> muq::SamplingAlgorithms::Diagnostics::ComputeRanks(std::vector<std::shared_ptr<SampleCollection>> const& collections,
                                                                   unsigned int                                      dim)
{

  // A vector of sample indices [chainInd, sampInd]
  std::vector<std::pair<unsigned int, unsigned int>> sampInds;

  for(unsigned int chainInd=0; chainInd<collections.size(); ++chainInd){
    for(unsigned int sampInd=0; sampInd<collections.at(chainInd)->size(); ++sampInd)
      sampInds.push_back(std::make_pair(chainInd,sampInd));
  }

  // Sort the vector of indices according to the value of the parameters
  auto compLambda = [&collections, dim](std::pair<unsigned int, unsigned int> const& p1, std::pair<unsigned int, unsigned int> const& p2) {
                      return collections.at(p1.first)->at(p1.second)->StateValue(dim) < collections.at(p2.first)->at(p2.second)->StateValue(dim);
                    };

  std::stable_sort(sampInds.begin(), sampInds.end(), compLambda);

  // Set up empty vectors for storing the ranks
  std::vector<Eigen::VectorXd> ranks(collections.size());
  for(unsigned int i=0; i<ranks.size(); ++i)
    ranks.at(i).resize(collections.at(i)->size());

  // Figure out the rank of each sample
  unsigned int rawRank = 0;
  double currVal, nextVal;
  unsigned int numRepeat, chainInd, sampInd;

  while(rawRank < sampInds.size()){
    std::tie(chainInd, sampInd) = sampInds.at(rawRank);
    currVal = collections.at(chainInd)->at(sampInd)->StateValue(dim);

    // Look ahead and find the next sample with a new value
    numRepeat = 1;
    for(numRepeat=1; numRepeat<sampInds.size()-rawRank; ++numRepeat){
      std::tie(chainInd, sampInd) = sampInds.at(rawRank+numRepeat);
      nextVal = collections.at(chainInd)->at(sampInd)->StateValue(dim);

      if(std::abs(currVal-nextVal)>1e-15){
        break;
      }
    }

    // Compute the average rank across all of the duplicates
    double avgRank = 0.5*(rawRank + rawRank+numRepeat-1);

    // Set the ranks to the average value
    for(int i=rawRank; i<rawRank+numRepeat; ++i){
      std::tie(chainInd, sampInd) = sampInds.at(i);
      ranks.at(chainInd)(sampInd) = avgRank;
    }

    // Update how many samples we've moved through
    rawRank += numRepeat;
  }

  return ranks;
}
