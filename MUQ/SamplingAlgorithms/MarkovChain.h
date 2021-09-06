#ifndef MarkovChain_H
#define MarkovChain_H

#include "MUQ/SamplingAlgorithms/SampleCollection.h"

namespace muq {
  namespace SamplingAlgorithms{

    /**
    @ingroup MCMC
    @class MarkovChain
    @brief A class for storing and working with the results of Markov chain Monte Carlo algorithms.
    @details The MarkovChain class is a child of SampleCollection where the sample
    weights correspond to the number of consecutive steps taking the same value,
    and the weights are unnormalized (i.e., do not sum to one).  This is a useful
    class for storing the chain produced by an MCMC algorithm without storing the
    duplicate points that result from rejected proposals.
    */
    class MarkovChain : public SampleCollection
    {
    public:

      MarkovChain() = default;

      virtual ~MarkovChain() = default;

      /** Computes the effective sample size of the Markov chain.  

        If method=="Wolff", the spectral method described in
            "Monte Carlo errors with less error" by Ulli Wolff is employed.
            This returns an ESS for each component of the chain.

        If method=="Batch" (default) The overlapping batch method (OBM) described in \cite Flegal2010 
            is used.  This method is also applied to each component independently,
            resulting in an ESS estimate for each component.

        If method=="MultiBatch",  The multivariate method of \cite Vats2019 is employed.  This 
            method takes into account the joint correlation of all components of the chain and 
            returns a single ESS.   This approach is preferred in high dimensional settings.
      */
      Eigen::VectorXd ESS(std::string const& method="Batch") const override{return ESS(-1,method);};
      virtual Eigen::VectorXd ESS(int blockDim) const override{return ESS(blockDim,"");};
      virtual Eigen::VectorXd ESS(int blockDim, std::string const& method) const override;


      Eigen::VectorXd WolffESS(int blockDim) const;
      
      static double SingleComponentWolffESS(Eigen::Ref<const Eigen::VectorXd> const& trace);

      virtual std::shared_ptr<SampleCollection> segment(unsigned int startInd, unsigned int length, unsigned int skipBy=1) const override;

    private:


      std::vector<std::unordered_map<std::string, boost::any> > meta;

    }; // class MarkovChain
  }
}

#endif
