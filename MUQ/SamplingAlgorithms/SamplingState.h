#ifndef SAMPLINGSTATE_H_
#define SAMPLINGSTATE_H_

#include <iostream>

#include <boost/any.hpp>

namespace muq {
  namespace SamplingAlgorithms {
    /// Each state is one sample generated by a sampling algorithm
    /**
       A state could include: the state variables, a weight, the target density evaluation
     */
    class SamplingState {
    public:

      SamplingState(boost::any const& state, double const weight);

      ~SamplingState();

      /// The state variables
      const boost::any state;

      /// The weight of this state
      const double weight;

    private:
      
    };
    
  } // namespace SamplingAlgoritms
} // namespace muq

#endif
