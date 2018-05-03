#ifndef SAMPLINGSTATE_H_
#define SAMPLINGSTATE_H_

#include <iostream>
#include <unordered_map>
#include <vector>

#include <boost/any.hpp>

namespace muq {
  namespace SamplingAlgorithms {
    /// Each state is one sample generated by a sampling algorithm
    /**
       A state could include: the state variables, a weight, the target density evaluation
     */
    class SamplingState {
    public:

      SamplingState(boost::any const& stateIn, double weight = 1.0);
      SamplingState(std::vector<boost::any> const& stateIn, double weight = 1.0);

      ~SamplingState();

      /// The state variables
      const std::vector<boost::any> state;

      /// The weight of this state
      double weight;

      /// Checks to see if the meta map contains a particular key
      bool HasMeta(std::string const& metaKey);

      /// The total number of parameters in the state, i.e., the sum of state[i].size()
      int TotalDim() const;

      /// A map containing extra information like the target density, run time, forward model output, etc...
      std::unordered_map<std::string, boost::any> meta;

    private:

    };

  } // namespace SamplingAlgoritms
} // namespace muq

#endif
