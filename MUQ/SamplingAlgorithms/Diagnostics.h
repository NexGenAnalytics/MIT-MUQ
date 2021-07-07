#ifndef DIAGNOSTICS_H_
#define DIAGNOSTICS_H_

#include <vector>

#include "MUQ/SamplingAlgorithms/SampleCollection.h"

#include <boost/property_tree/ptree.hpp>

namespace muq{
  namespace SamplingAlgorithms{

    namespace Diagnostics{
      /**
          @ingroup mcmcdiag
          @{
      */


      /** Computes the standard \f$\hat{R}\f$ diagnostic from \cite Gelman2013 or one of the modifications presented in \cite Vehtari2021

      Parameter Key | Type | Default Value | Description |
      ------------- | ------------- | ------------- | ------------- |
      "Split"  | boolean | True  | Whether the chains should be split in half as proposed by \cite Vehtari2021. |
      "Normalize"   | boolean | True  | If the parameters should be rank-transformed (i.e., normalized) before computing Rhat. |

      @param[in] collections A vector of SampleCollection variables returned by independent runs of an MCMC algorithm.  Assumes that all of the chains have the same length and were initialized with diffuse initial conditions.
      @returns A vector of \f$\hat{R}\f$ values for each component of the parameters.
      */
      Eigen::VectorXd Rhat(std::vector<std::shared_ptr<SampleCollection>> const& collections,
                           boost::property_tree::ptree options = boost::property_tree::ptree());

      /** For a set of scalar values \f$\{x_1,\ldots, x_S\}\f$, the rank of \f$x_i\f$ is the index of \f$x_i\f$ after sorting this set into a list that satisfies \f$x_i\leq x_{i+1}\f$.  We use \f$r_i\f$  to denote the rank of \f$x_i\f$ and adopt the convention that for repeated values (i.e., \f$x_{i}=x_{i+1}\f$), \f$r_i\f$ is given by the average rank of consecutive repeated values.

      This function computes and returns the values of \f$r_i\f$ when the initial set is given by the combined states of multiple SampleCollections.   The states are vector-valued, but this function operates only on a single component of the state vector.

      @param[in] collections Several SampleCollection instances used to define the set of values we want to rank
      @param[in] dim The component of the vector-valued states that we want to rank.
      @return A std::vector of Eigen::VectorXd containing the ranks of all samples in all sample collections.  The std::vector has the same size as the std::vector of collections.  Following \cite Vehtari2021, "Average rank for ties are used to conserve the number of unique values of discrete quantities."
      */
      std::vector<Eigen::VectorXd> ComputeRanks(std::vector<std::shared_ptr<SampleCollection>> const& collections,
                                   unsigned int                         dim);

      /**
      @}
      */
      
    } // namespace Diagnostics
  } // namespace SamplingAlgorithms
} // namespace muq

#endif // #ifndef DIAGNOSTICS_H_
