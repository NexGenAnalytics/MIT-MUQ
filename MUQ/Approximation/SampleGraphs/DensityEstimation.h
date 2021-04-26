#ifndef DENSITYESTIMATION_H_
#define DENSITYESTIMATION_H_

#include "MUQ/Approximation/SampleGraphs/SampleGraph.h"

namespace muq {
namespace Approximation {

/// Estimate the probability density function \f$\psi(x)\f$ given samples
class DensityEstimation : public SampleGraph {
public:

  /// Construct the sample graph by sampling a random variable from \f$\psi\f$
  /**
  @param[in] rv The random variable that we wish to sample
  @param[in] options Setup options
  */
  DensityEstimation(std::shared_ptr<muq::Modeling::RandomVariable> const& rv, boost::property_tree::ptree const& options);

  /// Construct the sample graph given samples from the underlying distribution \f$\psi\f$
  /**
  @param[in] samples Samples from the underlying distribution \f$\psi\f$
  @param[in] options Setup options
  */
  DensityEstimation(std::shared_ptr<muq::SamplingAlgorithms::SampleCollection> const& samples, boost::property_tree::ptree const& options);

  virtual ~DensityEstimation() = default;

  /// Estimate the density at each sample
  /**
  This function computes the bandwidth \f$r_i^2 = \frac{1}{k} \sum_{j=1}^{k} \| \boldsymbol{x}^{(i)}-\boldsymbol{x}^{(I(i,j))} \|^2\f$ and the kernel marix \f$\boldsymbol{K}_{\epsilon}\f$. We then estimate the density at each sample as
  \f{equation*}{
  \psi^{(i)}_{\epsilon} = \sum_{j=1}^{n} \frac{K_{\epsilon}^{(ij)}}{n (\pi \epsilon r_i^2)^{m/2}},
  \f}
  where \f$m\f$ is the manifold dimension.

  If <tt>tune</tt> is <tt>true</tt>, we must choose the optimal \f$\epsilon\f$ for the density estimation. Given the bandwidth exponent \f$l \in [l_{min}, l_{max}]\f$, we choose \f$\epsilon \in [2^{(l_{min})}, 2^{(l_{max})}]\f$. We do this by defining the parameters \f$\Sigma_{l} = \frac{1}{n^2} \sum_{i,j=1}^{n} K_{\epsilon_l}^{(ij)}\f$ and \f$\Sigma_{\epsilon_l}^{\prime} = (\log{(\Sigma_{\epsilon_{l+1}})}-\log{(\Sigma_{\epsilon_l})})/\log{(2)}\f$

  The optimal \f$\epsilon\f$ maximizes \f$\Sigma_{\epsilon}^{\prime}\f$ and the corresponding manifold dimension estimate is \f$m = 2 \max{(\Sigma_{\epsilon}^{\prime})}\f$.

  @param[in] tune <tt>true</tt> (default): Tune the bandwidth parameter values; <tt>false</tt>: use the stored parameters
  \return The density estimation at each sample \f$\psi^{(i)} \approx \psi(\boldsymbol{x}^{(i)})\f$
  */
  Eigen::VectorXd Estimate(bool const tune = true);

private:
};

} // Approximation
} // namespace muq

#endif
