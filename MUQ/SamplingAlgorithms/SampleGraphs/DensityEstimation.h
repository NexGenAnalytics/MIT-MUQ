#ifndef DENSITYESTIMATION_H_
#define DENSITYESTIMATION_H_

#include "MUQ/SamplingAlgorithms/SampleGraphs/SampleGraph.h"

namespace muq {
namespace SamplingAlgorithms {

/// Estimate the probability density function \f$\psi(x)\f$ given samples
/**
<B>Configuration Parameters:</B>
Parameter Key | Type | Default Value | Description |
------------- | ------------- | ------------- | ------------- |
"NumNearestNeighbors"   | <tt>std::size_t</tt> | <tt>25</tt> | The number of nearest neighbors used to compute the bandwidth parameter.   |
"ManifoldDimension"   | <tt>double</tt> | <tt>1.0</tt> | The manifold dimension (if this is not known, then we can estimate it).   |
"SparsityTolerance"   | <tt>double</tt> | <tt>0.1</tt> | The sparsity tolerance for kernel matrix construction (note this may be different than the sparsity tolerance for the optimization).   |
*/
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

  @param[in] epsilon The bandwidth parameter. If we use the automatic tuning, then this is the initial guess for the optimizer.
  @param[in] tune <tt>true</tt> (default): Tune the bandwidth parameter values; <tt>false</tt>: use the stored parameters
  \return The density estimation at each sample \f$\psi^{(i)} \approx \psi(\boldsymbol{x}^{(i)})\f$
  @param[in] estimate dimension  <tt>true</tt>: Use the computed dimension as the manifold dimension and update the stored dimension; <tt>false</tt> (defalt): use the stored manifold dimension
  \return The density estimation at each sample \f$\psi^{(i)} \approx \psi(\boldsymbol{x}^{(i)})\f$
  */
  Eigen::VectorXd EstimateDensity(double epsilon = 1.0, bool const tune = true, bool const tuneDimension = false) const;

protected:

  /// The manifold dimension (estimated if not known)
  mutable double manifoldDim;

  /// The sparsity tolerance for kernel matrix construction (note this may be different than the sparsity tolerance for the optimization).
  const double sparsityTol;

private:

  /// The number of nearest neighbors used to compute the bandwidth parameter
  const std::size_t numNearestNeighbors;
};

} // SamplingAlgorithms
} // namespace muq

#endif
