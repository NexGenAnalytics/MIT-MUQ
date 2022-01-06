#ifndef DENSITYESTIMATION_H_
#define DENSITYESTIMATION_H_

#include "MUQ/SamplingAlgorithms/SampleGraphs/SampleGraph.h"

namespace muq {
namespace SamplingAlgorithms {

/// Estimate the probability density function \f$\psi(x)\f$ given samples
/**
References:
- <a href="https://www.sciencedirect.com/science/article/pii/S1063520315000020">"Variable bandwidth diffusion kernels" by T. Berry & J. Harlim</a>
- <a href="https://www.sciencedirect.com/science/article/pii/S1063520317300982">"Data-driven spectral decomposition and forecasting of ergodic dynamical systems" by D. Giannakis</a>
- <a href="https://arxiv.org/abs/2104.15124">"Graph-theoretic algorithms for Kolmogorov operators: Approximating solutions and their gradients in elliptic and parabolic problems on manifolds" by A.D. Davis & D. Giannakis</a>

<B>Configuration Parameters:</B>
Parameter Key | Type | Default Value | Description |
------------- | ------------- | ------------- | ------------- |
"NumNearestNeighbors"   | <tt>std::size_t</tt> | <tt>25</tt> | The number of nearest neighbors used to compute the bandwidth parameter.   |
"ManifoldDimension"   | <tt>double</tt> | <tt>1.0</tt> | The manifold dimension (if this is not known, then we can estimate it).   |
"SparsityTolerance"   | <tt>double</tt> | <tt>0.1</tt> | The sparsity tolerance for kernel matrix construction (note this may be different than the sparsity tolerance for the optimization).   |
"TuneDimension"   | <tt>bool</tt> | <tt>false</tt> | Use the parameter tuning to tune the manifold dimension   |
*/
class DensityEstimation : public SampleGraph {
public:

  /// Construct the density estimation by sampling a random variable from \f$\psi\f$
  /**
  @param[in] rv The random variable that we wish to sample
  @param[in] options Setup options
  */
  DensityEstimation(std::shared_ptr<muq::Modeling::RandomVariable> const& rv, boost::property_tree::ptree const& options);

  /// Construct the density estimation given samples from the underlying distribution \f$\psi\f$
  /**
  @param[in] samples Samples from the underlying distribution \f$\psi\f$
  @param[in] options Setup options
  */
  DensityEstimation(std::shared_ptr<muq::SamplingAlgorithms::SampleCollection> const& samples, boost::property_tree::ptree const& options);

  /// Construct the density estimation given samples from the underlying distribution \f$\psi\f$
  /**
  @param[in] mat Each column is a sample from the underlying distribution \f$\psi\f$
  @param[in] options Setup options
  */
  DensityEstimation(Eigen::MatrixXd const& mat, boost::property_tree::ptree const& options);

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
  */
  Eigen::VectorXd EstimateDensity(double epsilon = std::numeric_limits<double>::quiet_NaN(), bool const tune = true) const;

  /// Tune the density bandwidth parameter (and the manifold dimension estimate if <tt>tuneDimension</tt> is <tt>true</tt>)
  void TuneDensityBandwidth() const;

  /// The density bandwidth parameter
  double DensityBandwidthParameter() const;

protected:

  /// The manifold dimension (estimated if not known)
  mutable double manifoldDim;

  /// The sparsity tolerance for kernel matrix construction (note this may be different than the sparsity tolerance for the optimization).
  const double sparsityTol;

  /// <tt>true</tt>: Use the computed dimension as the manifold dimension and update the stored dimension; <tt>false</tt> (defalt): use the stored manifold dimension
  const bool tuneDimension;

private:

  /// The bandwidth tuning parameter for the density estimation problem
  mutable double densityBandwidthParameter = 1.0;

  /// The number of nearest neighbors used to compute the bandwidth parameter
  const std::size_t numNearestNeighbors;

  const double powminDens;

  const double powmaxDens;
};

} // SamplingAlgorithms
} // namespace muq

#endif
