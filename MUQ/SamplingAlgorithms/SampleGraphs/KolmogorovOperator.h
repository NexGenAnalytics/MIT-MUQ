#ifndef KOLMOGOROVOPERATOR_H_
#define KOLMOGOROVOPERATOR_H_

#include "MUQ/SamplingAlgorithms/SampleGraphs/DensityEstimation.h"

namespace muq {
namespace SamplingAlgorithms {

/// Estimate an elliptic Kolmogorov operator of the form \f$\mathcal{L}_{\psi, c} f = \Delta f + c \nabla f \cdot \frac{\nabla \psi}{\psi}\f$
/**
Let \f$\Omega \subseteq \mathbb{R}^{m}\f$ be a smooth, \f$d\f$-dimensional Riemannian manifold without boundary---note that \f$m\f$ is the dimension of the ambient space and \f$d\f$ is the intrinsic manifold dimension. Suppose that \f$\Omega\f$ is equipped with a probability measure with a smooth density \f$\psi: \Omega \to \mathbb{R}^{+}\f$ relative to the volume form of \f$\Omega\f$. Given \f$c \geq 0\f$, let \f$\mathcal{H}_{\psi, c}(\Omega)\f$ be the Hilbert space of real-valued functions on \f$\Omega\f$ whose inner product is defined by \f$\psi\f$, \f$\langle f, g \rangle_{\psi,c} = \int_{\Omega} f(x) g(x) \psi^c(x) \, dx\f$. We restrict choices of \f$c\f$ such that \f$\int_{\Omega} \psi^{c}(x) \, d x < \infty\f$---clearly, \f$c=1\f$ is always valid since \f$\psi\f$ is a probability density function. The user supplies \f$c\f$, and we assume it is valid. Given a smooth function \f$f \in \mathcal{H}_{\psi, c}(\Omega)\f$, we focus an elliptic Kolmogorov operators of the form
\f{equation*}{
    \mathcal{L}_{\psi, c} f = \Delta f + c \nabla f \cdot \frac{\nabla \psi}{\psi},
    \label{eq:Kolmogorov}
\f}
where \f$\Delta\f$ and \f$\nabla\f$ denote the Laplacian and gradient operators, respectively, and the dot operator \f$\cdot\f$ denotes the Riemannian inner product between tangent vectors on \f$\Omega\f$.

<B>Configuration Parameters:</B>
Parameter Key | Type | Default Value | Description |
------------- | ------------- | ------------- | ------------- |
"VariableBandwidth"   | <tt>double</tt> | <tt>-0.5</tt> | The variable bandwidth parameter \f$\beta\f$---parameterizes the bandwidth function for the unnormalized kernel matrix.   |
"OperatorParameter"   | <tt>double</tt> | <tt>1.0</tt> | The parameter \f$c\f$---determines how much the density function is weighted in the Kolmogorov operator.   |
*/
class KolmogorovOperator : public DensityEstimation {
public:

  /// Construct the sample graph by sampling a random variable from \f$\psi\f$
  /**
  @param[in] rv The random variable that we wish to sample
  @param[in] options Setup options
  */
  KolmogorovOperator(std::shared_ptr<muq::Modeling::RandomVariable> const& rv, boost::property_tree::ptree const& options);

  /// Construct the sample graph given samples from the underlying distribution \f$\psi\f$
  /**
  @param[in] samples Samples from the underlying distribution \f$\psi\f$
  @param[in] options Setup options
  */
  KolmogorovOperator(std::shared_ptr<muq::SamplingAlgorithms::SampleCollection> const& samples, boost::property_tree::ptree const& options);

  virtual ~KolmogorovOperator() = default;

  /// Compute the discrete Kolmogorov operator
  /**
  @param[out] Lhat The symmetirc discrete Kolmogorov operator---this is related to the discrete Kolmogorov operator by a similarity transformation
  @param[out] similarity The diagonal of the matrix that defines the similarity transformation
  @param[in] epsilonOperator The bandwidth parameter for the operator estimation. If we use the automatic tuning, then this is the initial guess for the optimizer.
  @param[in] epsilonDensity The bandwidth parameter for the density estimation. If we use the automatic tuning, then this is the initial guess for the optimizer.
  @param[in] tune <tt>true</tt> (default): Tune the bandwidth parameter values; <tt>false</tt>: use the stored parameters
  */
  void DiscreteOperator(Eigen::SparseMatrix<double>& Lhat, Eigen::Ref<Eigen::VectorXd> similarity, double epsilonOperator = std::numeric_limits<double>::quiet_NaN(), double epsilonDensity = std::numeric_limits<double>::quiet_NaN(), bool const tune = true) const;
private:

private:

  /// The variable bandwidth parameter \f$\beta\f$---parameterizes the bandwidth function for the unnormalized kernel matrix.
  const double beta;

  /// The parameter \f$c\f$---determines how much the density function is weighted in the Kolmogorov operator.
  const double operatorParameter;

  /// The second variable bandwidth parameter \f$\alpha\f$---parameterizes the bandwidth function for the unnormalized kernel matrix.
  mutable double alpha;

  /// The bandwidth tuning parameter for the operator estimation problem
  mutable double operatorBandwidthParameter = 1.0e-4;

};

} // namespace SamplingAlgorithms
} // namespace muq

#endif
