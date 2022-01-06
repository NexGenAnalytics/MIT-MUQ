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

References:
- <a href="https://www.sciencedirect.com/science/article/pii/S1063520315000020">"Variable bandwidth diffusion kernels" by T. Berry & J. Harlim</a>
- <a href="https://www.sciencedirect.com/science/article/pii/S1063520317300982">"Data-driven spectral decomposition and forecasting of ergodic dynamical systems" by D. Giannakis</a>
- <a href="https://arxiv.org/abs/2104.15124">"Graph-theoretic algorithms for Kolmogorov operators: Approximating solutions and their gradients in elliptic and parabolic problems on manifolds" by A.D. Davis & D. Giannakis</a>

<B>Configuration Parameters:</B>
Parameter Key | Type | Default Value | Description |
------------- | ------------- | ------------- | ------------- |
"VariableBandwidth"   | <tt>double</tt> | <tt>-0.5</tt> | The variable bandwidth parameter \f$\beta\f$---parameterizes the bandwidth function for the unnormalized kernel matrix.   |
"OperatorParameter"   | <tt>double</tt> | <tt>1.0</tt> | The parameter \f$c\f$---determines how much the density function is weighted in the Kolmogorov operator.   |
"NumEigenpairs"   | <tt>std::size_t</tt> | <tt>5 log(n)</tt> | The number of eigenpairs to compute and store---defaults to \f$5 log{(n)}\f$, where \f$n\f$ is the number of samples.   |
*/
class KolmogorovOperator : public DensityEstimation {
public:

  /// Construct the Kolmogorov operator by sampling a random variable from \f$\psi\f$
  /**
  @param[in] rv The random variable that we wish to sample
  @param[in] options Setup options
  */
  KolmogorovOperator(std::shared_ptr<muq::Modeling::RandomVariable> const& rv, boost::property_tree::ptree const& options);

  /// Construct the Kolmogorov operator given samples from the underlying distribution \f$\psi\f$
  /**
  @param[in] samples Samples from the underlying distribution \f$\psi\f$
  @param[in] options Setup options
  */
  KolmogorovOperator(std::shared_ptr<muq::SamplingAlgorithms::SampleCollection> const& samples, boost::property_tree::ptree const& options);

  /// Construct the Kolmogorov operator given samples from the underlying distribution \f$\psi\f$
  /**
  @param[in] mat Each column is a sample from the underlying distribution \f$\psi\f$
  @param[in] options Setup options
  */
  KolmogorovOperator(Eigen::MatrixXd const& mat, boost::property_tree::ptree const& options);

  virtual ~KolmogorovOperator() = default;

  /// Compute the eigendecomposition of the discrete Kolmogorov operator
  /**
  @param[in] density The density \f$\psi\f$ evaluated at each sample
  @param[in] epsilonOperator The bandwidth parameter for the operator estimation. If we use the automatic tuning, then this is the initial guess for the optimizer.
  @param[in] tune <tt>true</tt> (default): Tune the bandwidth parameter values; <tt>false</tt>: use the stored parameters
  \return First: The similarity transformation between \f$\hat{L}\f$ and \f$L\f$, Second: The eigenvalues of the discrete Kolmogorov operator, Third: The eigenvectors of the discrete Kolmogorov operator
  */
  std::tuple<Eigen::VectorXd, Eigen::VectorXd, Eigen::MatrixXd> Eigendecomposition(Eigen::VectorXd const& density, double epsilonOperator = std::numeric_limits<double>::quiet_NaN(), bool const tune = true);

  /// Compute the discrete Kolmogorov operator
  /**
  @param[in] density The density \f$\psi\f$ evaluated at each sample
  @param[out] matrix A sparse matrix to store the discrete Kolmogorov operator \f$\hat{L}\f$ (which is related to \f$L\f$ by a similarity transformation)
  @param[in] epsilon The bandwidth parameter for the operator estimation. If we use the automatic tuning, then this is the initial guess for the optimizer.
  @param[in] tune <tt>true</tt> (default): Tune the bandwidth parameter values; <tt>false</tt>: use the stored parameters
  \return The similarity transformation between \f$\hat{L}\f$ and \f$L\f$
  */
  Eigen::VectorXd DiscreteOperator(Eigen::VectorXd const& density, Eigen::SparseMatrix<double>& matrix, double epsilon = std::numeric_limits<double>::quiet_NaN(), bool const tune = true) const;

  /// Compute the the inner product between gradient vector fields
  /**
  The computed matrix is actually \f$\hat{L}\f$, a symmetric matrix that is related to the discrete Kolmogorov operator by a similarity transformation \f$L = S \hat{L} S^{-1}\f$ (\f$S\f$ is diagonal).
  @param[in] similarity The similarity transformation between \f$\hat{L}\f$ and \f$L\f$
  @param[in] eigvals The eigenvalues of the discrete Kolmogorov operator
  @param[in] eigvecs The eigenvectors of the discrete Kolmogorov operator
  @param[in] func1 A function \f$u\f$ evaluated at each sample
  @param[in] func2 The \f$j^{th}\f$ column is a function \f$v_j\f$ evaluated at each sample
  \return The \f$j^{th}\f$ column is the field \f$\nabla u \cdot \nabla v_j\f$ evaluated at each sample
  */
  Eigen::MatrixXd GradientVectorInnerProduct(Eigen::VectorXd const& similarity, Eigen::VectorXd const& eigvals, Eigen::MatrixXd const& eigvecs, Eigen::VectorXd const& func1, Eigen::MatrixXd const& func2) const;

  /// Compute the the inner product between gradient vector fields
  /**
  @param[in] similarity The similarity transformation between \f$\hat{L}\f$ and \f$L\f$
  @param[in] eigvals The eigenvalues of the discrete Kolmogorov operator
  @param[in] eigvecs The eigenvectors of the discrete Kolmogorov operator
  @param[in] func1 A function \f$u\f$ evaluated at each sample
  \return Each row is the graident \f$\nabla u\f$ evaluated at a sample
  */
  Eigen::MatrixXd GradientVectorField(Eigen::VectorXd const& similarity, Eigen::VectorXd const& eigvals, Eigen::MatrixXd const& eigvecs, Eigen::VectorXd const& func1) const;

  /// Invert the Kolmogorov operator
  /**
  Solve the discrete version of \f$\mathcal{L}_{\psi,c} u = f\f$ with \f$\int_{\Omega} u \psi^c dx = 0\f$.
  @param[in] similarity The similarity transformation between \f$\hat{L}\f$ and \f$L\f$
  @param[in] eigvals The eigenvalues of the discrete Kolmogorov operator
  @param[in] eigvecs The eigenvectors of the discrete Kolmogorov operator
  @param[in] rhs A right hand side function \f$f\f$ evaluated at each sample
  \return The solution \f$u\f$ evaluated at each sample
  */
  Eigen::VectorXd KolmogorovProblemSolution(Eigen::VectorXd const& similarity, Eigen::VectorXd const& eigvals, Eigen::MatrixXd const& eigvecs, Eigen::VectorXd const& func) const;

  /// Tune the density bandwidth parameter (and the parameter \f$\alpha\f$ if <tt>tuneDimension</tt> is <tt>true</tt>)
  /**
  @param[in] density The probability density function (or an approximation of it) evaluated at each sample
  */
  void TuneOperatorBandwidth(Eigen::VectorXd const& density) const;

  /// The operator bandwidth parameter
  double OperatorBandwidthParameter() const;

  /// The parameter \f$c\f$---determines how much the density function is weighted in the Kolmogorov operator.
  const double operatorParameter;

  /// The variable bandwidth parameter \f$\beta\f$---parameterizes the bandwidth function for the unnormalized kernel matrix.
  const double beta;

private:

  /// The second variable bandwidth parameter \f$\alpha\f$---parameterizes the bandwidth function for the unnormalized kernel matrix.
  mutable double alpha;

  /// The number of eigenpairs to compute and store
  const std::size_t neigs;

  /// The bandwidth tuning parameter for the operator estimation problem
  mutable double operatorBandwidthParameter = 1.0e-1;

  const double powminOperator;

  const double powmaxOperator;
};

} // namespace SamplingAlgorithms
} // namespace muq

#endif
