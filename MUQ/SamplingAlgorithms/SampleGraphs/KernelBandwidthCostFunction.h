#ifndef KERNELBANDWIDTHCOSTFUNCTION_H_
#define KERNELBANDWIDTHCOSTFUNCTION_H_

#include "MUQ/Optimization/CostFunction.h"

#include "MUQ/SamplingAlgorithms/SampleGraphs/SampleGraph.h"

namespace muq {
namespace SamplingAlgorithms {

/// The cost function used to tune the bandwidth parameter
/**
References:
- <a href="https://www.sciencedirect.com/science/article/pii/S1063520315000020">"Variable bandwidth diffusion kernels" by T. Berry & J. Harlim</a>
- <a href="https://www.sciencedirect.com/science/article/pii/S1063520317300982">"Data-driven spectral decomposition and forecasting of ergodic dynamical systems" by D. Giannakis</a>
- <a href="https://arxiv.org/abs/2104.15124">"Graph-theoretic algorithms for Kolmogorov operators: Approximating solutions and their gradients in elliptic and parabolic problems on manifolds" by A.D. Davis & D. Giannakis</a>

<B>Configuration Parameters:</B>
Parameter Key | Type | Default Value | Description |
------------- | ------------- | ------------- | ------------- |
"StepSize"   | <tt>double</tt> | <tt>1.0</tt> | The step size parameter \f$\delta\f$.   |
"SparsityTolerance"   | <tt>double</tt> | <tt>0.1</tt> | The sparsity tolerance for the kernel matrix construction.   |
*/
class KernelBandwidthCostFunction : public muq::Optimization::CostFunction {
public:

  /**
  @parma[in] graph The graph that stores the samples used to construct the kernel matrix
  @param[in] bandwidth The bandwidth function evaluated (or approximated) at each sample
  @param[in] pt Options for the cost function
  */
  KernelBandwidthCostFunction(std::shared_ptr<const SampleGraph> const& graph, Eigen::VectorXd const& bandwidth, boost::property_tree::ptree const& pt);

  virtual ~KernelBandwidthCostFunction() = default;

  /// The cost associated with each bandwidth parameter
  /**
  @param[in] para The bandwidth parameter
  \return The cost associated with this bandwidth parameter
  */
  double BandwidthCost(double const para) const;

protected:

  /// The cost associated with each bandwidth parameter
  /**
  @param[in] input There is only one input, the (one dimensional) bandwidth parameter
  \return The cost associated with this bandwidth parameter
  */
  virtual double CostImpl(muq::Modeling::ref_vector<Eigen::VectorXd> const& input) override;

private:

  /// The graph that stores the samples
  std::weak_ptr<const SampleGraph> graph;

  /// The bandwidth function evaluated (or approximated) at each sample
  const Eigen::VectorXd& bandwidth;

  /// The stepsize parameter
  const double delta;

  /// The sparsity tolerance for the kernel matrix
  const double sparsityTol;
};

} // SamplingAlgorithms
} // namespace muq

#endif
