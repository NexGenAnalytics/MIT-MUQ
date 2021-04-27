#ifndef SAMPLEGRAPH_H_
#define SAMPLEGRAPH_H_

#include <boost/property_tree/ptree.hpp>

#include <Eigen/Sparse>

#include <nanoflann.hpp>

#include "MUQ/Modeling/Distributions/RandomVariable.h"

#include "MUQ/SamplingAlgorithms/SampleCollection.h"

namespace muq {
namespace Approximation {

/// Given samples from a distribution, create a graph that connects nearest neighbors
/**
Given \f$\{x_i\}_{i=1}^{n}\f$ samples from the distribution \f$\psi\f$, this class finds the
\f$k\f$ nearest neighbors \f$\{x_j\}_{j=1}^{k}\f$ to a given point \f$x\f$.

<B>Configuration Parameters:</B>
Parameter Key | Type | Default Value | Description |
------------- | ------------- | ------------- | ------------- |
"NumSamples"   | <tt>std::size_t</tt> | - | In the case where a random variable is passed to the constructor, we draw \f$n\f$ samples from the distribution.   |
"MaxLeaf"   | <tt>std::size_t</tt> | <tt>10</tt> | The maximum leaf size for the kd tree (nanoflann parameter). |
"Stride"   | <tt>std::size_t</tt> | <tt>1</tt> | Build \f$i \in [0, m-1]\f$ \f$k\f$-\f$d\f$ trees that ignore the first \f$i d\f$ samples, where \f$d = n/(i+1)\f$ (the stride parameter is number \f$i\f$). |
"NumThreads"   | <tt>std::size_t</tt> | <tt>1</tt> | The number of <tt>openMP</tt> threads available to this object. |
*/
class SampleGraph {
public:

  /// Construct the sample graph by sampling a random variable from \f$\psi\f$
  /**
  @param[in] rv The random variable that we wish to sample
  @param[in] options Setup options
  */
  SampleGraph(std::shared_ptr<muq::Modeling::RandomVariable> const& rv, boost::property_tree::ptree const& options);

  /// Construct the sample graph given samples from the underlying distribution \f$\psi\f$
  /**
  @param[in] samples Samples from the underlying distribution \f$\psi\f$
  @param[in] options Setup options
  */
  SampleGraph(std::shared_ptr<muq::SamplingAlgorithms::SampleCollection> const& samples, boost::property_tree::ptree const& options);

  virtual ~SampleGraph() = default;

  /// Construct the \f$k\f$-\f$d\f$ trees
  void BuildKDTrees() const;

  /// Construct the \f$k\f$-\f$d\f$ trees
  /**
  Before (re-)building the \f$k\f$-\f$d\f$ trees, reorder the points in the cloud so that the ones with the largest bandwidth come first. This makes seaching for the nearest neighbors based on the bandwidth radius more efficient.
  @param[in] bandwidth The bandwidth function at each sample
  */
  void BuildKDTrees(Eigen::VectorXd const& bandwidth) const;

  /// Get the \f$i^{th}\f$ sample
  /**
  @param[in] i We want this sample
  \return Get the \f$i^{th}\f$ sample
  */
  Eigen::VectorXd Point(std::size_t const i) const;

  /// The number of samples that make up the graph
  std::size_t NumSamples() const;

  /// The dimension of the sample state space
  /**
  All the samples are the same size, so just return the dimension of the first one.
  */
  std::size_t StateDimension() const;

  /// Find the neighbors within a specified squared radius
  /**
  @param[in] point We want the nearest neighbors to this point
  @param[in] radius2 Find all of the points with this this squared radius
  @param[out] neighbors Each component is first: the index of the nearest neighbor and second: the squared distance to the nearest neighbor
  @param[in] lag Ignore the first <tt>lag</tt> points (defaults to \f$0\f$)
  */
  void FindNeighbors(Eigen::VectorXd const& point, double const radius2, std::vector<std::pair<std::size_t, double> >& neighbors, std::size_t const& lag = 0) const;

  /// Find the closest \f$k\f$ neighbors
  /**
  @param[in] point We want the nearest neighbors to this point
  @param[in] k Find this many nearest neighbors
  @param[out] neighbors Each component is first: the index of the nearest neighbor and second: the squared distance to the nearest neighbor
  @param[in] lag Ignore the first <tt>lag</tt> points (defaults to \f$0\f$)
  */
  void FindNeighbors(Eigen::VectorXd const& point, std::size_t const k, std::vector<std::pair<std::size_t, double> >& neighbors, std::size_t const& lag = 0) const;

  /// Find the squared bandwidth parameter for each sample
  /**
  The squared bandwidth is defined as \f$b^2(x) = \sum_{j=1}^{k} \| x-x^{(I(i,j))} \|^2\f$ for each sample, where \f$I(x,j)\f$ is the index of the \f$j^{th}\f$ closest neighbor to \f$x\f$---note that \f$I(x,0)\f$ is the nearest neighbor.
  @param[in] x The point \f$x\f$
  @param[in] k The number of nearest neighbors (\f$k\f$)
  \return The squared bandwidth \f$b^2(x)\f$
  */
  double SquaredBandwidth(Eigen::VectorXd const& x, std::size_t const k) const;

  /// Compute the kernel matrix using the brute-force algorithm
  /**
  Compute the kernel using nested <tt>for</tt> loops. Since we are going to store a dense matrix, don't bother setting entries to zero if they are below a sparsity threshold. This method is impractical and should only be used for testing or toy examples.

  The \f$(i,j)\f$ entry of the kernel matrix is
  \f{equation*}{
  K_{ij} = \exp{\left( -\frac{\|x_i - x_j\|^2}{\epsilon^{2} b(x_i)b(x_j)} \right)}
  \f}
  @param[in] bandwidthParameter The bandwidth parameter \f$\epsilon\f$.
  @param[in] bandwidth The bandwidth function evaluated (or approximated) at each sample \f$b(x_i)\f$.
  @param[out] kernel The kernel matrix \f$K\f$
  \return The rowsum of the kernel matrix \f$\sum_{j=1}^{n} K_{ij}\f$
  */
  Eigen::VectorXd KernelMatrix(double bandwidthParameter, Eigen::VectorXd const& bandwidth, Eigen::Ref<Eigen::MatrixXd> kernel) const;

  /// Compute the kernel matrix using the \f$k\f$-\f$d\f$ trees
  /**
  The \f$(i,j)\f$ entry of the kernel matrix is
  \f{equation*}{
  K_{ij} = \exp{\left( -\frac{\|x_i - x_j\|^2}{\epsilon^{2} b(x_i)b(x_j)} \right)}
  \f}
  However, if the entry is less than the sparsity tolerance, we set it equal to zero.
  @param[in] sparsityTol The sparsity tolerance
  @param[in] bandwidthParameter The bandwidth parameter \f$\epsilon\f$.
  @param[in] bandwidth The bandwidth function evaluated (or approximated) at each sample \f$b(x_i)\f$.
  @param[out] kernel The kernel matrix \f$K\f$
  \return The rowsum of the kernel matrix \f$\sum_{j=1}^{n} K_{ij}\f$
  */
  Eigen::VectorXd KernelMatrix(double const sparsityTol, double bandwidthParameter, Eigen::VectorXd const& bandwidth, Eigen::SparseMatrix<double>& kernel) const;

private:

  /// Create a sample collection by sampling a random variable
  /**
  @param[in] rv The random variable that we wish to sample
  @param[in] n Sample the random variable \f$n\f$ times
  */
  static std::shared_ptr<muq::SamplingAlgorithms::SampleCollection> SampleRandomVariable(std::shared_ptr<muq::Modeling::RandomVariable> const& rv, std::size_t const n);

  /// Initialize the sample graph object
  /**
  @param[in] options Setup options
  */
  void Initialize(boost::property_tree::ptree const& options);

  /// Evaluate the kernel
  /**
  \f{equation*}{
  K(\epsilon, x_i, x_j) = \exp{\left( -\frac{\|x_i - x_j\|^2}{\epsilon} \right)}
  \f}
  @param[in] scale The scale parameter \f$\epsilon\f$
  @param[in] xi The point \f$x_i\f$
  @param[in] xj The point \f$x_j\f$
  \return The kernel evaluation
  */
  double Kernel(double const scale, Eigen::VectorXd const& xi, Eigen::VectorXd const& xj) const;

  /// A point cloud used by nanoflann to find the nearest neighbors
  class Cloud {
  public:
    /**
    @param[in] samples Samples from the underlying distribution \f$\psi\f$
    @param[in] indices The order of the samples (sorted by the bandwidth)
    @param[in] lag Ignore the first <tt>lag</tt> samples
    */
    Cloud(std::shared_ptr<const muq::SamplingAlgorithms::SampleCollection> const& samples, std::vector<std::size_t> const& indices, std::size_t const lag);

    virtual ~Cloud() = default;

    /// Get the number of samples (from \f$\psi\f$)
    /**
    \return The number of points in the sample collection (GraphLaplacian::PointCloud::samples)
    */
    std::size_t kdtree_get_point_count() const;

    /// Get the \f$i^{th}\f$ component of the \f$p^{th}\f$ point
    /**
    @param[in] p We want to access this particle number
    @param[in] i We want this index of the particle location
    \return The \f$i^{th}\f$ component of the \f$p^{th}\f$ point in GraphLaplacian::PointCloud::samples
    */
    double kdtree_get_pt(std::size_t const p, std::size_t const i) const;

    /// Optional bounding-box computation
    /**
    \return Return <tt>false</tt> to default to a standard bounding box computation loop
    */
    template<class BBOX>
    inline bool kdtree_get_bbox(BBOX& bb) const { return false; }

    /// Reorder the samples based on the bandwidth parameter
    /**
    We want the samples with the largest bandwidth to come first.
    @param[in] bandwidth The sample bandwidth
    */
    void ReorderSamples(Eigen::VectorXd const& bandwidth);

    /// Ignore the first <tt>lag</tt> samples
    const std::size_t lag;

  private:

    /// Samples from the distribution \f$\psi\f$
    std::shared_ptr<const muq::SamplingAlgorithms::SampleCollection> samples;

    /// The order of the samples in the \f$k\f$-\f$d\f$ tree
    const std::vector<std::size_t>& indices;
  };

  /// The order of the samples in the \f$k\f$-\f$d\f$ tree
  mutable std::vector<std::size_t> indices;

  /// The point clouds
	std::vector<Cloud> clouds;

  /// The nanoflann kd-tree type
  typedef nanoflann::KDTreeSingleIndexAdaptor<nanoflann::L2_Simple_Adaptor<double, Cloud>, Cloud> NanoflannKDTree;

  /// The nanoflann kd-trees
	std::vector<std::shared_ptr<NanoflannKDTree> > kdtrees;

  /// Samples from the distribution \f$\psi\f$
  std::shared_ptr<const muq::SamplingAlgorithms::SampleCollection> samples;

  /// The number of <tt>openMP</tt> threads available to this object.
  const std::size_t numThreads;
};

} // Approximation
} // namespace muq

#endif
