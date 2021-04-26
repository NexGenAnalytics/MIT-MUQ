#ifndef SAMPLEGRAPH_H_
#define SAMPLEGRAPH_H_

#include <boost/property_tree/ptree.hpp>

#include <nanoflann.hpp>

#include <MUQ/Modeling/Distributions/RandomVariable.h>

#include <MUQ/SamplingAlgorithms/SampleCollection.h>

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

  // Construct the \f$k\f$-\f$d\f$-trees
  void BuildKDTrees() const;

  /// Get the \f$i^{th}\f$ sample
  /**
  @param[in] i We want this sample
  \return Get the \f$i^{th}\f$ sample
  */
  Eigen::VectorXd Point(std::size_t const i) const;

  /// Get the \f$i^{th}\f$ sample
  /**
  @param[in] i We want this sample
  \return Get the \f$i^{th}\f$ sample
  */
  Eigen::Ref<Eigen::VectorXd> Point(std::size_t const i);

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
  void FindNeighbors(Eigen::Ref<const Eigen::VectorXd> const& point, double const radius2, std::vector<std::pair<std::size_t, double> >& neighbors, std::size_t const& lag = 0) const;

  /// Find the closest \f$k\f$ neighbors
  /**
  @param[in] point We want the nearest neighbors to this point
  @param[in] k Find this many nearest neighbors
  @param[out] neighbors Each component is first: the index of the nearest neighbor and second: the squared distance to the nearest neighbor
  @param[in] lag Ignore the first <tt>lag</tt> points (defaults to \f$0\f$)
  \return The average squared distnace from the input point to the nearest neighbors
  */
  double FindNeighbors(Eigen::Ref<const Eigen::VectorXd> const& point, std::size_t const k, std::vector<std::pair<std::size_t, double> >& neighbors, std::size_t const& lag = 0) const;

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

  /// A point cloud used by nanoflann to find the nearest neighbors
  class Cloud {
  public:
    /**
    @param[in] samples Samples from the underlying distribution \f$\psi\f$
    @param[in] lag Ignore the first <tt>lag</tt> samples
    */
    Cloud(std::shared_ptr<const muq::SamplingAlgorithms::SampleCollection> const& samples, std::size_t const lag);

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

    /// Ignore the first <tt>lag</tt> samples
    const std::size_t lag;

  private:

    /// Samples from the distribution \f$\psi\f$
    std::shared_ptr<const muq::SamplingAlgorithms::SampleCollection> samples;
  };

  /// The point clouds
	std::vector<Cloud> clouds;

  /// The nanoflann kd-tree type
  typedef nanoflann::KDTreeSingleIndexAdaptor<nanoflann::L2_Simple_Adaptor<double, Cloud>, Cloud> NanoflannKDTree;

  /// The nanoflann kd-trees
	std::vector<std::shared_ptr<NanoflannKDTree> > kdtrees;

  /// Samples from the distribution \f$\psi\f$
  std::shared_ptr<const muq::SamplingAlgorithms::SampleCollection> samples;
};

} // Approximation
} // namespace muq

#endif
