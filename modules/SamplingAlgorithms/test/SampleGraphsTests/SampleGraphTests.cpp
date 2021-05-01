#include <gtest/gtest.h>

#include <omp.h>

#include "MUQ/Modeling/Distributions/Gaussian.h"

#include "MUQ/SamplingAlgorithms/SampleGraphs/SampleGraph.h"

namespace pt = boost::property_tree;
using namespace muq::Modeling;
using namespace muq::SamplingAlgorithms;

class SampleGraphTests : public::testing::Test {
public:
  // Set up information to test the nearest neighbor construction
  virtual void SetUp() override {
    // create a standard Gaussian random variable
    rv = std::make_shared<Gaussian>(dim)->AsVariable();

    // set the options for the graph laplacian
    options.put("MaxLeaf", maxLeaf);
    options.put("NumSamples", n);
    options.put("Stride", 14);
  }

  /// Create the nearest neighbor searcher from samples
  /**
  \return The sample collection used to create the graph Laplacian
  */
  inline std::shared_ptr<SampleCollection> CreateFromSamples() {
    // add random samples into a sample collection
    auto samples = std::make_shared<SampleCollection>();
    assert(rv);
    for( std::size_t i=0; i<n; ++i ) { samples->Add(std::make_shared<SamplingState>(rv->Sample())); }

    // create the graph laplacian
    graph = std::make_shared<SampleGraph>(samples, options);

    // return the samples
    return samples;
  }

  /// Create the nearest neighbor searcher from a matrix of samples
  /**
  \return The sample matrix
  */
  inline Eigen::MatrixXd CreateFromMatrix() {
    // add random samples into a sample collection
    Eigen::MatrixXd samples(dim, n);
    assert(rv);
    for( std::size_t i=0; i<n; ++i ) { samples.col(i) = rv->Sample(); }

    // create the graph laplacian
    graph = std::make_shared<SampleGraph>(samples, options);

    // return the samples
    return samples;
  }

  /// Make sure everything is constructed correctly
  virtual void TearDown() override {
    EXPECT_TRUE(graph);

    EXPECT_EQ(graph->NumSamples(), n);
    EXPECT_EQ(graph->StateDimension(), dim);
  }

  /// The sparsity tolerance
  const double sparsityTol = 0.01;

  /// The dimension of state spaces
  const std::size_t dim = 1;

  /// The number of samples
  const std::size_t n = 2500;

  /// The max leaf size for the kd tree
  const std::size_t maxLeaf = 15;

  /// The options for the graph Laplacian
  pt::ptree options;

  /// The random variable that lets us sample from the underlying distribution
  std::shared_ptr<RandomVariable> rv;

  /// The sample graph
  std::shared_ptr<SampleGraph> graph;
};

TEST_F(SampleGraphTests, RandomVariableConstruction) {
  graph = std::make_shared<SampleGraph>(rv, options);

  // choose a random point
  const Eigen::VectorXd x = rv->Sample();

  // find the nearest neighbors to the point
  {
    const double radius = 0.1;
    std::vector<std::pair<std::size_t, double> > neighbors;
    graph->FindNeighbors(x, radius*radius, neighbors);
    for( const auto& neigh : neighbors ) {
      EXPECT_TRUE(neigh.first<n);
      EXPECT_TRUE(neigh.second<radius*radius);
    }
  }

  // find the nearest neighbors to the point with lag
  for( std::size_t lag=0; lag<n; lag+=10 ) {
    const double radius = 1.0;
    std::vector<std::pair<std::size_t, double> > neighbors;
    graph->FindNeighbors(x, radius*radius, neighbors, lag);
    for( const auto& neigh : neighbors ) {
      EXPECT_TRUE(neigh.first>=lag);
      EXPECT_TRUE(neigh.first<n);
      EXPECT_TRUE(neigh.second<radius*radius);
    }
  }

  // find the nearest neighbors to the point
  {
    const std::size_t k = 10;
    std::vector<std::pair<std::size_t, double> > neighbors;
    graph->FindNeighbors(x, k, neighbors);
    EXPECT_EQ(neighbors.size(), k);
    for( const auto& neigh : neighbors ) {
      EXPECT_TRUE(neigh.first<n);
    }
  }

  // find the nearest neighbors to the point with lag
  for( std::size_t lag=0; lag<n; lag+=10 ) {
    const std::size_t k = 10;
    std::vector<std::pair<std::size_t, double> > neighbors;
    graph->FindNeighbors(x, k, neighbors, lag);
    EXPECT_TRUE(neighbors.size()<=k);
    for( const auto& neigh : neighbors ) {
      EXPECT_TRUE(neigh.first>=lag);
      EXPECT_TRUE(neigh.first<n);
    }
  }
}

TEST_F(SampleGraphTests, SampleCollectionConstruction) {
  auto samples = CreateFromSamples();
  EXPECT_TRUE(samples);

  // check to make sure the samples match
  for( std::size_t i=0; i<n; ++i ) { EXPECT_NEAR((samples->at(i)->state[0]-graph->Point(i)).norm(), 0.0, 1.0e-10); }

  // choose a random point
  const Eigen::VectorXd x = rv->Sample();

  // find the nearest neighbors to the point
  {
    const double radius = 0.1;
    std::vector<std::pair<std::size_t, double> > neighbors;
    graph->FindNeighbors(x, radius*radius, neighbors);
    for( const auto& neigh : neighbors ) {
      EXPECT_TRUE(neigh.first<n);
      EXPECT_TRUE(neigh.second<radius*radius);
    }
  }

  // find the nearest neighbors to the point with lag
  for( std::size_t lag=0; lag<n; lag+=10 ) {
    const double radius = 1.0;
    std::vector<std::pair<std::size_t, double> > neighbors;
    graph->FindNeighbors(x, radius*radius, neighbors, lag);
    for( const auto& neigh : neighbors ) {
      EXPECT_TRUE(neigh.first>=lag);
      EXPECT_TRUE(neigh.first<n);
      EXPECT_TRUE(neigh.second<radius*radius);
    }
  }

  // find the nearest neighbors to the point
  {
    const std::size_t k = 10;
    std::vector<std::pair<std::size_t, double> > neighbors;
    graph->FindNeighbors(x, k, neighbors);
    EXPECT_EQ(neighbors.size(), k);
    for( const auto& neigh : neighbors ) {
      EXPECT_TRUE(neigh.first<n);
    }
  }

  // find the nearest neighbors to the point with lag
  for( std::size_t lag=0; lag<n; lag+=10 ) {
    const std::size_t k = 10;
    std::vector<std::pair<std::size_t, double> > neighbors;
    graph->FindNeighbors(x, k, neighbors, lag);
    EXPECT_TRUE(neighbors.size()<=k);
    for( const auto& neigh : neighbors ) {
      EXPECT_TRUE(neigh.first>=lag);
      EXPECT_TRUE(neigh.first<n);
    }
  }
}

TEST_F(SampleGraphTests, MatrixConstruction) {
  const Eigen::MatrixXd samples = CreateFromMatrix();

  // check to make sure the samples match
  for( std::size_t i=0; i<n; ++i ) { EXPECT_NEAR((samples.col(i)-graph->Point(i)).norm(), 0.0, 1.0e-10); }

  // choose a random point
  const Eigen::VectorXd x = rv->Sample();

  // find the nearest neighbors to the point
  {
    const double radius = 0.1;
    std::vector<std::pair<std::size_t, double> > neighbors;
    graph->FindNeighbors(x, radius*radius, neighbors);
    for( const auto& neigh : neighbors ) {
      EXPECT_TRUE(neigh.first<n);
      EXPECT_TRUE(neigh.second<radius*radius);
    }
  }

  // find the nearest neighbors to the point with lag
  for( std::size_t lag=0; lag<n; lag+=10 ) {
    const double radius = 1.0;
    std::vector<std::pair<std::size_t, double> > neighbors;
    graph->FindNeighbors(x, radius*radius, neighbors, lag);
    for( const auto& neigh : neighbors ) {
      EXPECT_TRUE(neigh.first>=lag);
      EXPECT_TRUE(neigh.first<n);
      EXPECT_TRUE(neigh.second<radius*radius);
    }
  }

  // find the nearest neighbors to the point
  {
    const std::size_t k = 10;
    std::vector<std::pair<std::size_t, double> > neighbors;
    graph->FindNeighbors(x, k, neighbors);
    EXPECT_EQ(neighbors.size(), k);
    for( const auto& neigh : neighbors ) {
      EXPECT_TRUE(neigh.first<n);
    }
  }

  // find the nearest neighbors to the point with lag
  for( std::size_t lag=0; lag<n; lag+=10 ) {
    const std::size_t k = 10;
    std::vector<std::pair<std::size_t, double> > neighbors;
    graph->FindNeighbors(x, k, neighbors, lag);
    EXPECT_TRUE(neighbors.size()<=k);
    for( const auto& neigh : neighbors ) {
      EXPECT_TRUE(neigh.first>=lag);
      EXPECT_TRUE(neigh.first<n);
    }
  }
}

TEST_F(SampleGraphTests, SquaredBandwidth) {
  graph = std::make_shared<SampleGraph>(rv, options);

  // compute the squared bandwidth and a random point using 10 nearest neighbors
  const Eigen::VectorXd point = rv->Sample();
  const double squaredBandwidth = graph->SquaredBandwidth(point, 10);

  double expected = 0.0;
  { // compute the expected squared bandwidth
    std::vector<std::pair<std::size_t, double> > neighbors;
    graph->FindNeighbors(point, (std::size_t)10, neighbors);
    for( const auto& it : neighbors ) { expected += it.second; }
  }
  EXPECT_NEAR(squaredBandwidth, expected, 1.0e-10);
}

TEST_F(SampleGraphTests, KernelMatrixComputation) {
  graph = std::make_shared<SampleGraph>(rv, options);

  // the bandwidth parameters
  const double eps = 0.5;
  const Eigen::VectorXd bandwidth = Eigen::VectorXd::Random(graph->NumSamples()).array().abs();

  // the sparsity tolerance
  const double delta = 1.0e-1;

  // compute the kernel matrix
  Eigen::MatrixXd kernel(graph->NumSamples(), graph->NumSamples());
  const Eigen::VectorXd rowsumDense = graph->KernelMatrix(eps, bandwidth, kernel);
  EXPECT_EQ(rowsumDense.size(), graph->NumSamples());
  Eigen::SparseMatrix<double> sparseKernel(graph->NumSamples(), graph->NumSamples());
  const Eigen::VectorXd rowsumSparse = graph->KernelMatrix(delta, eps, bandwidth, sparseKernel);
  EXPECT_EQ(rowsumSparse.size(), graph->NumSamples());

  for( std::size_t i=0; i<graph->NumSamples(); ++i ) {
    double rowsum = 0.0;
    double rowsumSparseExpected = 0.0;
    for( std::size_t j=0; j<graph->NumSamples(); ++j ) {
      const Eigen::VectorXd diff = graph->Point(i)-graph->Point(j);
      EXPECT_NEAR(kernel(i, j), std::exp(-diff.dot(diff)/(eps*bandwidth(i)*bandwidth(j))), 1.0e-10);
      EXPECT_NEAR(kernel(j, i), std::exp(-diff.dot(diff)/(eps*bandwidth(i)*bandwidth(j))), 1.0e-10);
      rowsum += kernel(i, j);

      if( kernel(i, j)>=delta ) {
        EXPECT_NEAR(sparseKernel.coeff(i, j), kernel(i, j), 1.0e-10);
        EXPECT_NEAR(sparseKernel.coeff(j, i), kernel(j, i), 1.0e-10);
        rowsumSparseExpected += kernel(i, j);
      } else {
        EXPECT_TRUE(kernel(i, j)<=delta);
        EXPECT_DOUBLE_EQ(sparseKernel.coeff(i, j), 0.0);
        EXPECT_DOUBLE_EQ(sparseKernel.coeff(j, i), 0.0);
      }
    }
    EXPECT_NEAR(rowsumDense(i), rowsum, 1.0e-10);
    EXPECT_NEAR(rowsumSparse(i), rowsumSparseExpected, 1.0e-10);
  }
}

TEST_F(SampleGraphTests, SharedMemoryKernelMatrixComputation) {
  options.put("NumThreads", omp_get_max_threads());
  graph = std::make_shared<SampleGraph>(rv, options);

  // the bandwidth parameters
  const double eps = 0.5;
  const Eigen::VectorXd bandwidth = Eigen::VectorXd::Random(graph->NumSamples()).array().abs();

  // the sparsity tolerance
  const double delta = 1.0e-1;

  // compute the kernel matrix
  Eigen::MatrixXd kernel(graph->NumSamples(), graph->NumSamples());
  const Eigen::VectorXd rowsumDense = graph->KernelMatrix(eps, bandwidth, kernel);
  EXPECT_EQ(rowsumDense.size(), graph->NumSamples());
  Eigen::SparseMatrix<double> sparseKernel(graph->NumSamples(), graph->NumSamples());
  const Eigen::VectorXd rowsumSparse = graph->KernelMatrix(delta, eps, bandwidth, sparseKernel);
  EXPECT_EQ(rowsumSparse.size(), graph->NumSamples());

  for( std::size_t i=0; i<graph->NumSamples(); ++i ) {
    double rowsum = 0.0;
    double rowsumSparseExpected = 0.0;
    for( std::size_t j=0; j<graph->NumSamples(); ++j ) {
      const Eigen::VectorXd diff = graph->Point(i)-graph->Point(j);
      EXPECT_NEAR(kernel(i, j), std::exp(-diff.dot(diff)/(eps*bandwidth(i)*bandwidth(j))), 1.0e-10);
      EXPECT_NEAR(kernel(j, i), std::exp(-diff.dot(diff)/(eps*bandwidth(i)*bandwidth(j))), 1.0e-10);
      rowsum += kernel(i, j);

      if( kernel(i, j)>=delta ) {
        EXPECT_NEAR(sparseKernel.coeff(i, j), kernel(i, j), 1.0e-10);
        EXPECT_NEAR(sparseKernel.coeff(j, i), kernel(j, i), 1.0e-10);
        rowsumSparseExpected += kernel(i, j);
      } else {
        EXPECT_TRUE(kernel(i, j)<=delta);
        EXPECT_DOUBLE_EQ(sparseKernel.coeff(i, j), 0.0);
        EXPECT_DOUBLE_EQ(sparseKernel.coeff(j, i), 0.0);
      }
    }
    EXPECT_NEAR(rowsumDense(i), rowsum, 1.0e-10);
    EXPECT_NEAR(rowsumSparse(i), rowsumSparseExpected, 1.0e-10);
  }
}


TEST_F(SampleGraphTests, BandwidthParameterTuning) {
  options.put("NumThreads", omp_get_max_threads());
  graph = std::make_shared<SampleGraph>(rv, options);

  // compute the squared bandwidth and a random point using 10 nearest neighbors
  Eigen::VectorXd bandwidth(graph->NumSamples());
  for( std::size_t i=0; i<graph->NumSamples(); ++i ) { bandwidth(i) = std::sqrt(graph->SquaredBandwidth(graph->Point(i), 25)); }

  // tune the bandwidth parameter
  double optPara, cost;
  std::tie(optPara, cost) = graph->TuneKernelBandwidth(bandwidth);
  EXPECT_NEAR(2.0*cost, dim, 1.0e-1);
}
