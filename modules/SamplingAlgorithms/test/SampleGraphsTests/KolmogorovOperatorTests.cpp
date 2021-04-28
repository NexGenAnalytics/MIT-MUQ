#include <gtest/gtest.h>

#include <omp.h>

#include "MUQ/Modeling/Distributions/Density.h"
#include "MUQ/Modeling/Distributions/Gaussian.h"

#include "MUQ/SamplingAlgorithms/SampleGraphs/KolmogorovOperator.h"

namespace pt = boost::property_tree;
using namespace muq::Modeling;
using namespace muq::SamplingAlgorithms;

class KolmogorovOperatorTests : public::testing::Test {
public:
  // Set up information to test the nearest neighbor construction
  virtual void SetUp() override {
    // create a standard Gaussian random variable
    rv = std::make_shared<Gaussian>(dim)->AsVariable();

    // set the options for the graph laplacian
    options.put("NumSamples", n);
    options.put("Stride", 14);
    options.put("ManifoldDimension", dim);
    options.put("NumThreads", omp_get_max_threads());
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
    kolmogorov = std::make_shared<KolmogorovOperator>(samples, options);

    // return the samples
    return samples;
  }

  /// Make sure everything is constructed correctly
  virtual void TearDown() override {
    EXPECT_TRUE(kolmogorov);

    EXPECT_EQ(kolmogorov->NumSamples(), n);
    EXPECT_EQ(kolmogorov->StateDimension(), dim);
  }

  /// The dimension of state spaces
  const std::size_t dim = 1;

  /// The number of samples
  const std::size_t n = 10000;

  /// The options for the graph Laplacian
  pt::ptree options;

  /// The random variable that lets us sample from the underlying distribution
  std::shared_ptr<RandomVariable> rv;

  /// The Kolmogorov operator
  std::shared_ptr<KolmogorovOperator> kolmogorov;
};

TEST_F(KolmogorovOperatorTests, RandomVariableConstruction) {
  kolmogorov = std::make_shared<KolmogorovOperator>(rv, options);

  // choose a random point
  const Eigen::VectorXd x = rv->Sample();

  // find the nearest neighbors to the point
  {
    const double radius = 0.1;
    std::vector<std::pair<std::size_t, double> > neighbors;
    kolmogorov->FindNeighbors(x, radius*radius, neighbors);
    for( const auto& neigh : neighbors ) {
      EXPECT_TRUE(neigh.first<n);
      EXPECT_TRUE(neigh.second<radius*radius);
    }
  }

  // find the nearest neighbors to the point with lag
  for( std::size_t lag=0; lag<n; lag+=10 ) {
    const double radius = 1.0;
    std::vector<std::pair<std::size_t, double> > neighbors;
    kolmogorov->FindNeighbors(x, radius*radius, neighbors, lag);
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
    kolmogorov->FindNeighbors(x, k, neighbors);
    EXPECT_EQ(neighbors.size(), k);
    for( const auto& neigh : neighbors ) {
      EXPECT_TRUE(neigh.first<n);
    }
  }

  // find the nearest neighbors to the point with lag
  for( std::size_t lag=0; lag<n; lag+=10 ) {
    const std::size_t k = 10;
    std::vector<std::pair<std::size_t, double> > neighbors;
    kolmogorov->FindNeighbors(x, k, neighbors, lag);
    EXPECT_TRUE(neighbors.size()<=k);
    for( const auto& neigh : neighbors ) {
      EXPECT_TRUE(neigh.first>=lag);
      EXPECT_TRUE(neigh.first<n);
    }
  }
}

TEST_F(KolmogorovOperatorTests, SampleCollectionConstruction) {
  auto samples = CreateFromSamples();
  EXPECT_TRUE(samples);

  // check to make sure the samples match
  for( std::size_t i=0; i<n; ++i ) { EXPECT_NEAR((samples->at(i)->state[0]-kolmogorov->Point(i)).norm(), 0.0, 1.0e-10); }

  // choose a random point
  const Eigen::VectorXd x = rv->Sample();

  // find the nearest neighbors to the point
  {
    const double radius = 0.1;
    std::vector<std::pair<std::size_t, double> > neighbors;
    kolmogorov->FindNeighbors(x, radius*radius, neighbors);
    for( const auto& neigh : neighbors ) {
      EXPECT_TRUE(neigh.first<n);
      EXPECT_TRUE(neigh.second<radius*radius);
    }
  }

  // find the nearest neighbors to the point with lag
  for( std::size_t lag=0; lag<n; lag+=10 ) {
    const double radius = 1.0;
    std::vector<std::pair<std::size_t, double> > neighbors;
    kolmogorov->FindNeighbors(x, radius*radius, neighbors, lag);
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
    kolmogorov->FindNeighbors(x, k, neighbors);
    EXPECT_EQ(neighbors.size(), k);
    for( const auto& neigh : neighbors ) {
      EXPECT_TRUE(neigh.first<n);
    }
  }

  // find the nearest neighbors to the point with lag
  for( std::size_t lag=0; lag<n; lag+=10 ) {
    const std::size_t k = 10;
    std::vector<std::pair<std::size_t, double> > neighbors;
    kolmogorov->FindNeighbors(x, k, neighbors, lag);
    EXPECT_TRUE(neighbors.size()<=k);
    for( const auto& neigh : neighbors ) {
      EXPECT_TRUE(neigh.first>=lag);
      EXPECT_TRUE(neigh.first<n);
    }
  }
}

TEST_F(KolmogorovOperatorTests, DiscreteOperatorComputation) {
  kolmogorov = std::make_shared<KolmogorovOperator>(rv, options);

  Eigen::SparseMatrix<double> Lhat;
  Eigen::VectorXd similarity(kolmogorov->NumSamples());
  kolmogorov->DiscreteOperator(Lhat, similarity);
  EXPECT_EQ(Lhat.rows(), kolmogorov->NumSamples());
  EXPECT_EQ(Lhat.cols(), kolmogorov->NumSamples());

  // make sure the matrix is symmetric
  for( std::size_t i=0; i<kolmogorov->NumSamples(); ++i ) {
    for( std::size_t j=0; j<kolmogorov->NumSamples(); ++j ) {
      EXPECT_NEAR(Lhat.coeff(i, j), Lhat.coeff(j, i), 1.0e-10);
    }
  }
}
