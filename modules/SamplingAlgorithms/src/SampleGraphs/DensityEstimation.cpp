#include "MUQ/SamplingAlgorithms/SampleGraphs/DensityEstimation.h"

namespace pt = boost::property_tree;
using namespace muq::Modeling;
using namespace muq::SamplingAlgorithms;

DensityEstimation::DensityEstimation(std::shared_ptr<RandomVariable> const& rv, pt::ptree const& options) : SampleGraph(rv, options),
manifoldDim(options.get<double>("ManifoldDimension", 1.0)),
numNearestNeighbors(options.get<std::size_t>("NumNearestNeighbors", 25)),
sparsityTol(options.get<double>("SparsityTolerance", 0.1))
{}

DensityEstimation::DensityEstimation(std::shared_ptr<SampleCollection> const& samples, pt::ptree const& options) : SampleGraph(samples, options),
manifoldDim(options.get<double>("ManifoldDimension", 1.0)),
numNearestNeighbors(options.get<std::size_t>("NumNearestNeighbors", 25)),
sparsityTol(options.get<double>("SparsityTolerance", 0.1))
{}

Eigen::VectorXd DensityEstimation::EstimateDensity(double epsilon, bool const tune, bool const tuneDimension) const {
  // compute the squared bandwidth and a random point using 10 nearest neighbors
  Eigen::VectorXd bandwidth(NumSamples());
  for( std::size_t i=0; i<NumSamples(); ++i ) { bandwidth(i) = std::sqrt(SquaredBandwidth(Point(i), numNearestNeighbors)); }

  // get the optimal bandwidth parameter
  if( tune ) {
    double dimensionEstimate;
    std::tie(epsilon, dimensionEstimate) = TuneKernelBandwidth(bandwidth, epsilon);
    if( tuneDimension ) { manifoldDim = 2.0*dimensionEstimate; }
  }

  // compute the kernel
  Eigen::SparseMatrix<double> kernel(NumSamples(), NumSamples());
  KernelMatrix(sparsityTol, epsilon, bandwidth, kernel, !tune);

  // comptue the normalization (store it in the bandwidth vector since we don't need it anymore)
  for( std::size_t i=0; i<NumSamples(); ++i ) { bandwidth(i) = NumSamples()*std::pow(M_PI*epsilon*epsilon*bandwidth(i)*bandwidth(i), manifoldDim/2.0); }

  // compute the density estimate
  return bandwidth.array().inverse().matrix().asDiagonal()*kernel*Eigen::VectorXd::Ones(NumSamples());
}
