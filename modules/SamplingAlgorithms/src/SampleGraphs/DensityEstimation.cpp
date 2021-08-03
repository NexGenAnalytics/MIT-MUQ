#include "MUQ/SamplingAlgorithms/SampleGraphs/DensityEstimation.h"

namespace pt = boost::property_tree;
using namespace muq::Modeling;
using namespace muq::SamplingAlgorithms;

DensityEstimation::DensityEstimation(std::shared_ptr<RandomVariable> const& rv, pt::ptree const& options) :
SampleGraph(rv, options),
manifoldDim(options.get<double>("ManifoldDimension", 1.0)),
numNearestNeighbors(options.get<std::size_t>("NumNearestNeighbors", 25)),
sparsityTol(options.get<double>("SparsityTolerance", 1.0e-1)),
tuneDimension(options.get<bool>("TuneDimension", false))
{}

DensityEstimation::DensityEstimation(std::shared_ptr<SampleCollection> const& samples, pt::ptree const& options) :
SampleGraph(samples, options),
manifoldDim(options.get<double>("ManifoldDimension", 1.0)),
numNearestNeighbors(options.get<std::size_t>("NumNearestNeighbors", 25)),
sparsityTol(options.get<double>("SparsityTolerance", 1.0e-1)),
tuneDimension(options.get<bool>("TuneDimension", false))
{}

DensityEstimation::DensityEstimation(Eigen::MatrixXd const& mat, pt::ptree const& options) :
SampleGraph(mat, options),
manifoldDim(options.get<double>("ManifoldDimension", 1.0)),
numNearestNeighbors(options.get<std::size_t>("NumNearestNeighbors", 25)),
sparsityTol(options.get<double>("SparsityTolerance", 1.0e-1)),
tuneDimension(options.get<bool>("TuneDimension", false))
{}

double DensityEstimation::DensityBandwidthParameter() const { return densityBandwidthParameter; }

void DensityEstimation::TuneDensityBandwidth() const {
  // compute the squared bandwidth and a random point using 10 nearest neighbors
  Eigen::VectorXd bandwidth(NumSamples());
  for( std::size_t i=0; i<NumSamples(); ++i ) { bandwidth(i) = std::sqrt(SquaredBandwidth(Point(i), numNearestNeighbors)); }

  double dimensionEstimate;
  std::tie(densityBandwidthParameter, dimensionEstimate) = TuneKernelBandwidth(bandwidth, densityBandwidthParameter);
  if( tuneDimension ) { manifoldDim = 2.0*dimensionEstimate; }
}

Eigen::VectorXd DensityEstimation::EstimateDensity(double epsilon, bool const tune) const {
  // compute the squared bandwidth and a random point using 10 nearest neighbors
  Eigen::VectorXd bandwidth(NumSamples());
  for( std::size_t i=0; i<NumSamples(); ++i ) { bandwidth(i) = std::sqrt(SquaredBandwidth(Point(i), numNearestNeighbors)); }

  // if not supplied, use the stored value
  if( std::isnan(epsilon) ) { epsilon = densityBandwidthParameter; }

  // get the optimal bandwidth parameter
  if( tune ) {
    double dimensionEstimate;
    std::tie(epsilon, dimensionEstimate) = TuneKernelBandwidth(bandwidth, epsilon);

    // update the stored values
    densityBandwidthParameter = epsilon;
    if( tuneDimension ) { manifoldDim = 2.0*dimensionEstimate; }
  }

  // compute the kernel
  Eigen::SparseMatrix<double> kernel(NumSamples(), NumSamples());
  assert(epsilon>1.0e-14);
  //assert(bandwidth.norm()>1.0e-14);
  KernelMatrix(sparsityTol, epsilon, bandwidth, kernel, !tune);

  // comptue the normalization (store it in the bandwidth vector since we don't need it anymore)
  for( std::size_t i=0; i<NumSamples(); ++i ) { bandwidth(i) = NumSamples()*std::pow(M_PI*epsilon*bandwidth(i)*bandwidth(i), manifoldDim/2.0); }

  // compute the density estimate
  return kernel*bandwidth.array().inverse().matrix();
}
