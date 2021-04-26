#include "MUQ/Approximation/SampleGraphs/DensityEstimation.h"

namespace pt = boost::property_tree;
using namespace muq::Modeling;
using namespace muq::SamplingAlgorithms;
using namespace muq::Approximation;

DensityEstimation::DensityEstimation(std::shared_ptr<RandomVariable> const& rv, pt::ptree const& options) : SampleGraph(rv, options) {}

DensityEstimation::DensityEstimation(std::shared_ptr<SampleCollection> const& samples, pt::ptree const& options) : SampleGraph(samples, options) {}

Eigen::VectorXd DensityEstimation::Estimate(bool const tune) {
  return Eigen::VectorXd();
}
