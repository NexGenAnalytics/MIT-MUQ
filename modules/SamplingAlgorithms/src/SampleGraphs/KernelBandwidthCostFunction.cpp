#include "MUQ/SamplingAlgorithms/SampleGraphs/KernelBandwidthCostFunction.h"

namespace pt = boost::property_tree;
using namespace muq::Optimization;
using namespace muq::SamplingAlgorithms;

KernelBandwidthCostFunction::KernelBandwidthCostFunction(std::shared_ptr<const SampleGraph> const& graph, Eigen::VectorXd const& bandwidth, pt::ptree const& pt) : CostFunction(Eigen::VectorXi::Ones(1)), graph(graph), bandwidth(bandwidth), delta(pt.get<double>("StepSize", 1.0)), sparsityTol(pt.get<double>("SparsityTolerance", 0.1)) {
  assert(delta>0.0);
}

double KernelBandwidthCostFunction::BandwidthCost(double const para) const {
  auto g = graph.lock();
  assert(g);

  double eps = std::sqrt(std::pow(2.0, para));
  Eigen::SparseMatrix<double> kernel(g->NumSamples(), g->NumSamples());
  const double chi = std::log2(g->KernelMatrix(sparsityTol, eps, bandwidth, kernel, false).sum());

  eps = std::sqrt(std::pow(2.0, para+delta));
  const double chiplus = std::log2(g->KernelMatrix(sparsityTol, eps, bandwidth, kernel, false).sum());

  return (chiplus-chi)/delta;
}

double KernelBandwidthCostFunction::CostImpl(muq::Modeling::ref_vector<Eigen::VectorXd> const& input) { return BandwidthCost(input[0] (0)); }
