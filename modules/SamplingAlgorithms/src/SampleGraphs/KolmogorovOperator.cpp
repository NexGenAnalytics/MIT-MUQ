#include "MUQ/SamplingAlgorithms/SampleGraphs/KolmogorovOperator.h"

#include <Spectra/SymEigsSolver.h>
#include <Spectra/MatOp/SparseGenMatProd.h>

namespace pt = boost::property_tree;
using namespace muq::Modeling;
using namespace muq::SamplingAlgorithms;

KolmogorovOperator::KolmogorovOperator(std::shared_ptr<RandomVariable> const& rv, pt::ptree const& options) :
DensityEstimation(rv, options),
beta(options.get<double>("VariableBandwidth", -0.5)),
operatorParameter(options.get<double>("OperatorParameter", 1.0))
{
  // the normalization parameter
  alpha = 1.0 + beta + (manifoldDim*beta - operatorParameter)/2.0;
}

KolmogorovOperator::KolmogorovOperator(std::shared_ptr<SampleCollection> const& samples, pt::ptree const& options) :
DensityEstimation(samples, options),
beta(options.get<double>("VariableBandwidth", -0.5)),
operatorParameter(options.get<double>("OperatorParameter", 1.0))
{
  // the normalization parameter
  alpha = 1.0 + beta + (manifoldDim*beta - operatorParameter)/2.0;
}

void KolmogorovOperator::DiscreteOperator(Eigen::SparseMatrix<double>& Lhat, Eigen::Ref<Eigen::VectorXd> similarity, double epsilonOperator, double epsilonDensity, bool const tune) const {
  // estimate the density function
  const Eigen::VectorXd psi = EstimateDensity(epsilonDensity, tune);

  // compute the bandwidth function
  Eigen::VectorXd scaledDens = psi.array().pow(beta);

  // if not supplied, use the stored value---the scaling means that the parameter the kernel requries is actually 2*eps
  if( std::isnan(epsilonOperator) ) { epsilonOperator = 2.0*operatorBandwidthParameter; }

  // get the optimal kolmogorov bandwidth parameter
  if( tune ) {
    double dimensionEstimate;
    std::tie(epsilonOperator, dimensionEstimate) = TuneKernelBandwidth(scaledDens, epsilonOperator);

    // update the stored values
    operatorBandwidthParameter = epsilonOperator/2.0;
    if( tuneDimension ) { alpha = 1.0 + beta + (manifoldDim*beta - operatorParameter)/2.0; }
  }

  // compute the unnormalized kernel and store the normalization in the bandwidth vector
  Lhat.resize(NumSamples(), NumSamples());
  similarity = (KernelMatrix(sparsityTol, epsilonOperator, scaledDens, Lhat, !tune).array()*scaledDens.array().pow(-manifoldDim)).pow(alpha);

  // loop through the non-zeros
  {
    std::vector<Eigen::Triplet<double> > entries;
    assert(Lhat.outerSize()==NumSamples());
    for( std::size_t k=0; k<Lhat.outerSize(); ++k ) {
      for( Eigen::SparseMatrix<double>::InnerIterator it(Lhat, k); it; ++it ) { entries.emplace_back(it.row(), it.col(), it.value()/(similarity(it.row())*similarity(it.col()))); }
    }
    Lhat.setFromTriplets(entries.begin(), entries.end());

    // recompute the normalization (in the bandwidth vector)
    similarity = Eigen::VectorXd::Zero(NumSamples());
    for( const auto& entry : entries ) { similarity(entry.row()) += entry.value(); }
  }

  similarity = (scaledDens.array()*similarity.array().sqrt()).matrix();
  Lhat = similarity.array().inverse().matrix().asDiagonal()*Lhat*similarity.array().inverse().matrix().asDiagonal();
  Lhat -= (scaledDens.array()*scaledDens.array()).matrix().asDiagonal();
  Lhat *= epsilonOperator*epsilonOperator/4.0;
}
