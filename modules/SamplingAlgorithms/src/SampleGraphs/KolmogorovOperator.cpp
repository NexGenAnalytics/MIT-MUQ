#include "MUQ/SamplingAlgorithms/SampleGraphs/KolmogorovOperator.h"

#include <Spectra/SymEigsSolver.h>
#include <Spectra/GenEigsSolver.h>
#include <Spectra/MatOp/SparseGenMatProd.h>
#include <Spectra/Util/SelectionRule.h>

namespace pt = boost::property_tree;
using namespace Spectra;
using namespace muq::Modeling;
using namespace muq::SamplingAlgorithms;

KolmogorovOperator::KolmogorovOperator(std::shared_ptr<RandomVariable> const& rv, pt::ptree const& options) :
DensityEstimation(rv, options),
beta(options.get<double>("VariableBandwidth", -0.5)),
operatorParameter(options.get<double>("OperatorParameter", 1.0)),
neigs(options.get<std::size_t>("NumEigenpairs", 5*std::log((double)NumSamples())))
{
  // the normalization parameter
  alpha = 1.0 + beta + (manifoldDim*beta - operatorParameter)/2.0;
}

KolmogorovOperator::KolmogorovOperator(std::shared_ptr<SampleCollection> const& samples, pt::ptree const& options) :
DensityEstimation(samples, options),
beta(options.get<double>("VariableBandwidth", -0.5)),
operatorParameter(options.get<double>("OperatorParameter", 1.0)),
neigs(options.get<std::size_t>("NumEigenpairs", 5*std::log((double)NumSamples())))
{
  // the normalization parameter
  alpha = 1.0 + beta + (manifoldDim*beta - operatorParameter)/2.0;
}

KolmogorovOperator::KolmogorovOperator(Eigen::MatrixXd const& mat, pt::ptree const& options) :
DensityEstimation(mat, options),
beta(options.get<double>("VariableBandwidth", -0.5)),
operatorParameter(options.get<double>("OperatorParameter", 1.0)),
neigs(options.get<std::size_t>("NumEigenpairs", 5*std::log((double)NumSamples())))
{
  // the normalization parameter
  alpha = 1.0 + beta + (manifoldDim*beta - operatorParameter)/2.0;
}

void KolmogorovOperator::DiscreteOperator(Eigen::VectorXd const& density, double epsilon, bool const tune) {
  // compute the bandwidth function
  const Eigen::VectorXd scaledDens = density.array().pow(beta);

  //std::cout << scaledDens.transpose() << std::endl;

  // if not supplied, use the stored value---the scaling means that the parameter the kernel requries is actually 2*eps
  epsilon = (std::isnan(epsilon)? 4.0*operatorBandwidthParameter : 4.0*epsilon);

  // get the optimal kolmogorov bandwidth parameter
  if( tune ) {
    std::cout << "TUNING DISCRETE OPERATOR" << std::endl;
    double dimensionEstimate;
    std::tie(epsilon, dimensionEstimate) = TuneKernelBandwidth(scaledDens, epsilon);
    std::cout << "operator dimension estimate: " << 2.0*dimensionEstimate << std::endl;

    // update the stored values
    operatorBandwidthParameter = epsilon/4.0;
    if( tuneDimension ) { alpha = 1.0 + beta + (manifoldDim*beta - operatorParameter)/2.0; }
  }

  std::cout << "epsilon: " << epsilon/4.0 << std::endl;
  //std::cout << "epsilon2: " << epsilon*epsilon/4.0 << std::endl;
  std::cout << "beta: " << beta << std::endl;
  std::cout << "alpha: " << alpha << std::endl;
  std::cout << "dim: " << manifoldDim << std::endl;

  // compute the unnormalized kernel and store the normalization in the bandwidth vector
  std::vector<Eigen::Triplet<double> > entries;
  Eigen::VectorXd rowsum = KernelMatrix(sparsityTol, epsilon, scaledDens, entries, !tune);
  similarity = rowsum.array()/scaledDens.array().pow(beta*manifoldDim);

  //std::cout << similarity.transpose() << std::endl;

  // loop through the non-zeros
  #pragma omp parallel num_threads(numThreads)
  for( auto& entry : entries ) { entry = Eigen::Triplet<double>(entry.row(), entry.col(), entry.value()/std::pow(similarity(entry.row())*similarity(entry.col()), alpha)); }

  Eigen::SparseMatrix<double> kernel(NumSamples(), NumSamples());
  //matrix.resize(NumSamples(), NumSamples());
  //matrix.setZero();
  kernel.setFromTriplets(entries.begin(), entries.end());

  // recompute the normalization (in the bandwidth vector)
  rowsum = Eigen::VectorXd::Zero(NumSamples());
  for( const auto& entry : entries ) { rowsum(entry.row()) += entry.value(); }

  similarity = scaledDens.array()*rowsum.array().sqrt();
  rowsum = similarity.array().inverse(); // hold the inverse for computational efficiency
  //matrix = (scaledDens.array()*scaledDens.array()).inverse().matrix().asDiagonal();

  //std::cout << rowsum.transpose() << std::endl;

  //std::cout << "BANDWIDTH: " << epsilon/4.0 << std::endl;

  //Eigen::MatrixXd mat(kernel);
  //std::cout << kernel << std::endl;

  //matrix = (rowsum.asDiagonal()*kernel*rowsum.asDiagonal()-matrix)/(epsilon/4.0);
  matrix = rowsum.asDiagonal()*kernel*rowsum.asDiagonal();
  matrix -= (scaledDens.array()*scaledDens.array()).inverse().matrix().asDiagonal();
  matrix /= epsilon/4.0;
  //matrix /= epsilon;
}

#include <Eigen/Eigenvalues>

std::tuple<Eigen::VectorXd, Eigen::VectorXd, Eigen::MatrixXd> KolmogorovOperator::Eigendecomposition(Eigen::VectorXd const& density, double epsilon, bool const tune) {
  //std::cout << density << std::endl;

  // compute and store the discrete Kolmogorov operator
  DiscreteOperator(density, epsilon, tune);



  //Spectra::SymEigsSolver<double, Spectra::SMALLEST_MAGN, Spectra::SparseGenMatProd<double> > eigsolverf;



  //std::cout << similarity.transpose() << std::endl;

  //std::cout << "neigs: " << neigs << std::endl;

  // compute the eigen decomposition
  Spectra::SparseGenMatProd<double> matrixvec(matrix);
  Spectra::SymEigsSolver<Spectra::SparseGenMatProd<double> > eigsolver(matrixvec, neigs, std::min(10*neigs+1, NumSamples()));
  eigsolver.init();
  const std::size_t nconvergedEigs = eigsolver.compute(Spectra::SortRule::SmallestMagn, 1000, 1.0e-10);
  std::cout << "Number of converged eigs: " << nconvergedEigs << std::endl;

  std::cout << "eigs: " << eigsolver.eigenvalues().transpose() << std::endl;

  std::cout << std::endl;
  //std::cout << eigsolver.eigenvectors() << std::endl;

  //Eigen::MatrixXd mat(matrix);
  //std::cout << mat.eigenvalues() << std::endl;

  return std::tuple<Eigen::VectorXd, Eigen::VectorXd, Eigen::MatrixXd>(similarity, eigsolver.eigenvalues(), eigsolver.eigenvectors());
}
