#include "MUQ/SamplingAlgorithms/SampleGraphs/KolmogorovOperator.h"

#include <Spectra/SymEigsSolver.h>
#include <Spectra/MatOp/SparseGenMatProd.h>
//#include <Spectra/GenEigsSolver.h>
//#include <Spectra/Util/SelectionRule.h>

namespace pt = boost::property_tree;
using namespace Spectra;
using namespace muq::Modeling;
using namespace muq::SamplingAlgorithms;

KolmogorovOperator::KolmogorovOperator(std::shared_ptr<RandomVariable> const& rv, pt::ptree const& options) :
DensityEstimation(rv, options),
operatorParameter(options.get<double>("OperatorParameter", 1.0)),
beta(options.get<double>("VariableBandwidth", (operatorParameter-2.0)/(manifoldDim+2.0))),
neigs(options.get<std::size_t>("NumEigenpairs", 5*std::log((double)NumSamples()))),
powminOperator(options.get<double>("OperatorMinLog2Bandwidth", -4)),
powmaxOperator(options.get<double>("OperatorMaxLog2Bandwidth", 0))
{
  // the normalization parameter
  alpha = 1.0 + beta + (manifoldDim*beta - operatorParameter)/2.0;
}

KolmogorovOperator::KolmogorovOperator(std::shared_ptr<SampleCollection> const& samples, pt::ptree const& options) :
DensityEstimation(samples, options),
operatorParameter(options.get<double>("OperatorParameter", 1.0)),
beta(options.get<double>("VariableBandwidth", (operatorParameter-2.0)/(manifoldDim+2.0))),
neigs(options.get<std::size_t>("NumEigenpairs", 5*std::log((double)NumSamples()))),
powminOperator(options.get<double>("OperatorMinLog2Bandwidth", -4)),
powmaxOperator(options.get<double>("OperatorMaxLog2Bandwidth", 0))
{
  // the normalization parameter
  alpha = 1.0 + beta + (manifoldDim*beta - operatorParameter)/2.0;
}

KolmogorovOperator::KolmogorovOperator(Eigen::MatrixXd const& mat, pt::ptree const& options) :
DensityEstimation(mat, options),
operatorParameter(options.get<double>("OperatorParameter", 1.0)),
beta(options.get<double>("VariableBandwidth", (operatorParameter-2.0)/(manifoldDim+2.0))),
neigs(options.get<std::size_t>("NumEigenpairs", 5*std::log((double)NumSamples()))),
powminOperator(options.get<double>("OperatorMinLog2Bandwidth", -4)),
powmaxOperator(options.get<double>("OperatorMaxLog2Bandwidth", 0))
{
  // the normalization parameter
  alpha = 1.0 + beta + (manifoldDim*beta - operatorParameter)/2.0;
}

double KolmogorovOperator::OperatorBandwidthParameter() const { return operatorBandwidthParameter; }

void KolmogorovOperator::TuneOperatorBandwidth(Eigen::VectorXd const& density) const {
  // compute the bandwidth function
  const Eigen::VectorXd scaledDens = density.array().pow(beta);

  double dimensionEstimate;
  std::tie(operatorBandwidthParameter, dimensionEstimate) = TuneKernelBandwidth(scaledDens, powminOperator, powmaxOperator, operatorBandwidthParameter);
  if( tuneDimension ) { alpha = 1.0 + beta + (manifoldDim*beta - operatorParameter)/2.0; }
}

Eigen::VectorXd KolmogorovOperator::DiscreteOperator(Eigen::VectorXd const& density, Eigen::SparseMatrix<double>& matrix, double epsilon, bool const tune) const {
  const std::size_t n = NumSamples();

  // compute the bandwidth function
  const Eigen::VectorXd scaledDens = density.array().pow(beta);

  const double scalem = 4.0;

  // if not supplied, use the stored value---the scaling means that the parameter the kernel requries is actually 2*eps
  epsilon = (std::isnan(epsilon)? scalem*operatorBandwidthParameter : scalem*epsilon);

  // get the optimal kolmogorov bandwidth parameter
  if( tune ) {
    double dimensionEstimate;
    std::tie(epsilon, dimensionEstimate) = TuneKernelBandwidth(scaledDens, powminOperator, powmaxOperator, epsilon);
    //epsilon = 0.8*epsilon;

    // update the stored values
    operatorBandwidthParameter = epsilon/scalem;
    if( tuneDimension ) { alpha = 1.0 + beta + (manifoldDim*beta - operatorParameter)/2.0; }
  }

  // compute the unnormalized kernel and store the normalization in the bandwidth vector
  std::vector<Eigen::Triplet<double> > entries;
  assert(epsilon>1.0e-14);
  Eigen::VectorXd rowsum = KernelMatrix(sparsityTol, epsilon, scaledDens, entries, !tune);
  Eigen::VectorXd similarity = (rowsum.array()*scaledDens.array().pow(-manifoldDim)).pow(-alpha);

  // loop through the non-zeros
  for( auto& entry : entries ) {
    entry = Eigen::Triplet<double>(entry.row(), entry.col(), entry.value()*similarity(entry.row())*similarity(entry.col()));
  }

  // recompute the normalization (in the bandwidth vector)
  rowsum = Eigen::VectorXd::Zero(NumSamples());
  for( const auto& entry : entries ) { rowsum(entry.row()) += entry.value(); }

  similarity = scaledDens.array()*rowsum.array().sqrt();

  epsilon = epsilon/scalem;

  for( auto& entry : entries ) { entry = Eigen::Triplet<double>(entry.row(), entry.col(), entry.value()/(epsilon*similarity(entry.row())*similarity(entry.col())) ); }

  for( std::size_t i=0; i<n; ++i ) { entries.emplace_back(i, i, -1.0/(epsilon*scaledDens(i)*scaledDens(i))); }

  matrix.resize(n, n);
  matrix.setFromTriplets(entries.begin(), entries.end());

  return similarity;
}

std::tuple<Eigen::VectorXd, Eigen::VectorXd, Eigen::MatrixXd> KolmogorovOperator::Eigendecomposition(Eigen::VectorXd const& density, double epsilon, bool const tune) {
  // compute and store the discrete Kolmogorov operator
  Eigen::SparseMatrix<double> matrix;
  const Eigen::VectorXd similarity = DiscreteOperator(density, matrix, epsilon, tune);

  // compute the eigen decomposition
  Spectra::SparseGenMatProd<double> matrixvec(matrix);
  Spectra::SymEigsSolver<Spectra::SparseGenMatProd<double> > eigsolver(matrixvec, neigs, std::min(10*neigs+1, NumSamples()));
  eigsolver.init();
  const std::size_t nconvergedEigs = eigsolver.compute(Spectra::SortRule::SmallestMagn, 1000, 1.0e-10);

  assert(similarity.size()==eigsolver.eigenvectors().rows());

  return std::tuple<Eigen::VectorXd, Eigen::VectorXd, Eigen::MatrixXd>(similarity, eigsolver.eigenvalues(), eigsolver.eigenvectors());
}

Eigen::MatrixXd KolmogorovOperator::GradientVectorField(Eigen::VectorXd const& similarity, Eigen::VectorXd const& eigvals, Eigen::MatrixXd const& eigvecs, Eigen::VectorXd const& func1) const {
  const std::size_t n = NumSamples();
  assert(n>0);
  const std::size_t dim = Point(0).size();

  Eigen::MatrixXd v(n, dim);
  for( std::size_t i=0; i<n; ++i ) { v.row(i) = Point(i); }

  return GradientVectorInnerProduct(similarity, eigvals, eigvecs, func1, v);
}

Eigen::MatrixXd KolmogorovOperator::GradientVectorInnerProduct(Eigen::VectorXd const& similarity, Eigen::VectorXd const& eigvals, Eigen::MatrixXd const& eigvecs, Eigen::VectorXd const& func1, Eigen::MatrixXd const& func2) const {
  assert(eigvecs.rows()==func1.size());
  assert(eigvecs.rows()==func2.rows());
  assert(eigvecs.cols()==eigvals.size());
  assert(eigvals.size()==neigs);

  const Eigen::MatrixXd eigvecsRight = similarity.array().inverse().matrix().asDiagonal()*eigvecs;
  const Eigen::MatrixXd eigvecsLeft = similarity.asDiagonal()*eigvecs;

  // compute the coefficients for the functions
  const Eigen::VectorXd ucoeff = eigvecsLeft.transpose()*func1;
  const Eigen::MatrixXd vcoeff = eigvecsLeft.transpose()*func2;

  // compute the coordinate coefficients
  const std::size_t n = NumSamples();
  Eigen::MatrixXd gradient = Eigen::MatrixXd::Zero(n, func2.cols());
  for( std::size_t j=0; j<neigs; ++j ) {
    for( std::size_t k=j; k<neigs; ++k ) {
      const Eigen::VectorXd phijk = eigvecsRight.col(j).array()*eigvecsRight.col(k).array();
      assert(phijk.size()==n);
      for( std::size_t l=0; l<neigs; ++l ) {
        const double Cjkl = phijk.dot(eigvecsLeft.col(l))/2.0;
        const double scale = Cjkl*(eigvals(l)-eigvals(k)-eigvals(j));
        const double scalej = ucoeff(j)*scale;

        for( std::size_t col=0; col<gradient.cols(); ++col ) { gradient.col(col) += (scalej*vcoeff(k, col))*eigvecsRight.col(l); }

        if( j!=k ) {
          const double scalek = ucoeff(k)*scale;
          for( std::size_t col=0; col<gradient.cols(); ++col ) { gradient.col(col) += (scalek*vcoeff(j, col))*eigvecsRight.col(l); }
        }
      }
    }
  }

  return gradient;
}

Eigen::VectorXd KolmogorovOperator::KolmogorovProblemSolution(Eigen::VectorXd const& similarity, Eigen::VectorXd const& eigvals, Eigen::MatrixXd const& eigvecs, Eigen::VectorXd const& func) const {
  assert(eigvecs.rows()==func.size());
  assert(eigvecs.cols()==eigvals.size());
  assert(eigvals.size()==neigs);

  // compute the coefficients for the functions
  const Eigen::VectorXd coeff = eigvecs.transpose()*similarity.asDiagonal()*func;
  assert(coeff.size()==neigs);

  // the coefficients of the solution
  Eigen::VectorXd solnCoeff = Eigen::VectorXd::Zero(neigs);
  for( std::size_t i=0; i<neigs; ++i ) {
    if( std::abs(eigvals(i))>1.0e-10 ) { solnCoeff(i) = coeff(i)/eigvals(i); }
  }

  return similarity.array().inverse().matrix().asDiagonal()*eigvecs*solnCoeff;
}
