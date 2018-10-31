#include "MUQ/Approximation/TransportMaps/KLSamplesCost.h"

using namespace muq::Modeling;
using namespace muq::Optimization;
using namespace muq::Approximation;

KLSamplesCost::KLSamplesCost(Eigen::MatrixXd const& vand, Eigen::MatrixXd const& deriv) : CostFunction(Eigen::VectorXi::Constant(1, vand.cols())), vand(vand), deriv(deriv) {}

double KLSamplesCost::CostImpl(ref_vector<Eigen::VectorXd> const& input) {
  const Eigen::VectorXd& c = input[0];

  // reference variables
  const Eigen::VectorXd r = vand*c;

  /*std::cout << vand << std::endl;
  std::cout << std::endl << std::endl;
  std::cout << deriv << std::endl;*/

  // the derivative
  /*Eigen::VectorXd dSdz = (deriv*c).array().log();
  for( unsigned int i=0; i<dSdz.size(); ++i ) {
    if( std::isnan(dSdz(i)) ) { return RAND_MAX; }
  }
  return (0.5*r.dot(r)-dSdz.sum())/(double)vand.rows();*/

  // the cost function
  return (0.5*r.dot(r)-(deriv*c).array().abs().log().sum())/(double)vand.rows();
}

void KLSamplesCost::GradientImpl(unsigned int const inputDimWrt, ref_vector<Eigen::VectorXd> const& input, Eigen::VectorXd const& sensitivity) {
  const Eigen::VectorXd& c = input[0];

  // reference variables
  gradient = vand.transpose()*vand*c;

  const Eigen::VectorXd dSdz = deriv*c;

  for( unsigned int i=0; i<c.size(); ++i ) {
    for( unsigned int j=0; j<deriv.rows(); ++j ) {
      gradient(i) -= std::fabs(deriv(j,i))/dSdz(i);
    }
  }

  gradient /= (double)vand.rows();
}

KLSamplesConstraint::KLSamplesConstraint(Eigen::MatrixXd const& deriv) : ModPiece(Eigen::VectorXi::Constant(1, deriv.cols()), Eigen::VectorXi::Constant(1, deriv.rows())), deriv(deriv) {}

void KLSamplesConstraint::EvaluateImpl(ref_vector<Eigen::VectorXd> const& input) {
  const Eigen::VectorXd& c = input[0];

  std::cout << "EVALUATING THE CONSTRAINT" << std::endl;

  outputs.resize(1);
  outputs[0] = deriv*c;;
}

void KLSamplesConstraint::JacobianImpl(unsigned int const outputDimWrt, unsigned int const inputDimWrt, muq::Modeling::ref_vector<Eigen::VectorXd> const& input) {
  std::cout << "CHECK DIMENSIONS" << std::endl;
  jacobian = -deriv;
}
