#include "MUQ/Approximation/TransportMaps/KLSamplesCost.h"

using namespace muq::Modeling;
using namespace muq::Optimization;
using namespace muq::Approximation;

KLSamplesCost::KLSamplesCost(Eigen::MatrixXd const& vand, Eigen::MatrixXd const& deriv) : CostFunction(Eigen::VectorXi::Constant(1, vand.cols())), vand(vand), deriv(deriv) {}

double KLSamplesCost::CostImpl(ref_vector<Eigen::VectorXd> const& input) {
  const Eigen::VectorXd& c = input[0];

  // reference variables
  const Eigen::VectorXd r = vand*c;

  // the cost function
  return (0.5*r.dot(r)-(deriv*c).array().abs().log().sum())/(double)vand.rows();
}

void KLSamplesCost::GradientImpl(unsigned int const inputDimWrt, ref_vector<Eigen::VectorXd> const& input, Eigen::VectorXd const& sensitivity) {
  const Eigen::VectorXd& c = input[0];

  Eigen::VectorXd fgamma = vand*c;
  Eigen::VectorXd ggamma = deriv*c;

  gradient = (fgamma.transpose()*vand - ggamma.array().inverse().matrix().transpose() * deriv).transpose();
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
