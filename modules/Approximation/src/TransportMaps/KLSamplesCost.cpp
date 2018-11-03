#include "MUQ/Approximation/TransportMaps/KLSamplesCost.h"

using namespace muq::Modeling;
using namespace muq::Optimization;
using namespace muq::Approximation;

KLSamplesCost::KLSamplesCost(Eigen::MatrixXd const& vand,
                             Eigen::MatrixXd const& deriv) :
                               CostFunction(vand.cols()),
                               vand(vand),
                               deriv(deriv),
                               numSamps(vand.rows()) {

 assert(vand.cols()==deriv.cols());
 assert(vand.rows()==deriv.rows());
}

void KLSamplesCost::SetPoint(Eigen::VectorXd const& evalPt) {

  x = evalPt;
  vandApp = vand*x;
  derivApp = deriv*x;
}

double KLSamplesCost::Cost() {
  return (0.5*vandApp.squaredNorm()-derivApp.array().abs().log().sum()) / numSamps;
}


Eigen::VectorXd const& KLSamplesCost::Gradient() {
  gradient = (vandApp.transpose()*vand - derivApp.array().inverse().matrix().transpose() * deriv).transpose();
  gradient /= numSamps;
  return gradient;
}

Eigen::MatrixXd KLSamplesCost::Hessian() {
  return (vand.transpose() * vand + deriv.transpose() * (derivApp.array().square().inverse().matrix().asDiagonal() * deriv))/numSamps;
};

Eigen::VectorXd KLSamplesCost::ApplyHessian(Eigen::VectorXd const& vec) {
   return (vand.transpose() * vand * vec + deriv.transpose()*derivApp.array().square().inverse().matrix().asDiagonal() * deriv*vec)/numSamps;
};

KLSamplesConstraint::KLSamplesConstraint(Eigen::MatrixXd const& deriv) :
                      ModPiece(Eigen::VectorXi::Constant(1, deriv.cols()),
                               Eigen::VectorXi::Constant(1, deriv.rows())),
                      deriv(deriv) {}

void KLSamplesConstraint::EvaluateImpl(ref_vector<Eigen::VectorXd> const& input) {
  const Eigen::VectorXd& c = input[0];

  outputs.resize(1);
  outputs[0] = deriv*c;
}

void KLSamplesConstraint::JacobianImpl(unsigned int const outputDimWrt, unsigned int const inputDimWrt, muq::Modeling::ref_vector<Eigen::VectorXd> const& input) {
  jacobian = -deriv;
}
