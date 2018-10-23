#include "MUQ/Approximation/TransportMaps/PolynomialMap.h"

using namespace muq::Approximation;

PolynomialMap::PolynomialMap(std::vector<std::shared_ptr<BasisExpansion> > const& expansions) : ConditionableMap(expansions.size()), expansions(expansions) {
  // make sure the output of each expansion is size 1
  for( unsigned int i=0; i<expansions.size(); ++i ) { assert(expansions[i]->outputSizes(0)==1); }
}

Eigen::VectorXd PolynomialMap::EvaluateForward(Eigen::VectorXd const& x) const {
  assert(x.size()==expansions.size());

  // evaluate each expansion
  Eigen::VectorXd result(expansions.size());
  for( unsigned int i=0; i<expansions.size(); ++i ) { result(i) = boost::any_cast<Eigen::VectorXd>(expansions[i]->Evaluate((Eigen::VectorXd)x.head(i+1)) [0]) (0); }

  // return the result
  return result;
}

Eigen::VectorXd PolynomialMap::EvaluateInverse(Eigen::VectorXd const& refPt, Eigen::VectorXd const& tgtPt0) const {
  assert(refPt.size()==expansions.size());
  assert(tgtPt0.size()==expansions.size());

  // set initial guess
  Eigen::VectorXd result = tgtPt0;

  for( unsigned int i=0; i<expansions.size(); ++i ) {
    // initialize counter and a large error
    unsigned int cnt = 0;
    double eval = 1.0;

    // Newton's method
    while( std::fabs(eval)>tol && cnt++<maxit ) {
      eval = boost::any_cast<Eigen::VectorXd>(expansions[i]->Evaluate((Eigen::VectorXd)result.head(i+1)) [0]) (0) - refPt(i);
      const double deriv = expansions[i]->Derivative(i, (Eigen::VectorXd)result.head(i+1)) (0);
      assert(std::fabs(deriv)>tol);
      result(i) -= eval/deriv;
    }
  }

  // return the result
  return result;
}

double PolynomialMap::LogDeterminant(Eigen::VectorXd const& evalPt) const {
  assert(evalPt.size()==expansions.size());

  double logdet = 1.0;
  for( unsigned int i=0; i<expansions.size(); ++i ) {
    logdet += std::log(std::fabs(expansions[i]->Derivative(i, (Eigen::VectorXd)evalPt.head(i+1)) (0)));
  }

  return logdet;
}
