#include "MUQ/Approximation/TransportMaps/PolynomialMap.h"

using namespace muq::Approximation;

PolynomialMap::PolynomialMap(std::vector<std::shared_ptr<BasisExpansion> > const& expansions) : ConditionableMap(expansions.size()), expansions(expansions) {
  // make sure the output of each expansion is size 1
  for( unsigned int i=0; i<expansions.size(); ++i ) { assert(expansions[i]->outputSizes(0)==1); }
}

Eigen::VectorXd PolynomialMap::EvaluateForward(Eigen::VectorXd const& x) const {
  // evaluate each expansion
  Eigen::VectorXd result(expansions.size());
  for( unsigned int i=0; i<expansions.size(); ++i ) { result(i) = boost::any_cast<Eigen::VectorXd>(expansions[i]->Evaluate((Eigen::VectorXd)x.head(i+1)) [0]) (0); }

  // return the result
  return result;
}

Eigen::VectorXd PolynomialMap::EvaluateInverse(Eigen::VectorXd const& refPt, Eigen::VectorXd const& tgtPt0) const {
  // evaluate each expansion
  Eigen::VectorXd result = tgtPt0;
  for( unsigned int i=0; i<expansions.size(); ++i ) {
    std::cout << std::endl;
    std::cout << expansions[i]->Jacobian(0, 0, (Eigen::VectorXd)result.head(i+1)) << std::endl;
    std::cout << std::endl;
    //std::cout << expansions[i]->Derivative(i, (Eigen::VectorXd)result.head(i+1)).transpose() << std::endl;
    std::cout << std::endl;
    unsigned int cnt = 0;
    double eval = 1.0;
    while( std::fabs(eval)>tol && cnt++<100 ) {
      eval = boost::any_cast<Eigen::VectorXd>(expansions[i]->Evaluate((Eigen::VectorXd)result.head(i+1)) [0]) (0) - refPt(i);
      const double deriv = expansions[i]->Jacobian(0, 0, (Eigen::VectorXd)result.head(i+1)) (i);
      assert(std::fabs(deriv)>tol);
      result(i) -= eval/deriv;
      std::cout << eval << std::endl;
    }
    std::cout << std::endl << std::endl;
    //result(i) = boost::any_cast<Eigen::VectorXd>(expansions[i]->Evaluate((Eigen::VectorXd)x.head(i+1)) [0]) (0);
  }

  // return the result
  return result;
}
