#include "MUQ/Approximation/TransportMaps/PolynomialMap.h"

#include "MUQ/Approximation/Polynomials/Monomial.h"
#include "MUQ/Approximation/Polynomials/OrthogonalPolynomial.h"

using namespace muq::Approximation;

PolynomialMap::PolynomialMap(std::vector<std::shared_ptr<BasisExpansion> > const& expansions, PolynomialMap::InverseMethod const& invMethod) : ConditionableMap(expansions.size()), expansions(expansions), invMethod(invMethod) {
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
    if( invMethod==InverseMethod::Newton ) {
      NewtonsMethod(result, refPt(i), i);
    } else {
      PolynomialRootFinding(result, refPt(i), i, invMethod);
    }
  }

  // return the result
  return result;
}

void PolynomialMap::PolynomialRootFinding(Eigen::VectorXd& result, double const refPt, unsigned int const component, PolynomialMap::InverseMethod const& invMethod) const {
  const Eigen::MatrixXd& coeff = expansions[component]->GetCoeffs();
  const std::vector<std::shared_ptr<IndexedScalarBasis> >& basisComps = expansions[component]->BasisComponents();
  Eigen::VectorXd polyCoeff = Eigen::VectorXd::Zero(expansions[component]->Multis()->GetMaxOrders() (component)+1);
  for( unsigned int term=0; term<expansions[component]->NumTerms(); ++term ) {
    double scale = coeff(term);
    int p = 0;

    // get the coefficients for non costant terms
    for( auto it=expansions[component]->Multis()->at(term)->GetNzBegin(); it!=expansions[component]->Multis()->at(term)->GetNzEnd(); ++it ) {
      if( it->first==component ) { // get the coeffs for this basis
        p = it->second;
      } else {
        assert(it->first<component);
        scale *= basisComps[it->first]->BasisEvaluate(it->second, result(it->first));
      }
    }

    polyCoeff(p) += scale;
  }

  // scale the constant coeff. by the reference point
  polyCoeff(0) -= refPt/basisComps[component]->BasisEvaluate(0, 0.0);

  // choose the closest root to the input point
  auto basis = std::dynamic_pointer_cast<OrthogonalPolynomial>(basisComps[component]);
  assert(basis);
  const Eigen::VectorXd& roots = invMethod==InverseMethod::Comrade? basis->GetRootsComrade(polyCoeff) : basis->GetRootsSturm(polyCoeff, tol);
  assert(roots.size()>0); // make sure that we have at least one root
  const Eigen::VectorXd diff = (roots-Eigen::VectorXd::Constant(roots.size(), result(component))).array().abs();
  int rt;
  diff.minCoeff(&rt);
  result(component) = roots(rt);
}

void PolynomialMap::NewtonsMethod(Eigen::VectorXd& result, double const refPt, unsigned int const component) const {
  // initialize counter and a large error
  unsigned int cnt = 0;
  double eval = 1.0;

  // Newton's method
  while( std::fabs(eval)>tol && cnt++<maxit ) {
    eval = boost::any_cast<Eigen::VectorXd>(expansions[component]->Evaluate((Eigen::VectorXd)result.head(component+1)) [0]) (0) - refPt;
    const double deriv = expansions[component]->Derivative(component, (Eigen::VectorXd)result.head(component+1)) (0);
    assert(std::fabs(deriv)>tol);
    result(component) -= eval/deriv;
  }
}

double PolynomialMap::LogDeterminant(Eigen::VectorXd const& evalPt) const {
  assert(evalPt.size()==expansions.size());

  double logdet = 0.0;
  for( unsigned int i=0; i<expansions.size(); ++i ) {
    logdet += std::log(std::fabs(expansions[i]->Derivative(i, (Eigen::VectorXd)evalPt.head(i+1)) (0)));
  }

  return logdet;
}
