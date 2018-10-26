#include "MUQ/Approximation/Polynomials/Monomial.h"

using namespace muq::Approximation;

Monomial::Monomial() : Polynomial() {}

Monomial::~Monomial() {}

double Monomial::DerivativeEvaluate(int const polyOrder, int const derivOrder, double const x) const {

    if((derivOrder > polyOrder) || (polyOrder==0))
        return 0.0;

    double c = 1.0;
    for(int k=polyOrder; k>polyOrder-derivOrder; --k)
        c *= k;

    return c*std::pow(x, polyOrder-derivOrder);

}

double Monomial::BasisEvaluate(int const order, double const x) const {
    return std::pow(x, order);
}

Eigen::VectorXd Monomial::GetMonomialCoeffs(unsigned int polyOrder) const{
  Eigen::VectorXd coeffs = Eigen::VectorXd::Zero(polyOrder+1);
  coeffs(polyOrder) = 1.0;
  return coeffs;
}

REGISTER_SCALARBASIS_FAMILY(Monomial)
