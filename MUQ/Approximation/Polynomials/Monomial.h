#ifndef MONOMIAL_H_
#define MONOMIAL_H_

#include "MUQ/Approximation/Polynomials/Polynomial.h"

namespace muq {
  namespace Approximation {

    /** @ingroup Polynomials
        @class Monomial
        @brief Family of monomial polynomials, i.e. ()\f$1\f$, \f$x\f$, \f$x^2\f$, ect. ...)
        @details This is a simple polynomial basis but could cause conditioning problems in some cases ...
     */
    class Monomial : public Polynomial {
    public:

      Monomial();

      virtual ~Monomial();

      virtual double BasisEvaluate(int const order, double const x) const override;

      virtual double DerivativeEvaluate(int const polyOrder, int const derivOrder, double const x) const override;

      virtual Eigen::VectorXd GetMonomialCoeffs(unsigned int polyOrder) const override;

    protected:
      /// Implement \f$a_k(x)\f$
      /**
	     @param[in] k The order of the polynomial
       */
      virtual double ak(unsigned int k) const override{return 1.0;};

      /// Implement \f$b_k(x)\f$
      /**
	       @param[in] k The order of the polynomial
       */
      virtual double bk(unsigned int k) const override{return 0.0;};

      /// Implement \f$c_k(x)\f$
      /**
	       @param[in] k The order of the polynomial
       */
      virtual double ck(unsigned int k) const override{return 0.0;};

      /// Implement \f$\phi_0(x)\f$
      /**
	       @param[in] x The point where w are evaluating the polynomial
       */
      virtual double phi0(double x) const override{return 1.0;};

      /// Implement \f$\phi_1(x)\f$
      /**
	       @param[in] x The point where w are evaluating the polynomial
       */
      virtual double phi1(double x) const override{return x;};

    };
  } // namespace Approximation
} // namespace muq

#endif
