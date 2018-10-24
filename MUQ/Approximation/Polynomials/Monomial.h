#ifndef MONOMIAL_H_
#define MONOMIAL_H_

#include "MUQ/Approximation/Polynomials/IndexedScalarBasis.h"

namespace muq {
  namespace Approximation {

    /** @ingroup Polynomials
        @class Monomial
        @brief Family of monomial polynomials, i.e. ()\f$1\f$, \f$x\f$, \f$x^2\f$, ect. ...)
        @details This is a simple polynomial basis but could cause conditioning problems in some cases ...
     */
    class Monomial : public IndexedScalarBasis {
    public:

      Monomial();

      virtual ~Monomial();

      virtual double BasisEvaluate(int const order, double const x) const override;

      virtual double DerivativeEvaluate(int const polyOrder, int const derivOrder, double const x) const override;

      /** Given a vector of N coefficients, this function evaluates a monomial
          expansion of order N-1.
          @param[in] P The polynomial coefficients (poly order increases with index)
          @param[in] x Value where we want to evaluate the expansion
      */
      static double MonomialEvaluate(Eigen::VectorXd const& P, double x);

      /** Use a Sturm sequence and bisection method to compute the roots of a monomial expansion.
          @param[in] P The polynomial coefficients.
          @param[in] tol The desired accuracy of root estimates.
      */
      static Eigen::VectorXd MonomialRoots(Eigen::VectorXd const& Pin, double tol);


    protected:
      static void MonomialDivision(Eigen::VectorXd const& A,
                                   Eigen::VectorXd const& B,
                                   Eigen::VectorXd      & Q,
                                   Eigen::VectorXd      & R);

    };
  } // namespace Approximation
} // namespace muq

#endif
