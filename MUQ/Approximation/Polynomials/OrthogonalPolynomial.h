#ifndef ORTHOGONALPOLYNOMIAL_H_
#define ORTHOGONALPOLYNOMIAL_H_

#include <functional>
#include <string>

#include "MUQ/Approximation/Polynomials/Polynomial.h"

namespace muq {
  namespace Approximation {

    class GaussQuadrature;

    /**
    @ingroup Polynomials
    @class OrthogonalPolynomial
    @brief A 1D orthogonal polynomial
    @details In general, we use an recursive formula to evaluate a \f$d^{th}\f$ degree
        polynomial using a three term recurrence:
       \f{eqnarray*}{
       p_0(x) &=& \phi_0(x) \\
       p_1(x) &=& \phi_1(x) \\
       p_{k}(x) & = & (a_k x + b_k) p_{k-1}(x) - c_k p_{k-2}(x)
       \f}
       Subclasses specialize for particular polynomials (e.g., Hermite) by
      implementing functions for \f$a_k\f$, \f$b_k\f$, \f$c_k\f$,
      \f$\phi_0(x)\f$, and \f$\phi_1(x)\f$.

       The BasisEvaluate function Uses the Clenshaw algorithm from: http://en.wikipedia.org/wiki/Clenshaw_algorithm.
     */
    class OrthogonalPolynomial : public Polynomial {

      friend class GaussQuadrature;

    public:

      enum RootMethod {
          Sturm,
          Comrade
      };

      /// Create a polynomial
      OrthogonalPolynomial() : Polynomial(){};

      virtual ~OrthogonalPolynomial() = default;

      /**
         @brief Returns the normalization constant for the polynomial of order \f$p\f$.
         @details Many polynomial families are orthogonal with respect to some measure.  This means that for two polynomials \f$P_n(x)\f$ and \f$P_m(x)\f$ in the same family,
\f[
\int P_n(x) P_m(x) w(x) dx = a_n \delta_{nm},
\f]
where \f$\delta_{mn}\f$ is the Kronecker delta function, \f$w(x)\f$ is a weighting function tied to the polynomial family (e.g., the Askey scheme), and \f$a_n\in\mathbb{R}\f$ is a scalar normalization constant.  This function computes the normalization constant \f$a_n\f$.  Polynomial families that are not orthogonal do not have a normalizing constant \f$a_n\f$ and will throw an exception if this function is called.
         @param[in] polyOrder The order of the polynomial
         @return The normalization constant for a particular polynomial family and its weighting function \f$w(x)\f$.
       */
      virtual double Normalization(unsigned int polyOrder) const;

    };
  } // namespace Approximation
} // namespace muq


#endif
