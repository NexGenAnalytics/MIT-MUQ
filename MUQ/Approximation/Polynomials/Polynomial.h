#ifndef POLYNOMIAL_H_
#define POLYNOMIAL_H_

#include <functional>
#include <string>

#include "MUQ/Approximation/Polynomials/IndexedScalarBasis.h"

namespace muq {
  namespace Approximation {

    class GaussQuadrature;


    class Polynomial : public IndexedScalarBasis {

      friend class GaussQuadrature;

    public:

      enum RootMethod {
          Sturm,
          Comrade
      };

      /// Create a polynomial
      Polynomial() : IndexedScalarBasis(){};

      virtual ~Polynomial() = default;

      /// Evaluate the specific polynomial type (must be implemented by the child)
      /**
	 Inputs:
	 <ol>
	 <li> The order of the polynomial (unsigned int)
	 <li> The point where we are evaluating the polynomial
	 </ol>
	 \return The polynomial value
       */
      virtual double BasisEvaluate(int const order, double const x) const override;

      virtual Eigen::VectorXd EvaluateAllTerms(int    const maxOrder,
                                               double const x) const override;


      /** An orthogonal polynomial of order \f$n\f$ can be represented as
          the sum of one or more monomials.  This function returns the coefficients
          in that monomial expansion.  The expansion takes the form:
          \f$
          \phi_n(x) = \sum_{i=0}^n c_i x^i
          \f$
          @param[in] polyOrder Order of the orthogonal polynomial in question.
          @return Coefficients in the expansion.
      */
      virtual Eigen::VectorXd GetMonomialCoeffs(unsigned int polyOrder) const;

      /** Computes the roots of a 1d polynomial expansion. */
      virtual Eigen::VectorXd GetRoots(Eigen::VectorXd const& coeffs,
                                       RootMethod             method,
                                       double                 tol=1e-10) const;

      /** Computes the roots of a 1d polynomial expansion by forming the Comrade matrix. */
      virtual Eigen::VectorXd GetRootsComrade(Eigen::VectorXd const& coeffs) const;

      /** Computes the roots of a 1d polynomial using a Sturm sequence to identify intervals
          containing the roots and then a bisection solver.
      */
      virtual Eigen::VectorXd GetRootsSturm(Eigen::VectorXd const& coeffs, double tol) const;

    private:

      /// Implement \f$a_k(x)\f$
      /**
	     @param[in] k The order of the polynomial
       */
      virtual double ak(unsigned int k) const = 0;

      /// Implement \f$b_k(x)\f$
      /**
	       @param[in] k The order of the polynomial
       */
      virtual double bk(unsigned int k) const = 0;

      /// Implement \f$c_k(x)\f$
      /**
	       @param[in] k The order of the polynomial
       */
      virtual double ck(unsigned int k) const = 0;

      /// Implement \f$\phi_0(x)\f$
      /**
	       @param[in] x The point where w are evaluating the polynomial
       */
      virtual double phi0(double x) const = 0;

      /// Implement \f$\phi_1(x)\f$
      /**
	       @param[in] x The point where w are evaluating the polynomial
       */
      virtual double phi1(double x) const = 0;
    };
  } // namespace Approximation
} // namespace muq


#endif
