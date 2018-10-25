#ifndef STURMSOLVER_H
#define STURMSOLVER_H

#include <Eigen/Core>
#include <vector>

namespace muq {
  namespace Approximation {

    /** Computes the roots of a one dimensional monomial expansion using a Sturm
    sequence to identify intervals containing the roots and then a bisection
    method to find the roots themselves.
    */
    class SturmSolver {

    public:
      SturmSolver();

      /** Construct the solver and call "Compute" */
      SturmSolver(Eigen::VectorXd const& coeffs, double tol=1e-10);

      /**
      Compute the roots (up to the specified tolerance) of a monomial polynomial expansion.
      */
      void Compute(Eigen::VectorXd const& coeffs, double tol=1e-10);

      /** Returns the roots of the polynomial (only valid after Compute has been called)
      */
      Eigen::VectorXd const& Roots() const{return roots;};

      /** Returns the number of computed roots.  Only valid after Compute has been called.
      Will return -1 if the roots have not yet been computed.
      */
      int const& NumRoots() const{return numRealRoots;};

      /** Given a vector of N coefficients, this function evaluates a monomial
          expansion of order N-1.
          @param[in] P The polynomial coefficients (poly order increases with index)
          @param[in] x Value where we want to evaluate the expansion
      */
      static double MonomialEvaluate(Eigen::VectorXd const& P, double x);


    private:

      std::vector<Eigen::VectorXd> BuildSturmSeq(Eigen::VectorXd const& P) const;

      std::pair<int, std::pair<int,int>> ComputeNumRoots(std::vector<Eigen::VectorXd> const& sturmSeq) const;

      Eigen::VectorXd BoundRoots(std::vector<Eigen::VectorXd> const& sturmSeq,
                                 std::pair<int,int>           const& numSignChanges,
                                 double                              tol) const;

      Eigen::VectorXd FindRoots(std::vector<Eigen::VectorXd> const& sturmSeq,
                                Eigen::VectorXd              const& bounds,
                                double                              tol);

      static void MonomialDivision(Eigen::VectorXd const& A,
                                   Eigen::VectorXd const& B,
                                   Eigen::VectorXd      & Q,
                                   Eigen::VectorXd      & R);

      int numRealRoots;

      Eigen::VectorXd roots;

    };
  }
}


#endif // #ifndef STURMSOLVER_H
