#include "MUQ/Approximation/Polynomials/Polynomial.h"

#include "MUQ/Approximation/Polynomials/SturmSolver.h"

#include <cstdlib>
#include <typeinfo>
#include <memory>
#include <cxxabi.h>

using namespace muq::Modeling;
using namespace muq::Approximation;

double Polynomial::BasisEvaluate(int const order, double const x) const {

    if(order==0){
        return phi0(x);
    }else if(order==1){
        return phi1(x);
    }else{

      // "Downward" Clenshaw algorithm  http://mathworld.wolfram.com/ClenshawRecurrenceFormula.html
      double yk2 = 0.0;
      double yk1 = 0.0;
      double yk = 1.0;
      double alpha, beta;

      for( int k=order-1; k>=0; k-- ) {
        yk2 = yk1;
        yk1 = yk;

        alpha = ak(k+1)*x + bk(k+1);
        beta = -ck(k+2);
        yk = alpha*yk1 + beta*yk2;
      }
      beta = -ck(2);
      return yk1*phi1(x) + beta * phi0(x)*yk2;
    }
}

Eigen::VectorXd Polynomial::EvaluateAllTerms(int    const maxOrder,
                                                       double const x) const {

  Eigen::VectorXd output(maxOrder+1);
  output(0) = phi0(x);
  if(maxOrder>0)
    output(1) = phi1(x);

  for(int i=2; i<=maxOrder; ++i)
    output(i) = (ak(i)*x + bk(i))*output(i-1) - ck(i)*output(i-2);

  return output;
}

Eigen::VectorXd Polynomial::GetMonomialCoeffs(unsigned int polyOrder) const {

  Eigen::VectorXd oldCoeffs, oldOldCoeffs;

  oldOldCoeffs = Eigen::VectorXd::Zero(polyOrder+1);
  oldOldCoeffs(0) = phi0(0);

  if(polyOrder>=1){
    oldCoeffs = Eigen::VectorXd::Zero(polyOrder+1);
    double slope = (phi1(0.9)-phi1(0.1))/0.8;
    double intercept = phi1(0.5) - slope*0.5;
    oldCoeffs(0) = intercept;
    oldCoeffs(1) = slope;
  }

  if(polyOrder==0){
    return oldOldCoeffs;
  }else if(polyOrder==1){
    return oldCoeffs;
  }else{
    Eigen::VectorXd monoCoeffs = Eigen::VectorXd::Zero(polyOrder+1);

    for(int p=2; p<=polyOrder; ++p){
      monoCoeffs = Eigen::VectorXd::Zero(polyOrder+1);
      monoCoeffs.tail(polyOrder) = ak(p)*oldCoeffs.head(polyOrder);
      monoCoeffs += bk(p)*oldCoeffs;
      monoCoeffs -= ck(p)*oldOldCoeffs;

      oldOldCoeffs = oldCoeffs;
      oldCoeffs = monoCoeffs;
    }

    return monoCoeffs;
  }

}


Eigen::VectorXd Polynomial::GetRoots(Eigen::VectorXd const& coeffs,
                                     RootMethod             method,
                                     double                 tol) const {

  if(method == RootMethod::Sturm){
    return GetRootsSturm(coeffs, tol);
  }else if(method == RootMethod::Comrade){
    return GetRootsComrade(coeffs);
  }else{
    assert(false);
    return Eigen::VectorXd();
  }
}

Eigen::VectorXd Polynomial::GetRootsSturm(Eigen::VectorXd const& coeffs, double tol) const
{
  Eigen::VectorXd monoCoeffs = Eigen::VectorXd::Zero(coeffs.size());
  for(int order=0; order<coeffs.size(); ++order)
    monoCoeffs.head(order+1) += coeffs(order)*GetMonomialCoeffs(order);

  return SturmSolver(monoCoeffs,tol).Roots();
}


Eigen::VectorXd Polynomial::GetRootsComrade(Eigen::VectorXd const& coeffs) const
{
  const int N = coeffs.size()-1;

  // If we're constant, we assume there are no roots (could be an infinite number though)
  if(N==0){
    return Eigen::VectorXd();

  // If linear, compute the roots
  }else if(N==1){
    if(coeffs(1)==0){
      return Eigen::VectorXd();
    }else{
      return -coeffs(0)/(coeffs(1)*BasisEvaluate(1,1.0)) * Eigen::VectorXd::Ones(1);
    }
  }

  // Normalize the coefficients
  Eigen::VectorXd normCoeffs = coeffs/coeffs(N);

  // Use the three term recurrence relationship to build the Comrade matrix
  Eigen::MatrixXd C = Eigen::MatrixXd::Zero(N,N);

  C(0,0) = -bk(1)/ak(1);
  C(0,1) = 1.0/ak(1);
  for(int i=1; i<N-1; ++i){
    C(i,i-1) = ck(i+1)/ak(i+1);
    C(i,i) = -bk(i+1)/ak(i+1);
    C(i,i+1) = 1.0/ak(i+1);
  }

  for(int i=0; i<N; ++i)
    C(N-1,i) = -normCoeffs(i)/ak(N);

  C(N-1,N-2) += ck(N)/ak(N);
  C(N-1,N-1) -= bk(N)/ak(N);

  // get the real eigenvalues
  const auto& eigens = C.eigenvalues();
  const auto isReal = eigens.imag().array().abs() < 1.0e-14;
  Eigen::VectorXd eigs(isReal.count());

  int cnt = 0;
  for( unsigned int i=0; i<eigens.size(); ++i ) {
    if( isReal(i) )
      eigs(cnt++) = eigens(i).real();
  }

  if(eigs.size()>1)
    std::sort(&eigs[0], &eigs[eigs.size()-1]);

  return eigs;
}
