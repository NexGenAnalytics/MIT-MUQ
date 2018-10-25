#include "MUQ/Approximation/Polynomials/OrthogonalPolynomial.h"

#include "MUQ/Approximation/Polynomials/Monomial.h"

#include "MUQ/Utilities/Exceptions.h"

#include <cstdlib>
#include <typeinfo>
#include <memory>
#include <cxxabi.h>

using namespace muq::Modeling;
using namespace muq::Approximation;


double OrthogonalPolynomial::Normalization(unsigned int polyOrder) const {

    std::string rawName = typeid(*this).name();

    int status = -4; // some arbitrary value to eliminate the compiler warning

    // enable c++11 by passing the flag -std=c++11 to g++
    std::unique_ptr<char, void(*)(void*)> res {
        abi::__cxa_demangle(rawName.c_str(), NULL, NULL, &status),
        std::free
    };

    std::string className = (status==0) ? res.get() : rawName;

    throw muq::NotImplementedError("The Normalization function has not been implemented for the class \"" + className + "\".  Is this polynomial family orthogonal?");

    return std::numeric_limits<double>::quiet_NaN();
};

double OrthogonalPolynomial::BasisEvaluate(int const order, double const x) const {

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

Eigen::VectorXd OrthogonalPolynomial::EvaluateAllTerms(int    const maxOrder,
                                                       double const x) const {

  Eigen::VectorXd output(maxOrder+1);
  output(0) = phi0(x);
  if(maxOrder>0)
    output(1) = phi1(x);

  for(int i=2; i<=maxOrder; ++i)
    output(i) = (ak(i)*x + bk(i))*output(i-1) - ck(i)*output(i-2);

  return output;
}

Eigen::VectorXd OrthogonalPolynomial::GetMonomialCoeffs(unsigned int polyOrder) const {

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


Eigen::VectorXd OrthogonalPolynomial::GetRoots(Eigen::VectorXd const& coeffs, std::string const& method) const {

  if(method == "Sturm"){
    return GetRootsSturm(coeffs, 1e-10);
  }else if(method=="Comrade"){
    return GetRootsComrade(coeffs);
  }else{
    throw std::invalid_argument("Invalid option to OrthogonalPolynomial::GetRoots.  Valid methods are \"Sturm\" and \"Comrade\".");
  }
}

Eigen::VectorXd OrthogonalPolynomial::GetRootsSturm(Eigen::VectorXd const& coeffs, double tol) const
{
  Eigen::VectorXd monoCoeffs = Eigen::VectorXd::Zero(coeffs.size());
  for(int order=0; order<coeffs.size(); ++order)
    monoCoeffs.head(order+1) += coeffs(order)*GetMonomialCoeffs(order);

  return Monomial::MonomialRoots(monoCoeffs,tol);
}


Eigen::VectorXd OrthogonalPolynomial::GetRootsComrade(Eigen::VectorXd const& coeffs) const
{
  const int N = coeffs.size()-1;

  // Normalize the coefficients
  Eigen::VectorXd normCoeffs = coeffs/coeffs(N);

  // Use the three term recurrence relationship to build the Comrade matrix
  Eigen::MatrixXd C = Eigen::MatrixXd::Zero(N,N);

  C(0,0) = -bk(1)/ak(1);
  if( N>1 ) { C(0,1) = 1.0/ak(1); }
  for(int i=1; i<N-1; ++i){
    C(i,i-1) = ck(i+1)/ak(i+1);
    C(i,i) = -bk(i+1)/ak(i+1);
    C(i,i+1) = 1.0/ak(i+1);
  }

  for(int i=0; i<N; ++i)
    C(N-1,i) = -normCoeffs(i)/ak(N);

  if( N>1 ) { C(N-1,N-2) += ck(N)/ak(N); }
  C(N-1,N-1) -= bk(N)/ak(N);

  // get the real eigenvalues
  const auto& eigens = C.eigenvalues();
  auto ims = eigens.imag().array().abs()<1.0e-14;
  Eigen::VectorXd eigs(ims.count());
  int cnt = 0;
  for( unsigned int i=0; i<eigens.size(); ++i ) {
    if( ims(i) ) { eigs(cnt++) = eigens(i).real(); }
  }

  return eigs;

}
