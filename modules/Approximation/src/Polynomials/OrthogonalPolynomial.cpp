#include "MUQ/Approximation/Polynomials/OrthogonalPolynomial.h"

#include "MUQ/Approximation/Polynomials/SturmSolver.h"

#include "MUQ/Utilities/Exceptions.h"

#include <cstdlib>
#include <typeinfo>
#include <memory>
#include <cxxabi.h>

using namespace muq::Modeling;
using namespace muq::Approximation;


std::shared_ptr<OrthogonalPolynomial> OrthogonalPolynomial::Construct(std::string const& polyName)
{
  return std::dynamic_pointer_cast<OrthogonalPolynomial>(IndexedScalarBasis::Construct(polyName));
}


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
