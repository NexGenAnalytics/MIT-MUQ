#include "MUQ/Approximation/Polynomials/Polynomial.h"

#include "MUQ/Utilities/Exceptions.h"

#include <cstdlib>
#include <typeinfo>
#include <memory>
#include <cxxabi.h>

using namespace muq::Modeling;
using namespace muq::Approximation;

Polynomial::Polynomial() :
  WorkPiece(std::vector<std::string>({typeid(unsigned int).name(), typeid(double).name()}), // input types (order, point)
	    std::vector<std::string>({typeid(double).name()})) // output times (polynomial evaluation)
{}

Polynomial::~Polynomial() {}

void Polynomial::EvaluateImpl(ref_vector<boost::any> const& inputs) {
  // extract the inputs
  const unsigned int order = boost::any_cast<unsigned int>(inputs[0]);
  const double x = boost::any_cast<double>(inputs[1]);

  // evaluate the polynomial
  outputs.resize(1);
  outputs[0] = PolynomialEvaluate(order, x);
}


double Polynomial::Normalization(unsigned int polyOrder) const {

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

double Polynomial::PolynomialEvaluate(int const order, double const x) const {

    if(order==0){
        return phi0(x);
    }else if(order==1){
        return phi1(x);
    }else{

        double bkp2 = 0.0;
        double bkp1 = 0.0;
        double bk = 1.0;
        
        for( int k=order-1; k>=0; k-- ) {
            // increment
            bkp2 = bkp1;
            bkp1 = bk;
            
            // compute new bk
            bk = -alpha(k, x)*bkp1 - beta(k+1, x)*bkp2;
        }
        
        return bk*phi0(x) + bkp1*(phi1(x)+alpha(0, x)*phi0(x));
    }
}

std::shared_ptr<Polynomial> Polynomial::Construct(std::string const& polyName){

  auto map = GetPolynomialMap();
  auto it = map->find(polyName);
  if(it == map->end()){
      std::string message = "The polynomial family, \"" + polyName + "\" has not been registered with the polynomial factory.  Does the class exist?\n  Registered families include:\n";
      for(auto iter = map->begin(); iter!=map->end(); ++iter)
          message += "    " + iter->first + "\n";
      message += "\n";
      
      throw muq::NotRegisteredError(message);
      return nullptr;
  }else{    
      return it->second();
  }
}

std::shared_ptr<Polynomial::PolynomialMapType> Polynomial::GetPolynomialMap() {
  // define a static map from type to constructor
  static std::shared_ptr<PolynomialMapType> map;

  if( !map ) { // if the map has not yet been created ...
    // ... create the map
    map = std::make_shared<PolynomialMapType>();
  }

  return map;
}

