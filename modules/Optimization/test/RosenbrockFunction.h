#include "MUQ/Optimization/CostFunction.h"

class RosenbrockFunction : public muq::Optimization::CostFunction {
public: 
  inline RosenbrockFunction() : CostFunction(Eigen::Vector2i(2, 1)) {}

  virtual inline ~RosenbrockFunction() {}

 private:

  inline virtual double CostImpl(muq::Modeling::ref_vector<Eigen::VectorXd> const& input) override {
    const Eigen::VectorXd& xc = input[0];
    const Eigen::VectorXd& a = input[1];

    return (1.0-xc(0))*(1.0-xc(0))+a(0)*(xc(1)-xc(0)*xc(0))*(xc(1)-xc(0)*xc(0));
  }

  virtual void GradientImpl(unsigned int const inputDimWrt, muq::Modeling::ref_vector<Eigen::VectorXd> const& input, Eigen::VectorXd const& sensitivity) override {
    assert(inputDimWrt==0);
    
    const Eigen::VectorXd& xc = input[0];
    const Eigen::VectorXd& a = input[1];

    gradient = Eigen::Vector2d::Constant(2, std::numeric_limits<double>::quiet_NaN());
    gradient(0) = -4.0*a(0)*(xc(1)-xc(0)*xc(0))*xc(0)-2.0*(1.0-xc(0));
    gradient(1) = 2.0*a(0)*(xc(1)-xc(0)*xc(0));

    gradient.array()*sensitivity.array();
  }
};
