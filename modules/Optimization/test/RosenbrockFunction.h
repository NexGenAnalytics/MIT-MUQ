#include "MUQ/Optimization/CostFunction.h"

class RosenbrockFunction : public muq::Optimization::CostFunction {
public:
  RosenbrockFunction() : CostFunction(2) {}

  virtual ~RosenbrockFunction() = default;

  virtual Eigen::MatrixXd Hessian() override {

    Eigen::MatrixXd hess(2,2);
    hess(0,0) = 2.0 - 4.0*a*(x(1)-3.0*x(0)*x(0));
    hess(1,0) = -4.0*a*x(0);
    hess(0,1) = hess(1,0);
    hess(1,1) = 2.0*a;

    return hess;
  }

  virtual double Cost() override {
    return (1.0-x(0))*(1.0-x(0))+a*(x(1)-x(0)*x(0))*(x(1)-x(0)*x(0));
  }

  virtual Eigen::VectorXd const& Gradient() override {

    gradient = Eigen::Vector2d::Constant(2, std::numeric_limits<double>::quiet_NaN());
    gradient(0) = -4.0*a*(x(1)-x(0)*x(0))*x(0)-2.0*(1.0-x(0));
    gradient(1) = 2.0*a*(x(1)-x(0)*x(0));

    return gradient;
  }

 private:

  const double a = 5.0;

};
