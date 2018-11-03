#include "MUQ/Optimization/CostFunction.h"

class RosenbrockFunction : public muq::Optimization::CostFunction {
public:
  RosenbrockFunction() : CostFunction(2*Eigen::VectorXi::Ones(1)) {}

  virtual ~RosenbrockFunction() = default;

  virtual Eigen::MatrixXd Hessian(unsigned int const inputDimWrt,
                                  std::vector<Eigen::VectorXd> const& input) override {
    assert(inputDimWrt==0);
    const Eigen::VectorXd& xc = input[0];

    Eigen::MatrixXd hess(2,2);
    hess(0,0) = 2.0 - 4.0*a*(xc(1)-3.0*xc(0)*xc(0));
    hess(1,0) = -4.0*a*xc(0);
    hess(0,1) = hess(1,0);
    hess(1,1) = 2.0*a;

    return hess;
  }

 private:

  const double a = 5.0;

  virtual double CostImpl(muq::Modeling::ref_vector<Eigen::VectorXd> const& input) override {
    const Eigen::VectorXd& xc = input.at(0);

    return (1.0-xc(0))*(1.0-xc(0))+a*(xc(1)-xc(0)*xc(0))*(xc(1)-xc(0)*xc(0));

  }

  virtual void GradientImpl(unsigned int const inputDimWrt,
                    muq::Modeling::ref_vector<Eigen::VectorXd> const& input,
                    Eigen::VectorXd const& sensitivity) override {

    assert(inputDimWrt==0);

    const Eigen::VectorXd& xc = input[0];

    gradient = Eigen::Vector2d::Constant(2, std::numeric_limits<double>::quiet_NaN());
    gradient(0) = -4.0*a*(xc(1)-xc(0)*xc(0))*xc(0)-2.0*(1.0-xc(0));
    gradient(1) = 2.0*a*(xc(1)-xc(0)*xc(0));

    gradient *= sensitivity(0);
  }

};
