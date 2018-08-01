#include <gtest/gtest.h>

#include "RosenbrockFunction.h"

using namespace muq::Modeling;
using namespace muq::Optimization;

TEST(CostFunctionTests, RosenbrockCost) {
  // the Rosenbrock cost function
  auto rosen = std::make_shared<RosenbrockFunction>();

  // choose a random point
  const Eigen::VectorXd x = Eigen::Vector2d::Random();

  // the true value
  const double cst = (1.0-x(0))*(1.0-x(0))+100.0*(x(1)-x(0)*x(0))*(x(1)-x(0)*x(0));

  // check the cost evaluations
  EXPECT_DOUBLE_EQ(cst, rosen->Evaluate(x).at(0) (0));
  EXPECT_DOUBLE_EQ(cst, rosen->Cost(ref_vector<Eigen::VectorXd>(1, x)));
  EXPECT_DOUBLE_EQ(cst, rosen->Cost(x));

  // the true gradient
  const Eigen::Vector2d grad_true(-400.0*(x(1)-x(0)*x(0))*x(0)-2.0*(1.0-x(0)), 200.0*(x(1)-x(0)*x(0)));

  // compute the gradient
  const Eigen::VectorXd& grad_test0 = rosen->Gradient(0, x, (Eigen::VectorXd)Eigen::VectorXd::Ones(2));
  const Eigen::VectorXd& grad_test1 = rosen->Gradient(0, ref_vector<Eigen::VectorXd>(1, x), (Eigen::VectorXd)Eigen::VectorXd::Ones(2));

  EXPECT_DOUBLE_EQ((grad_true-grad_test0).norm(), 0.0);
  EXPECT_DOUBLE_EQ((grad_true-grad_test1).norm(), 0.0);
}
