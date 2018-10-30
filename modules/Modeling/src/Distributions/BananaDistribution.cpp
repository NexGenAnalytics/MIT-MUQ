#include "MUQ/Modeling/Distributions/BananaDistribution.h"

using namespace muq::Modeling;


BananaDistribution::BananaDistribution(double aIn, double bIn) : Distribution(2),
                                                                 a(aIn),
                                                                 b(bIn),
                                                                 stdNorm(std::make_shared<Gaussian>(2)){}

double BananaDistribution::LogDensityImpl(ref_vector<Eigen::VectorXd> const& inputs) {
  Eigen::VectorXd const& x = inputs.at(0).get();
  Eigen::VectorXd r(2);
  r(0) = x(0)/a;
  r(1) = (x(1) - b*(x(0)*x(0)+a*a))*a;

  return stdNorm->LogDensity(r);
}

Eigen::VectorXd BananaDistribution::GradLogDensity(unsigned int wrt,
                                                   ref_vector<Eigen::VectorXd> const& inputs) {

  assert(wrt==0);

  Eigen::VectorXd const& x = inputs.at(0).get();
  Eigen::VectorXd r(2);
  r(0) = x(0)/a;
  r(1) = (x(1) - b*(x(0)*x(0)+a*a))*a;

  Eigen::MatrixXd jac(2,2);
  jac << 1.0/a, 0,
         -2.0*b*a*x(0), a;

  return jac.transpose()*stdNorm->GradLogDensity(0,r);
};


Eigen::VectorXd BananaDistribution::SampleImpl(ref_vector<Eigen::VectorXd> const& inputs) {

  Eigen::VectorXd r = stdNorm->Sample();

  Eigen::VectorXd x(2);
  x << a*r(0),
       r(1)/a + b*a*a*(r(0)*r(0) + 1.0);

  return x;
}
