#include <gtest/gtest.h>

#include <Eigen/Core>

#include "MUQ/Utilities/RandomGenerator.h"
#include "MUQ/Modeling/Distributions/Density.h"
#include "MUQ/Modeling/Distributions/BananaDistribution.h"

using namespace muq::Utilities;
using namespace muq::Modeling;

TEST(BananaDistributionTests, Sample) {

  double a = 0.5;
  double b = 1.2;

  std::shared_ptr<Distribution> dist = std::make_shared<BananaDistribution>(a,b);

  const unsigned int numSamps = 1000;
  Eigen::MatrixXd samps(2,numSamps);

  for(int k=0; k<numSamps; ++k)
    samps.col(k) = dist->Sample();

  Eigen::VectorXd sampMean = samps.rowwise().mean();
  EXPECT_NEAR(0.0,sampMean(0),2.0/sqrt(numSamps));

}

TEST(BananaDistributionTests, Density) {

  double a = 0.5;
  double b = 1.2;

  std::shared_ptr<Distribution> dist = std::make_shared<BananaDistribution>(a,b);
  auto stdNorm = std::make_shared<Gaussian>(2);

  Eigen::VectorXd x = Eigen::VectorXd::Random(2);
  double logDens = dist->LogDensity(x);

  Eigen::VectorXd r(2);
  r(0) = x(0)/a;
  r(1) = a*(x(1) - b*(x(0)*x(0) + a*a));

  double trueLogDens = stdNorm->LogDensity(r);

  EXPECT_DOUBLE_EQ(trueLogDens,logDens);

  Eigen::VectorXd grad = dist->GradLogDensity(0,x);
  auto xvec = std::vector<Eigen::VectorXd>(1,x);
  Eigen::VectorXd fdgrad = dist->AsDensity()->JacobianByFD(0,0,xvec).transpose();

  EXPECT_NEAR(fdgrad(0), grad(0), 1e-4);
  EXPECT_NEAR(fdgrad(1), grad(1), 1e-4);
}
