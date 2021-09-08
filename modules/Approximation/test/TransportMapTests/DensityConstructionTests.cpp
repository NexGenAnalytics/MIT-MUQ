#include "MUQ/Approximation/TransportMaps/PolynomialMap.h"

#include "MUQ/Modeling/Distributions/BananaDistribution.h"

#include <boost/property_tree/ptree.hpp>
#include <gtest/gtest.h>

using namespace muq::Approximation;
using namespace muq::Modeling;

TEST(DensityConstructionTests, BananaPolynomial) {

  // First, construct the BananaDensity
  double a = 1.0;
  double b = 1.0;
  std::shared_ptr<Density> nanner = std::make_shared<BananaDistribution>(a,b)->AsDensity();

  // Now, set up the options
  boost::property_tree::ptree options;
  options.put("IndexSetType","TotalOrder");
  options.put("Order", 3);
  options.put("OptBlock", "Optimizer");
  options.put("Optimizer.Algorithm", "LBFGS");

  // Construct the map
  std::shared_ptr<PolynomialMap> map = PolynomialMap::FromDensity(nanner, options);

  // Test the map-induced density vs the true density
  Eigen::VectorXd x = Eigen::VectorXd::Random(2);
  std::shared_ptr<Density> stdNorm = std::make_shared<Gaussian>(2)->AsDensity();

  double mapLogDens = stdNorm->LogDensity(map->EvaluateForward(x)) + map->LogDeterminant(x);
  double trueLogDens = nanner->LogDensity(x);

  EXPECT_NEAR(trueLogDens, mapLogDens, 1e-10);
}
