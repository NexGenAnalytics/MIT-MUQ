#include <gtest/gtest.h>

#include "MUQ/Approximation/PolynomialChaos/AdaptiveSmolyakPCE.h"
#include "MUQ/Approximation/Quadrature/GaussQuadrature.h"
#include "MUQ/Approximation/Polynomials/Legendre.h"

#include "MUQ/Modeling/CwiseOperators/CwiseUnaryOperator.h"

#include "MUQ/Utilities/MultiIndices/MultiIndexFactory.h"
#include "MUQ/Utilities/RandomGenerator.h"

#include <Eigen/Core>

using namespace muq::Approximation;
using namespace muq::Utilities;
using namespace muq::Modeling;


TEST(PolynomialChaos, AdaptiveSmolyak_Cos2d) {

  unsigned int dim = 2;

  // Define the model
  auto model = std::make_shared<CosOperator>(dim);

  // Define the quadrature rules
  auto poly1d = std::make_shared<Legendre>();
  auto quad1d = std::make_shared<GaussQuadrature>(poly1d);

  AdaptiveSmolyakPCE smolyPCE(model, {quad1d, quad1d}, {poly1d,poly1d});

  unsigned int maxOrder = 5;
  auto multiSet = MultiIndexFactory::CreateTotalOrder(dim, maxOrder);

  std::shared_ptr<PolynomialChaosExpansion> pce = smolyPCE.Compute(multiSet);

  Eigen::VectorXd testPt = 2.0*RandomGenerator::GetUniform(dim) - Eigen::VectorXd::Ones(dim);

  Eigen::VectorXd pceOut = pce->Evaluate(testPt).at(0);
  Eigen::VectorXd trueOut = model->Evaluate(testPt).at(0);
  for(unsigned int d=0; d<dim; ++d)
    EXPECT_NEAR(trueOut(d),pceOut(d), 1e-4);

}
