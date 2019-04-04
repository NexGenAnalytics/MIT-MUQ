#include <gtest/gtest.h>

#include "MUQ/Approximation/Quadrature/AdaptiveSmolyakQuadrature.h"
#include "MUQ/Approximation/Quadrature/SmolyakQuadrature.h"
#include "MUQ/Approximation/Quadrature/GaussQuadrature.h"
#include "MUQ/Approximation/Polynomials/Legendre.h"

#include "MUQ/Modeling/CwiseOperators/CwiseUnaryOperator.h"
#include "MUQ/Utilities/MultiIndices/MultiIndexFactory.h"

#include <Eigen/Core>

using namespace muq::Approximation;
using namespace muq::Utilities;
using namespace muq::Modeling;


TEST(Quadrature, AdaptiveSmolyak_GaussQuad) {

  unsigned int dim = 2;

  // Define the model
  auto model = std::make_shared<CosOperator>(dim);

  // Define the quadrature rules
  auto poly1d = std::make_shared<Legendre>();
  auto quad1d = std::make_shared<GaussQuadrature>(poly1d);

  AdaptiveSmolyakQuadrature smolyQuad(model, {quad1d, quad1d});

  unsigned int maxOrder = 5;
  auto multiSet = MultiIndexFactory::CreateTotalOrder(dim, maxOrder);

  Eigen::VectorXd res = smolyQuad.Compute(multiSet);

  EXPECT_NEAR(4.0*std::sin(1.0), res(0), 1e-10);
  EXPECT_NEAR(4.0*std::sin(1.0), res(1), 1e-10);
}
