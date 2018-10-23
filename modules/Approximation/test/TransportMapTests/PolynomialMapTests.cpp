#include <gtest/gtest.h>

#include "MUQ/Utilities/MultiIndices/MultiIndexFactory.h"

#include "MUQ/Approximation/Polynomials/Legendre.h"

#include "MUQ/Approximation/TransportMaps/PolynomialMap.h"

using namespace muq::Utilities;
using namespace muq::Approximation;

class PolynomialMapTests : public::testing::Test {
public:

  PolynomialMapTests() {
    // create the basis expension for each component of the map
    expansion.resize(dim);
    auto legendre = std::make_shared<Legendre>();
    for( unsigned int i=1; i<=dim; ++i ) {
      // create the bases, multi index, and coefficents
      std::vector<std::shared_ptr<IndexedScalarBasis> > bases(i, legendre);
      std::shared_ptr<MultiIndexSet> multis = MultiIndexFactory::CreateTotalOrder(i, 3);
      Eigen::MatrixXd coeffs = Eigen::MatrixXd::Ones(1, multis->Size());

      // the basis expansion for component i
      expansion[i-1] = std::make_shared<BasisExpansion>(bases, multis, coeffs);
    }

    map = std::make_shared<PolynomialMap>(expansion);
    EXPECT_TRUE(map->inputSizes.size()==1);
    EXPECT_TRUE(map->outputSizes.size()==1);
    EXPECT_TRUE(map->inputSizes(0)==dim);
    EXPECT_TRUE(map->outputSizes(0)==dim);
  }

  virtual ~PolynomialMapTests() = default;

protected:

  /// The dimension of the transport map
  const unsigned int dim = 5;

  /// Polynomial basis
  std::vector<std::shared_ptr<BasisExpansion> > expansion;

  /// Polynomial transport map
  std::shared_ptr<PolynomialMap> map;
};

TEST_F(PolynomialMapTests, ForwardEvaluation) {
  // choose a random point to evaluate the function
  const Eigen::VectorXd xpnt = Eigen::VectorXd::Random(dim);

  // evaluate the transport map
  const Eigen::VectorXd rpnt = map->EvaluateForward(xpnt);
  EXPECT_EQ(rpnt.size(), dim);
  for( unsigned int i=0; i<dim; ++i ) {
    EXPECT_DOUBLE_EQ(rpnt(i), boost::any_cast<Eigen::VectorXd>(expansion[i]->Evaluate((Eigen::VectorXd)xpnt.head(i+1)) [0]) (0));
  }

  // check the evaluate function
  const Eigen::VectorXd rpnteval = map->Evaluate(xpnt) [0];
  EXPECT_EQ(rpnteval.size(), dim);
  EXPECT_DOUBLE_EQ((rpnteval-rpnt).norm(), 0.0);
}

TEST_F(PolynomialMapTests, InverseEvaluation) {
  // choose a random point to evaluate the function
  const Eigen::VectorXd rpnt = Eigen::VectorXd::Random(dim);

  // an initial guess
  const Eigen::VectorXd xpnt0 = Eigen::VectorXd::Random(dim);

  // evaluate the transport map
  const Eigen::VectorXd xpnt = map->EvaluateInverse(rpnt, xpnt0);
}
