#include <gtest/gtest.h>

#include "MUQ/Utilities/MultiIndices/MultiIndexFactory.h"
#include "MUQ/Utilities/RandomGenerator.h"

#include "MUQ/Approximation/Polynomials/Legendre.h"
#include "MUQ/Approximation/Polynomials/Monomial.h"

#include "MUQ/Approximation/TransportMaps/PolynomialMap.h"

#include <boost/property_tree/ptree.hpp>

using namespace muq::Utilities;
using namespace muq::Approximation;

class PolynomialMapTests : public::testing::Test {
public:

  PolynomialMapTests() {
    // create the basis expension for each component of the map
    expansion.resize(dim);
    auto legendre = std::make_shared<Legendre>();
    std::vector<std::shared_ptr<MultiIndexSet>> multis = MultiIndexFactory::CreateTriTotalOrder(dim, 3);

    // create the bases, multi index, and coefficents
    std::vector<std::shared_ptr<IndexedScalarBasis> > bases(dim, legendre);

    for( unsigned int i=0; i<dim; ++i ) {

      Eigen::MatrixXd coeffs = Eigen::MatrixXd::Random(1, multis.at(i)->Size());
      // the basis expansion for component i
      expansion.at(i) = std::make_shared<BasisExpansion>(bases, multis.at(i), coeffs);
    }
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
  map = std::make_shared<PolynomialMap>(expansion);
  EXPECT_TRUE(map->inputSizes.size()==1);
  EXPECT_TRUE(map->outputSizes.size()==1);
  EXPECT_TRUE(map->inputSizes(0)==dim);
  EXPECT_TRUE(map->outputSizes(0)==dim);

  // choose a random point to evaluate the function
  const Eigen::VectorXd xpnt = Eigen::VectorXd::Random(dim);

  // evaluate the transport map
  const Eigen::VectorXd rpnt = map->EvaluateForward(xpnt);
  EXPECT_EQ(rpnt.size(), dim);
  for( unsigned int i=0; i<dim; ++i ) {
    EXPECT_DOUBLE_EQ(rpnt(i), expansion[i]->Evaluate( xpnt ).at(0)(0) );
  }

  // check the evaluate function
  const Eigen::VectorXd rpnteval = map->Evaluate(xpnt) [0];
  EXPECT_EQ(rpnteval.size(), dim);
  EXPECT_DOUBLE_EQ((rpnteval-rpnt).norm(), 0.0);
}

TEST_F(PolynomialMapTests, NewtonInverseEvaluation) {
  map = std::make_shared<PolynomialMap>(expansion, PolynomialMap::InverseMethod::Newton);
  EXPECT_TRUE(map->inputSizes.size()==1);
  EXPECT_TRUE(map->outputSizes.size()==1);
  EXPECT_TRUE(map->inputSizes(0)==dim);
  EXPECT_TRUE(map->outputSizes(0)==dim);

  // choose a random point to evaluate the function
  const Eigen::VectorXd rpnt = Eigen::VectorXd::Random(dim);

  // an initial guess
  const Eigen::VectorXd xpnt0 = Eigen::VectorXd::Zero(dim);

  // evaluate the transport map
  const Eigen::VectorXd xpnt = map->EvaluateInverse(rpnt, xpnt0);
  EXPECT_EQ(xpnt.size(), dim);

  // the inverse should be equal to the reference
  const Eigen::VectorXd result = map->EvaluateForward(xpnt);
  EXPECT_EQ(result.size(), dim);
  EXPECT_NEAR((result-rpnt).norm(), 0.0, 1.0e-10);
}

TEST_F(PolynomialMapTests, SturmInverseEvaluation) {
  map = std::make_shared<PolynomialMap>(expansion, PolynomialMap::InverseMethod::Sturm);
  EXPECT_TRUE(map->inputSizes.size()==1);
  EXPECT_TRUE(map->outputSizes.size()==1);
  EXPECT_TRUE(map->inputSizes(0)==dim);
  EXPECT_TRUE(map->outputSizes(0)==dim);

  // choose a random point to evaluate the function
  const Eigen::VectorXd rpnt = Eigen::VectorXd::Random(dim);

  // an initial guess
  const Eigen::VectorXd xpnt0 = Eigen::VectorXd::Zero(dim);

  // evaluate the transport map
  const Eigen::VectorXd xpnt = map->EvaluateInverse(rpnt, xpnt0);
  EXPECT_EQ(xpnt.size(), dim);

  // the inverse should be equal to the reference
  const Eigen::VectorXd result = map->EvaluateForward(xpnt);
  EXPECT_EQ(result.size(), dim);
  EXPECT_NEAR((result-rpnt).norm(), 0.0, 1.0e-10);
}

TEST_F(PolynomialMapTests, ComradeInverseEvaluation) {
  map = std::make_shared<PolynomialMap>(expansion, PolynomialMap::InverseMethod::Comrade); // this is also the default
  EXPECT_TRUE(map->inputSizes.size()==1);
  EXPECT_TRUE(map->outputSizes.size()==1);
  EXPECT_TRUE(map->inputSizes(0)==dim);
  EXPECT_TRUE(map->outputSizes(0)==dim);

  // choose a random point to evaluate the function
  const Eigen::VectorXd rpnt = Eigen::VectorXd::Random(dim);

  // an initial guess
  const Eigen::VectorXd xpnt0 = Eigen::VectorXd::Zero(dim);

  // evaluate the transport map
  const Eigen::VectorXd xpnt = map->EvaluateInverse(rpnt, xpnt0);
  EXPECT_EQ(xpnt.size(), dim);

  // the inverse should be equal to the reference
  const Eigen::VectorXd result = map->EvaluateForward(xpnt);
  EXPECT_EQ(result.size(), dim);
  EXPECT_NEAR((result-rpnt).norm(), 0.0, 1.0e-10);
}

TEST_F(PolynomialMapTests, LogDeterminate) {
  map = std::make_shared<PolynomialMap>(expansion);
  EXPECT_TRUE(map->inputSizes.size()==1);
  EXPECT_TRUE(map->outputSizes.size()==1);
  EXPECT_TRUE(map->inputSizes(0)==dim);
  EXPECT_TRUE(map->outputSizes(0)==dim);

  // an initial guess
  const Eigen::VectorXd xpnt = Eigen::VectorXd::Random(dim);

  // evaluate the log determinate
  const double logdet = map->LogDeterminant(xpnt);
  EXPECT_TRUE(!std::isnan(logdet));
}

TEST_F(PolynomialMapTests, CreateIdentity) {

  boost::property_tree::ptree options;
  options.put("IndexSetType","TotalOrder");
  options.put("Order",3);

  const int dim = 4;
  auto id = PolynomialMap::BuildIdentity(dim, options);

  Eigen::VectorXd randInp = RandomGenerator::GetNormal(dim);

  Eigen::VectorXd output = id->EvaluateForward(randInp);

  for(int i=0; i<dim; ++i)
    EXPECT_DOUBLE_EQ(randInp(i), output(i));

  output = id->EvaluateInverse(randInp, randInp);

  for(int i=0; i<dim; ++i)
    EXPECT_DOUBLE_EQ(randInp(i), output(i));
}
