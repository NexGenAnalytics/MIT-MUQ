#include <gtest/gtest.h>

#include "MUQ/Utilities/MultiIndices/MultiIndexFactory.h"

#include "MUQ/Modeling/Distributions/Gaussian.h"

#include "MUQ/Approximation/Polynomials/Legendre.h"
#include "MUQ/Approximation/Polynomials/Monomial.h"

#include "MUQ/Approximation/TransportMaps/PolynomialMap.h"

using namespace muq::Utilities;
using namespace muq::Modeling;
using namespace muq::Approximation;

class PolynomialMapTests : public::testing::Test {
public:

  PolynomialMapTests() {

    // create the basis expansion for each component of the map
    expansion.resize(dim);
    auto legendre = std::make_shared<Legendre>();

    for( unsigned int i=1; i<=dim; ++i ) {
      // create the basis, multi index, and coefficents
      std::vector<std::shared_ptr<IndexedScalarBasis> > basis(i, legendre);
      std::shared_ptr<MultiIndexSet> multis = MultiIndexFactory::CreateTotalOrder(i, 3);
      Eigen::MatrixXd coeffs = Eigen::MatrixXd::Random(1, multis->Size());

      // the basis expansion for component i
      expansion[i-1] = std::make_shared<BasisExpansion>(basis, multis, coeffs);
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
    EXPECT_DOUBLE_EQ(rpnt(i), boost::any_cast<Eigen::VectorXd>(expansion[i]->Evaluate((Eigen::VectorXd)xpnt.head(i+1)) [0]) (0));
  }

  // check the evaluate function
  const Eigen::VectorXd rpnteval = map->Evaluate(xpnt) [0];
  EXPECT_EQ(rpnteval.size(), dim);
  EXPECT_DOUBLE_EQ((rpnteval-rpnt).norm(), 0.0);

}

TEST_F(PolynomialMapTests, NewtonInverseEvaluation) {
  map = std::make_shared<PolynomialMap>(expansion, PolynomialMap::Newton);
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
  map = std::make_shared<PolynomialMap>(expansion, PolynomialMap::Sturm);
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
  map = std::make_shared<PolynomialMap>(expansion, PolynomialMap::Comrade); // this is also the default
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
  // create the basis expension for each component of the map
  expansion.resize(dim);
  auto legendre = std::make_shared<Legendre>();
  for( unsigned int i=1; i<=dim; ++i ) {
    // create the bases, multi index, and coefficents
    std::vector<std::shared_ptr<IndexedScalarBasis> > bases(i, legendre);
    std::shared_ptr<MultiIndexSet> multis = MultiIndexFactory::CreateTotalOrder(i, 1);
    Eigen::MatrixXd coeffs = Eigen::MatrixXd::Random(1, multis->Size());

    // the basis expansion for component i
    expansion[i-1] = std::make_shared<BasisExpansion>(bases, multis, coeffs);
  }

  map = std::make_shared<PolynomialMap>(expansion);
  EXPECT_TRUE(map->inputSizes.size()==1);
  EXPECT_TRUE(map->outputSizes.size()==1);
  EXPECT_TRUE(map->inputSizes(0)==dim);
  EXPECT_TRUE(map->outputSizes(0)==dim);

  // an initial guess
  const Eigen::VectorXd rpnt0 = Eigen::VectorXd::Random(dim);

  // evaluate the log determinate
  const double logdet0 = map->LogDeterminant(rpnt0);

  // an initial guess
  const Eigen::VectorXd rpnt1 = Eigen::VectorXd::Random(dim);

  // evaluate the log determinate
  const double logdet1 = map->LogDeterminant(rpnt1);

  // linear determinates should be the same
  EXPECT_DOUBLE_EQ(logdet0, logdet1);
}

TEST_F(PolynomialMapTests, ScaledDensity) {

  Eigen::VectorXd c = Eigen::VectorXd::Zero(dim);
  Eigen::MatrixXd L = Eigen::MatrixXd::Zero(dim, dim);

  // create the basis, multi index, and coefficents
  auto mono = std::make_shared<Monomial>();

  for (int i=1; i<=dim; ++i) {
    std::vector<std::shared_ptr<IndexedScalarBasis> > basis(i, mono);
    std::shared_ptr<MultiIndexSet> multis = MultiIndexFactory::CreateTotalOrder(i,1);
    Eigen::MatrixXd coeffs = Eigen::MatrixXd::Random(1, multis->Size());

    c(i-1) = coeffs(0,0);

    for (int j=0; j<multis->Size()-1; ++j) {
      L(i-1,i-1-j) = coeffs(0,j+1);
    }

    expansion[i-1] = std::make_shared<BasisExpansion>(basis, multis, coeffs);
  }

  map = std::make_shared<PolynomialMap>(expansion);

  // reference density
  auto stdnormal = std::make_shared<Gaussian>(dim);

  // make sure we get what we expect
  const Eigen::VectorXd rpnt = stdnormal->Sample();
  const Eigen::VectorXd xpnt = map->EvaluateForward(rpnt);
  const Eigen::VectorXd xexpect = L*rpnt + c;
  EXPECT_NEAR((xpnt-xexpect).norm(), 0.0, 1.0e-10);

  // target density
  auto target = std::make_shared<Gaussian>(c, L*L.transpose());

  // the target density should be defined by the reference denisty and the inverse transport map; note: \pi_X(x) = \pi_X(T(r)) = \pi_R(r) |\det{T^{-1}(r)}|
  EXPECT_NEAR(target->LogDensity(xpnt), stdnormal->LogDensity(rpnt) - map->LogDeterminant(rpnt), 1.0e-10);

}


TEST_F(PolynomialMapTests, NonlinearDensity) {

  double a = 2.0;
  double b = 0.5;
  
  std::vector<std::shared_ptr<BasisExpansion> > expansions(2);
  
  // create the basis, multi index, and coefficents
  auto poly = std::make_shared<Legendre>();

  // Create first map T(x_1)
  std::vector<std::shared_ptr<IndexedScalarBasis> > basis0(1, poly);
  std::shared_ptr<MultiIndexSet> multis = MultiIndexFactory::CreateTotalOrder(1,1);
  Eigen::MatrixXd coeffs = Eigen::MatrixXd::Zero(1, multis->Size());
  
  coeffs(0,1) = a;
  expansions[0] = std::make_shared<BasisExpansion>(basis0, multis, coeffs);


  // Create second map T(x_1, x_2)
  std::vector<std::shared_ptr<IndexedScalarBasis> > basis1(2, poly);
  multis = MultiIndexFactory::CreateTotalOrder(2,2);
  coeffs = Eigen::MatrixXd::Zero(1, multis->Size());

  coeffs(0,0) = a*a*b+1.0/3.0*a*a*b;
  coeffs(0,1) = 1.0/a;
  coeffs(0,5) = 2.0/3.0*a*a*b;

  expansions[1] = std::make_shared<BasisExpansion>(basis1, multis, coeffs);

  // Create Polynomial map
  map = std::make_shared<PolynomialMap>(expansions, PolynomialMap::Newton);

  // Test push forward evaluation
  auto stdnormal = std::make_shared<Gaussian>(2);

  const Eigen::VectorXd rpnt = stdnormal->Sample();
  const Eigen::VectorXd xpnt = map->EvaluateForward(rpnt);

  const Eigen::Vector2d xexpect(a*rpnt(0),
                                1.0/a*rpnt(1)+a*a*b*(rpnt(0)*rpnt(0)+1.0));

  EXPECT_NEAR((xpnt-xexpect).norm(), 0.0, 1.0e-10);

  // Test pull back evaluation and density
  const Eigen::VectorXd rpnt0 = Eigen::VectorXd::Zero(2);
  const Eigen::VectorXd rmap = map->EvaluateInverse(xpnt, rpnt);

  const Eigen::VectorXd rexpect =
    Eigen::Vector2d(1.0/a*xpnt(0),
                    a*xpnt(1) - a*b*(xpnt(0)*xpnt(0)+a*a));

  EXPECT_NEAR((rmap-rexpect).norm(), 0.0, 1.0e-10);
  EXPECT_NEAR(stdnormal->LogDensity(rexpect),
              stdnormal->LogDensity(rmap), 1.0e-10);
  EXPECT_NEAR(map->LogDeterminant(rmap), 0.0, 1.0e-10);
  
}
