#include <gtest/gtest.h>

#include "MUQ/Utilities/MultiIndices/MultiIndexFactory.h"

#include "MUQ/Approximation/Polynomials/Legendre.h"

#include "MUQ/Approximation/TransportMaps/PolynomialMap.h"

using namespace muq::Utilities;
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
      Eigen::MatrixXd coeffs = Eigen::MatrixXd::Ones(1, multis->Size());

      // the basis expansion for component i
      expansion[i-1] = std::make_shared<BasisExpansion>(basis, multis, coeffs);
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
  const Eigen::VectorXd xpnt0 = Eigen::VectorXd::Zero(dim);

  // evaluate the transport map
  const Eigen::VectorXd xpnt = map->EvaluateInverse(rpnt, xpnt0);
  EXPECT_EQ(xpnt.size(), dim);

  // the inverse should be equal to the reference
  const Eigen::VectorXd result = map->EvaluateForward(xpnt);
  EXPECT_EQ(result.size(), dim);
  EXPECT_NEAR((result-rpnt).norm(), 0.0, 1.0e-14);
}


TEST_F(PolynomialMapTests, LogDeterminate) {
  // an initial guess
  const Eigen::VectorXd xpnt = Eigen::VectorXd::Random(dim);

  // evaluate the log determinate
  const double logdet = map->LogDeterminant(xpnt);
  EXPECT_TRUE(logdet>0.0);
}


TEST_F(PolynomialMapTests, Densities) {
  
  Eigen::VectorXd c = Eigen::VectorXd::Zero(dim);
  Eigen::MatrixXd L = Eigen::MatrixXd::Zero(dim, dim);

  // create the basis, multi index, and coefficents
  auto legendre = std::make_shared<Legendre>();

  for (int i=1; i<=dim; ++i) {

    std::vector<std::shared_ptr<IndexedScalarBasis> > basis(i, legendre);
    std::shared_ptr<MultiIndexSet> multis = MultiIndexFactory::CreateTotalOrder(i,1);
    Eigen::MatrixXd coeffs = Eigen::MatrixXd::Random(1, multis->Size()).cwiseAbs();

    c(i-1) = coeffs(0,0);
    
    for (int j=0; j<multis->Size()-1; ++j) {
      L(i-1,j) = coeffs(0,j+1);
    }

    expansion[i-1] = std::make_shared<BasisExpansion>(basis, multis, coeffs);

  }

  auto densitymap = std::make_shared<PolynomialMap>(expansion);

  // choose a random point to evaluate the function
  Eigen::VectorXd rpnt = Eigen::VectorXd::Random(dim);
  Eigen::VectorXd rmu = Eigen::VectorXd::Zero(dim);  
  Eigen::MatrixXd rcov = Eigen::MatrixXd::Zero(dim, dim);
  for (int i=0; i<dim; ++i) 
    rcov(i,i) = 1.0;

  const Eigen::VectorXd xpnt = densitymap->EvaluateForward(rpnt);

  //  const Eigen::MatrixXd L = L_inv.inverse();
  const Eigen::VectorXd xmu = L*rmu + c;
  const Eigen::MatrixXd xcov = L*rcov*L.transpose();

  std::cout << "rpnt: " << rpnt << "\n" << std::endl;
  std::cout << "rmu: " << rpnt << "\n" << std::endl;
  std::cout << "xpnt: " << xpnt << "\n" << std::endl;
  std::cout << "xmu: " << xmu << "\n" << std::endl;

  std::cout << exp(map->LogDeterminant(rpnt)) << std::endl;

  double  A = 1.0/std::sqrt((2.0*M_PI*xcov).determinant());
  std::cout << A << std::endl;
  std::cout << A*exp(-0.5*(xpnt-xmu).transpose()*xcov.inverse()*(xpnt-xmu)) << std::endl;

  A = 1.0/std::sqrt((2.0*M_PI*rcov).determinant());
  std::cout << exp(map->LogDeterminant(rpnt))*A*exp(-0.5*(rpnt-rmu).transpose()*rcov.inverse()*(rpnt-rmu)) << std::endl;

  
}
    
  
  
  
