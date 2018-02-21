#include "MUQ/Approximation/Polynomials/BasisExpansion.h"

#include "MUQ/Utilities/MultiIndices/MultiIndexFactory.h"

using namespace muq::Approximation;
using namespace muq::Utilities;

BasisExpansion::BasisExpansion(std::vector<std::shared_ptr<IndexedScalarBasis>> const& basisCompsIn) :
                              BasisExpansion(basisCompsIn,
                                             MultiIndexFactory::CreateTotalOrder(basisCompsIn.size(),0))
{
}

BasisExpansion::BasisExpansion(std::vector<std::shared_ptr<IndexedScalarBasis>> const& basisCompsIn,
                               std::shared_ptr<muq::Utilities::MultiIndexSet>          multisIn) :
                               BasisExpansion(basisCompsIn,
                                              multisIn,
                                              Eigen::MatrixXd::Zero(1,multisIn->Size()))
{
};

BasisExpansion::BasisExpansion(std::vector<std::shared_ptr<IndexedScalarBasis>> const& basisCompsIn,
                               std::shared_ptr<muq::Utilities::MultiIndexSet>          multisIn,
                               Eigen::MatrixXd                                  const& coeffsIn) :
                               basisComps(basisCompsIn),
                               multis(multisIn),
                               coeffs(coeffsIn)
{
  assert(basisComps.size() == multis->GetMultiLength());
  assert(multis->Size() == coeffs.cols());
}

Eigen::VectorXd BasisExpansion::GetAllTerms(Eigen::VectorXd const& x) const{

  // Get the maximum orders
  Eigen::VectorXi maxOrders = multis->GetMaxOrders();

  // Evaluate each dimension up to the maximum order
  std::vector<std::vector<double>> uniEvals(basisComps.size());
  assert(uniEvals.size() == maxOrders.size());

  for(int i=0; i<uniEvals.size(); ++i){
    uniEvals.at(i).resize(maxOrders(i)+1);
    for(int j=0; j<=maxOrders(i); ++j){
      uniEvals.at(i).at(j) = basisComps.at(i)->BasisEvaluate(j, x(i));
    }
  }

  // Now that we have all the univariate terms evaluated, evaluate the expansion
  Eigen::VectorXd allTerms = Eigen::VectorXd::Ones(multis->Size());
  for(int i=0; i<multis->Size(); ++i){

    for(auto it = multis->at(i)->GetNzBegin(); it != multis->at(i)->GetNzEnd(); ++it)
      allTerms(i) *= uniEvals.at(it->first).at(it->second);
  }

  return allTerms;
}

void BasisExpansion::EvaluateImpl(muq::Modeling::ref_vector<boost::any> const& inputs) {

  // If there are two inputs, then the second input will reset he coefficients
  if(inputs.size()>1){
      Eigen::MatrixXd const& newCoeffs = *boost::any_cast<Eigen::MatrixXd>(&inputs.at(1).get());
      assert(newCoeffs.rows() == coeffs.rows());
      assert(newCoeffs.cols() == coeffs.cols());

      coeffs = newCoeffs;
  }else if(inputs.size()==0){
    throw std::logic_error("Could not evaluate BasisExpansion because no input point was provided.  BasisExpansion::EvaluateImpl requires at least 1 input.");
  }

  // Extract the point where we want to evaluate the expansion

  Eigen::VectorXd const& x = *boost::any_cast<Eigen::VectorXd>(&inputs.at(0).get());

  // Compute the output
  outputs.resize(1);
  outputs.at(0) = (coeffs*GetAllTerms(x)).eval();

}
