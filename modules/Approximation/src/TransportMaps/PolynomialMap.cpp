#include "MUQ/Approximation/TransportMaps/PolynomialMap.h"

#include "MUQ/Approximation/Polynomials/Monomial.h"
#include "MUQ/Approximation/Polynomials/OrthogonalPolynomial.h"

#include "MUQ/Utilities/MultiIndices/MultiIndexSet.h"
#include "MUQ/Utilities/MultiIndices/MultiIndexFactory.h"

#include <boost/property_tree/ptree.hpp>

using namespace muq::Approximation;
using namespace muq::Utilities;

PolynomialMap::PolynomialMap(std::vector<std::shared_ptr<BasisExpansion> > const& expansions, PolynomialMap::InverseMethod const& invMethod) : ConditionableMap(expansions.size()), expansions(expansions), invMethod(invMethod) {
  // make sure the output of each expansion is size 1
  for( unsigned int i=0; i<expansions.size(); ++i ) { assert(expansions[i]->outputSizes(0)==1); }
}

Eigen::VectorXd PolynomialMap::EvaluateForward(Eigen::VectorXd const& x) const {
  assert(x.size()==expansions.size());

  // evaluate each expansion
  Eigen::VectorXd result(expansions.size());
  for( unsigned int i=0; i<expansions.size(); ++i )
     result(i) = expansions[i]->Evaluate( x ).at(0)(0);

  // return the result
  return result;
}

Eigen::VectorXd PolynomialMap::EvaluateInverse(Eigen::VectorXd const& refPt, Eigen::VectorXd const& tgtPt0) const {
  assert(refPt.size()==expansions.size());
  assert(tgtPt0.size()==expansions.size());

  // set initial guess
  Eigen::VectorXd result = tgtPt0;

  for( unsigned int i=0; i<expansions.size(); ++i ) {
    if( invMethod==InverseMethod::Newton ) {
      NewtonsMethod(result, refPt(i), i);
    } else {
      PolynomialRootFinding(result, refPt(i), i, invMethod);
    }
  }

  // return the result
  return result;
}

void PolynomialMap::PolynomialRootFinding(Eigen::VectorXd& result, double const refPt, unsigned int const component, PolynomialMap::InverseMethod const& invMethod) const {

  const Eigen::MatrixXd& coeff = expansions[component]->GetCoeffs();
  const std::vector<std::shared_ptr<IndexedScalarBasis> >& basisComps = expansions[component]->BasisComponents();
  Eigen::VectorXd polyCoeff = Eigen::VectorXd::Zero(expansions[component]->Multis()->GetMaxOrders() (component)+1);
  for( unsigned int term=0; term<expansions[component]->NumTerms(); ++term ) {
    double scale = coeff(term);
    int p = 0;

    // get the coefficients for non costant terms
    for( auto it=expansions[component]->Multis()->at(term)->GetNzBegin(); it!=expansions[component]->Multis()->at(term)->GetNzEnd(); ++it ) {
      if( it->first==component ) { // get the coeffs for this basis
        p = it->second;
      } else {
        assert(it->first<component);
        scale *= basisComps[it->first]->BasisEvaluate(it->second, result(it->first));
      }
    }

    polyCoeff(p) += scale;
  }

  // scale the constant coeff. by the reference point
  polyCoeff(0) -= refPt/basisComps[component]->BasisEvaluate(0, 0.0);

  // choose the closest root to the input point
  auto basis = std::dynamic_pointer_cast<OrthogonalPolynomial>(basisComps[component]);
  assert(basis);
  Eigen::VectorXd roots = invMethod==InverseMethod::Comrade? basis->GetRootsComrade(polyCoeff) : basis->GetRootsSturm(polyCoeff, tol);
  assert(roots.size()>0); // make sure that we have at least one root
  const Eigen::VectorXd diff = (roots-Eigen::VectorXd::Constant(roots.size(), result(component))).array().abs();
  int rt;
  diff.minCoeff(&rt);

  result(component) = roots(rt);
}

void PolynomialMap::NewtonsMethod(Eigen::VectorXd& result, double const refPt, unsigned int const component) const {
  // initialize counter and a large error
  unsigned int cnt = 0;
  double eval = 1.0;

  // Newton's method
  while( std::fabs(eval)>tol && cnt++<maxit ) {
    eval = expansions[component]->Evaluate( result ).at(0)(0) - refPt;
    const double deriv = expansions[component]->Derivative(component, result) (0);
    assert(std::fabs(deriv)>tol);
    result(component) -= eval/deriv;
  }
}

double PolynomialMap::LogDeterminant(Eigen::VectorXd const& evalPt) const {
  assert(evalPt.size()==expansions.size());

  double logdet = 1.0;
  for( unsigned int i=0; i<expansions.size(); ++i ) {
    logdet += std::log(std::fabs(expansions[i]->Derivative(i, evalPt )(0)));
  }

  return logdet;
}


std::shared_ptr<PolynomialMap> PolynomialMap::BuildIdentity(unsigned int dim,
                                                            boost::property_tree::ptree& options) {

  const std::string setType = options.get("IndexSetType","TotalOrder");
  const unsigned int order = options.get<int>("Order");

  // Set up the multiIndices
  std::vector<std::shared_ptr<MultiIndexSet>> multis;
  if(setType == "TotalOrder"){
    multis = MultiIndexFactory::CreateTriTotalOrder(dim, order);

  }else if(setType=="Hyperbolic"){
    const double q = options.get("NormScale", 0.5);
    multis = MultiIndexFactory::CreateTriHyperbolic(dim, order, q);

  } else {
    std::stringstream msg;
    msg << "In PolynomialMap::BuildIdentity, invalid value for option \"IndexSetType\".";
    msg << "Received \"" << setType << "\", but valid options are either \"TotalOrder\" or \"Hyperbolic\".";
    throw std::invalid_argument(msg.str());
  }

  // Set up the basis functions
  const std::string basisType = options.get("BasisType", "ProbabilistHermite");
  std::shared_ptr<IndexedScalarBasis> basis = IndexedScalarBasis::Construct(basisType);

  // Construct the expansions by putting together the basis and multiindices
  std::vector<std::shared_ptr<IndexedScalarBasis>> basisVec = std::vector<std::shared_ptr<IndexedScalarBasis>>(dim,basis);
  std::vector<std::shared_ptr<BasisExpansion>> expansions(dim);

  for(int d=0; d<dim; ++d){
    expansions.at(d) = std::make_shared<BasisExpansion>(basisVec, multis.at(d));

    // find the linear component
    std::shared_ptr<MultiIndex> idMulti = std::make_shared<MultiIndex>(dim);
    idMulti->SetValue(d,1);
    int idInd = expansions.at(d)->Multis()->MultiToIndex(idMulti);
    expansions.at(d)->SetCoeff(0,idInd,1.0/basis->BasisEvaluate(1,1.0));
  }

  return std::make_shared<PolynomialMap>(expansions, InverseMethod::Comrade);
}

std::shared_ptr<PolynomialMap> PolynomialMap::BuildFromSamples(Eigen::MatrixXd const& samps,
                                                               boost::property_tree::ptree& options) {

  const unsigned int mapDim = samps.rows();
  const unsigned int numSamps = samps.cols();

  // Initialize all the dimensions to  with the Identity map.
  std::shared_ptr<PolynomialMap> map = BuildIdentity(mapDim, options);

  // Loop over each dimension of the map -- the map coefficients can be computed independently for each output
  for(int d=0; d<mapDim; ++d){

    // Set up Vandermonde and Derivative matrices

    // Set up optimization problem

    // Call Optimizer to compute the coefficients

  }

  // return the map
  return map;

}
