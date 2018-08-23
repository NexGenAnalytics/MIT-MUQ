#include "MUQ/SamplingAlgorithms/ExpensiveSamplingProblem.h"

#include "MUQ/Utilities/RandomGenerator.h"

#include "MUQ/SamplingAlgorithms/SamplingState.h"

namespace pt = boost::property_tree;
using namespace muq::Utilities;
using namespace muq::Approximation;
using namespace muq::SamplingAlgorithms;

ExpensiveSamplingProblem::ExpensiveSamplingProblem(std::shared_ptr<muq::Modeling::ModPiece> target, pt::ptree const& pt) : SamplingProblem(target) {
  // create the local regressor
  reg = std::make_shared<LocalRegression>(target, pt.get_child(pt.get<std::string>("RegressionOptions")));

  beta = std::pair<double, double>(pt.get<double>("BetaScale", 0.0), -pt.get<double>("BetaExponent", RAND_MAX));
  assert(beta.second<0.0);
}

double ExpensiveSamplingProblem::LogDensity(unsigned int const step, std::shared_ptr<SamplingState> state) {
  RefineSurrogate(step);
  
  return target->Evaluate(state->state).at(0)(0);
}

void ExpensiveSamplingProblem::RefineSurrogate(unsigned int const step) {
  // BETA refinement
  if( RandomGenerator::GetUniform()<beta.first*std::pow((double)step, beta.second) ) {
    std::cout << "BETA refinement" << std::endl;
  }
  
  // check to see if we should increment the level
}
