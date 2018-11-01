#include "MUQ/Optimization/Optimizer.h"
#include "MUQ/Utilities/Exceptions.h"

namespace pt = boost::property_tree;
using namespace muq::Modeling;
using namespace muq::Optimization;


Optimizer::Optimizer(std::shared_ptr<CostFunction> cost,
                     pt::ptree const& pt) :
  WorkPiece(cost->InputTypes(),
            cost->numInputs,
            std::vector<std::string>({typeid(Eigen::VectorXd).name(), typeid(double).name()})),
  ftol_rel(pt.get<double>("Ftol.AbsoluteTolerance", 1.0e-8)),
  ftol_abs(pt.get<double>("Ftol.RelativeTolerance", 1.0e-8)),
  xtol_rel(pt.get<double>("Xtol.AbsoluteTolerance", 0.0)),
  xtol_abs(pt.get<double>("Xtol.RelativeTolerance", 0.0)),
  constraint_tol(pt.get<double>("ConstraintTolerance", 1.0e-8)),
  maxEvals(pt.get<unsigned int>("MaxEvaluations", 100)) {}


Optimizer::~Optimizer() {}

void Optimizer::AddInequalityConstraint(std::vector<std::shared_ptr<ModPiece>> const& ineq) {
  ineqConstraints.insert(ineqConstraints.end(), ineq.begin(), ineq.end());
}

void Optimizer::AddInequalityConstraint(std::shared_ptr<ModPiece> const& ineq) {
  ineqConstraints.push_back(ineq);
}

void Optimizer::ClearInequalityConstraint() {
  ineqConstraints.clear();
}

void Optimizer::AddEqualityConstraint(std::vector<std::shared_ptr<ModPiece>> const& eq) {
  eqConstraints.insert(eqConstraints.end(), eq.begin(), eq.end());
}

void Optimizer::AddEqualityConstraint(std::shared_ptr<ModPiece> const& eq) {
  eqConstraints.push_back(eq);
}

void Optimizer::ClearEqualityConstraint() {
  eqConstraints.clear();
}


void Optimizer::ListMethods(std::string prefix)
{
  auto map = GetOptimizerMap();
  for(auto pair : *map)
    std::cout << prefix << pair.first << std::endl;
}

std::shared_ptr<Optimizer> Optimizer::Construct(std::shared_ptr<CostFunction> cost,
                                                boost::property_tree::ptree const& options) {

  std::string method = options.get<std::string>("Algorithm");
  Optimizer::OptimizerMap const& map = *GetOptimizerMap();
  auto iter = map.find(method);
  if(iter != map.end()) {
    return iter->second(cost,options);
  }else{
    std::stringstream msg;
    msg << "Invalid \"Algorithm\" passed to Optimizer::Construct.  A value of \"" << method << "\" was used, but valid options are:\n";
    for(auto& part : map)
      msg << "  " << part.first << "\n";
    throw std::invalid_argument(msg.str());
  }
}

std::shared_ptr<Optimizer::OptimizerMap> Optimizer::GetOptimizerMap()
{
  static std::shared_ptr<Optimizer::OptimizerMap> map;

  if( !map )
    map = std::make_shared<Optimizer::OptimizerMap>();

  return map;
}
