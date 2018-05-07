#include "MUQ/SamplingAlgorithms/MHKernel.h"
#include "MUQ/SamplingAlgorithms/MCMCProposal.h"

#include "MUQ/Utilities/AnyHelpers.h"
#include "MUQ/Utilities/RandomGenerator.h"

namespace pt = boost::property_tree;
using namespace muq::Utilities;
using namespace muq::Modeling;
using namespace muq::SamplingAlgorithms;

REGISTER_TRANSITION_KERNEL(MHKernel)

MHKernel::MHKernel(pt::ptree const& pt, std::shared_ptr<AbstractSamplingProblem> problem) : TransitionKernel(pt, problem)
{
  // Extract the proposal parts from the ptree
  std::string proposalName = pt.get<std::string>("Proposal");

  boost::property_tree::ptree subTree = pt.get_child(proposalName);
  subTree.put("BlockIndex", blockInd);

  // Construct the proposal
  proposal = MCMCProposal::Construct(subTree, problem);
  assert(proposal);
}



MHKernel::~MHKernel() {}

void MHKernel::PostStep(unsigned int const t, std::vector<std::shared_ptr<SamplingState>> const& state){
  proposal->Adapt(t,state);
}

std::vector<std::shared_ptr<SamplingState>> MHKernel::Step(std::shared_ptr<SamplingState> prevState){

  assert(proposal);

  // propose a new point
  boost::any propAny = proposal->Sample(prevState);
  std::shared_ptr<SamplingState> prop = AnyCast(propAny);

  // compute acceptance probability
  double propTarget;
  double currentTarget;

  if( prevState->HasMeta("LogTarget") ){
    currentTarget = AnyCast( prevState->meta["LogTarget"]);
  }else{
    currentTarget = problem->LogDensity(prevState);
    prevState->meta["LogTarget"] = currentTarget;
  }

  propTarget = problem->LogDensity(prop);
  prop->meta["LogTarget"] = propTarget;

  // Aceptance probability
  const double forwardPropDens = proposal->LogDensity(prevState, prop);
  const double backPropDens = proposal->LogDensity(prop, prevState);
  const double alpha = std::exp(propTarget + forwardPropDens - currentTarget - backPropDens);

  // accept/reject
  if( RandomGenerator::GetUniform()<alpha ) {
    return std::vector<std::shared_ptr<SamplingState>>(1,prop);
  } else {
    prevState->weight += 1.0;
    return std::vector<std::shared_ptr<SamplingState>>(1,prevState);
  }
}
