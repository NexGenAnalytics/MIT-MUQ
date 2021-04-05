#include <iostream>
#include <string>

#include <boost/property_tree/ptree.hpp>

#include <MUQ/Modeling/OneStepCachePiece.h>
#include <MUQ/Modeling/ModPiece.h>
#include <MUQ/Modeling/ModGraphPiece.h>

#include <MUQ/SamplingAlgorithms/SamplingProblem.h>
#include <MUQ/SamplingAlgorithms/SingleChainMCMC.h>
#include <MUQ/SamplingAlgorithms/MCMCFactory.h>
#include "MUQ/SamplingAlgorithms/CrankNicolsonProposal.h"
#include "MUQ/SamplingAlgorithms/ParallelizableMIComponentFactory.h"
#include "MUQ/SamplingAlgorithms/SamplingProblem.h"
#include "MUQ/SamplingAlgorithms/SubsamplingMIProposal.h"
#include "MUQ/SamplingAlgorithms/AMProposal.h"

#include "HTTPModPiece.h"
#include "ExaHyPEModelGraph.h"

namespace pt = boost::property_tree;

using namespace muq::Utilities;
using namespace muq::Modeling;
using namespace muq::SamplingAlgorithms;


class MyInterpolation : public MIInterpolation {
public:
  std::shared_ptr<SamplingState> Interpolate (std::shared_ptr<SamplingState> const& coarseProposal, std::shared_ptr<SamplingState> const& fineProposal) override {
    return std::make_shared<SamplingState>(coarseProposal->state);
  }
};

class MyMIComponentFactory : public ParallelizableMIComponentFactory {
public:
  MyMIComponentFactory(pt::ptree pt, std::shared_ptr<parcer::Communicator> globalComm) : _pt(pt), _globalComm(globalComm) {}

  virtual std::shared_ptr<MCMCProposal> Proposal (std::shared_ptr<MultiIndex> const& index, std::shared_ptr<AbstractSamplingProblem> const& samplingProblem) override {
    pt::ptree pt;
    pt.put("BlockIndex",0);

    /*Eigen::VectorXd mu = Eigen::VectorXd::Zero(NUM_PARAM);
    Eigen::MatrixXd cov = Eigen::MatrixXd::Identity(NUM_PARAM,NUM_PARAM);

    auto prior = std::make_shared<Gaussian>(mu, cov);

    return std::make_shared<CrankNicolsonProposal>(pt, samplingProblem, prior);*/

    pt.put("AdaptStart", 100);
    pt.put("AdaptSteps", 100);
    pt.put("InitialVariance", 100);
    return std::make_shared<AMProposal>(pt, samplingProblem);

    /*auto mu = Eigen::VectorXd::Zero(NUM_PARAM);
    Eigen::MatrixXd cov = Eigen::MatrixXd::Identity(NUM_PARAM, NUM_PARAM);
    cov *= 0.5;

    auto prior = std::make_shared<Gaussian>(mu, cov);

    return std::make_shared<MHProposal>(pt, samplingProblem, prior);*/
  }

  void SetComm(std::shared_ptr<parcer::Communicator> const& comm) override {
    _comm = comm;
  }

  virtual std::shared_ptr<MultiIndex> FinestIndex() override {
    auto index = std::make_shared<MultiIndex>(1);
    index->SetValue(0, 0);
    return index;
  }

  virtual std::shared_ptr<MCMCProposal> CoarseProposal (std::shared_ptr<MultiIndex> const& fineIndex,
                                                        std::shared_ptr<MultiIndex> const& coarseIndex,
                                                        std::shared_ptr<AbstractSamplingProblem> const& coarseProblem,
                                                        std::shared_ptr<SingleChainMCMC> const& coarseChain) override {
    pt::ptree ptProposal = _pt;
    ptProposal.put("BlockIndex",0);
    return std::make_shared<SubsamplingMIProposal>(ptProposal, coarseProblem, coarseIndex, coarseChain);
  }

  virtual std::shared_ptr<AbstractSamplingProblem> SamplingProblem (std::shared_ptr<MultiIndex> const& index) override {
    // TODO!!
    //return std::make_shared<MySamplingProblem>(_comm,_globalComm,index);
  }

  virtual std::shared_ptr<MIInterpolation> Interpolation (std::shared_ptr<MultiIndex> const& index) override {
    return std::make_shared<MyInterpolation>();
  }

  virtual Eigen::VectorXd StartingPoint (std::shared_ptr<MultiIndex> const& index) override {
    //Starting guess: zero
    Eigen::VectorXd start = Eigen::VectorXd::Ones(2);
    start(0) = 10.0;
    start(1) = 10.0;
    return start;
  }

private:
  std::shared_ptr<parcer::Communicator> _comm;
  std::shared_ptr<parcer::Communicator> _globalComm;
  pt::ptree _pt;

};

int main(int argc, char* argv[])
{
  if (argc <= 1) {
    std::cerr << "Call with: ./binary HOSTNAME:PORT [BEARER_TOKEN]" << std::endl;
    exit(-1);
  }
  std::string host(argv[1]);

  std::string bearer_token = "";
  if (argc > 2) {
    bearer_token = std::string(argv[2]);
  } else {
    std::cout << "No bearer token was passed. Proceeding without token.";
  }
  int level = 0;
  auto graph = createWorkGraph(host, bearer_token, level);
  auto posterior = graph->CreateModPiece("likelihood");
  auto qoi = graph->CreateModPiece("qoi");

  auto samplingProblem = std::make_shared<SamplingProblem>(posterior, qoi);

  pt::ptree pt;
  pt.put("NumSamples", 100); // number of MCMC steps
  pt.put("BurnIn", 0);
  pt.put("PrintLevel",3);
  pt.put("KernelList", "Kernel1"); // Name of block that defines the transition kernel
  pt.put("Kernel1.Method","MHKernel");  // Name of the transition kernel class
  pt.put("Kernel1.Proposal", "MyProposal"); // Name of block defining the proposal distribution
  pt.put("Kernel1.MyProposal.Method", "AMProposal"); // Name of proposal class
  pt.put("Kernel1.MyProposal.InitialVariance", 500);
  pt.put("Kernel1.MyProposal.AdaptSteps", 100);
  pt.put("Kernel1.MyProposal.AdaptStart", 100);
  pt.put("Kernel1.MyProposal.AdaptEnd", 1000);

  Eigen::VectorXd startPt = Eigen::VectorXd::Zero(2);
  auto mcmc = MCMCFactory::CreateSingleChain(pt, samplingProblem);

  mcmc->Run(startPt);

  std::cout << "Sample mu: " << std::endl << mcmc->GetSamples()->Mean() << std::endl;
  std::cout << "Sample sigma: " << std::endl << mcmc->GetSamples()->Variance() << std::endl;

  mcmc->GetSamples()->WriteToFile("mcmc-samples.hdf5");
  mcmc->GetQOIs()->WriteToFile("mcmc-qois.hdf5");

  return 0;
}