#include <iostream>
#include <string>

#include <boost/property_tree/ptree.hpp>

#include <MUQ/Modeling/OneStepCachePiece.h>
#include <MUQ/Modeling/ModPiece.h>
#include <MUQ/Modeling/ModGraphPiece.h>

#include <MUQ/SamplingAlgorithms/ExpensiveSamplingProblem.h>
#include <MUQ/SamplingAlgorithms/MHKernel.h>
#include <MUQ/SamplingAlgorithms/SamplingProblem.h>
#include <MUQ/SamplingAlgorithms/SingleChainMCMC.h>
#include <MUQ/SamplingAlgorithms/MCMCFactory.h>
#include "MUQ/SamplingAlgorithms/CrankNicolsonProposal.h"
#include "MUQ/SamplingAlgorithms/ParallelizableMIComponentFactory.h"
#include "MUQ/SamplingAlgorithms/SamplingProblem.h"
#include "MUQ/SamplingAlgorithms/SubsamplingMIProposal.h"
#include "MUQ/SamplingAlgorithms/AMProposal.h"
#include "MUQ/SamplingAlgorithms/ParallelFixedSamplesMIMCMC.h"
#include <MUQ/Modeling/HTTPModel/HTTPModPiece.h>

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

class MyMIComponentFactory : public MIComponentFactory {
public:
  MyMIComponentFactory(pt::ptree pt, std::string host)
   : _pt(pt),
     host(host) {}

  ~MyMIComponentFactory() {
    std::cout << "==== qoi evals" << std::endl;
    for (auto qoi : qoi_list)
      std::cout << qoi->GetNumCalls() << std::endl;
    std::cout << "==== posterior evals" << std::endl;
    for (auto posterior : posterior_list)
      std::cout << posterior->GetNumCalls() << std::endl;
    std::cout << "==== indices" << std::endl;
    for (auto index : index_list)
      std::cout << *index << std::endl;
  }

  virtual std::shared_ptr<MCMCProposal> Proposal (std::shared_ptr<MultiIndex> const& index, std::shared_ptr<AbstractSamplingProblem> const& samplingProblem) override {
    pt::ptree pt;
    pt.put("BlockIndex",0);

    /*Eigen::VectorXd mu = Eigen::VectorXd::Zero(NUM_PARAM);
    Eigen::MatrixXd cov = Eigen::MatrixXd::Identity(NUM_PARAM,NUM_PARAM);

    auto prior = std::make_shared<Gaussian>(mu, cov);

    return std::make_shared<CrankNicolsonProposal>(pt, samplingProblem, prior);*/

    pt.put("AdaptStart", 100);
    pt.put("AdaptSteps", 100);
    pt.put("AdaptEnd", 1000);
    pt.put("InitialVariance", 10);
    //return std::make_shared<AMProposal>(pt, samplingProblem);

    pt.put("ProposalVariance",10);
    return std::make_shared<MHProposal>(pt, samplingProblem);

    /*auto mu = Eigen::VectorXd::Zero(NUM_PARAM);
    Eigen::MatrixXd cov = Eigen::MatrixXd::Identity(NUM_PARAM, NUM_PARAM);
    cov *= 0.5;

    auto prior = std::make_shared<Gaussian>(mu, cov);

    return std::make_shared<MHProposal>(pt, samplingProblem, prior);*/
  }

  /*void SetComm(std::shared_ptr<parcer::Communicator> const& comm) override {
    _comm = comm;
  }*/

  virtual std::shared_ptr<MultiIndex> FinestIndex() override {
    auto index = std::make_shared<MultiIndex>(1);
    index->SetValue(0, 1);
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
    //auto graph = createWorkGraph(host, bearer_token, index->GetValue(0));
    //auto posterior = graph->CreateModPiece("likelihood");
    //auto qoi = graph->CreateModPiece("qoi");

    // Set up model
    json config;
    config["level"] = index->GetValue(0);
    auto posterior = std::make_shared<HTTPModPiece>(host, config);

    posterior_list.push_back(posterior);
    //qoi_list.push_back(qoi);
    index_list.push_back(index);

    Eigen::VectorXd centroid = Eigen::VectorXd::Ones(2);
    centroid(0) = 30.0;
    centroid(1) = 30.0;

    centroid = Eigen::VectorXd::Zero(1);

    return std::make_shared<ExpensiveSamplingProblem>(posterior, centroid, _pt);
    //return std::make_shared<ExpensiveSamplingProblem>(posterior, qoi, centroid, _pt);
  }

  virtual std::shared_ptr<MIInterpolation> Interpolation (std::shared_ptr<MultiIndex> const& index) override {
    return std::make_shared<MyInterpolation>();
  }

  virtual Eigen::VectorXd StartingPoint (std::shared_ptr<MultiIndex> const& index) override {
    //Starting guess: zero
    Eigen::VectorXd start = Eigen::VectorXd::Ones(2);
    start(0) = 30.0;
    start(1) = 30.0;

    start = Eigen::VectorXd::Zero(1);

    return start;
  }

private:
  //std::shared_ptr<parcer::Communicator> _comm;
  //std::shared_ptr<parcer::Communicator> _globalComm;
  pt::ptree _pt;
  std::string host;
  std::vector<std::shared_ptr<HTTPModPiece>> posterior_list;
  std::vector<std::shared_ptr<ModGraphPiece>> qoi_list;
  std::vector<std::shared_ptr<MultiIndex>> index_list;
};

