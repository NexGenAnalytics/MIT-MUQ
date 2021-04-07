#include <iostream>
#include <string>

#include <boost/property_tree/ptree.hpp>

#include <MUQ/Modeling/OneStepCachePiece.h>
#include <MUQ/Modeling/ModPiece.h>
#include <MUQ/Modeling/ModGraphPiece.h>

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
  MyMIComponentFactory(pt::ptree pt, std::shared_ptr<parcer::Communicator> globalComm, std::string host, std::string bearer_token)
   : _pt(pt), _globalComm(globalComm), host(host), bearer_token(bearer_token) {}

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
    auto graph = createWorkGraph(host, bearer_token, index->GetValue(0));
    auto posterior = graph->CreateModPiece("likelihood");
    auto qoi = graph->CreateModPiece("qoi");

    return std::make_shared<muq::SamplingAlgorithms::SamplingProblem>(posterior, qoi);
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
  std::string host, bearer_token;
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

  MPI_Init(&argc, &argv);
  auto comm = std::make_shared<parcer::Communicator>(MPI_COMM_WORLD);

  std::time_t result = std::time(nullptr);
  std::string timestamp = std::asctime(std::localtime(&result));
  comm->Bcast(timestamp, 0);
  auto tracer = std::make_shared<OTF2Tracer>("trace",timestamp);


{ // Inverse UQ


  pt::ptree pt;

  pt.put("verbosity", 1); // show some output
  pt.put("MCMC.BurnIn", 50);
  pt.put("NumSamples_0", 1e3);
  pt.put("NumSamples_1", 1e2);
  pt.put("MLMCMC.Scheduling", true);
  pt.put("MLMCMC.Subsampling_0", 49);
  pt.put("MLMCMC.Subsampling_1", 0);

  auto componentFactory = std::make_shared<MyMIComponentFactory>(pt, comm, host, bearer_token);


  StaticLoadBalancingMIMCMC parallelMIMCMC (pt, componentFactory, std::make_shared<RoundRobinStaticLoadBalancer>(), comm, tracer);
  if (comm->GetRank() == 0) {
    parallelMIMCMC.Run();
    Eigen::VectorXd meanQOI = parallelMIMCMC.MeanQOI();
    std::cout << "mean QOI: " << meanQOI.transpose() << std::endl;
  }

  if (comm->GetRank() == 0)
    remove("FullParallelMLMCMC.hdf5");

  parallelMIMCMC.WriteToFile("FullParallelMLMCMC.hdf5");
  parallelMIMCMC.Finalize();
}

  tracer->write();


  MPI_Finalize();
  return 0;
}