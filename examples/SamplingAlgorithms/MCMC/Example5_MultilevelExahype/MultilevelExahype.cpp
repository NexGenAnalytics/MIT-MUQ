#include "MUQ/SamplingAlgorithms/SLMCMC.h"
#include "MUQ/SamplingAlgorithms/GreedyMLMCMC.h"
#include "MUQ/SamplingAlgorithms/MIMCMC.h"

#include "MUQ/Modeling/Distributions/Gaussian.h"
#include "MUQ/Modeling/Distributions/Density.h"

#include "MUQ/SamplingAlgorithms/MHKernel.h"
#include "MUQ/SamplingAlgorithms/MHProposal.h"
#include "MUQ/SamplingAlgorithms/CrankNicolsonProposal.h"
#include "MUQ/SamplingAlgorithms/SamplingProblem.h"
#include "MUQ/SamplingAlgorithms/SubsamplingMIProposal.h"

#include "MUQ/SamplingAlgorithms/ParallelizableMIComponentFactory.h"
#include "MUQ/SamplingAlgorithms/ParallelMIComponentFactory.h"
#include "MUQ/SamplingAlgorithms/ParallelFixedSamplesMIMCMC.h"

#include "MUQ/SamplingAlgorithms/DummyKernel.h"

#include <boost/property_tree/ptree.hpp>

#include <initandsoon.h>

#include "calculateLikelihood.hh"

namespace pt = boost::property_tree;
using namespace muq::Modeling;
using namespace muq::SamplingAlgorithms;
using namespace muq::Utilities;

#include <fstream>

//TODO: read in
const int NUM_PARAM = 2;
int count;
const static Eigen::IOFormat CSVFormat(Eigen::FullPrecision, 0, ", ", ";\n", "", "", "", ";\n");

class MySamplingProblem : public AbstractSamplingProblem {
public:
  MySamplingProblem(std::shared_ptr<parcer::Communicator> comm,std::shared_ptr<MultiIndex> index)
  : AbstractSamplingProblem(Eigen::VectorXi::Constant(1,NUM_PARAM), Eigen::VectorXi::Constant(1,NUM_PARAM))
{
    this->comm = comm;
    this->index = index;
    muq::setCommunicator(comm->GetMPICommunicator());

}

  virtual ~MySamplingProblem(){
  }

  virtual double LogDensity(unsigned int const t, std::shared_ptr<SamplingState> state, AbstractSamplingProblem::SampleType type) override {
    comm->Barrier();
    lastState = state;

    //TODO pass directly to exahype
    //std::ofstream file("Input/parameters.csv");
    //file << state->state[0].format(CSVFormat);
    //file.close();
    unsigned int idx= index->GetValue(0); //TODO check if passing index works

    //Write parameters into ExaHyPE readable format
    std::vector<double> param(state->state[0].size());
    for(int i = 0; i < param.size(); i++){
        param[i] = state->state[0][i];
    }

    //Create some debug output
    std::cout << "parameter:" << state->state[0].transpose() << std::endl;

    std::ofstream ost;
    ost.open("parameters.log", std::ios::app);
    ost << param[0] << ", " << param[1] << std::endl;
    ost.close();

    std::cout << "Sample number: " << count++ << std::endl;
    std::cout << "Parameter " << param[0] << " " << param[1] << std::endl;

    //Discard stupid parameters
    if (param[0] >1.0 || param[0] < -0.1 || param[1]>1.0 || param[1]<-0.1)//reject parameters outside domain
        return -24;

    //run forward model
    int level = 4;
    auto output = muq::run_exahype(param,level);

    comm->Barrier();
    double sigma = 1.0;
    return calculateLikelihood(output);// - 0.5/(sigma*sigma)*state->state[0].squaredNorm();
  };

  virtual std::shared_ptr<SamplingState> QOI() override {
    assert (lastState != nullptr);
    return std::make_shared<SamplingState>(lastState->state, 1.0);
  }

private:
  std::shared_ptr<SamplingState> lastState = nullptr;
  std::shared_ptr<parcer::Communicator> comm;
  std::shared_ptr<MultiIndex> index;
};


class MyInterpolation : public MIInterpolation {
public:
  std::shared_ptr<SamplingState> Interpolate (std::shared_ptr<SamplingState> coarseProposal, std::shared_ptr<SamplingState> fineProposal) {
    return std::make_shared<SamplingState>(coarseProposal->state);
  }
};

class MyMIComponentFactory : public ParallelizableMIComponentFactory {
public:
  virtual std::shared_ptr<MCMCProposal> Proposal (std::shared_ptr<MultiIndex> index, std::shared_ptr<AbstractSamplingProblem> samplingProblem) override {
    pt::ptree pt;
    pt.put("BlockIndex",0);

    /*Eigen::VectorXd mu = Eigen::VectorXd::Zero(NUM_PARAM);
    Eigen::MatrixXd cov = Eigen::MatrixXd::Identity(NUM_PARAM,NUM_PARAM);

    auto prior = std::make_shared<Gaussian>(mu, cov);

    return std::make_shared<CrankNicolsonProposal>(pt, samplingProblem, prior);*/

    auto mu = Eigen::VectorXd::Zero(NUM_PARAM);
    Eigen::MatrixXd cov = Eigen::MatrixXd::Identity(NUM_PARAM, NUM_PARAM);
    cov *= 0.005;

    auto prior = std::make_shared<Gaussian>(mu, cov);

    return std::make_shared<MHProposal>(pt, samplingProblem, prior);
  }

  void SetComm(std::shared_ptr<parcer::Communicator> comm) override {
    _comm = comm;
  }

  virtual std::shared_ptr<MultiIndex> FinestIndex() override {
    auto index = std::make_shared<MultiIndex>(1);
    index->SetValue(0, 0);
    return index;
  }

  virtual std::shared_ptr<MCMCProposal> CoarseProposal (std::shared_ptr<MultiIndex> index,
                                                        std::shared_ptr<AbstractSamplingProblem> coarseProblem,
                                                           std::shared_ptr<SingleChainMCMC> coarseChain) override {
    pt::ptree ptProposal;
    ptProposal.put("BlockIndex",0);
    int subsampling = 5;
    ptProposal.put("subsampling", subsampling);
    return std::make_shared<SubsamplingMIProposal>(ptProposal, coarseProblem, coarseChain);
  }

  virtual std::shared_ptr<AbstractSamplingProblem> SamplingProblem (std::shared_ptr<MultiIndex> index) override {

    return std::make_shared<MySamplingProblem>(_comm,index);
  }

  virtual std::shared_ptr<MIInterpolation> Interpolation (std::shared_ptr<MultiIndex> index) override {
    return std::make_shared<MyInterpolation>();
  }

  virtual Eigen::VectorXd StartingPoint (std::shared_ptr<MultiIndex> index) override {
    //Starting guess: zero
      Eigen::VectorXd start = Eigen::VectorXd::Ones(NUM_PARAM);
    start(0) = .6;
    start(1) = .6;
    return start;
  }

private:
  std::shared_ptr<parcer::Communicator> _comm;

};

class MCSampleProposal : public MCMCProposal {
public:
  MCSampleProposal(boost::property_tree::ptree       const& pt,
                   std::shared_ptr<AbstractSamplingProblem> prob,
                   std::shared_ptr<Distribution> dist
                  )
   : MCMCProposal(pt, prob),
     dist(dist)
  {}

  std::shared_ptr<SamplingState> Sample(std::shared_ptr<SamplingState> currentState) override {
    return std::make_shared<SamplingState>(dist->Sample());
  }

  double LogDensity(std::shared_ptr<SamplingState> currState,
                    std::shared_ptr<SamplingState> propState) override {
    return 0.0;
  }


private:
  std::shared_ptr<Distribution> dist;
};


int main(int argc, char** argv){
  int initThreadProvidedThreadLevelSupport;
  //bool result = MPI_Init_thread( &argc, &argv, MPI_THREAD_MULTIPLE, &initThreadProvidedThreadLevelSupport );

  muq::init(argc,argv);
  count = 0;
/*{ // Forward UQ
  pt::ptree pt;
  pt.put("BlockIndex",0);
  pt.put("NumSamples",1);

  auto componentFactory = std::make_shared<MyMIComponentFactory>();
  auto finestIndex = componentFactory->FinestIndex();
  auto problem = componentFactory->SamplingProblem(finestIndex);

  Eigen::VectorXd mu = Eigen::VectorXd::Zero(NUM_PARAM);
  Eigen::MatrixXd cov = Eigen::MatrixXd::Identity(NUM_PARAM,NUM_PARAM);
  auto dist = std::make_shared<Gaussian>(mu, cov);

  auto proposal = std::make_shared<MCSampleProposal>(pt, problem, dist);

  std::vector<std::shared_ptr<TransitionKernel>> kernels(1);
  kernels[0] = std::make_shared<DummyKernel>(pt,problem,proposal); // Accept all proposals into chain

  auto startingPoint = componentFactory->StartingPoint(finestIndex);
  auto chain = std::make_shared<SingleChainMCMC>(pt,kernels,startingPoint);

  chain->Run();
}*/

{ // Inverse UQ
  auto localFactory = std::make_shared<MyMIComponentFactory>();

  pt::ptree pt;

  pt.put("NumSamples", 100); // number of samples for single level
  pt.put("NumInitialSamples", 3); //ignore// number of initial samples for greedy MLMCMC
  pt.put("GreedyTargetVariance", 0.05); //ignore// estimator variance to be achieved by greedy algorithm
  pt.put("verbosity", 1); // show some output
  pt.put("BurnIn", 10);

  /*std::cout << std::endl << "*************** greedy multillevel chain" << std::endl << std::endl;

  GreedyMLMCMC greedymlmcmc (pt, componentFactory);
  greedymlmcmc.Run();
  std::cout << "mean QOI: " << greedymlmcmc.MeanQOI().transpose() << std::endl;
  greedymlmcmc.Draw(false);*/

  auto comm = std::make_shared<parcer::Communicator>();
  localFactory->SetComm(comm);
  auto componentFactory = std::make_shared<ParallelMIComponentFactory>(comm, comm, localFactory);

  if (comm->GetRank() == 0) {
  std::cout << std::endl << "*************** single chain reference" << std::endl << std::endl;

  SLMCMC slmcmc (pt, componentFactory);
  slmcmc.Run();


  std::cout << "mean Param: " << slmcmc.MeanParameter().transpose() << std::endl;
  std::cout << "mean QOI: " << slmcmc.MeanQOI().transpose() << std::endl;

  //std::cout << "variance QOI: " << slmcmc.VarianceQOI().transpose() << std::endl;

  //Write mean to file
  /*std::ofstream file("Input/parameters.csv");
  file << (slmcmc.MeanQOI()).format(CSVFormat);
  file.close();*/
  }
  componentFactory->finalize();

  //plot mean
  //std::vector<double> param = {slmcmc.MeanQOI()(0), slmcmc.MeanQOI()(1)};
  //muq::run_exahype(param);
  //calculateLikelihood();

}

  muq::finalize();
  MPI_Finalize();
  return 0;
}
