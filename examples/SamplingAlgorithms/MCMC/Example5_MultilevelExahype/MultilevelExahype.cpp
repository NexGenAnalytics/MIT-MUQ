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

#include "MUQ/SamplingAlgorithms/MIComponentFactory.h"
#include "MUQ/SamplingAlgorithms/DummyKernel.h"

#include <boost/property_tree/ptree.hpp>

namespace pt = boost::property_tree;
using namespace muq::Modeling;
using namespace muq::SamplingAlgorithms;
using namespace muq::Utilities;

#include <fstream>

const int NUM_PARAM = 25;


const static Eigen::IOFormat CSVFormat(Eigen::FullPrecision, 0, ", ", ";\n", "", "", "", ";\n");

class MySamplingProblem : public AbstractSamplingProblem {
public:
  MySamplingProblem()
  : AbstractSamplingProblem(Eigen::VectorXi::Constant(1,NUM_PARAM), Eigen::VectorXi::Constant(1,NUM_PARAM))
     {}

  virtual ~MySamplingProblem() = default;


  virtual double LogDensity(unsigned int const t, std::shared_ptr<SamplingState> state, AbstractSamplingProblem::SampleType type) override {
    lastState = state;

    std::ofstream file("Input/parameters.csv");
    file << state->state[0].format(CSVFormat);
    file.close();

    // run it
    system("./ExaHyPE-SWE ../SWE_MC.exahype");

    std::cout << "ref:" << state->state[0] << std::endl;



    /*const int dim_measurements = 2;
    Eigen::VectorXd measurements(dim_measurements);
    std::ifstream infile("heightmeasurement.csv");
    std::string line;
    int row = 0;
    while (std::getline(infile, line))
    {
      if (line.back() == ';') {
        measurements[row] = std::stod(line);
        std::cout << "number: " << measurements[row] << std::endl;
        row++;
      }
    }
    assert(row == dim_measurements);


    return target->Evaluate(state->state).at(0)(0);*/
    return 0.5;
  };

  virtual std::shared_ptr<SamplingState> QOI() override {
    assert (lastState != nullptr);
    return std::make_shared<SamplingState>(lastState->state, 1.0);
  }

private:
  std::shared_ptr<SamplingState> lastState = nullptr;

};


class MyInterpolation : public MIInterpolation {
public:
  std::shared_ptr<SamplingState> Interpolate (std::shared_ptr<SamplingState> coarseProposal, std::shared_ptr<SamplingState> fineProposal) {
    return std::make_shared<SamplingState>(coarseProposal->state);
  }
};

class MyMIComponentFactory : public MIComponentFactory {
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
    cov *= 0.01;

    auto prior = std::make_shared<Gaussian>(mu, cov);

    return std::make_shared<MHProposal>(pt, samplingProblem, prior);
  }

  virtual std::shared_ptr<MultiIndex> FinestIndex() override {
    auto index = std::make_shared<MultiIndex>(1);
    index->SetValue(0, 3);
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

    return std::make_shared<MySamplingProblem>();
  }

  virtual std::shared_ptr<MIInterpolation> Interpolation (std::shared_ptr<MultiIndex> index) override {
    return std::make_shared<MyInterpolation>();
  }

  virtual Eigen::VectorXd StartingPoint (std::shared_ptr<MultiIndex> index) override {
    /*Eigen::VectorXd mu(2);
    mu << 1.0, 2.0;*/
    return Eigen::VectorXd::Zero(NUM_PARAM);
  }

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


int main(){

{ // Forward UQ
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
}

{ // Inverse UQ
  auto componentFactory = std::make_shared<MyMIComponentFactory>();

  pt::ptree pt;

  pt.put("NumSamples", 1e1); // number of samples for single level
  pt.put("NumInitialSamples", 1e3); // number of initial samples for greedy MLMCMC
  pt.put("GreedyTargetVariance", 0.05); // estimator variance to be achieved by greedy algorithm
  pt.put("verbosity", 1); // show some output

  /*std::cout << std::endl << "*************** greedy multillevel chain" << std::endl << std::endl;

  GreedyMLMCMC greedymlmcmc (pt, componentFactory);
  greedymlmcmc.Run();
  std::cout << "mean QOI: " << greedymlmcmc.MeanQOI().transpose() << std::endl;
  greedymlmcmc.Draw(false);*/

  std::cout << std::endl << "*************** single chain reference" << std::endl << std::endl;

  SLMCMC slmcmc (pt, componentFactory);
  slmcmc.Run();
  std::cout << "mean QOI: " << slmcmc.MeanQOI().transpose() << std::endl;
}
  return 0;
}
