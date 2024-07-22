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
#include "MUQ/SamplingAlgorithms/MultiIndexEstimator.h"

#include "MUQ/SamplingAlgorithms/MIComponentFactory.h"

#include <boost/property_tree/ptree.hpp>

#include <gtest/gtest.h>

namespace pt = boost::property_tree;
using namespace muq::Modeling;
using namespace muq::SamplingAlgorithms;
using namespace muq::Utilities;


class MySamplingProblem : public AbstractSamplingProblem {
public:
  MySamplingProblem(std::shared_ptr<muq::Modeling::ModPiece> targetIn)
   : AbstractSamplingProblem(Eigen::VectorXi::Constant(1,2), Eigen::VectorXi::Constant(1,2)),
     target(targetIn){
      std::cout << "      @ Constructing MySamplingProblem with ModPiece " << targetIn->Name() << std::endl;
      }

  virtual ~MySamplingProblem() = default;


  virtual double LogDensity(std::shared_ptr<SamplingState> const& state) override {
    // trigger PDE run here, and perform log density evaluation
    // but here just simply copy the data
    std::cout << "MySamplingAlgorithm: LogDensity: state: " << std::endl;
    for (auto const& s: state->state) {
      std::cout << "-----" << std::endl << s << std::endl;
    }
    lastState = state;
    return target->Evaluate(state->state).at(0)(0);
  };

  virtual std::shared_ptr<SamplingState> QOI() override {
    std::cout << "MySamplingAlgorithm: QOI() " << std::endl;
    assert (lastState != nullptr);
    return std::make_shared<SamplingState>(lastState->state, 1.0);
  }

private:
  std::shared_ptr<SamplingState> lastState = nullptr;

  std::shared_ptr<muq::Modeling::ModPiece> target;

};


class MyMLInterpolation : public MIInterpolation {
public:
  std::shared_ptr<SamplingState> Interpolate (std::shared_ptr<SamplingState> const& coarseProposal, std::shared_ptr<SamplingState> const& fineProposal) {
    // in this case the same params dim is the same for all levels
    // so the state proposed at fine level is the same as the one at coarse
    return std::make_shared<SamplingState>(coarseProposal->state);
  }
};

class MyMLComponentFactory : public MIComponentFactory {
public:
  MyMLComponentFactory(pt::ptree pt)
   : pt(pt) {
    std::cout << "Constructing empty MyMLComponentFactory" << std::endl;
   }

  virtual std::shared_ptr<MCMCProposal> Proposal (std::shared_ptr<MultiIndex> const& index, std::shared_ptr<AbstractSamplingProblem> const& samplingProblem) override {
    // proposal for any given level
    // in this case we have the same problem for every level because the same uncertain param dim
    // if index is 0, then i need to provide proposal for that level
    std::cout << "    - MyMLComponentFactory: Proposal(): index: -" << index->ToString() << std::endl;
    pt::ptree pt;
    pt.put("BlockIndex",0);

    Eigen::VectorXd mu(2);
    mu << 1.0, 2.0;
    Eigen::MatrixXd cov(2,2);
    cov << 0.7, 0.6,
    0.6, 1.0;
    cov *= 20.0;

    auto prior = std::make_shared<Gaussian>(mu, cov);
    std::cout << "    -----------------------------------------------" << std::endl;
    return std::make_shared<CrankNicolsonProposal>(pt, samplingProblem, prior);
  }

  virtual std::shared_ptr<MultiIndex> FinestIndex() override {
    // @UNDERSTOOD
    std::cout << "    - MyMLComponentFactory: FinestIndex() -" << std::endl;
    auto index = std::make_shared<MultiIndex>(1); // 1 = dimension of the hierarchy, for ML this is 1 alway
    index->SetValue(0, 3);// 0,1,2,3 defines the levels, where 3 is the finest
    std::cout << "    ---------------------------------------" << std::endl;
    return index;
  }

  virtual std::shared_ptr<MCMCProposal> CoarseProposal (std::shared_ptr<MultiIndex> const& fineIndex,
                                                        std::shared_ptr<MultiIndex> const& coarseIndex,
                                                        std::shared_ptr<AbstractSamplingProblem> const& coarseProblem,
                                                        std::shared_ptr<SingleChainMCMC> const& coarseChain) override {
    // @DONOTTOUCH Linus says leave it as it is.
    // this basically drives the coarse chain by a bit and then take the final result
    std::cout << "    - MyMLComponentFactory: CoarseProposal(): -" << std::endl
              << "fineIndex: " << fineIndex->ToString() << std::endl
              << "coarseIndex: " << coarseIndex->ToString() << std::endl
              << "coarseChain status: " << std::endl;
    coarseChain->PrintStatus("");
    pt::ptree ptProposal = pt;
    ptProposal.put("BlockIndex",0);
    std::cout << "    --------------------------------------------" << std::endl;
    return std::make_shared<SubsamplingMIProposal>(ptProposal, coarseProblem, coarseIndex, coarseChain);
  }

  virtual std::shared_ptr<AbstractSamplingProblem> SamplingProblem (std::shared_ptr<MultiIndex> const& index) override {
    // return the "log posterior" to use at a given level
    std::cout << "    - MyMLComponentFactory: SamplingProblem(): -" << std::endl;
    Eigen::VectorXd mu(2);
    mu << 1.0, 2.0;
    Eigen::MatrixXd cov(2,2);
    cov << 0.7, 0.6,
           0.6, 1.0;

    // shifting mean and broaden to emulate
    if (index->GetValue(0) == 0) {
      std::cout << "       index->GetValue(0) == 0" << std::endl;
      mu *= 0.8;
      cov *= 2.0;
    } else if (index->GetValue(0) == 1) {
      std::cout << "       index->GetValue(0) == 1" << std::endl;
      mu *= 0.9;
      cov *= 1.5;
    } else if (index->GetValue(0) == 2) {
      std::cout << "       index->GetValue(0) == 2" << std::endl;
      mu *= 0.99;
      cov *= 1.1;
    } else if (index->GetValue(0) == 3) {
      std::cout << "       index->GetValue(0) == 3" << std::endl;
      mu *= 1.0;
      cov *= 1.0;
    } else {
      std::cerr << "Sampling problem not defined!" << std::endl;
      assert (false);
    }

    auto coarseTargetDensity = std::make_shared<Gaussian>(mu, cov)->AsDensity();
    auto res = std::make_shared<MySamplingProblem>(coarseTargetDensity);
    std::cout << "    --------------------------------------------" << std::endl;
    return res;
  }

  virtual std::shared_ptr<MIInterpolation> Interpolation (std::shared_ptr<MultiIndex> const& index) override {
    // only need if the parameter dim changes across levels
    // for example, level 0 has 10 KL coeffs, and level 1 has 56 KL coeffs, etc
    // always betweeen two subsquenet levels
    std::cout << "MyMLComponentFactory: Interpolation(): index: " << index->ToString() << std::endl;
    return std::make_shared<MyMLInterpolation>();
  }

  virtual Eigen::VectorXd StartingPoint (std::shared_ptr<MultiIndex> const& index) override {
    std::cout << "MyMLComponentFactory: StartingPoint(): index: " << index->ToString() << std::endl;
    Eigen::VectorXd mu(2);
    mu << 1.0, 2.0;
    return mu;
  }
  pt::ptree pt;
};

TEST(MLMCMCTest, GreedyMLMCMC)
{

  pt::ptree pt;

  pt.put("NumSamples", 1e4); // number of samples for single level
  pt.put("NumInitialSamples", 1e3); // number of initial samples for greedy MLMCMC
  pt.put("GreedyTargetVariance", 0.05); // estimator variance to be achieved by greedy algorithm
  // take 5 steps on the coarse before using for level 1
  pt.put("MLMCMC.Subsampling_0", 5); // estimator variance to be achieved by greedy algorithm
  // take 5 steps on the coarse before using for level 2
  pt.put("MLMCMC.Subsampling_1", 3); // estimator variance to be achieved by greedy algorithm
  // take 5 steps on the coarse before using for level 3
  pt.put("MLMCMC.Subsampling_2", 1); // estimator variance to be achieved by greedy algorithm
  pt.put("MLMCMC.Subsampling_3", 0); // estimator variance to be achieved by greedy algorithm

  auto componentFactory = std::make_shared<MyMLComponentFactory>(pt);

  std::cout << "Creating greedymlmcmc" <<  std::endl;
  // greedy chooses automatically how many samples to draw
  GreedyMLMCMC greedymlmcmc (pt, componentFactory);
  // this calls the GreedyMLMCMC constructor, that uses the component Factory to get the FinestIndex; it also constructs our samplingProblem with the finest index and gets its numBlocksQOI and tests if this is > 0 to set the GreedyMLMCMCs useQOI boolean.
  // after the initialization list, it loops through all levels 0->levels (included) and creates a multiIndex of 1 value that is the level (boxHighestIndex); with that and the componentFactory it creates a MIMCMCBox.
  // MIMCMCBox constructor:
  // 1) Set up root index sampling
  // amont other things, of interest, usage of the componentFactory:
  //  - constructs 1 samplingProblem
  // - goes through 2 loops based on different indices and uses the componentFactory in different ways


  std::cout << "Running greedymlmcmc" <<  std::endl;
  greedymlmcmc.Run();
  std::cout << "Setting Draw to false for greedymlmcmc" <<  std::endl;
  greedymlmcmc.Draw(false);

  // questions:
  // - understand each method
  // - what are the inputs to the algo

  Eigen::VectorXd trueMu(2);
  trueMu << 1.0, 2.0;
  Eigen::MatrixXd trueCov(2,2);
  trueCov << 0.7, 0.6,
             0.6, 1.0;

  auto params = greedymlmcmc.GetSamples();
  Eigen::VectorXd mean = params->Mean();
  Eigen::VectorXd mcse = params->StandardError();
  EXPECT_NEAR(trueMu(0), mean(0), 3.*mcse(0));
  EXPECT_NEAR(trueMu(1), mean(1), 3.0*mcse(1));

  Eigen::VectorXd variance = params->Variance();
  EXPECT_NEAR(trueCov(0,0), variance(0), 5.0*mcse(0));
  EXPECT_NEAR(trueCov(1,1), variance(1), 5.0*mcse(1));

  Eigen::VectorXd skewness = params->Skewness();
  EXPECT_NEAR(0.0, skewness(0), 0.5);
  EXPECT_NEAR(0.0, skewness(0), 0.5);

  Eigen::MatrixXd covariance = params->Covariance();
  EXPECT_NEAR(trueCov(0,0), covariance(0,0), 0.2);
  EXPECT_NEAR(trueCov(1,1), covariance(1,1), 0.2);
  EXPECT_NEAR(trueCov(0,1), covariance(0,1), 0.2);
  EXPECT_NEAR(trueCov(1,0), covariance(1,0), 0.2);


  auto qois = greedymlmcmc.GetQOIs();
  mean = qois->Mean();
  EXPECT_NEAR(trueMu(0), mean(0), 0.3);
  EXPECT_NEAR(trueMu(1), mean(1), 0.3);

  variance = qois->Variance();
  EXPECT_NEAR(trueCov(0,0), variance(0), 0.3);
  EXPECT_NEAR(trueCov(1,1), variance(1), 0.3);

  skewness = qois->Skewness();
  EXPECT_NEAR(0.0, skewness(0), 0.5);
  EXPECT_NEAR(0.0, skewness(0), 0.5);
}
