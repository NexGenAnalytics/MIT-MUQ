#include "MUQ/SamplingAlgorithms/SLMCMC.h"

namespace muq {
  namespace SamplingAlgorithms {

    SLMCMC::SLMCMC (pt::ptree pt, std::shared_ptr<MIComponentFactory> componentFactory)
     : SamplingAlgorithm(nullptr),
       componentFactory(componentFactory)
    {
      std::cout << "constr" << std::endl;
      auto index = componentFactory->FinestIndex();
      std::cout << "constr2" << std::endl;

      auto problem = componentFactory->SamplingProblem(index);
      auto proposal = componentFactory->Proposal(index, problem);

      pt::ptree ptBlockID;
      ptBlockID.put("BlockIndex",0);
      std::vector<std::shared_ptr<TransitionKernel>> kernels(1);
      kernels[0] = std::make_shared<MHKernel>(ptBlockID,problem,proposal);

      Eigen::VectorXd startingPoint = componentFactory->StartingPoint(index);

      coarse_chain = std::make_shared<SingleChainMCMC>(pt,kernels,startingPoint);
    }

    std::shared_ptr<SampleCollection> SLMCMC::GetSamples() const {
      return nullptr;
    }
    std::shared_ptr<SampleCollection> SLMCMC::GetQOIs() const {
      return nullptr;
    }

    std::shared_ptr<SampleCollection> SLMCMC::RunImpl() {
      return coarse_chain->Run();
    }

    Eigen::VectorXd SLMCMC::MeanQOI() {
      return coarse_chain->GetQOIs()->Mean();
    }

    Eigen::VectorXd SLMCMC::MeanParameter() {
      auto samps = coarse_chain->GetSamples();
      return samps->Mean();
    }

  }
}
