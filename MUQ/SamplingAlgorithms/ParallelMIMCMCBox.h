#ifndef PARALLELMIMCMCBOX_H_
#define PARALLELMIMCMCBOX_H_

#if MUQ_HAS_MPI

#if !MUQ_HAS_PARCER
#error
#endif

#include <parcer/Communicator.h>


namespace muq {
  namespace SamplingAlgorithms {

    class ParallelMIMCMCBox {
    public:

        ParallelMIMCMCBox(std::shared_ptr<MIComponentFactory> componentFactory, std::shared_ptr<MultiIndex> boxHighestIndex, std::shared_ptr<parcer::Communicator> global_comm, std::shared_ptr<PhonebookClient> phonebookClient)
       : componentFactory(componentFactory),
         boxHighestIndex(boxHighestIndex)
      {
        pt::ptree ptChains;
        ptChains.put("NumSamples", 1e4); // number of MCMC steps expected, so we'll pass it in
        pt::ptree ptBlockID;
        ptBlockID.put("BlockIndex",0);


        boxLowestIndex = MultiIndex::Copy(boxHighestIndex);
        --(*boxLowestIndex);

        std::shared_ptr<MultiIndex> boxSize = std::make_shared<MultiIndex>(*boxHighestIndex - *boxLowestIndex);

        // Set up Multiindex box
        boxIndices = MultiIndexFactory::CreateFullTensor(boxSize->GetVector());
        boxChains.resize(boxIndices->Size());

        std::shared_ptr<SingleChainMCMC> coarse_chain = nullptr;
        std::shared_ptr<AbstractSamplingProblem> coarse_problem = nullptr;

        for (uint i = 0; i < boxIndices->Size(); i++) {
          std::shared_ptr<MultiIndex> boxIndex = (*boxIndices)[i];

          if (boxIndex->Max() == 0) { // We're the root node of the MI box
            if (boxLowestIndex->Max() == 0) { // and also the absolute root node, so run the coarsest chain by ourselves
              coarse_problem = componentFactory->SamplingProblem(boxLowestIndex);
              auto proposal_coarse = componentFactory->Proposal(boxLowestIndex, coarse_problem);

              std::vector<std::shared_ptr<TransitionKernel>> coarse_kernels(1);
              coarse_kernels[0] = std::make_shared<MHKernel>(ptBlockID,coarse_problem,proposal_coarse);

              Eigen::VectorXd startPtCoarse = componentFactory->StartingPoint(boxLowestIndex);

              coarse_chain = std::make_shared<SingleChainMCMC>(ptChains,coarse_kernels,startPtCoarse);
              boxChains[boxIndices->MultiToIndex(boxIndex)] = coarse_chain;

            } else { // or we have to request proposals from the next coarser chain
              std::shared_ptr<MultiIndex> remoteIndex = MultiIndex::Copy(boxLowestIndex);
              int maxCoeffId;
              int maxEntry = remoteIndex->GetVector().maxCoeff(&maxCoeffId);
              if (maxEntry <= 0) {
                std::cerr << "whoopsie" << std::endl;
                exit(45);
              }

              remoteIndex->SetValue(maxCoeffId, maxEntry - 1);

              pt::ptree pt;
              pt.put("BlockIndex",0);
              pt.put("NumSamples",1);

              coarse_problem = componentFactory->SamplingProblem(remoteIndex); // TODO: Try to avoid this

              auto coarse_proposal = std::make_shared<RemoteMIProposal>(pt, coarse_problem, global_comm, remoteIndex, phonebookClient);
              auto startingPoint = componentFactory->StartingPoint(boxLowestIndex);
              auto problem = componentFactory->SamplingProblem(boxLowestIndex);
              auto proposal = componentFactory->Proposal(boxLowestIndex, problem);
              auto proposalInterpolation = componentFactory->Interpolation(boxLowestIndex);

              std::vector<std::shared_ptr<TransitionKernel>> kernels(1);
              kernels[0] = std::make_shared<MIKernel>(ptBlockID,problem,coarse_problem,proposal,coarse_proposal,proposalInterpolation,nullptr);

              auto startingState = std::make_shared<SamplingState>(startingPoint);
              startingState->meta["coarseSample"] = std::make_shared<SamplingState>(componentFactory->StartingPoint(remoteIndex));

              coarse_chain = std::make_shared<SingleChainMCMC>(pt,kernels,startingState);

              boxChains[boxIndices->MultiToIndex(boxIndex)] = coarse_chain;
              coarse_problem = problem;
            }
          } else { // All in all, you're just another node in the MI box (We don't need no communication... We don't need no remote control...) - Pink Floyd

            std::shared_ptr<MultiIndex> index = std::make_shared<MultiIndex>(*boxLowestIndex + *boxIndex);

            auto problem = componentFactory->SamplingProblem(index);
            auto proposal = componentFactory->Proposal(index, problem);
            auto coarse_proposal = componentFactory->CoarseProposal(index, coarse_problem, coarse_chain);
            //auto coarse_proposal = componentFactory->CoarseProposal(index, coarse_problem, coarse_chain);
            auto proposalInterpolation = componentFactory->Interpolation(index);
            auto startingPoint = componentFactory->StartingPoint(index);

            std::vector<std::shared_ptr<TransitionKernel>> kernels(1);
            kernels[0] = std::make_shared<MIKernel>(ptBlockID,problem,coarse_problem,proposal,coarse_proposal,proposalInterpolation,coarse_chain);
            //kernels[0] = std::make_shared<MIKernel>(ptBlockID,problem,coarse_problem,proposal,coarse_proposal,proposalInterpolation,coarse_chain);

            //auto startingState = std::make_shared<SamplingState>(startingPoint);
            //startingState->meta["coarseSample"] = std::make_shared<SamplingState>(componentFactory->StartingPoint(index));

            auto chain = std::make_shared<SingleChainMCMC>(ptChains,kernels,startingPoint);

            boxChains[boxIndices->MultiToIndex(boxIndex)] = chain;

            if (boxIndex->Max() == 1)
              finestProblem = problem;
          }

        }
      }

      std::shared_ptr<AbstractSamplingProblem> GetFinestProblem() {
        return finestProblem;
      }

      std::shared_ptr<SingleChainMCMC> FinestChain() {
        std::shared_ptr<MultiIndex> boxSize = std::make_shared<MultiIndex>(*boxHighestIndex - *boxLowestIndex);
        return boxChains[boxIndices->MultiToIndex(boxSize)];
      }

      std::shared_ptr<SingleChainMCMC> GetChain(int index) {
        //std::shared_ptr<MultiIndex> boxSize = std::make_shared<MultiIndex>(*boxHighestIndex - *boxLowestIndex);
        return boxChains[index];
      }

      int NumChains() {
        //std::shared_ptr<MultiIndex> boxSize = std::make_shared<MultiIndex>(*boxHighestIndex - *boxLowestIndex);
        return boxIndices->Size();
      }

      void Sample() {
        for (uint i = 0; i < boxIndices->Size(); i++) {
          std::shared_ptr<MultiIndex> boxIndex = (*boxIndices)[i];
          auto chain = boxChains[boxIndices->MultiToIndex(boxIndex)];
          chain->Sample();
        }
      }

      Eigen::VectorXd MeanQOI() {
        Eigen::VectorXd sampMean = Eigen::VectorXd::Zero(GetFinestProblem()->blockSizesQOI.sum());

        for (uint i = 0; i < boxIndices->Size(); i++) {
          std::shared_ptr<MultiIndex> boxIndex = (*boxIndices)[i];
          auto chain = boxChains[boxIndices->MultiToIndex(boxIndex)];
          auto samps = chain->GetQOIs();

          std::shared_ptr<MultiIndex> index = std::make_shared<MultiIndex>(*boxLowestIndex + *boxIndex);
          auto indexDiffFromTop = std::make_shared<MultiIndex>(*boxHighestIndex - *index);

          if (indexDiffFromTop->Sum() % 2 == 0) {
            sampMean += samps->Mean();
          } else {
            sampMean -= samps->Mean();
          }
        }
        return sampMean;
      }

    private:


      std::shared_ptr<MIComponentFactory> componentFactory;
      std::shared_ptr<MultiIndex> boxHighestIndex;
      std::shared_ptr<MultiIndex> boxLowestIndex;
      std::shared_ptr<MultiIndexSet> boxIndices;
      std::vector<std::shared_ptr<SingleChainMCMC>> boxChains;
      std::vector<std::shared_ptr<SingleChainMCMC>> tailChains;
      std::shared_ptr<AbstractSamplingProblem> finestProblem;
    };

  }
}

#endif

#endif
