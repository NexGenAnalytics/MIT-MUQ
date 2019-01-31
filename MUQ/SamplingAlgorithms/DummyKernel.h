#ifndef DUMMYKERNEL_H_
#define DUMMYKERNEL_H_

#include "MUQ/SamplingAlgorithms/TransitionKernel.h"

#include "MUQ/SamplingAlgorithms/MCMCProposal.h"

#include "MUQ/SamplingAlgorithms/SingleChainMCMC.h"

namespace muq {
  namespace SamplingAlgorithms {

    /** @brief MCMC kernel for Multiindex methods.
        @details This kernel combines a coarse proposal from a coarse chain
        with a fine one, as needed for MIMCMC.
     */
    class DummyKernel : public TransitionKernel {
    public:

      DummyKernel(boost::property_tree::ptree const& pt,
               std::shared_ptr<AbstractSamplingProblem> problem,
               std::shared_ptr<MCMCProposal> proposal);

      ~DummyKernel();

      virtual std::shared_ptr<MCMCProposal> Proposal(){return proposal;};

      virtual void PostStep(unsigned int const t, std::vector<std::shared_ptr<SamplingState>> const& state) override;

      virtual void PrintStatus(std::string prefix) const override;

      virtual std::vector<std::shared_ptr<SamplingState>> Step(unsigned int const t, std::shared_ptr<SamplingState> prevState) override;

    protected:
      std::shared_ptr<MCMCProposal> proposal;

      unsigned int numCalls = 0;
    };
  } // namespace SamplingAlgorithms
} // namespace muq

#endif
