#ifndef PARALLELABSTRACTSAMPLINGPROBLEM_H_
#define PARALLELABSTRACTSAMPLINGPROBLEM_H_

#if MUQ_HAS_MPI

#if !MUQ_HAS_PARCER
#error
#endif

#include <parcer/Communicator.h>
#include "MUQ/SamplingAlgorithms/AbstractSamplingProblem.h"
#include "MUQ/SamplingAlgorithms/ParallelFlags.h"

namespace muq {
  namespace SamplingAlgorithms {

    class ParallelAbstractSamplingProblem : public AbstractSamplingProblem
    {
    public:

      ParallelAbstractSamplingProblem(std::shared_ptr<parcer::Communicator> comm, int id, std::shared_ptr<AbstractSamplingProblem> abstractSamplingProblem)
       : AbstractSamplingProblem(abstractSamplingProblem->blockSizes, abstractSamplingProblem->blockSizesQOI)
      {
        this->comm = comm;
        this->id = id;
        this->abstractSamplingProblem = abstractSamplingProblem;
      }

      virtual ~ParallelAbstractSamplingProblem() = default;

      virtual double LogDensity(unsigned int const t, std::shared_ptr<SamplingState> state, AbstractSamplingProblem::SampleType type) override {

        if (comm->GetRank() == 0) {
          for (int dest = 1; dest < comm->GetSize(); dest++) {
            comm->Send(ControlFlag::LOGDENSITY, dest, WorkgroupTag);
            comm->Send(id, dest, WorkgroupTag);
            comm->Send(state->state[0], dest, WorkgroupTag);
          }
        }

        return abstractSamplingProblem->LogDensity(t, state, type);
      }

      virtual std::shared_ptr<SamplingState> QOI() override {
        return abstractSamplingProblem->QOI();
      }

      std::shared_ptr<AbstractSamplingProblem> GetSequentialProblem() {
        return abstractSamplingProblem;
      }

    private:
      int id;
      std::shared_ptr<parcer::Communicator> comm;
      std::shared_ptr<AbstractSamplingProblem> abstractSamplingProblem;
    };

  }
}

#endif

#endif
