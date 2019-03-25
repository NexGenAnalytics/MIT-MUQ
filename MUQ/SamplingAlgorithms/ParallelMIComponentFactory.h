#ifndef PARALLELMICOMPONENTFACTORY_H_
#define PARALLELMICOMPONENTFACTORY_H_

#if MUQ_HAS_MPI

#if !MUQ_HAS_PARCER
#error
#endif


#include "MUQ/SamplingAlgorithms/ParallelAbstractSamplingProblem.h"

namespace muq {
  namespace SamplingAlgorithms {

    class ParallelMIComponentFactory : public MIComponentFactory {

    public:

      ParallelMIComponentFactory (std::shared_ptr<parcer::Communicator> comm, std::shared_ptr<parcer::Communicator> global_comm, std::shared_ptr<MIComponentFactory> componentFactory)
       : comm(comm), global_comm(global_comm), componentFactory(componentFactory)
      {

        if (comm->GetRank() != 0) {
          //Dune::Timer timer_idle;
          //Dune::Timer timer_full;
          while (true) {
            std::cout << "Rank " << comm->GetRank() << " waiting..." << std::endl;
            //timer_idle.start();
            ControlFlag command = comm->Recv<ControlFlag>(0, WorkgroupTag);
            //timer_idle.stop();
            std::cout << "Rank " << comm->GetRank() << " received " << command << std::endl;
            if (command == ControlFlag::FINALIZE) {
              comm->Barrier();
              std::cout << "Rank " << comm->GetRank() << " went through barrier" << std::endl;
              break;
            }
            if (command == ControlFlag::INIT_PROBLEM) {
              auto index = std::make_shared<MultiIndex>(1, comm->Recv<int>(0, WorkgroupTag));
              int id = comm->Recv<int>(0, WorkgroupTag);
              std::cout << "Rank " << comm->GetRank() << " building model index " << *index << std::endl;
              samplingProblems[id] = componentFactory->SamplingProblem(index);//std::make_shared<MySamplingProblem>(index, comm, id, measurements);
            }
            else if (command == ControlFlag::LOGDENSITY) {
              int id = comm->Recv<int>(0, WorkgroupTag);
              auto state = std::make_shared<SamplingState>(comm->Recv<Eigen::VectorXd>(0, WorkgroupTag));
              samplingProblems[id]->LogDensity(0.0, state, AbstractSamplingProblem::SampleType::Proposed); // TODO: Correct arguments!
            }
            else if (command == ControlFlag::TEST) {
              int id = comm->Recv<int>(0, WorkgroupTag);
              auto state = std::make_shared<SamplingState>(comm->Recv<Eigen::VectorXd>(0, WorkgroupTag));

              double density;
              //std::dynamic_pointer_cast<MySamplingProblem>(samplingProblems[id])->test(state, "test", density);
              std::cerr << "Not implemented!!!" << std::endl;
            }
            else if (command == ControlFlag::QOI) {
              int id = comm->Recv<int>(0, WorkgroupTag);

              samplingProblems[id]->QOI();
            } else {
              std::cerr << "Unexpected command!" << std::endl;
              exit(43);
            }
          }
          //std::cout << "Worker " << global_comm->GetRank() << " idle time:\t" << timer_idle.elapsed() << " of:\t" << timer_full.elapsed() << std::endl;
        }
      }

      void finalize() {
        if (comm->GetRank() != 0)
          return;
        std::cout << "Send finalize" << std::endl;
        for (int dest = 1; dest < comm->GetSize(); dest++)
          comm->Send(ControlFlag::FINALIZE, dest, WorkgroupTag);
        comm->Barrier();
        std::cout << "Finalized" << std::endl;
      }

      virtual std::shared_ptr<MCMCProposal> Proposal (std::shared_ptr<MultiIndex> index, std::shared_ptr<AbstractSamplingProblem> samplingProblem) override {
        return componentFactory->Proposal(index, samplingProblem);
      }

      virtual std::shared_ptr<MultiIndex> FinestIndex() override {
        return componentFactory->FinestIndex();
      }

      virtual std::shared_ptr<MCMCProposal> CoarseProposal (std::shared_ptr<MultiIndex> index,
                                                            std::shared_ptr<AbstractSamplingProblem> coarseProblem,
                                                            std::shared_ptr<SingleChainMCMC> coarseChain) override {
        return componentFactory->CoarseProposal(index, coarseProblem, coarseChain);
      }

      virtual std::shared_ptr<AbstractSamplingProblem> SamplingProblem (std::shared_ptr<MultiIndex> index) override {
        //int idcnt = index->GetValue(0);
        idcnt++;
        if (comm->GetRank() == 0) {
          std::cout << "Rank " << comm->GetRank() << " requesting model " << *index << std::endl;
          for (int dest = 1; dest < comm->GetSize(); dest++) {
            comm->Send(ControlFlag::INIT_PROBLEM, dest, WorkgroupTag);
            comm->Send(index->GetValue(0), dest, WorkgroupTag);
            comm->Send(idcnt, dest, WorkgroupTag);
          }
        }
        return std::make_shared<ParallelAbstractSamplingProblem>(comm, idcnt, componentFactory->SamplingProblem(index));
      }

      virtual std::shared_ptr<MIInterpolation> Interpolation (std::shared_ptr<MultiIndex> index) override {
        return componentFactory->Interpolation(index);
      }

      virtual Eigen::VectorXd StartingPoint (std::shared_ptr<MultiIndex> index) override {
        return componentFactory->StartingPoint(index);
      }

    private:
      int idcnt = 0;
      std::shared_ptr<parcer::Communicator> comm;
      std::shared_ptr<parcer::Communicator> global_comm;
      std::shared_ptr<MIComponentFactory> componentFactory;

      std::map<int, std::shared_ptr<AbstractSamplingProblem>> samplingProblems;

    };

  }
}

#endif

#endif
