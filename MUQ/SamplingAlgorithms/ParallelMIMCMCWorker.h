#ifndef PARALLELMIMCMCWORKER_H_
#define PARALLELMIMCMCWORKER_H_

#if MUQ_HAS_MPI

#if !MUQ_HAS_PARCER
#error
#endif


#include "MUQ/SamplingAlgorithms/MarkovChain.h"
#include "MUQ/SamplingAlgorithms/DistributedCollection.h"
#include "MUQ/SamplingAlgorithms/ParallelFlags.h"
#include "MUQ/SamplingAlgorithms/ParallelMIComponentFactory.h"
#include "MUQ/SamplingAlgorithms/ParallelMIMCMCBox.h"
#include "MUQ/Utilities/Cereal/MultiIndexSerializer.h"

namespace muq {
  namespace SamplingAlgorithms {

    class CollectorClient {
    public:
      CollectorClient(std::shared_ptr<parcer::Communicator> comm, std::vector<int> subgroup, std::shared_ptr<MultiIndex> modelindex)
       : comm(comm), subgroup(subgroup), boxHighestIndex(modelindex) {
         for (int dest : subgroup) {
           comm->Send(ControlFlag::ASSIGN_COLLECTOR, dest, ControlTag);
           comm->Send(subgroup, dest, ControlTag);
           comm->Send(*boxHighestIndex, dest, ControlTag);
         }

         // Set up Multiindex box
         boxLowestIndex = MultiIndex::Copy(boxHighestIndex);
         --(*boxLowestIndex);
         std::shared_ptr<MultiIndex> boxSize = std::make_shared<MultiIndex>(*boxHighestIndex - *boxLowestIndex);
         boxIndices = MultiIndexFactory::CreateFullTensor(boxSize->GetVector());
      }

      void CollectSamples (int numSamples) {
        sampling = true;
        for (int dest : subgroup) {
          comm->Send(ControlFlag::SAMPLE_BOX, dest, ControlTag);
          comm->Send(numSamples, dest, ControlTag);
        }
      }

      void ComputeMeans() {
        computingMeans = true;
        for (int dest : subgroup) {
          comm->Send(ControlFlag::MEANS, dest, ControlTag);
        }
      }

      bool Receive (ControlFlag command, const MPI_Status& status) {
        if (status.MPI_SOURCE != subgroup[0])
          return false;

        if (command == ControlFlag::SAMPLE_BOX_DONE) {
          sampling = false;
        } else if (command == ControlFlag::MEANS_DONE) {

          for (uint i = 0; i < boxIndices->Size(); i++) {
            std::shared_ptr<MultiIndex> boxIndex = (*boxIndices)[i];

            Eigen::VectorXd chainSampleMean = comm->Recv<Eigen::VectorXd>(status.MPI_SOURCE, ControlTag);
            Eigen::VectorXd chainQOIMean = comm->Recv<Eigen::VectorXd>(status.MPI_SOURCE, ControlTag);

            std::shared_ptr<MultiIndex> index = std::make_shared<MultiIndex>(*boxLowestIndex + *boxIndex);
            auto indexDiffFromTop = std::make_shared<MultiIndex>(*boxHighestIndex - *index);

            if (i == 0) {
              if (indexDiffFromTop->Sum() % 2 == 0) {
                boxQOIMean = chainQOIMean;
              } else {
                boxQOIMean = -chainQOIMean;
              }
            } else {
              if (indexDiffFromTop->Sum() % 2 == 0) {
                boxQOIMean += chainQOIMean;
              } else {
                boxQOIMean -= chainQOIMean;
              }
            }
          }
          computingMeans = false;
        } else {
          std::cerr << "Unexpected command!" << std::endl;
          exit(43);
        }

        return true;
      }

      Eigen::VectorXd GetQOIMean() {
        return boxQOIMean;
      }

      void Unassign() {
        for (int dest : subgroup) {
          comm->Send(ControlFlag::UNASSIGN, dest, ControlTag);
        }
      }

      bool IsSampling() {
        return sampling;
      }
      bool IsComputingMeans() {
        return computingMeans;
      }

    private:
      bool sampling = false;
      bool computingMeans = false;
      std::shared_ptr<MultiIndexSet> boxIndices;
      Eigen::VectorXd boxQOIMean;

      std::shared_ptr<MultiIndex> boxHighestIndex;
      std::shared_ptr<MultiIndex> boxLowestIndex;

      std::shared_ptr<parcer::Communicator> comm;
      std::vector<int> subgroup;
      std::map<std::shared_ptr<MultiIndex>, Eigen::VectorXd, MultiPtrComp> means;

    };

    class WorkerClient {
    public:
      WorkerClient(std::shared_ptr<parcer::Communicator> comm, std::shared_ptr<PhonebookClient> phonebookClient, int RootRank)
       : comm(comm), phonebookClient(phonebookClient) {
      }


      void assignGroup (std::vector<int> subgroup, std::shared_ptr<MultiIndex> modelindex,
                        std::shared_ptr<Eigen::VectorXd> measurements) {
        for (int dest : subgroup) {
          comm->Send(ControlFlag::ASSIGN, dest, ControlTag);
          comm->Send(subgroup, dest, ControlTag);
          comm->Send(*modelindex, dest, ControlTag);
          comm->Send(*measurements, dest, ControlTag);
        }
        phonebookClient->Register(modelindex, subgroup[0]);
      }

      void UnassignGroup (std::shared_ptr<MultiIndex> modelIndex, int groupRootRank) {
        std::cout << "Sending unassign to " << groupRootRank << std::endl;
        phonebookClient->UnRegister(modelIndex, groupRootRank);
        comm->Ssend(ControlFlag::UNASSIGN, groupRootRank, ControlTag);
      }

      void UnassignAll() {
        std::shared_ptr<MultiIndex> largest = nullptr;
        do {
          largest = phonebookClient->LargestIndex();
          std::cout << "Unassigning model " << *largest << std::endl;
          std::vector<int> ranks = phonebookClient->GetWorkgroups(largest);
          for (int rank : ranks) {
            UnassignGroup(largest, rank);
          }
        } while (largest->Max() != 0);
      }

      void Finalize() {
        for (int dest = 1; dest < comm->GetSize(); dest++)
          comm->Send(ControlFlag::QUIT, dest, ControlTag);
      }

    private:
      std::shared_ptr<parcer::Communicator> comm;
      std::shared_ptr<PhonebookClient> phonebookClient;
    };

    template <typename COMPONENT_FACTORY>
    class WorkerServer {
    public:
      WorkerServer(boost::property_tree::ptree const& pt, std::shared_ptr<parcer::Communicator> comm, std::shared_ptr<PhonebookClient> phonebookClient, int RootRank) {

        while (true) {
          ControlFlag command = comm->Recv<ControlFlag>(RootRank, ControlTag);
          if (command == ControlFlag::ASSIGN) {
            std::vector<int> subgroup_proc = comm->Recv<std::vector<int>>(0, ControlTag);
            auto samplingProblemIndex = std::make_shared<MultiIndex>(comm->Recv<MultiIndex>(0, ControlTag));
            auto measurements = std::make_shared<Eigen::VectorXd>(comm->Recv<Eigen::VectorXd>(0, ControlTag));

            MPI_Group world_group;
            MPI_Comm_group (MPI_COMM_WORLD, &world_group);
            MPI_Group subgroup;
            MPI_Group_incl (world_group, subgroup_proc.size(), &subgroup_proc[0], &subgroup);

            MPI_Comm subcomm_raw;
            MPI_Comm_create_group(MPI_COMM_WORLD, subgroup, ControlTag, &subcomm_raw);
            auto subcomm = std::make_shared<parcer::Communicator>(subcomm_raw);

            auto componentFactory = std::make_shared<COMPONENT_FACTORY>(subcomm, comm, measurements);
            auto parallelComponentFactory = std::make_shared<ParallelMIComponentFactory>(subcomm, comm, componentFactory);

            if (subcomm->GetRank() == 0) {
              std::cout << "Subgroup root is global " << comm->GetRank() << std::endl;
              auto finestProblem = parallelComponentFactory->SamplingProblem(parallelComponentFactory->FinestIndex());

              auto box = std::make_shared<ParallelMIMCMCBox>(parallelComponentFactory, samplingProblemIndex, comm, phonebookClient);

              // Burn in coarsest chains
              if (samplingProblemIndex->Max() == 0) {
                std::cout << "Burning in" << std::endl;
                for (int i = 0; i < configuration.get<int>("MCMC.burnin"); i++) // TODO: Proper solution for burnin
                  box->Sample();
                std::cout << "Burned in" << std::endl;
              }

              for (int i = 0; i < 5; i++) // TODO: Really subsampling on every level? Maybe subsample when requesting samples?
                box->Sample();
              phonebookClient->SetWorkerReady(samplingProblemIndex, comm->GetRank());

              //Dune::Timer timer_idle;
              //Dune::Timer timer_full;
              while (true) {
                MPI_Status status;
                //timer_idle.start();
                command = comm->Recv<ControlFlag>(MPI_ANY_SOURCE, ControlTag, &status);
                //timer_idle.stop();
                if (command == ControlFlag::UNASSIGN) {
                  break;
                } else if (command == ControlFlag::SAMPLE) {
                  //std::cout << "Rank " << comm->GetRank() << " sending sample to rank " << status.MPI_SOURCE << std::endl;
                  auto sampleCollection = box->FinestChain()->GetSamples(); // TODO: last() function for collection? // TODO: Do not store chains here
                  auto latestSample = sampleCollection->at(sampleCollection->size()-1);
                  // TODO: Send "full" sample via parcer?
                  comm->Send<Eigen::VectorXd>(latestSample->state[0], status.MPI_SOURCE, ControlTag);
                  if (latestSample->HasMeta("QOI")) {
                    std::shared_ptr<SamplingState> qoi = AnyCast(latestSample->meta["QOI"]);
                    comm->Send<Eigen::VectorXd>(qoi->state[0], status.MPI_SOURCE, ControlTag);
                  } else {
                    std::cerr << "No QOI!" << std::endl;
                    exit(47);
                  }


                  for (int i = 0; i < 5; i++)
                    box->Sample();
                  phonebookClient->SetWorkerReady(samplingProblemIndex, comm->GetRank());
                } else if (command == ControlFlag::SAMPLE_BOX) {
                  //std::cout << "Rank " << comm->GetRank() << " sending sample to rank " << status.MPI_SOURCE << std::endl;
                  for (int i = 0; i < box->NumChains(); i++) {
                    auto sampleCollection = box->GetChain(i)->GetSamples(); // TODO: last() function for collection? // TODO: Do not store chains here
                    auto latestSample = sampleCollection->at(sampleCollection->size()-1);
                    // TODO: Send "full" sample via parcer?
                    comm->Send<Eigen::VectorXd>(latestSample->state[0], status.MPI_SOURCE, ControlTag);
                    if (latestSample->HasMeta("QOI")) {
                      std::shared_ptr<SamplingState> qoi = AnyCast(latestSample->meta["QOI"]);
                      comm->Send<Eigen::VectorXd>(qoi->state[0], status.MPI_SOURCE, ControlTag);
                    } else {
                      std::cerr << "No QOI!" << std::endl;
                      exit(47);
                    }
                  }

                  for (int i = 0; i < 5; i++)
                    box->Sample();
                  phonebookClient->SetWorkerReady(samplingProblemIndex, comm->GetRank());
                } else {
                  std::cerr << "Unexpected command!" << std::endl;
                  exit(43);
                }

              }
              //std::cout << "Worker Controller " << comm->GetRank() << " idle time:\t" << timer_idle.elapsed() << " of:\t" << timer_full.elapsed() << std::endl;

            }

            parallelComponentFactory->finalize();
            std::cout << "Rank " << comm->GetRank() << " finalized" << std::endl;
          } else if (command == ControlFlag::ASSIGN_COLLECTOR) {
            std::vector<int> subgroup_proc = comm->Recv<std::vector<int>>(RootRank, ControlTag);
            auto boxHighestIndex = std::make_shared<MultiIndex>(comm->Recv<MultiIndex>(RootRank, ControlTag));

            // Set up subcommunicator
            MPI_Group world_group;
            MPI_Comm_group (MPI_COMM_WORLD, &world_group);
            MPI_Group subgroup;
            MPI_Group_incl (world_group, subgroup_proc.size(), &subgroup_proc[0], &subgroup);

            MPI_Comm subcomm_raw;
            MPI_Comm_create_group(MPI_COMM_WORLD, subgroup, ControlTag, &subcomm_raw);
            auto subcomm = std::make_shared<parcer::Communicator>(subcomm_raw);


            // Set up Multiindex box
            auto boxLowestIndex = MultiIndex::Copy(boxHighestIndex);
            --(*boxLowestIndex);
            std::shared_ptr<MultiIndex> boxSize = std::make_shared<MultiIndex>(*boxHighestIndex - *boxLowestIndex);
            std::shared_ptr<MultiIndexSet> boxIndices = MultiIndexFactory::CreateFullTensor(boxSize->GetVector());

            std::vector<std::shared_ptr<DistributedCollection>> sampleCollections(boxIndices->Size());
            std::vector<std::shared_ptr<DistributedCollection>> qoiCollections(boxIndices->Size());
            for (uint i = 0; i < boxIndices->Size(); i++) {
              auto sampleCollection = std::make_shared<MarkovChain>();
              sampleCollections[i] = std::make_shared<DistributedCollection>(sampleCollection, subcomm);
              auto qoiCollection = std::make_shared<MarkovChain>();
              qoiCollections[i] = std::make_shared<DistributedCollection>(qoiCollection, subcomm);
            }


            //Dune::Timer timer_idle;
            //Dune::Timer timer_full;
            while (true) {
              MPI_Status status;
              //timer_idle.start();
              command = comm->Recv<ControlFlag>(MPI_ANY_SOURCE, ControlTag, &status);
              //timer_idle.stop();
              if (command == ControlFlag::UNASSIGN) {
                break;
              } else if (command == ControlFlag::SAMPLE_BOX) {

                int numSamples = comm->Recv<int>(0, ControlTag);

            		//std::cout << "Rank " << comm->GetRank() << " requesting sample from rank " << remoteRank << std::endl;

                for (int i = 0; i < numSamples; i++) {
                  int remoteRank = phonebookClient->Query(boxHighestIndex);
                  comm->Send(ControlFlag::SAMPLE_BOX, remoteRank, ControlTag); // TODO: Receive sample in one piece?
                  for (uint i = 0; i < boxIndices->Size(); i++) {
                    //std::shared_ptr<MultiIndex> boxIndex = (*boxIndices)[i];
                    sampleCollections[i]->Add(std::make_shared<SamplingState>(comm->Recv<Eigen::VectorXd>(remoteRank, ControlTag)));
                    qoiCollections[i]->Add(std::make_shared<SamplingState>(comm->Recv<Eigen::VectorXd>(remoteRank, ControlTag)));
                		//Eigen::VectorXd remoteQOI = comm->Recv<Eigen::VectorXd>(remoteRank, ControlTag);
                  }
                  std::cout << "Collector " << comm->GetRank() << " received sample box " << *boxHighestIndex << std::endl;
                }

                if (subcomm->GetRank() == 0)
                  comm->Send(ControlFlag::SAMPLE_BOX_DONE, RootRank, ControlTag); // TODO: Receive sample in one piece?
                //  box->Sample();
              } else if (command == ControlFlag::MEANS) {
                std::list<Eigen::VectorXd> sampleMeans;
                std::list<Eigen::VectorXd> qoiMeans;
                for (uint i = 0; i < boxIndices->Size(); i++) {
                  sampleMeans.push_back(sampleCollections[i]->GlobalMean());
                  qoiMeans.push_back(qoiCollections[i]->GlobalMean());
                }
                if (subcomm->GetRank() == 0) {
                  comm->Send(ControlFlag::MEANS_DONE, RootRank, ControlTag);
                  auto qoiMean = qoiMeans.begin();
                  for (auto sampleMean = sampleMeans.begin(); sampleMean != sampleMeans.end(); sampleMean++) {
                    comm->Send(*sampleMean, RootRank, ControlTag);
                    comm->Send(*qoiMean, RootRank, ControlTag);
                    qoiMean++;
                  }
                }
              } else {
                std::cerr << "Unexpected command!" << std::endl;
                exit(43);
              }

            }
            //std::cout << "Collector " << comm->GetRank() << " idle time:\t" << timer_idle.elapsed() << " of:\t" << timer_full.elapsed() << std::endl;


          } else if (command == ControlFlag::QUIT) {
            std::cout << "Rank " << comm->GetRank() << " quit" << std::endl;
            break;
          } else {
            std::cerr << "Unexpected command!" << std::endl;
            exit(42);
          }
        }

      }

    };

  }
}

#endif

#endif
