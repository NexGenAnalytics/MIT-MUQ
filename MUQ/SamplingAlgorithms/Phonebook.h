#ifndef PHONEBOOK_H_
#define PHONEBOOK_H_

#if MUQ_HAS_MPI

#if !MUQ_HAS_PARCER
#error
#endif

#include <parcer/Communicator.h>
#include "MUQ/SamplingAlgorithms/ParallelFlags.h"

namespace muq {
  namespace SamplingAlgorithms {

    class PhonebookServer {
    public:

      PhonebookServer(std::shared_ptr<parcer::Communicator> comm)
        : comm(comm)
      {
      }

      void Run() {
        //Dune::Timer timer_idle;
        //Dune::Timer timer_full;

        std::vector< std::tuple<std::shared_ptr<MultiIndex>, int> > requests;

        while (true) {
          MPI_Status status;
          //timer_idle.start();
          ControlFlag command = comm->Recv<ControlFlag>(MPI_ANY_SOURCE, ControlTag, &status);
          //timer_idle.stop();
          if (command == ControlFlag::GET_WORKGROUP) {
            auto index = std::make_shared<MultiIndex>(1, comm->Recv<int>(status.MPI_SOURCE, ControlTag));

            requests.push_back(std::tie(index, status.MPI_SOURCE));
          } else if (command == ControlFlag::SET_WORKGROUP) {
            auto index = std::make_shared<MultiIndex>(1, comm->Recv<int>(status.MPI_SOURCE, ControlTag));
            int rank = comm->Recv<int>(status.MPI_SOURCE, ControlTag, &status);

            if (!phonebook.count(index)) {
              phonebook[index] = WorkerList();
            }
            phonebook[index].AddWorker(rank);
            std::cout << "Phonebook entry for " << *index << " set" << std::endl;
          } else if (command == ControlFlag::UNSET_WORKGROUP) {
            auto index = std::make_shared<MultiIndex>(1, comm->Recv<int>(status.MPI_SOURCE, ControlTag));
            int rank = comm->Recv<int>(status.MPI_SOURCE, ControlTag, &status);

            if (!phonebook.count(index)) {
              std::cerr << "unsetting nonexistent entry!" << std::endl;
            }
            phonebook[index].RemoveWorker(rank);
            if (phonebook[index].NumWorkers() == 0) {
              std::cout << "erasing index " << *index << std::endl;
              phonebook.erase(index);
            }
            std::cout << "Phonebook entry for " << *index << " unset" << std::endl;
          } else if (command == ControlFlag::GET_WORKGROUPS) {
            auto index = std::make_shared<MultiIndex>(1, comm->Recv<int>(status.MPI_SOURCE, ControlTag));

            if (!phonebook.count(index)) {
              std::cerr << "getting workers for nonexistent model!" << std::endl;
            }
            comm->Send(phonebook[index].GetWorkers(), status.MPI_SOURCE, ControlTag);
          } else if (command == ControlFlag::GET_LARGEST_INDEX) {
            if (phonebook.empty()) {
              comm->Send(-1, status.MPI_SOURCE, ControlTag);
              std::cout << "Sent empty largest index" << std::endl;
            } else {
              comm->Send(phonebook.rbegin()->first->GetValue(0), status.MPI_SOURCE, ControlTag);
              std::cout << "Sent largest index " << phonebook.rbegin()->first->GetValue(0) << std::endl;
            }
          } else if (command == ControlFlag::SET_WORKER_READY) {
            auto index = std::make_shared<MultiIndex>(1, comm->Recv<int>(status.MPI_SOURCE, ControlTag));
            if (!phonebook.count(index)) {
              std::cerr << "setting ready for nonexistent model!" << std::endl;
              continue;
            }
            int rank = comm->Recv<int>(status.MPI_SOURCE, ControlTag, &status);
            phonebook[index].SetWorkerReady(rank);
          } else if (command == ControlFlag::QUIT) {
            std::cout << "Rank " << comm->GetRank() << " quit" << std::endl;
            break;
          }

          for (auto request_iter = requests.begin(); request_iter != requests.end();) {
            std::shared_ptr<MultiIndex> index = std::get<0>(*request_iter);
            int sender = std::get<1>(*request_iter);

            if (!phonebook.count(index)) {
              std::cerr << "checking request for nonexistent model!" << std::endl;
            }
            int workerRank = phonebook[index].NextWorker();
            if (workerRank == -1) {
              request_iter++;
            } else {
              comm->Send(workerRank, sender, ControlTag);
              request_iter = requests.erase(request_iter);
            }
          }

        }
        //std::cout << "Phonebook " << comm->GetRank() << " idle time:\t" << timer_idle.elapsed() << " of:\t" << timer_full.elapsed() << std::endl;

      }

    private:

      class WorkerList {
      public:
        int NextWorker() {
          //assert(workers.size() > 0);
          if (workers.size() == 0)
            return -1;

          if (workers_ready.size() == 0)
            return -1;

          int worker = workers_ready.back();
          workers_ready.pop_back();
          return worker;
        }
        void AddWorker(int worker) {
          workers.push_back(worker);
        }
        void RemoveWorker(int worker) {
          workers_ready.erase(std::remove(workers_ready.begin(), workers_ready.end(), worker), workers_ready.end());
          workers.erase(std::remove(workers.begin(), workers.end(), worker), workers.end());
          std::cout << "Removed worker " << worker << ", " << NumWorkers() << " remaining" << std::endl;
        }
        int NumWorkers() {
          return workers.size();
        }
        void SetWorkerReady(int worker) {
          workers_ready.push_back(worker);
        }
        std::vector<int> GetWorkers() {
          return workers;
        }
      private:
        std::vector<int> workers;
        std::vector<int> workers_ready;
      };

      std::map<std::shared_ptr<MultiIndex>, WorkerList, MultiPtrComp> phonebook;
      std::shared_ptr<parcer::Communicator> comm;
    };


    class RoundRobinPhonebookServer {
    public:

      RoundRobinPhonebookServer(std::shared_ptr<parcer::Communicator> comm)
        : comm(comm)
      {
      }

      void Run() {
        //Dune::Timer timer_idle;
        //Dune::Timer timer_full;

        while (true) {
          MPI_Status status;
          //timer_idle.start();
          ControlFlag command = comm->Recv<ControlFlag>(MPI_ANY_SOURCE, ControlTag, &status);
          //timer_idle.stop();
          if (command == ControlFlag::GET_WORKGROUP) {
            auto index = std::make_shared<MultiIndex>(1, comm->Recv<int>(status.MPI_SOURCE, ControlTag));
            int workerRank = phonebook[index].NextWorker();
            comm->Send(workerRank, status.MPI_SOURCE, ControlTag);
            //std::cout << "Phonebook rank " << workerRank << " for model " << *index << " sent to rank " << status.MPI_SOURCE << std::endl;
          } else if (command == ControlFlag::SET_WORKGROUP) {
            auto index = std::make_shared<MultiIndex>(1, comm->Recv<int>(status.MPI_SOURCE, ControlTag));
            int rank = comm->Recv<int>(status.MPI_SOURCE, ControlTag, &status);

            if (!phonebook.count(index)) {
              phonebook[index] = WorkerList();
            }
            phonebook[index].AddWorker(rank);
            std::cout << "Phonebook entry for " << *index << " set" << std::endl;
          } else if (command == ControlFlag::UNSET_WORKGROUP) {
            auto index = std::make_shared<MultiIndex>(1, comm->Recv<int>(status.MPI_SOURCE, ControlTag));
            int rank = comm->Recv<int>(status.MPI_SOURCE, ControlTag, &status);

            if (!phonebook.count(index)) {
              std::cerr << "unsetting nonexistent entry!" << std::endl;
            }
            phonebook[index].RemoveWorker(rank);
            if (phonebook[index].NumWorkers() == 0) {
              std::cout << "erasing index " << *index << std::endl;
              phonebook.erase(phonebook.find(index));
            }
            std::cout << "Phonebook entry for " << *index << " unset" << std::endl;
          } else if (command == ControlFlag::GET_LARGEST_INDEX) {
            if (phonebook.empty()) {
              comm->Send(-1, status.MPI_SOURCE, ControlTag);
              std::cout << "Sent empty largest index" << std::endl;
            } else
              comm->Send(phonebook.rbegin()->first->GetValue(0), status.MPI_SOURCE, ControlTag);
          } else if (command == ControlFlag::SET_WORKER_READY) {
            // Just a dummy receive to be compatible with other phonebooks
            comm->Recv<int>(status.MPI_SOURCE, ControlTag);
            comm->Recv<int>(status.MPI_SOURCE, ControlTag, &status);
          } else if (command == ControlFlag::QUIT) {
            std::cout << "Rank " << comm->GetRank() << " quit" << std::endl;
            break;
          }
        }
        //std::cout << "Phonebook " << comm->GetRank() << " idle time:\t" << timer_idle.elapsed() << " of:\t" << timer_full.elapsed() << std::endl;

      }

    private:

      class WorkerList {
      public:
        int NextWorker() {
          assert(workers.size() > 0);

          if (currentWorker >= workers.size())
            currentWorker = 0;
          int worker = workers[currentWorker];
          currentWorker++;
          return worker;
        }
        void AddWorker(int worker) {
          workers.push_back(worker);
        }
        void RemoveWorker(int worker) {
          workers.erase(std::remove(workers.begin(), workers.end(), worker), workers.end());
        }
        int NumWorkers() {
          return workers.size();
        }
      private:
        std::vector<int> workers;
        int currentWorker;
      };

      std::map<std::shared_ptr<MultiIndex>, WorkerList, MultiPtrComp> phonebook;
      std::shared_ptr<parcer::Communicator> comm;
    };


    class PhonebookClient {

    public:

      PhonebookClient(std::shared_ptr<parcer::Communicator> comm, int phonebookRank)
        : comm(comm), phonebookRank(phonebookRank)
      {
      }

      int Query(std::shared_ptr<MultiIndex> remoteIndex) {
        comm->Send(ControlFlag::GET_WORKGROUP, phonebookRank, ControlTag);
        comm->Send(remoteIndex->GetValue(0), phonebookRank, ControlTag);
        return comm->Recv<int>(phonebookRank, ControlTag);
      }

      std::vector<int> GetWorkgroups(std::shared_ptr<MultiIndex> modelIndex) {
        comm->Send(ControlFlag::GET_WORKGROUPS, phonebookRank, ControlTag);
        comm->Send(modelIndex->GetValue(0), phonebookRank, ControlTag);
        return comm->Recv<std::vector<int>>(phonebookRank, ControlTag);
      }

      std::shared_ptr<MultiIndex> LargestIndex() {
        comm->Send(ControlFlag::GET_LARGEST_INDEX, phonebookRank, ControlTag);
        return std::make_shared<MultiIndex>(1, comm->Recv<int>(phonebookRank, ControlTag));
      }

      void Register(std::shared_ptr<MultiIndex> modelIndex, int rank) {
        comm->Send(ControlFlag::SET_WORKGROUP, phonebookRank, ControlTag);
        comm->Send(modelIndex->GetValue(0), phonebookRank, ControlTag);
        comm->Ssend(rank, phonebookRank, ControlTag);
      }

      void UnRegister(std::shared_ptr<MultiIndex> modelIndex, int rank) {
        comm->Send(ControlFlag::UNSET_WORKGROUP, phonebookRank, ControlTag);
        comm->Send(modelIndex->GetValue(0), phonebookRank, ControlTag);
        comm->Ssend(rank, phonebookRank, ControlTag);
      }

      void SetWorkerReady(std::shared_ptr<MultiIndex> modelIndex, int rank) {
        comm->Send(ControlFlag::SET_WORKER_READY, phonebookRank, ControlTag);
        comm->Send(modelIndex->GetValue(0), phonebookRank, ControlTag);
        comm->Send(rank, phonebookRank, ControlTag);
      }

    private:
      std::shared_ptr<parcer::Communicator> comm;
      int phonebookRank;
    };
  }
}

#endif

#endif
