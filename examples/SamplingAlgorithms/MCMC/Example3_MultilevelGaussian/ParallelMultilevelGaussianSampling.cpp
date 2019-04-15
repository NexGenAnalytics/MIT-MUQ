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

#include <boost/property_tree/ptree.hpp>

namespace pt = boost::property_tree;
using namespace muq::Modeling;
using namespace muq::SamplingAlgorithms;
using namespace muq::Utilities;

#include "Problem.h"

#include "MUQ/SamplingAlgorithms/ParallelMIMCMCWorker.h"



int main(int argc, char **argv){

  MPI_Init(&argc, &argv);

  pt::ptree pt;
  pt.put("MCMC.burnin", 10); // number of samples for single level


  auto comm = std::make_shared<parcer::Communicator>();

  std::cout << "Rank: " <<  comm->GetRank() << std::endl;

  int rootRank = 0;
  int phonebookRank = 1;

  if (comm->GetRank() == rootRank) {

    auto comm_self = std::make_shared<parcer::Communicator>(MPI_COMM_SELF);

    std::shared_ptr<Eigen::VectorXd> measurements = std::make_shared<Eigen::VectorXd>(0);
    auto componentFactory = std::make_shared<MyMIComponentFactory>(comm_self, comm, measurements);


    auto phonebookClient = std::make_shared<PhonebookClient>(comm, phonebookRank);
    WorkerClient workerClient(comm, phonebookClient, rootRank);

    std::vector<CollectorClient> collectorClients;
    auto indices = MultiIndexFactory::CreateFullTensor(componentFactory->FinestIndex()->GetVector());

    int rank = 2;

    for (int i = 0; i < indices->Size(); i++) {
      std::shared_ptr<MultiIndex> index = (*indices)[i];
      collectorClients.push_back(CollectorClient(comm, {rank}, index));
      rank++;
    }

    for (CollectorClient& client : collectorClients) {
      client.CollectSamples(5);
    }

    /*collector_0.CollectSamples(10000);
    collector_1.CollectSamples(10000);
    collector_2.CollectSamples(1000);
    collector_3.CollectSamples(1000);*/
    //collector_2.CollectSamples(10);

    // Assign work groups
    int specialRanks = rank;
    //int specialRanks = 5;

    int modelID = 0;
    for (int rank = specialRanks; rank < comm->GetSize(); rank++) {

      //if (rank == comm->GetSize()-1)
      //  modelID = 2;
      workerClient.assignGroup({rank}, std::make_shared<MultiIndex>(1, modelID), measurements);

      modelID++;
      if (modelID > 3)
        modelID = 0;
    }

    //Dune::Timer timer_idle;
    //Dune::Timer timer_full;
    while (true) {
      MPI_Status status;
      //timer_idle.start();
      ControlFlag command = comm->Recv<ControlFlag>(MPI_ANY_SOURCE, ControlTag, &status);
      //timer_idle.stop();

      for (CollectorClient& client : collectorClients) {
        if (client.Receive(command, status))
          break;
      }

      /*if (!command_handled) {
       s td::cerr << "Unexpected command!" << std::*endl;
       exit(43);
    }*/

      bool isSampling = false;
      for (CollectorClient& client : collectorClients) {
        isSampling = isSampling || client.IsSampling();
      }
      if (!isSampling)
        break;
    }
    std::cout << "Sampling completed" << std::endl;


    for (CollectorClient& client : collectorClients) {
      client.ComputeMeans();
    }

    while (true) {
      MPI_Status status;
      //timer_idle.start();
      ControlFlag command = comm->Recv<ControlFlag>(MPI_ANY_SOURCE, ControlTag, &status);
      //timer_idle.stop();

      for (CollectorClient& client : collectorClients) {
        if (client.Receive(command, status))
          break;
      }

      /*if (!command_handled) {
       std::cerr << "Unexpected command!" << std::*endl;
       exit(43);
      }*/


      bool isComputingMeans = false;
      for (CollectorClient& client : collectorClients) {
        isComputingMeans = isComputingMeans || client.IsComputingMeans();
      }
      if (!isComputingMeans)
        break;
    }
    std::cout << "Computing means completed" << std::endl;
    //std::cout << "Root process " << comm->GetRank() << " idle time:\t" << timer_idle.elapsed() << " of:\t" << timer_full.elapsed() << std::endl;





    Eigen::VectorXd mean_box = collectorClients[0].GetQOIMean();
    mean_box.setZero();
    Eigen::VectorXd mean = mean_box;

    for (CollectorClient& client : collectorClients) {
      mean_box = client.GetQOIMean();
      mean += mean_box;
      std::cout << "Mean level:\t" << mean_box.transpose() << " adding up to:\t" << mean.transpose() << std::endl;
    }



    std::cout << "Starting unassign sequence" << std::endl;
    //for (int i = collectorClients.size() - 1; i >= 0; i--)
    //  collectorClients[i].Unassign();
    for (CollectorClient& client : collectorClients) {
      client.Unassign();
    }
    workerClient.UnassignAll();
    std::cout << "Finished unassign sequence" << std::endl;

    workerClient.Finalize();
    std::cout << "Rank " << comm->GetRank() << " quit" << std::endl;


  } else if (comm->GetRank() == phonebookRank) {
    PhonebookServer phonebook(comm);
    phonebook.Run();
  } else {
    auto phonebookClient = std::make_shared<PhonebookClient>(comm, phonebookRank);
    WorkerServer<MyMIComponentFactory> worker(pt, comm, phonebookClient, rootRank);
  }

  MPI_Finalize();
}
