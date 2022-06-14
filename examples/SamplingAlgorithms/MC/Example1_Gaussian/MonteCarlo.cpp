#include "MUQ/Modeling/Distributions/Gaussian.h"
#include "MUQ/Modeling/Distributions/UniformBox.h"
#include "MUQ/Modeling/Distributions/Density.h"
#include "MUQ/Modeling/UMBridge/UMBridgeModPiece.h"
#include "MUQ/Modeling/OneStepCachePiece.h"

#include "MUQ/SamplingAlgorithms/DummyKernel.h"
#include "MUQ/SamplingAlgorithms/MHProposal.h"
#include "MUQ/SamplingAlgorithms/MHKernel.h"
#include "MUQ/SamplingAlgorithms/SamplingProblem.h"
#include "MUQ/SamplingAlgorithms/SingleChainMCMC.h"
#include "MUQ/SamplingAlgorithms/MCMCFactory.h"

#include <boost/property_tree/ptree.hpp>
#include <thread>

namespace pt = boost::property_tree;
using namespace muq::Modeling;
using namespace muq::SamplingAlgorithms;
using namespace muq::Utilities;

#include "MCSampleProposal.h"

void thread_draw_samples(int index, int num_samples, std::string url) {
  json config;
  //config["level"] = index->GetValue(0);
  //int port = 4243 + index;
  //std::string url = "http://localhost:" + std::to_string(port);
  auto model = std::make_shared<UMBridgeModPiece>(url, config);
  auto model_cached = std::make_shared<OneStepCachePiece>(model);

  /*Eigen::MatrixXd box_bounds(3,2);
  box_bounds << 1.0, 1.05,
         1.0, 1.05,
         1.0, 1.05;*/
  Eigen::MatrixXd box_bounds(1,2);
  box_bounds << 1.0, 1.05;
  auto uniform_box_dist = std::make_shared<UniformBox>(box_bounds);


  auto problem = std::make_shared<SamplingProblem>(model_cached, model_cached);

  pt::ptree pt;
  pt.put("NumSamples", num_samples); // number of MCMC steps
  pt.put("BurnIn", 0);
  pt.put("PrintLevel",3);

  auto proposalDensity = uniform_box_dist->AsDensity();
  auto proposal = std::make_shared<MCSampleProposal>(pt, problem, proposalDensity);
  std::vector<std::shared_ptr<TransitionKernel> > kernels = {std::make_shared<DummyKernel>(pt, problem, proposal)};
  auto mc = std::make_shared<SingleChainMCMC>(pt, kernels);

  //Eigen::VectorXd startPt(3);
  Eigen::VectorXd startPt(1);
  mc->Run(startPt);

  std::shared_ptr<SampleCollection> samps = mc->GetQOIs();

  Eigen::VectorXd sampMean = samps->Mean();
  std::cout << "\nSample Mean = \n" << sampMean.transpose() << std::endl;
}

int main(int argc, char** argv){

  int num_threads = std::stoi(argv[1]);
  int num_samples = std::stoi(argv[2]);
  std::string url = argv[3];

  std::cout << "Computing " << num_samples << " samples on " << num_threads << " threads with model " << url << std::endl;

  std::vector<std::shared_ptr<std::thread>> threads(0);
  for (int thread_index = 0; thread_index < num_threads; thread_index++) {

    int thread_num_samples = num_samples / (num_threads - thread_index);
    num_samples -= thread_num_samples;

    auto thread = std::make_shared<std::thread>(thread_draw_samples, thread_index, thread_num_samples, url);
    threads.push_back(thread);
  }

  std::cout << "Samples not queued: " << num_samples << std::endl;

  for (auto thread : threads) {
    thread->join();
  }
  /*Eigen::VectorXd sampVar = samps->Variance();
  std::cout << "\nSample Variance = \n" << sampVar.transpose() << std::endl;

  Eigen::MatrixXd sampCov = samps->Covariance();
  std::cout << "\nSample Covariance = \n" << sampCov << std::endl;*/

  //Eigen::VectorXd sampMom3 = samps->CentralMoment(3);
  //std::cout << "\nSample Third Moment = \n" << sampMom3 << std::endl << std::endl;

  return 0;
}
