#include "MUQ/SamplingAlgorithms/SingleChainMCMC.h"

#include <chrono>

#include "MUQ/Utilities/AnyHelpers.h"
#include "MUQ/Utilities/StringUtilities.h"

#include "MUQ/SamplingAlgorithms/MarkovChain.h"
#include "MUQ/SamplingAlgorithms/ExpensiveSamplingProblem.h"

namespace pt = boost::property_tree;
using namespace muq::Utilities;
using namespace muq::SamplingAlgorithms;

SingleChainMCMC::SingleChainMCMC(pt::ptree pt, std::shared_ptr<AbstractSamplingProblem> problem, std::vector<Eigen::VectorXd> x0) :
  SamplingAlgorithm(std::make_shared<MarkovChain>()),
  printLevel(pt.get("PrintLevel",3)),
  prevState(std::make_shared<SamplingState>(x0))
{
  SetUp(pt, problem);
}

#if MUQ_HAS_PARCER
SingleChainMCMC::SingleChainMCMC(pt::ptree pt, std::shared_ptr<AbstractSamplingProblem> problem, std::vector<Eigen::VectorXd> x0, std::shared_ptr<parcer::Communicator> comm) :
  SamplingAlgorithm(std::make_shared<MarkovChain>(), comm),
  printLevel(pt.get("PrintLevel",3)),
  prevState(std::make_shared<SamplingState>(x0))
{
  SetUp(pt, problem);
}
#endif

SingleChainMCMC::SingleChainMCMC(boost::property_tree::ptree pt,
                //std::shared_ptr<AbstractSamplingProblem> problem,
                std::vector<std::shared_ptr<TransitionKernel>> kernelsIn,
                std::vector<Eigen::VectorXd> x0) :
                SamplingAlgorithm(std::make_shared<MarkovChain>()),
                printLevel(pt.get("PrintLevel",3)),
                kernels(kernelsIn),
                prevState(std::make_shared<SamplingState>(x0))
{
  // TODO: clean this up, maybe somehow merge with SetUp(..)
  numSamps = pt.get<unsigned int>("NumSamples");
  burnIn = pt.get("BurnIn",0);
  if(burnIn==0) {
    samples->Add(prevState);
    prevState = samples->at(0); // Add() copies state, so need to retrieve copy here
  }

  scheduler = std::make_shared<ThinScheduler>(pt);
  schedulerQOI = std::make_shared<ThinScheduler>(pt);
}

SingleChainMCMC::SingleChainMCMC(boost::property_tree::ptree pt,
                std::vector<std::shared_ptr<TransitionKernel>> kernelsIn,
                std::shared_ptr<SamplingState> x0) :
                SamplingAlgorithm(std::make_shared<MarkovChain>()),
                printLevel(pt.get("PrintLevel",3)),
                kernels(kernelsIn),
                prevState(x0)
{
  // TODO: clean this up, maybe somehow merge with SetUp(..)
  numSamps = pt.get<unsigned int>("NumSamples");
  burnIn = pt.get("BurnIn",0);
  if(burnIn==0) {
    samples->Add(prevState);
    prevState = samples->at(0); // Add() copies state, so need to retrieve copy here
  }

  scheduler = std::make_shared<ThinScheduler>(pt);
  schedulerQOI = std::make_shared<ThinScheduler>(pt);
}

void SingleChainMCMC::SetUp(pt::ptree pt, std::shared_ptr<AbstractSamplingProblem> problem) {
  numSamps = pt.get<unsigned int>("NumSamples");
  burnIn = pt.get("BurnIn",0);
  if(burnIn==0) {
    samples->Add(prevState);
    prevState = samples->at(0); // Add() copies state, so need to retrieve copy here
  }

  std::string kernelString = pt.get<std::string>("KernelList");

  std::vector<std::string> kernelNames = StringUtilities::Split(kernelString, ',');

  unsigned int numBlocks = kernelNames.size();
  kernels.resize(numBlocks);

  scheduler = std::make_shared<ThinScheduler>(pt);
  schedulerQOI = std::make_shared<ThinScheduler>(pt);

  // Add the block id to a child tree and construct a kernel for each block
  for(int i=0; i<kernels.size(); ++i) {
    boost::property_tree::ptree subTree = pt.get_child(kernelNames.at(i));
    subTree.put("BlockIndex",i);
    if( std::dynamic_pointer_cast<ExpensiveSamplingProblem>(problem) ) { subTree.put("ReevaluateAcceptedDensity", true); }
    kernels.at(i) = TransitionKernel::Construct(subTree, problem);

#if MUQ_HAS_PARCER
    kernels.at(i)->SetCommunicator(comm);
#endif
  }
}

void SingleChainMCMC::PrintStatus(std::string prefix) {
  PrintStatus(prefix, numSamps+1);
}


void SingleChainMCMC::PrintStatus(std::string prefix, unsigned int currInd) const
{
  std::cout << prefix << int(std::floor(double((currInd - 1) * 100) / double(numSamps))) << "% Complete" << std::endl;
  if(printLevel>1){
    for(int blockInd=0; blockInd<kernels.size(); ++blockInd){
      std::cout << prefix << "  Block " << blockInd << ":\n";
      kernels.at(blockInd)->PrintStatus(prefix + "    ");
    }
  }
}

std::shared_ptr<SampleCollection> SingleChainMCMC::RunImpl()
{

  // What is the next iteration that we want to print at
  const unsigned int printIncr = std::floor(numSamps / double(10));
  unsigned int nextPrintInd = printIncr;

  // Run until we've received the desired number of samples
  if(printLevel>0)
    std::cout << "Starting single chain MCMC sampler..." << std::endl;

  while(sampNum < numSamps)
  {
    // Should we print
    if(sampNum > nextPrintInd){
      if(printLevel>0){
        PrintStatus("  ", sampNum);
      }
      nextPrintInd += printIncr;
    }

    Sample();
  }


  if(printLevel>0){
    PrintStatus("  ", numSamps+1);
    std::cout << "Completed in " << totalTime << " seconds." << std::endl;
  }

  return samples;
}

void SingleChainMCMC::Sample() {
  auto startTime = std::chrono::high_resolution_clock::now();

  std::vector<std::shared_ptr<SamplingState>> newStates;

  // Loop through each parameter block
  for(int blockInd=0; blockInd<kernels.size(); ++blockInd){
    // kernel prestep
    kernels.at(blockInd)->PreStep(sampNum, prevState);

    // use the kernel to get the next state(s)
    newStates = kernels.at(blockInd)->Step(sampNum, prevState);

    // kernel post-processing
    kernels.at(blockInd)->PostStep(sampNum, newStates);

    // add the new states to the SampleCollection (this also increments sampNum)
    prevState = SaveSamples(newStates, lastSavedState, sampNum);
  }

  auto endTime = std::chrono::high_resolution_clock::now();
  totalTime += std::chrono::duration<double>(endTime - startTime).count();
}

std::shared_ptr<SamplingState> SingleChainMCMC::SaveSamples(std::vector<std::shared_ptr<SamplingState> > const& newStates, std::shared_ptr<SamplingState>& lastSavedState, unsigned int& sampNum) const {
  for( auto it : newStates ) {
    // save the sample, if we want to
    if( ShouldSave(sampNum) ) {
      samples->Add(it);
      if (it->HasMeta("QOI")) {
        std::shared_ptr<SamplingState> qoi = AnyCast(it->meta["QOI"]);
        QOIs->Add(qoi);
      }
    }

    // increment the number of samples and break of we hit the max. number
    if( ++sampNum>=numSamps ) { return it; }
  }

  return newStates.back();
}

bool SingleChainMCMC::ShouldSave(unsigned int const sampNum) const { return sampNum>=burnIn && scheduler->ShouldSave(sampNum); }
