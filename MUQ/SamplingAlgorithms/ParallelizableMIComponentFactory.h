#ifndef PARALLELIZABLEMICOMPONENTFACTORY_H_
#define PARALLELIZABLEMICOMPONENTFACTORY_H_

#if MUQ_HAS_MPI

#if !MUQ_HAS_PARCER
#error
#endif


namespace muq {
  namespace SamplingAlgorithms {

    class ParallelizableMIComponentFactory : public MIComponentFactory {
    public:
      virtual void SetComm(std::shared_ptr<parcer::Communicator> comm) = 0;
    };

  }
}

#endif

#endif
