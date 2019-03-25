#ifndef SAMPLINGALGORITHM_H_
#define SAMPLINGALGORITHM_H_

#include "MUQ/config.h"

#if MUQ_HAS_PARCER
#include <parcer/Communicator.h>
#endif

#include "MUQ/Modeling/WorkPiece.h"

#include "MUQ/SamplingAlgorithms/AbstractSamplingProblem.h"
#include "MUQ/SamplingAlgorithms/SampleCollection.h"

/**
@defgroup SamplingAlgorithms

*/

namespace muq {
  namespace SamplingAlgorithms {

    class SamplingAlgorithm {//} : public muq::Modeling::WorkPiece {
    public:

      SamplingAlgorithm(std::shared_ptr<SampleCollection> samplesIn,
                        std::shared_ptr<SampleCollection> QOIsIn);

      SamplingAlgorithm(std::shared_ptr<SampleCollection> samplesIn);

#if MUQ_HAS_PARCER
      SamplingAlgorithm(std::shared_ptr<SampleCollection> samplesIn, std::shared_ptr<parcer::Communicator> comm);
#endif

      virtual ~SamplingAlgorithm() = default;

      virtual std::shared_ptr<SampleCollection> GetSamples() const;
      virtual std::shared_ptr<SampleCollection> GetQOIs() const;

      virtual std::shared_ptr<SampleCollection> Run();

#if MUQ_HAS_PARCER
      std::shared_ptr<parcer::Communicator> GetCommunicator() const;
#endif

    protected:

      virtual std::shared_ptr<SampleCollection> RunImpl() = 0;

      /**
	 Inputs:
	 <ol>
	 <li> Parameters for the algorithm
	 <li> The muq::SamplingAlgorithms::SamplingProblem that evaluates/samples the target distribution
	 </ol>
	 @param[in] inputs Inputs to the algorithm
       */
      //virtual void EvaluateImpl(muq::Modeling::ref_vector<boost::any> const& inputs) override;

      std::shared_ptr<SampleCollection> samples;

      std::shared_ptr<SampleCollection> QOIs;

#if MUQ_HAS_PARCER
      std::shared_ptr<parcer::Communicator> comm;// = std::make_shared<parcer::Communicator>();
#endif

    private:
      
    };
  } // namespace SamplingAlgorithms
} // namespace muq

#endif
