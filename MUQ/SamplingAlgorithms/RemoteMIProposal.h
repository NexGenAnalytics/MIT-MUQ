#ifndef REMOTEMIPROPOSAL_H_
#define REMOTEMIPROPOSAL_H_

#if MUQ_HAS_MPI

#if !MUQ_HAS_PARCER
#error
#endif


#include <boost/property_tree/ptree.hpp>
#include "MUQ/SamplingAlgorithms/Phonebook.h"

#include <boost/property_tree/ptree.hpp>
namespace pt = boost::property_tree;

namespace muq {
  namespace SamplingAlgorithms {


		class RemoteMIProposal : public MCMCProposal {
		public:
			RemoteMIProposal (pt::ptree const& pt, std::shared_ptr<AbstractSamplingProblem> prob, std::shared_ptr<parcer::Communicator> comm, std::shared_ptr<MultiIndex> remoteIndex, std::shared_ptr<PhonebookClient> phonebookClient)
				: MCMCProposal(pt,prob),
				subsampling(pt.get("subsampling",1)),
				comm(comm),
		    remoteIndex(remoteIndex),
				phonebookClient(phonebookClient)
			{
			}

			std::shared_ptr<SamplingState> Sample(std::shared_ptr<SamplingState> currentState) {

				int remoteRank = phonebookClient->Query(remoteIndex);

				comm->Send(ControlFlag::SAMPLE, remoteRank, ControlTag);
				Eigen::VectorXd remoteState = comm->Recv<Eigen::VectorXd>(remoteRank, ControlTag);
				Eigen::VectorXd remoteQOI = comm->Recv<Eigen::VectorXd>(remoteRank, ControlTag);

				auto proposal = std::make_shared<SamplingState>(remoteState); // TODO: Support non-QOI samples!!
				proposal->meta["QOI"] = std::make_shared<SamplingState>(remoteQOI);

				return proposal;
			}

			double LogDensity(std::shared_ptr<SamplingState> currState,
			                  std::shared_ptr<SamplingState> propState) {
				return 0;
			}

		private:
			int sampleID = 0;
			int sampleWeight = 0;
			const int subsampling;
			std::shared_ptr<parcer::Communicator> comm;
		  std::shared_ptr<MultiIndex> remoteIndex;
			std::shared_ptr<PhonebookClient> phonebookClient;
		};
	}
}

#endif

#endif
