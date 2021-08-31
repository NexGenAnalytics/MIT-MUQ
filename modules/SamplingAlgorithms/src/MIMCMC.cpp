#include "MUQ/SamplingAlgorithms/MIMCMC.h"
#include "spdlog/spdlog.h"

namespace muq {
  namespace SamplingAlgorithms {

    MIMCMC::MIMCMC (pt::ptree pt, std::shared_ptr<MIComponentFactory> componentFactory)
    : pt(pt),
      componentFactory(componentFactory)
    {
      gridIndices = MultiIndexFactory::CreateFullTensor(componentFactory->FinestIndex()->GetVector());

      for (int i = 0; i < gridIndices->Size(); i++) {
        std::shared_ptr<MultiIndex> boxHighestIndex = (*gridIndices)[i];
        auto box = std::make_shared<MIMCMCBox>(componentFactory, boxHighestIndex);
        boxes.push_back(box);
      }
    }

    std::shared_ptr<MIMCMCBox> MIMCMC::GetBox(std::shared_ptr<MultiIndex> index) {
      for (std::shared_ptr<MIMCMCBox> box : boxes) {
        if (*(box->GetHighestIndex()) == *index)
          return box;
      }
      return nullptr;
    }

    std::shared_ptr<MultiIndexEstimator> MIMCMC::GetSamples() const {
      return std::make_shared<MultiIndexEstimator>(boxes);
    }
    std::shared_ptr<MultiIndexEstimator> MIMCMC::GetQOIs() const {
      return std::make_shared<MultiIndexEstimator>(boxes,true);
    }

    std::shared_ptr<MultiIndexEstimator> MIMCMC::Run() {
      for (auto box : boxes) {
        assert(box);
        int numSamples = pt.get<int>("NumSamples" + multiindexToConfigString(box->GetHighestIndex()));
        for (int samp = 0; samp < numSamples; samp++) {
          box->Sample();
        }
      }

      return GetSamples();
    }

    std::shared_ptr<MIMCMCBox> MIMCMC::GetMIMCMCBox(std::shared_ptr<MultiIndex> index) {
      for (auto box : boxes) {
        if (*(box->GetHighestIndex()) == *index)
          return box;
      }
      return nullptr;
    }

    void MIMCMC::WriteToFile(std::string filename) {
      for (auto box : boxes) {
        box->WriteToFile(filename);
      }
    }


    std::shared_ptr<MultiIndexSet> MIMCMC::GetIndices() {
      return gridIndices;
    }


    std::string MIMCMC::multiindexToConfigString (std::shared_ptr<MultiIndex> index) {
      std::stringstream strs;
      for (int i = 0; i < index->GetLength(); i++) {
        strs << "_" << index->GetValue(i);
      }
      return strs.str();
    }

    void MIMCMC::Draw(bool drawSamples) {
      std::ofstream graphfile;
      graphfile.open ("graph");
      graphfile << "digraph {" << std::endl;
      graphfile << "nodesep=1.2;" << std::endl;
      graphfile << "splines=false;" << std::endl;
      for (auto box : boxes) {
        box->Draw(graphfile, drawSamples);
      }
      graphfile << "}" << std::endl;
      graphfile.close();
    }

  }
}
