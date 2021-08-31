#ifndef GreedyMLMCMC_H
#define GreedyMLMCMC_H

#include <boost/property_tree/ptree.hpp>

#include "MUQ/SamplingAlgorithms/MIMCMCBox.h"
#include "MUQ/SamplingAlgorithms/MIComponentFactory.h"
#include "MUQ/SamplingAlgorithms/SamplingAlgorithm.h"
#include "MUQ/SamplingAlgorithms/MultiIndexEstimator.h"

namespace pt = boost::property_tree;

namespace muq {
  namespace SamplingAlgorithms {

    /** @brief Greedy Multilevel MCMC method.
        @details A Multilevel MCMC method choosing
        the number of samples adaptively at runtime,
        estimating the most profitable level from
        statistical information on samples.
     */
    class GreedyMLMCMC {
    public:
      GreedyMLMCMC (pt::ptree pt, std::shared_ptr<MIComponentFactory> componentFactory);

      virtual std::shared_ptr<MultiIndexEstimator> GetSamples() const;
      virtual std::shared_ptr<MultiIndexEstimator> GetQOIs() const;

      Eigen::VectorXd MeanQOI();

      void Draw(bool drawSamples = true);

      std::shared_ptr<MIMCMCBox> GetBox(int index);
      std::vector<std::shared_ptr<MIMCMCBox>> GetBoxes();

      void WriteToFile(std::string filename);

      virtual std::shared_ptr<MultiIndexEstimator> Run();

    protected:
      
    private:
      std::shared_ptr<MIComponentFactory> componentFactory;
      const int numInitialSamples;
      const double e;
      const double beta;
      const int levels;
      int verbosity;
      std::vector<std::shared_ptr<MIMCMCBox>> boxes;
    };

  }
}

#endif
