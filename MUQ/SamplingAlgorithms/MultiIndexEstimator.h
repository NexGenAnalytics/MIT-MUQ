#ifndef MULTIINDEXESTIMATOR_H
#define MULTIINDEXESTIMATOR_H

#include "MUQ/SamplingAlgorithms/SampleEstimator.h"
#include "MUQ/SamplingAlgorithms/MIMCMCBox.h"

namespace muq{
namespace SamplingAlgorithms{

    /** @class MultiIndexEstimator
        @ingroup MCMC
        @brief Class for estimating expectations using multi-index samples from MC or MCMC.
    */
    class MultiIndexEstimator : public SampleEstimator
    {
    public:
        MultiIndexEstimator(std::vector<std::shared_ptr<MIMCMCBox>> const& boxesIn);

        virtual ~MultiIndexEstimator() = default;

        virtual unsigned int BlockSize(int blockInd) const override;

        virtual unsigned int NumBlocks() const override;

        virtual Eigen::VectorXd ExpectedValue(std::shared_ptr<muq::Modeling::ModPiece> const& f,
                                              std::vector<std::string> const& metains = std::vector<std::string>()) const override;

        virtual Eigen::MatrixXd Covariance(Eigen::VectorXd const& mean, 
                                           int                    blockInd=-1) const override;

    private:
        const Eigen::VectorXi blockSizes;

        std::vector<std::shared_ptr<MIMCMCBox>> boxes;
    };
}
}



#endif 