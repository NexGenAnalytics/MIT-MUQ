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

        /** Construct the multiindex estimator using MIMCMC boxes.  These boxes are typically constructed by 
            a MIMCMC methods such as the GreedyMLMCMC or MIMCMC classes.

            @param[in] boxesIn "Boxes" holding the differences between chains at different indices
            @param[in] useQoisIn (optional) Whether this estimator should use the QOIs in the chains or
                               the parameters themselves.  Defaults to false, which implies the parameters
                               will be used in the estimates.
         */ 
        MultiIndexEstimator(std::vector<std::shared_ptr<MIMCMCBox>> const& boxesIn, 
                            bool                                           useQoisIn = false);

        virtual ~MultiIndexEstimator() = default;

        virtual unsigned int BlockSize(int blockInd) const override;

        virtual unsigned int NumBlocks() const override;

        virtual Eigen::VectorXd ExpectedValue(std::shared_ptr<muq::Modeling::ModPiece> const& f,
                                              std::vector<std::string> const& metains = std::vector<std::string>()) const override;

        virtual Eigen::MatrixXd Covariance(int blockInd=-1) const override { return SampleEstimator::Covariance(blockInd);};
        virtual Eigen::MatrixXd Covariance(Eigen::VectorXd const& mean, 
                                           int                    blockInd=-1) const override;

    private:
        const Eigen::VectorXi blockSizes;
        const bool useQois;

        std::vector<std::shared_ptr<MIMCMCBox>> boxes;
        
    };
}
}



#endif 