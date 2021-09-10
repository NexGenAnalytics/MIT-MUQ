#include "MUQ/SamplingAlgorithms/MultiIndexEstimator.h"
#include "MUQ/SamplingAlgorithms/MarkovChain.h"

using namespace muq::SamplingAlgorithms;


MultiIndexEstimator::MultiIndexEstimator(std::vector<std::shared_ptr<MIMCMCBox>> const& boxesIn,
                                         bool                                           useQoisIn) : boxes(boxesIn),
                                                                                                     blockSizes(useQoisIn ? boxesIn.at(0)->GetFinestProblem()->blockSizesQOI : boxesIn.at(0)->GetFinestProblem()->blockSizes),
                                                                                                     useQois(useQoisIn)
{

}

unsigned int MultiIndexEstimator::BlockSize(int blockInd) const
{   
    if(blockInd<0){
        return blockSizes.sum();
    }else{
        return blockSizes(blockInd);
    }
}


unsigned int MultiIndexEstimator::NumBlocks() const
{
    return blockSizes.size();
}

Eigen::VectorXd MultiIndexEstimator::ExpectedValue(std::shared_ptr<muq::Modeling::ModPiece> const& f,
                                                   std::vector<std::string>                 const& metains) const
{
    assert(f->outputSizes.size()==1);

    Eigen::VectorXd telescopingSum = Eigen::VectorXd::Zero(f->outputSizes(0));
   
    // Add up the telescoping series of MI boxes
    for (auto& box : boxes) {

        // Compute the expected difference for one term in the telescoping series
        Eigen::VectorXd diffMean = Eigen::VectorXd::Zero(f->outputSizes(0));

        auto boxIndices = box->GetBoxIndices();
        for (int i = 0; i < boxIndices->Size(); i++) {

            std::shared_ptr<MultiIndex> boxIndex = (*boxIndices)[i];
            auto chain = box->GetChain(boxIndex);
            std::shared_ptr<MarkovChain> samps;
            if(useQois){
                samps = chain->GetQOIs();
            }else{
                samps = chain->GetSamples();
            } 

            MultiIndex index = *(box->GetLowestIndex()) + *boxIndex;
            MultiIndex indexDiffFromTop = *(box->GetHighestIndex()) - index;

            if (indexDiffFromTop.Sum() % 2 == 0) {
                diffMean += samps->ExpectedValue(f,metains);
            } else {
                diffMean -= samps->ExpectedValue(f,metains);
            }
        }

        // Add this term of the series to the running total      
        telescopingSum += diffMean;
    }

    return telescopingSum;
} 

Eigen::MatrixXd MultiIndexEstimator::Covariance(Eigen::VectorXd const& mean, 
                                                int                    blockInd) const
{
    Eigen::MatrixXd telescopingSum = Eigen::MatrixXd::Zero(BlockSize(blockInd),BlockSize(blockInd));

    // Add up the telescoping series of MI boxes
    for (auto& box : boxes) {

        // Compute the expected difference for one term in the telescoping series
        Eigen::MatrixXd diffMean = Eigen::MatrixXd::Zero(BlockSize(blockInd),BlockSize(blockInd));

        auto boxIndices = box->GetBoxIndices();
        for (int i = 0; i < boxIndices->Size(); i++) {

            std::shared_ptr<MultiIndex> boxIndex = (*boxIndices)[i];
            auto chain = box->GetChain(boxIndex);
            std::shared_ptr<MarkovChain> samps;
            if(useQois){
                samps = chain->GetQOIs();
            }else{
                samps = chain->GetSamples();
            } 

            MultiIndex index = *(box->GetLowestIndex()) + *boxIndex;
            MultiIndex indexDiffFromTop = *(box->GetHighestIndex()) - index;

            if (indexDiffFromTop.Sum() % 2 == 0) {
                diffMean += samps->Covariance(mean);
            } else {
                diffMean -= samps->Covariance(mean);
            }
        }

        // Add this term of the series to the running total      
        telescopingSum += diffMean;
    }

    return 0.5*(telescopingSum + telescopingSum.transpose());
}