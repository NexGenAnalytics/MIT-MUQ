CreateTestGroup(
    SAMPLING_ALGORITHM

    SampleCollectionTests.cpp
    MarkovChainTests.cpp
    MIMCMCTests.cpp
    MLMCMCTests.cpp

    GMHKernelTests.cpp

    #MonteCarloTests.cpp
    ImportanceSamplingTests.cpp
    MCMCTests.cpp
    DelayedRejectionTests.cpp
    DILITests.cpp
    CrankNicolsonTests.cpp
    InfMALATests.cpp
    InverseGammaProposalTests.cpp
)

CreateTestGroup(
  SAMPLING_PROBLEMS

  SamplingProblemTests.cpp
)

CreateTestGroup(
  FANCY_SAMPLING_PROBLEMS

  ExpensiveSamplingProblemTests.cpp
)

CreateTestGroup(
  SAMPLING_ALGORITHMS_GRAPH_BASED_ALGORITHMS

  SampleGraphsTests/SampleGraphTests.cpp
  SampleGraphsTests/DensityEstimationTests.cpp
  SampleGraphsTests/KolmogorovOperatorTests.cpp
)

CreateParallelTestGroup(
    SAMPLING_ALGORITHM

    SamplingStateParallelTests.cpp

    DistributedCollectionTests.cpp

    GMHKernelParallelTests.cpp

    SingleChainMCMCParallelTests.cpp

    ExpensiveSamplingProblemParallelTests.cpp

    ImportanceSamplingParallelTests.cpp
)
