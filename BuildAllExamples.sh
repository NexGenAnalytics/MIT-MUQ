#!/bin/bash

exit_code=0

summary=""

directories=(
"Approximation/GaussianProcess_CO2"
"Approximation/MonotoneRegression"

"Modeling/FlowEquation/cpp"
"Modeling/UMBridge"
"Modeling/CustomModPiece/cpp"
"Modeling/MemoryTest/cpp"
"Modeling/WorkGraphs/cpp"

"SamplingAlgorithms/MC/Example1_Gaussian"

"SamplingAlgorithms/MCMC/BasicMetropolisInGibbs/cpp"
"SamplingAlgorithms/MCMC/CustomGaussianProposal/cpp"
"SamplingAlgorithms/MCMC/EllipticInference/cpp"
"SamplingAlgorithms/MCMC/MultilevelMCMC_FlowModel/cpp"
"SamplingAlgorithms/MCMC/Example1_Gaussian/cpp"
"SamplingAlgorithms/MCMC/Example2_GaussianInverseGamma/cpp"
"SamplingAlgorithms/MCMC/Example3_MultilevelGaussian/cpp"
"SamplingAlgorithms/MCMC/Example4_MultiindexGaussian/cpp"
"SamplingAlgorithms/MCMC/Example5_MALASampling"

"Utilities/HDF5/BlockOperations"
"Utilities/HDF5/SimpleReadWrite"
)

binaries=(
    "Modeling/FlowEquation/cpp/build/FlowEquation"
    "Modeling/CustomModPiece/cpp/build/BasicModPiece"
    "Modeling/MemoryTest/cpp/build/MemoryUsage"
    "Modeling/WorkGraphs/cpp/build/SimpleWorkGraph"
    "Modeling/WorkGraphs/cpp/build/SplitSum"

    "SamplingAlgorithms/MC/Example1_Gaussian/build/MonteCarlo"
    "SamplingAlgorithms/MC/Example1_Gaussian/build/MultilevelMonteCarlo"

    "SamplingAlgorithms/MCMC/BasicMetropolisInGibbs/cpp/build/MetropolisInGibbs"
    "SamplingAlgorithms/MCMC/CustomGaussianProposal/cpp/build/CustomProposal"
    "SamplingAlgorithms/MCMC/EllipticInference/cpp/build/InvariantSampling"
    "SamplingAlgorithms/MCMC/Example1_Gaussian/cpp/build/GaussianSampling"
    "SamplingAlgorithms/MCMC/Example2_GaussianInverseGamma/cpp/build/GaussianGammaSampling"
    "SamplingAlgorithms/MCMC/Example3_MultilevelGaussian/cpp/build/BasicMultilevel"
    "SamplingAlgorithms/MCMC/Example3_MultilevelGaussian/cpp/build/AdvancedMultilevel"
    "SamplingAlgorithms/MCMC/Example4_MultiindexGaussian/cpp/build/MultiindexGaussianSampling"
    "./SamplingAlgorithms/MCMC/Example5_MALASampling/build/malaSampling"

    "Utilities/HDF5/BlockOperations/build/BlockOperations"
    "Utilities/HDF5/SimpleReadWrite/build/SimpleReadWrite"
)

for dir in "${directories[@]}"
do
    echo "======================"
    echo "Building example in examples/$dir"
    echo ""

    mkdir -p examples/$dir/build
    cmake -S examples/$dir -B examples/$dir/build -DMUQ_DIR=$PWD/build/install/lib/cmake/MUQ
    cmake_exit=$?
    if [ $cmake_exit -eq 0 ]
    then
        echo "CMake successful for example examples/$dir"
    else
        echo "CMake failed for example examples/$dir"
        summary="$summary\nCMake FAILED for example examples/$dir"
    fi
    exit_code=$(($exit_code + $cmake_exit)) # Simply add exit codes; will be zero if all tests successful

    make -C examples/$dir/build
    make_exit=$?
    if [ $make_exit -eq 0 ]
    then
        echo "Make successful for example examples/$dir"
        summary="$summary\nBuild successful for example examples/$dir"
    else
        echo "Make failed for example examples/$dir"
        summary="$summary\nMake FAILED for example examples/$dir"
    fi
    exit_code=$(($exit_code + $make_exit)) # Simply add exit codes; will be zero if all tests successful
done

echo "======================"

for bin in "${binaries[@]}"
do
    echo "======================"
    echo "Running example $bin"
    echo ""

    examples/$bin
    run_exit=$?
    if [ $run_exit -eq 0 ]
    then
        echo "Run successful for example $bin"
        summary="$summary\nRun successful for example $bin"
    else
        echo "Run failed for example $bin"
        summary="$summary\nRun FAILED for example $bin"
    fi
    exit_code=$(($exit_code + $run_exit)) # Simply add exit codes; will be zero if all tests successful
done

echo "======================"

echo -e $summary

echo ""

if [ $exit_code -eq 0 ]
then
    echo "All examples successful."
else
    echo "Some examples FAILED!"
fi

exit $exit_code
