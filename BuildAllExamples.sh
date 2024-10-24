#!/bin/bash

exit_code=0

if [ -z "${MUQ_DIR}" ]; then
    echo "MUQ_DIR is unset or set to the empty string."
    echo "You need to set it to point to a valid installation path before running this script."
    exit 21
fi

ARRAY=( 
    "Utilities/HDF5/BlockOperations:BlockOperations"
    "Utilities/HDF5/SimpleReadWrite:SimpleReadWrite"

    "Approximation/GaussianProcess_CO2:GaussianProcess_CO2_exe"
    "Approximation/MonotoneRegression:MonotoneRegression"

    "Modeling/MemoryTest/cpp:MemoryUsage"
    "Modeling/FlowEquation/cpp:FlowEquation"
    "Modeling/CustomModPiece/cpp:BasicModPiece"
    "Modeling/WorkGraphs/cpp:SimpleWorkGraph"
    "Modeling/WorkGraphs/cpp:SplitSum"
    "Modeling/UMBridge:void" # for this only build as it was already in master branch

    "SamplingAlgorithms/MC/Example1_Gaussian:MonteCarlo"
    "SamplingAlgorithms/MC/Example1_Gaussian:MultilevelMonteCarlo"
    "SamplingAlgorithms/MC/Example1_Gaussian:ParallelMultilevelMonteCarlo"

    "SamplingAlgorithms/MCMC/BasicMetropolisInGibbs/cpp/:MetropolisInGibbs"
    "SamplingAlgorithms/MCMC/CustomGaussianProposal/cpp/:CustomProposal"
    "SamplingAlgorithms/MCMC/EllipticInference/cpp/:InvariantSampling"
    "SamplingAlgorithms/MCMC/Example1_Gaussian/cpp/:GaussianSampling"
    "SamplingAlgorithms/MCMC/Example2_GaussianInverseGamma/cpp/:GaussianGammaSampling"
    "SamplingAlgorithms/MCMC/Example3_MultilevelGaussian/cpp/:BasicMultilevel"
    "SamplingAlgorithms/MCMC/Example3_MultilevelGaussian/cpp/:AdvancedMultilevel"
    "SamplingAlgorithms/MCMC/Example4_MultiindexGaussian/cpp/:MultiindexGaussianSampling"
    "SamplingAlgorithms/MCMC/Example5_MALASampling:malaSampling"
)


for it in "${ARRAY[@]}" ; do
    exampleDir=${it%%:*}
    executableName=${it##*:}

    exampleSourceDir=$PWD/examples/${exampleDir}
    exampleBuildDir=${exampleSourceDir}/build

    echo "======================================"
    echo "Building ${exampleSourceDir}"
    echo ""

    if [ -d ${exampleBuildDir} ]; then
        rm -rf ${exampleBuildDir}/*
    else
        mkdir -p ${exampleBuildDir}
    fi

    cmake -S ${exampleSourceDir} -B ${exampleBuildDir}
    cmake_exit=$?
    if [ $cmake_exit -eq 0 ]
    then
        echo "CMake successful for example ${exampleSourceDir}"
    else
        echo "CMake failed for example ${exampleSourceDir}"
        summary="$summary\nCMake FAILED for example ${exampleSourceDir}"
    fi
    # Simply add exit codes; will be zero if all tests successful
    exit_code=$(($exit_code + $cmake_exit)) 

    make -C ${exampleBuildDir}
    make_exit=$?
    if [ $make_exit -eq 0 ]
    then
        echo "Make successful for ${exampleSourceDir}"
        summary="$summary\nBuild successful for ${exampleSourceDir}"
    else
        echo "Make failed for ${exampleSourceDir}"
        summary="$summary\nMake FAILED for ${exampleSourceDir}"
    fi
    # Simply add exit codes; will be zero if all tests successful
    exit_code=$(($exit_code + $make_exit)) 

    # in some cases there might be that the executable is not present
    if [ -f ${exampleBuildDir}/${executableName} ]; then
        if [ "$executableName" != "void" ]; then
            cd ${exampleBuildDir}
            ./${executableName}
            run_exit=$?
            exit_code=$(($exit_code + $run_exit)) 
            cd -
        fi
    fi
done

echo -e $summary
echo ""

if [ $exit_code -eq 0 ]
then
    echo "All examples successful."
else
    echo "Some examples FAILED!"
fi

exit $exit_code
