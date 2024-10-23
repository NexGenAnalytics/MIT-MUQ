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

    # in some cases there might be that the executable is not present, 
    # for example for MPI tests like ParallelMultilevelMonteCarlo
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










# summary=""

# directories=(
# # "Approximation/GaussianProcess_CO2"
# # "Approximation/MonotoneRegression"

# # "Modeling/FlowEquation/cpp"
# # "Modeling/UMBridge"
# # "Modeling/CustomModPiece/cpp"
# # "Modeling/MemoryTest/cpp"
# # "Modeling/WorkGraphs/cpp"

# # "SamplingAlgorithms/MC/Example1_Gaussian"

# # "SamplingAlgorithms/MCMC/BasicMetropolisInGibbs/cpp"
# # "SamplingAlgorithms/MCMC/CustomGaussianProposal/cpp"
# # "SamplingAlgorithms/MCMC/EllipticInference/cpp"
# # "SamplingAlgorithms/MCMC/MultilevelMCMC_FlowModel/cpp"
# # "SamplingAlgorithms/MCMC/Example1_Gaussian/cpp"
# # "SamplingAlgorithms/MCMC/Example2_GaussianInverseGamma/cpp"
# # "SamplingAlgorithms/MCMC/Example3_MultilevelGaussian/cpp"
# # "SamplingAlgorithms/MCMC/Example4_MultiindexGaussian/cpp"
# # "SamplingAlgorithms/MCMC/Example5_MALASampling"

# "Utilities/HDF5/BlockOperations"
# "Utilities/HDF5/SimpleReadWrite"
# )

# binaries=(
#     # "Modeling/FlowEquation/cpp/build/FlowEquation"
#     # "Modeling/CustomModPiece/cpp/build/BasicModPiece"
#     # "Modeling/MemoryTest/cpp/build/MemoryUsage"
#     # "Modeling/WorkGraphs/cpp/build/SimpleWorkGraph"
#     # "Modeling/WorkGraphs/cpp/build/SplitSum"

#     # "SamplingAlgorithms/MC/Example1_Gaussian/build/MonteCarlo"
#     # "SamplingAlgorithms/MC/Example1_Gaussian/build/MultilevelMonteCarlo"

#     # "SamplingAlgorithms/MCMC/BasicMetropolisInGibbs/cpp/build/MetropolisInGibbs"
#     # "SamplingAlgorithms/MCMC/CustomGaussianProposal/cpp/build/CustomProposal"
#     # "SamplingAlgorithms/MCMC/EllipticInference/cpp/build/InvariantSampling"
#     # "SamplingAlgorithms/MCMC/Example1_Gaussian/cpp/build/GaussianSampling"
#     # "SamplingAlgorithms/MCMC/Example2_GaussianInverseGamma/cpp/build/GaussianGammaSampling"
#     # "SamplingAlgorithms/MCMC/Example3_MultilevelGaussian/cpp/build/BasicMultilevel"
#     # "SamplingAlgorithms/MCMC/Example3_MultilevelGaussian/cpp/build/AdvancedMultilevel"
#     # "SamplingAlgorithms/MCMC/Example4_MultiindexGaussian/cpp/build/MultiindexGaussianSampling"
#     # "./SamplingAlgorithms/MCMC/Example5_MALASampling/build/malaSampling"

#     "Utilities/HDF5/BlockOperations/build/BlockOperations"
#     "Utilities/HDF5/SimpleReadWrite/build/SimpleReadWrite"
# )

# if [[ $flag == 1 ]]; then 

#     for dir in "${directories[@]}"
#     do
#         echo "======================"
#         echo "Building example in examples/$dir"
#         echo ""

#         mkdir -p examples/$dir/build
#         rm -rf examples/$dir/build/*
#         cmake -S examples/$dir -B examples/$dir/build
#         cmake_exit=$?
#         if [ $cmake_exit -eq 0 ]
#         then
#             echo "CMake successful for example examples/$dir"
#         else
#             echo "CMake failed for example examples/$dir"
#             summary="$summary\nCMake FAILED for example examples/$dir"
#         fi
#         # Simply add exit codes; will be zero if all tests successful
#         exit_code=$(($exit_code + $cmake_exit)) 

#         make -C examples/$dir/build
#         make_exit=$?
#         if [ $make_exit -eq 0 ]
#         then
#             echo "Make successful for example examples/$dir"
#             summary="$summary\nBuild successful for example examples/$dir"
#         else
#             echo "Make failed for example examples/$dir"
#             summary="$summary\nMake FAILED for example examples/$dir"
#         fi
#         # Simply add exit codes; will be zero if all tests successful
#         exit_code=$(($exit_code + $make_exit)) 

#     done
#     echo "======================"
# fi

# if [[ $flag == 2 ]]; then 

#     for bin in "${binaries[@]}"
#     do
#         echo "======================"
#         echo "Running example $bin"
#         echo ""

#         examples/$bin
#         run_exit=$?
#         if [ $run_exit -eq 0 ]
#         then
#             echo "Run successful for example $bin"
#             summary="$summary\nRun successful for example $bin"
#         else
#             echo "Run failed for example $bin"
#             summary="$summary\nRun FAILED for example $bin"
#         fi
#         exit_code=$(($exit_code + $run_exit)) # Simply add exit codes; will be zero if all tests successful
#     done
#     echo "======================"

# fi

# echo -e $summary

# echo ""

# if [ $exit_code -eq 0 ]
# then
#     echo "All examples successful."
# else
#     echo "Some examples FAILED!"
# fi

# exit $exit_code
