add_subdirectory(Approximation/GaussianProcess_CO2)
add_subdirectory(Approximation/MonotoneRegression)

add_subdirectory(Modeling/DockerModPiece/cpp)

add_subdirectory(SamplingAlgorithms/MC/Example1_Gaussian)

add_subdirectory(SamplingAlgorithms/MCMC/Example1_Gaussian/cpp)
add_subdirectory(SamplingAlgorithms/MCMC/Example2_GaussianInverseGamma)
add_subdirectory(SamplingAlgorithms/MCMC/Example3_MultilevelGaussian/cpp)
add_subdirectory(SamplingAlgorithms/MCMC/Example4_MultiindexGaussian/cpp)

add_subdirectory(Utilities/HDF5/BlockOperations)
add_subdirectory(Utilities/HDF5/SimpleReadWrite)

add_custom_target(examples
                  DEPENDS GaussianProcess_CO2_exe
                          MonotoneRegression

                          MonteCarlo
                          MultilevelMonteCarlo

                          GaussianSampling
                          GaussianGammaSampling
                          MultilevelGaussianSampling
                          MultiindexGaussianSampling

                          BlockOperations
                          SimpleReadWrite

                          DockerModPiece
                          DockerModPieceClient
                          )
if (MUQ_HAS_PARCER)
  add_dependencies(examples
                   ParallelMultilevelMonteCarlo

                   ModelParallelMultilevelGaussianSampling
                   FullParallelMultilevelGaussianSampling
                   SubsamplingTestMultilevelGaussianSampling
                   FullParallelMultiindexGaussianSampling)
endif()
