cmake_minimum_required (VERSION 3.10)

project(MultilevelGaussian)

set (CMAKE_CXX_STANDARD 17)

find_package(MUQ REQUIRED)

add_executable(BasicMultilevel MultilevelMCMC_BasicInterface.cpp)
target_link_libraries(BasicMultilevel muq::dev::muqSamplingAlgorithms)

add_executable(AdvancedMultilevel MultilevelMCMC_AdvancedInterface.cpp)
target_link_libraries(AdvancedMultilevel muq::dev::muqSamplingAlgorithms)

if (MUQ_HAS_PARCER)
  add_executable(ModelParallelMultilevelGaussianSampling ModelParallelMultilevelGaussianSampling.cpp)
  target_link_libraries(ModelParallelMultilevelGaussianSampling muq::dev::muqSamplingAlgorithms)

  add_executable(FullParallelMultilevelGaussianSampling FullParallelMultilevelGaussianSampling.cpp)
  target_link_libraries(FullParallelMultilevelGaussianSampling muq::dev::muqSamplingAlgorithms)

  add_executable(SubsamplingTestMultilevelGaussianSampling SubsamplingTestMultilevelGaussianSampling.cpp)
  target_link_libraries(SubsamplingTestMultilevelGaussianSampling muq::dev::muqSamplingAlgorithms)
endif()
