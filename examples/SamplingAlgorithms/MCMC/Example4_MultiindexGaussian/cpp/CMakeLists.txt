cmake_minimum_required (VERSION 3.10)

project(MultiindexGaussian)

set (CMAKE_CXX_STANDARD 17)

find_package(MUQ REQUIRED)

add_executable(MultiindexGaussianSampling MultiindexGaussianSampling.cpp)
target_link_libraries(MultiindexGaussianSampling muq::dev::muqSamplingAlgorithms)

if (MUQ_HAS_PARCER)
  add_executable(FullParallelMultiindexGaussianSampling FullParallelMultiindexGaussianSampling.cpp)
    target_link_libraries(FullParallelMultiindexGaussianSampling muq::dev::muqSamplingAlgorithms)
endif()
