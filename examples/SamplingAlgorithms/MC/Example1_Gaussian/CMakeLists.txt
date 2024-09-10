cmake_minimum_required (VERSION 3.10)

project(MCGaussian)

find_package(MUQ REQUIRED)

add_executable(MonteCarlo MonteCarlo.cpp)
target_link_libraries(MonteCarlo muq::dev::muqSamplingAlgorithms)

add_executable(MultilevelMonteCarlo MultilevelMonteCarlo.cpp)
target_link_libraries(MultilevelMonteCarlo muq::dev::muqSamplingAlgorithms)

if (MUQ_HAS_PARCER)
  add_executable(ParallelMultilevelMonteCarlo ParallelMultilevelMonteCarlo.cpp)
  target_link_libraries(ParallelMultilevelMonteCarlo muq::dev::muqSamplingAlgorithms)
endif()
