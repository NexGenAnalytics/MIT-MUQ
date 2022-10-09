cmake_minimum_required (VERSION 3.10)

project(MultilevelGaussian)

set (CMAKE_CXX_STANDARD 17)

find_package(MUQ REQUIRED)
include_directories(${MUQ_INCLUDE_DIRS})
link_libraries(${MUQ_LIBRARIES} ${MUQ_LINK_LIBRARIES})

add_executable(BasicMultilevel MultilevelMCMC_BasicInterface.cpp)

add_executable(AdvancedMultilevel MultilevelMCMC_AdvancedInterface.cpp)


if (MUQ_HAS_PARCER)
  add_executable(ModelParallelMultilevelGaussianSampling ModelParallelMultilevelGaussianSampling.cpp)

  add_executable(FullParallelMultilevelGaussianSampling FullParallelMultilevelGaussianSampling.cpp)

  add_executable(SubsamplingTestMultilevelGaussianSampling SubsamplingTestMultilevelGaussianSampling.cpp)
endif()
