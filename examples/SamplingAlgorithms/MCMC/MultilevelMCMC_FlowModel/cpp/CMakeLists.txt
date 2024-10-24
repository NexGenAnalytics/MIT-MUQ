cmake_minimum_required (VERSION 3.10)

project(MultilevelMCMC_FLowModel)

set (CMAKE_CXX_STANDARD 17)

find_package(MUQ REQUIRED)

set(SOURCE
    MultilevelMCMC_FLowModel.cpp
    FlowModelComponents.cpp
)

add_executable(MultilevelMCMC_FLowModel ${SOURCE})
target_link_libraries(MultilevelMCMC_FLowModel muq::muqSamplingAlgorithms)