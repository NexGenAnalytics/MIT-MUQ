cmake_minimum_required (VERSION 3.10)

project(MultilevelMCMC_FLowModel)

set (CMAKE_CXX_STANDARD 17)

find_package(MUQ REQUIRED)
include_directories(${MUQ_INCLUDE_DIRS})
link_libraries(${MUQ_LIBRARIES} ${MUQ_LINK_LIBRARIES})

set(SOURCE
    MultilevelMCMC_FLowModel.cpp
    FlowModelComponents.cpp
)

add_executable(MultilevelMCMC_FLowModel ${SOURCE})
