cmake_minimum_required (VERSION 3.10)

project(MCGaussian)

find_package(MUQ REQUIRED)
include_directories(${MUQ_INCLUDE_DIRS})
link_libraries(${MUQ_LIBRARIES} ${MUQ_LINK_LIBRARIES})

add_executable(MonteCarlo MonteCarlo.cpp)

add_executable(MultilevelMonteCarlo MultilevelMonteCarlo.cpp)

if (MUQ_HAS_PARCER)
  add_executable(ParallelMultilevelMonteCarlo ParallelMultilevelMonteCarlo.cpp)
endif()
