cmake_minimum_required (VERSION 3.10)

project(Gaussian)

set(CMAKE_CXX_STANDARD 17)

find_package(MUQ REQUIRED)

add_executable(CustomProposal CustomProposal.cpp)
target_link_libraries(CustomProposal muq::muqSamplingAlgorithms)