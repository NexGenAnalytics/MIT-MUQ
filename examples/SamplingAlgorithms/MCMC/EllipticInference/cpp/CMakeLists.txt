cmake_minimum_required (VERSION 3.10)

project(EllipticInference)

set (CMAKE_CXX_STANDARD 17)

find_package(MUQ REQUIRED)

add_executable(InvariantSampling InvariantSampling.cpp)
target_link_libraries(InvariantSampling muq::dev::muqSamplingAlgorithms)
