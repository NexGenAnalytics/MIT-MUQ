cmake_minimum_required (VERSION 3.10)

project(Gaussian)

set (CMAKE_CXX_STANDARD 17)

find_package(MUQ REQUIRED)

add_executable(GaussianSampling GaussianSampling.cpp)
target_link_libraries(GaussianSampling muq::dev::muqSamplingAlgorithms)
