cmake_minimum_required(VERSION 3.10)

project(GaussianInverseGamma)
set(CMAKE_CXX_STANDARD 17)

find_package(MUQ REQUIRED)

add_executable(GaussianGammaSampling GaussianGammaSampling.cpp)
target_link_libraries(GaussianGammaSampling muq::muqSamplingAlgorithms)
