cmake_minimum_required (VERSION 3.10)

project(MALASampler)

find_package(MUQ REQUIRED)

add_executable(malaSampling MALASampler.cpp)
target_link_libraries(malaSampling muq::dev::muqSamplingAlgorithms)
