cmake_minimum_required (VERSION 3.10)

project(SimpleWorkGraph)

set (CMAKE_CXX_STANDARD 17)

find_package(MUQ REQUIRED)

add_executable(SimpleWorkGraph SimpleWorkGraph.cpp)
target_link_libraries(SimpleWorkGraph muq::muqSamplingAlgorithms)

add_executable(SplitSum SplitSumWorkGraph.cpp)
target_link_libraries(SplitSum muq::muqSamplingAlgorithms)