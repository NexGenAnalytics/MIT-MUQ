cmake_minimum_required (VERSION 3.10)

project(SimpleWorkGraph)

set (CMAKE_CXX_STANDARD 17)

find_package(MUQ REQUIRED)
include_directories(${MUQ_INCLUDE_DIRS})
link_libraries(${MUQ_LIBRARIES} ${MUQ_LINK_LIBRARIES})

add_executable(SimpleWorkGraph SimpleWorkGraph.cpp)

add_executable(SplitSum SplitSumWorkGraph.cpp)
