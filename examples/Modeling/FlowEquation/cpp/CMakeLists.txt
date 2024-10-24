cmake_minimum_required (VERSION 3.10)

project(FlowEquation)

set (CMAKE_CXX_STANDARD 17)

find_package(MUQ REQUIRED)

add_executable(FlowEquation FlowEquation.cpp)
target_link_libraries(FlowEquation muq::muqModeling)
