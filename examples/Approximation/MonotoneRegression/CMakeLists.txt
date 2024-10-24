cmake_minimum_required (VERSION 3.10)

project(MonotoneRegression)

set (CMAKE_CXX_STANDARD 17)

find_package(MUQ REQUIRED)

add_executable(MonotoneRegression MonotoneRegression.cpp)
target_link_libraries(MonotoneRegression muq::muqApproximation)