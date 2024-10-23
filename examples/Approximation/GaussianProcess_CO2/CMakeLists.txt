cmake_minimum_required (VERSION 3.10)

project(GaussianProcess_CO2)

set (CMAKE_CXX_STANDARD 17)

find_package(MUQ REQUIRED)

add_executable(GaussianProcess_CO2_exe GaussianProcess_CO2.cpp)
target_link_libraries(GaussianProcess_CO2_exe muq::muqApproximation)
