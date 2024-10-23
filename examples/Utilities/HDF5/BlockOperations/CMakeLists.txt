cmake_minimum_required(VERSION 3.10)

project(BlockOperations)

find_package(MUQ REQUIRED)

add_executable(BlockOperations BlockOperations.cpp)
target_link_libraries(BlockOperations muq::muqUtilities)