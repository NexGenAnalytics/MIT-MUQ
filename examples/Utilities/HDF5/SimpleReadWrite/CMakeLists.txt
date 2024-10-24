cmake_minimum_required(VERSION 3.10)

project(SimpleReadWrite)

set(CMAKE_CXX_STANDARD 17)

find_package(MUQ REQUIRED)

add_executable(SimpleReadWrite SimpleReadWrite.cpp)
target_link_libraries(SimpleReadWrite muq::muqUtilities)