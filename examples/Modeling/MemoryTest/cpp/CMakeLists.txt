cmake_minimum_required (VERSION 3.10)

project(MemoryTest)

set (CMAKE_CXX_STANDARD 17)

find_package(MUQ REQUIRED)

add_executable(MemoryUsage MemoryUsage.cpp)

target_link_libraries(MemoryUsage muq::dev::muqModeling)
