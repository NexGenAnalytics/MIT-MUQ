cmake_minimum_required (VERSION 3.10)

project(MALASampler)

find_package(MUQ REQUIRED)
include_directories(${MUQ_INCLUDE_DIRS})
link_libraries(${MUQ_LIBRARIES} ${MUQ_LINK_LIBRARIES})

add_executable(malaSampling MALASampler.cpp)
