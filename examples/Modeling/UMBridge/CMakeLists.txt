cmake_minimum_required (VERSION 3.10)

project(UMBridge)

set (CMAKE_CXX_STANDARD 17)

find_package(MUQ REQUIRED)
include_directories(${MUQ_INCLUDE_DIRS})
link_libraries(${MUQ_LIBRARIES} ${MUQ_LINK_LIBRARIES})

add_executable(UMBridgeModelClient UMBridgeModelClient.cpp)

add_executable(UMBridgeModelServer UMBridgeModelServer.cpp)
