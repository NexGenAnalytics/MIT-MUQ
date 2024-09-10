cmake_minimum_required (VERSION 3.10)

project(UMBridge)

set (CMAKE_CXX_STANDARD 17)

find_package(MUQ REQUIRED)

add_executable(UMBridgeModelClient UMBridgeModelClient.cpp)
target_link_libraries(UMBridgeModelClient muq::dev::muqModeling)

add_executable(UMBridgeModelServer UMBridgeModelServer.cpp)
target_link_libraries(UMBridgeModelServer muq::dev::muqModeling)
