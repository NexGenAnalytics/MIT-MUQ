cmake_minimum_required (VERSION 3.10)

project(CustomModPiece)

set (CMAKE_CXX_STANDARD 17)

find_package(MUQ REQUIRED)

add_executable(BasicModPiece BasicModPiece.cpp)

target_link_libraries(BasicModPiece muq::muqModeling)
