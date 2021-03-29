cmake_minimum_required (VERSION 2.8)

project(DockerModPiece)

if(CMAKE_SOURCE_DIR STREQUAL CMAKE_CURRENT_SOURCE_DIR) # Only look for MUQ when building example as standalone project
  find_package(MUQ REQUIRED)

  include_directories(${MUQ_INCLUDE_DIRS})

  set(CMAKE_CXX_FLAGS ${MUQ_CXX_FLAGS})
  set(CMAKE_CXX_COMPILER ${MUQ_CXX_COMPILER})
endif()


add_executable(DockerModPiece DockerModPiece.cpp)
target_link_libraries(DockerModPiece ${MUQ_LIBRARIES} ${MUQ_LINK_LIBRARIES})

add_executable(DockerModPieceClient DockerModPieceClient.cpp)
target_link_libraries(DockerModPieceClient ${MUQ_LIBRARIES} ${MUQ_LINK_LIBRARIES})


add_executable(HTTPModPieceServer HTTPModPieceServer.cpp)
target_link_libraries(HTTPModPieceServer ${MUQ_LIBRARIES} ${MUQ_LINK_LIBRARIES})

add_executable(HTTPModPieceClientTest HTTPModPieceClientTest.cpp)
target_link_libraries(HTTPModPieceClientTest ${MUQ_LIBRARIES} ${MUQ_LINK_LIBRARIES} ssl crypto)

add_executable(HTTPModPieceClientMCMC HTTPModPieceClientMCMC.cpp)
target_link_libraries(HTTPModPieceClientMCMC ${MUQ_LIBRARIES} ${MUQ_LINK_LIBRARIES} ssl crypto)
