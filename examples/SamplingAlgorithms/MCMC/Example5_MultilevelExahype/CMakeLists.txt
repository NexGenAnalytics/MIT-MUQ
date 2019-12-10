cmake_minimum_required (VERSION 2.8)

find_package(MUQ_PRIVATE REQUIRED)

include_directories(${MUQ_INCLUDE_DIRS})

set(CMAKE_CXX_FLAGS ${MUQ_CXX_FLAGS})
set(CMAKE_CXX_COMPILER ${MUQ_CXX_COMPILER})

message("${MUQ_INCLUDE_DIRS}")
message("${MUQ_LIBRARIES}")

file(READ "${EXAHYPE_PATH}/exahype_link" EXAHYPE_LINKS)
string(STRIP ${EXAHYPE_LINKS} EXAHYPE_LINKS)

add_executable(MultilevelExahype MultilevelExahype.cpp)
target_link_libraries(MultilevelExahype ${MUQ_LIBRARIES} ${MUQ_LINK_LIBRARIES})
target_include_directories(MultilevelExahype PRIVATE ${EXAHYPE_PATH})
set_target_properties(MultilevelExahype PROPERTIES LINK_FLAGS ${EXAHYPE_LINKS})
