cmake_minimum_required (VERSION 2.8)

find_package(MUQ_PRIVATE REQUIRED)

include_directories(${MUQ_INCLUDE_DIRS})

set(CMAKE_CXX_FLAGS ${MUQ_CXX_FLAGS})
set(CMAKE_CXX_COMPILER ${MUQ_CXX_COMPILER})

message("${MUQ_INCLUDE_DIRS}")
message("${MUQ_LIBRARIES}")

add_executable(MultilevelGaussianSampling MultilevelGaussianSampling.cpp)
target_link_libraries(MultilevelGaussianSampling ${MUQ_LIBRARIES} ${MUQ_LINK_LIBRARIES})

add_executable(ModelParallelMultilevelGaussianSampling ModelParallelMultilevelGaussianSampling.cpp)
target_link_libraries(ModelParallelMultilevelGaussianSampling ${MUQ_LIBRARIES} ${MUQ_LINK_LIBRARIES})

add_executable(FullParallelMultilevelGaussianSampling FullParallelMultilevelGaussianSampling.cpp)
target_link_libraries(FullParallelMultilevelGaussianSampling ${MUQ_LIBRARIES} ${MUQ_LINK_LIBRARIES})

add_executable(SubsamplingTestMultilevelGaussianSampling SubsamplingTestMultilevelGaussianSampling.cpp)
target_link_libraries(SubsamplingTestMultilevelGaussianSampling ${MUQ_LIBRARIES} ${MUQ_LINK_LIBRARIES})

