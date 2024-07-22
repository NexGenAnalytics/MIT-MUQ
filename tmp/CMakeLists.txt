cmake_minimum_required(VERSION 3.29)

project(test-proj LANGUAGES CXX)

#find_package(HDF5 REQUIRED COMPONENTS C CXX HL)
#find_package(SUNDIALS REQUIRED)
#find_package(NLopt REQUIRED)
#find_package(Boost REQUIRED COMPONENTS system filesystem graph)

#if(HDF5_FOUND)
#    message(STATUS "HDF5 found")
#else()
#    message(FATAL_ERROR "HDF5 not found")
#endif()

#if(SUNDIALS_FOUND)
#    message(STATUS "SUNDIALS found")
#else()
#    message(FATAL_ERROR "SUNDIALS not found")
#endif()
#
#if(NLopt_FOUND)
#    message(STATUS "NLopt found")
#else()
#    message(FATAL_ERROR "NLopt not found")
#endif()
#
#if(Boost_FOUND)
#    message(STATUS "Boost found")
#else()
#    message(FATAL_ERROR "Boost not found")
#endif()

include(FetchContent)
FetchContent_Declare(HDF5
        GIT_REPOSITORY https://github.com/HDFGroup/hdf5.git
        GIT_TAG hdf5-1_8_19
)
FetchContent_MakeAvailable(HDF5)
#install(
#        TARGETS HDF5
#        EXPORT ${PROJECT_NAME}Targets
#)

add_executable(test test.cpp)
target_link_libraries(test PRIVATE
        hdf5::hdf5
        hdf5::hdf5_cpp
        hdf5::hdf5_hl
        hdf5::hdf5_hl_cpp
#        SUNDIALS::cvodes
#        SUNDIALS::idas
#        SUNDIALS::kinsol
#        SUNDIALS::nvecserial
#        NLopt::nlopt
#        Boost::system
#        Boost::filesystem
#        Boost::graph
)
