
file(MAKE_DIRECTORY ${CMAKE_INSTALL_PREFIX}/muq_external/)
file(MAKE_DIRECTORY ${CMAKE_INSTALL_PREFIX}/muq_external/include)
file(MAKE_DIRECTORY ${CMAKE_INSTALL_PREFIX}/muq_external/lib)
include_directories(${CMAKE_CURRENT_SOURCE_DIR}/external/include)

# OPENMP
IF(MUQ_USE_OPENMP)
  find_package(OpenMP)
  if (OpenMP_CXX_FOUND)
    list(APPEND MUQ_LINK_LIBS OpenMP::OpenMP_CXX)

    CHECK_CXX_COMPILER_FLAG("-pthread" HAS_PTHREAD)
    if(HAS_PTHREAD)
      set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -pthread")
    endif()

    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -ldl")
  else()
    message(WARNING "The flag MUQ_USE_OPENMP is ON, but cmake cannot find OpenMP for c++. OpenMP will not be used.")
  endif()
ENDIF(MUQ_USE_OPENMP)

# HDF5
list (FIND MUQ_REQUIRES HDF5 dindex)
if (${dindex} GREATER -1)
    find_package(HDF5 REQUIRED COMPONENTS C CXX HL)
    LIST(APPEND MUQ_LINK_LIBS hdf5::hdf5 hdf5::hdf5_cpp hdf5::hdf5_hl)
endif()

# NLOPT
list (FIND MUQ_REQUIRES NLOPT dindex)
if (${dindex} GREATER -1)
    find_package(NLopt REQUIRED)
    LIST(APPEND MUQ_LINK_LIBS NLopt::nlopt)
endif()

# SUNDIALS
set(MUQ_HAS_SUNDIALS 0) # needed for preprocessor directives in the MUQ source code
list (FIND MUQ_REQUIRES SUNDIALS dindex)
if (${dindex} GREATER -1)
    find_package(SUNDIALS 5.5.0...<6.0.0 REQUIRED)
    set(MUQ_HAS_SUNDIALS 1)
    LIST(APPEND MUQ_LINK_LIBS
        SUNDIALS::cvodes  SUNDIALS::idas SUNDIALS::kinsol SUNDIALS::nvecserial)
endif()

# EIGEN3
list (FIND MUQ_REQUIRES EIGEN3 dindex)
if (${dindex} GREATER -1)
    find_package(Eigen3 3.3 REQUIRED NO_MODULE)
    LIST(APPEND MUQ_LINK_LIBS Eigen3::Eigen)
endif()

# NANOFLANN
list (FIND MUQ_REQUIRES NANOFLANN dindex)
if (${dindex} GREATER -1)
    find_package(nanoflann REQUIRED)
    LIST(APPEND MUQ_LINK_LIBS nanoflann::nanoflann)
endif()

# STANMATH
list (FIND MUQ_REQUIRES STANMATH dindex)
if (${dindex} GREATER -1)
    if(EXISTS "${stanmath_SRC_DIR}/stan/math.hpp")
        include_directories(${stanmath_SRC_DIR})
        LIST(APPEND MUQ_EXTERNAL_INCLUDES ${stanmath_SRC_DIR})
    else()
        message(FATAL_ERROR "stanmath directory provided doesn't contain stan/math.hpp")
    endif()
endif()

# BOOST
list (FIND MUQ_REQUIRES BOOST dindex)
if (${dindex} GREATER -1)
    set(BOOST_MIN_VERSION "1.56.0")
    find_package(Boost ${BOOST_MIN_VERSION} REQUIRED COMPONENTS system filesystem graph regex)
    LIST(APPEND MUQ_LINK_LIBS
        Boost::system Boost::filesystem Boost::graph Boost::regex)
endif()

# ###############################################################
# ##### LOOK FOR Parallel Sampling Algorithm dependencies  ######
# ###############################################################

set(MUQ_HAS_PARCER 0) # needed for preprocessor directives
set(MUQ_HAS_OTF2 0)   # needed for preprocessor directives
set(MUQ_HAS_MPI 0)
if(MUQ_USE_MPI)
    # PARCER
    find_package(PARCER REQUIRED)
    include_directories(${PARCER_INCLUDE_DIRS})
    LIST(APPEND MUQ_EXTERNAL_INCLUDES ${PARCER_INCLUDE_DIRS})
    LIST(APPEND MUQ_LINK_LIBS ${PARCER_LIBRARIES})
    LIST(APPEND MUQ_LINK_LIBS_STATIC ${PARCER_LIBRARIES_STATIC})
    set(MUQ_HAS_PARCER 1)

    # otf2
    set(OTF2_LIBRARIES ${otf2_DIR}/lib/${CMAKE_STATIC_LIBRARY_PREFIX}otf2${CMAKE_STATIC_LIBRARY_SUFFIX})
    set(OTF2_LIBRARY ${OTF2_LIBRARIES})
    set(OTF2_INCLUDE_DIRS ${otf2_DIR}/include)
    include_directories(${OTF2_INCLUDE_DIRS})
    LIST(APPEND MUQ_EXTERNAL_INCLUDES ${OTF2_INCLUDE_DIRS})
    LIST(APPEND MUQ_LINK_LIBS ${OTF2_LIBRARIES})
    LIST(APPEND MUQ_LINK_LIBS_STATIC ${OTF2_LIBRARIES_STATIC})
    set(MUQ_HAS_OTF2 1)

    # spdlog
    find_package(spdlog REQUIRED)
    LIST(APPEND MUQ_LINK_LIBS spdlog::spdlog)

    # MPI
    find_package(MPI REQUIRED)
    list(APPEND MUQ_LINK_LIBS MPI::MPI_CXX)
    set(MUQ_HAS_MPI 1)
endif()


# REMOVE DUPLICATES
list( REMOVE_DUPLICATES MUQ_EXTERNAL_INCLUDES)
set(MUQ_EXTERNAL_INCLUDE_DIRS ${MUQ_EXTERNAL_INCLUDES}
    CACHE INTERNAL "List of external include directories for MUQ.")


########################################
##### LOOK FOR GTEST              ######
########################################
IF(MUQ_USE_GTEST)
  set(GTEST_VERSION "v1.14.0")
  set(BUILD_GMOCK   OFF)
  set(INSTALL_GTEST OFF)
  message(STATUS "GTest not found, fetching version ${GTEST_VERSION}")

  list(APPEND CMAKE_MESSAGE_INDENT "[GTest] ")
  include(FetchContent)
  FetchContent_Declare(
    googletest
    DOWNLOAD_EXTRACT_TIMESTAMP FALSE
    URL https://github.com/google/googletest/archive/refs/tags/${GTEST_VERSION}.tar.gz
    URL_HASH MD5=c8340a482851ef6a3fe618a082304cfc
  )
  FetchContent_MakeAvailable(googletest)
  list(POP_BACK CMAKE_MESSAGE_INDENT)

ELSE(MUQ_USE_GTEST)
    message(STATUS “MUQ_USE_GTEST is OFF.  Turning off tests.”)
    set(MUQ_BUILD_TESTS OFF)
    set(MUQ_NEEDS_GTEST OFF)
ENDIF(MUQ_USE_GTEST)


########################################
##### LOOK FOR PYTHON             ######
########################################
list (FIND MUQ_REQUIRES PYTHON dindex)
if (${dindex} GREATER -1)
    set(MUQ_NEEDS_PYTHON ON)

    if(MUQ_USE_PYTHON)
        set(PYBIND11_CPP_STANDARD -std=c++11)
        FIND_PACKAGE(pybind11)
        if(NOT pybind11_FOUND)
            message(STATUS "Falling back to internal pybind11 version")
            add_subdirectory(${CMAKE_SOURCE_DIR}/external/pybind11)
            include_directories(${CMAKE_SOURCE_DIR}/external/pybind11/include)
        endif()

        message("PYTHON_SITE_PACKAGES = ${PYTHON_SITE_PACKAGES}")

    endif()
else()
    set(MUQ_NEEDS_PYTHON OFF)
    set(MUQ_USE_PYTHON OFF)
endif()

########################################
##### LOOK FOR DOLFIN/Fenics      ######
########################################
list (FIND MUQ_REQUIRES DOLFIN dindex)
message(${MUQ_REQUIRES})
if (${dindex} GREATER -1)
    set(MUQ_NEEDS_DOLFIN ON)

    if(MUQ_USE_DOLFIN AND NOT MUQ_USE_PYTHON)
        message(WARNING "Requested compilation with Fenics/Dolfin, but not Python.  Building the Fenics/Dolfin bindings requires building MUQ with Python support.")
        set(MUQ_USE_DOLFIN OFF)
    endif()

    if(MUQ_USE_DOLFIN)
        find_package(DOLFIN)
        if(DOLFIN_FOUND)
            if (EXISTS ${DOLFIN_USE_FILE})
                include(${DOLFIN_USE_FILE})

                # Default build type (can be overridden by user)
                if (NOT CMAKE_BUILD_TYPE)
                    set(CMAKE_BUILD_TYPE "RelWithDebInfo" CACHE STRING
                        "Choose the type of build, options are: Debug MinSizeRel Release RelWithDebInfo." FORCE)
                endif()
            else()
                # Compiler definitions
                add_definitions(${DOLFIN_CXX_DEFINITIONS})

                # Include directories
                include_directories(${DOLFIN_INCLUDE_DIRS})
                include_directories(SYSTEM ${DOLFIN_3RD_PARTY_INCLUDE_DIRS})
            endif()
        else()
            set(MUQ_USE_DOLFIN OFF)
        endif()

    endif()
else()
    set(MUQ_NEEDS_DOLFIN OFF)
    set(MUQ_USE_DOLFIN OFF)
endif()


########################################
##### REMOVE DUPLICATE INCLUDES   ######
########################################
if(MUQ_LINK_LIBS)
  list( REMOVE_DUPLICATES MUQ_LINK_LIBS)
endif()
