
########################################
##### LOOK FOR OPENMP             ######
########################################
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


########################################
##### LOOK FOR MPI                ######
########################################
set(MUQ_HAS_MPI 0)
if(MUQ_USE_MPI)
  find_package(MPI REQUIRED)
  list(APPEND MUQ_LINK_LIBS MPI::MPI_CXX)
  # list(APPEND MUQ_LINK_LIBS ${MPI_CXX_LIBRARIES} ${MPI_CXX_LINK_FLAGS})
  # list(APPEND MUQ_EXTERNAL_INCLUDES ${MPI_CXX_INCLUDE_DIRS})
  # include_directories(${MPI_CXX_INCLUDE_DIRS})
  # link_directories(${MPI_CXX_LIBRARIES})
  set(MUQ_HAS_MPI 1)
endif(MUQ_USE_MPI)


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
