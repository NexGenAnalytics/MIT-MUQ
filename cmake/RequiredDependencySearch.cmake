
file(MAKE_DIRECTORY ${CMAKE_INSTALL_PREFIX}/muq_external/)
file(MAKE_DIRECTORY ${CMAKE_INSTALL_PREFIX}/muq_external/include)
file(MAKE_DIRECTORY ${CMAKE_INSTALL_PREFIX}/muq_external/lib)
include_directories(${CMAKE_CURRENT_SOURCE_DIR}/external/include)

########################################################
##### REQUIRED DEPENDENCIES ######
########################################################

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
    find_package(Boost ${BOOST_MIN_VERSION} COMPONENTS system filesystem graph regex)
    LIST(APPEND MUQ_LINK_LIBS 
        Boost::system Boost::filesystem Boost::graph Boost::regex)
endif()

# ###############################################################
# ##### LOOK FOR Parallel Sampling Algorithm dependencies  ######
# ###############################################################

set(MUQ_HAS_PARCER 0) # needed for preprocessor directives 
set(MUQ_HAS_OTF2 0)   # needed for preprocessor directives 
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
endif()

########################################
##### REMOVE DUPLICATE INCLUDES   ######
########################################

list( REMOVE_DUPLICATES MUQ_EXTERNAL_INCLUDES)
set(MUQ_EXTERNAL_INCLUDE_DIRS ${MUQ_EXTERNAL_INCLUDES} 
    CACHE INTERNAL "List of external include directories for MUQ.")

