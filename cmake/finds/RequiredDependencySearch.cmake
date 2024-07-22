# define a macro to look for a package and install a local copy if we can't find it
macro (GetDependency name)
        message(STATUS "===== GetDependency for ${name} =====")
        # check to see if this dependency is required by any group
        list (FIND MUQ_REQUIRES ${name} dindex)
        if (${dindex} GREATER -1)
	          set(MUQ_NEEDS_${name} ON)

            if(NOT DEFINED MUQ_FORCE_INTERNAL_${name})
                set(MUQ_FORCE_INTERNAL_${name} OFF)
            endif()

            if(${MUQ_FORCE_INTERNAL_${name}})

                set(USE_INTERNAL_${name} 1)

            else()

                find_package(${name})
                if(${name}_FOUND)
                    message(STATUS "===== Found ${name} =====")
                    # check to make sure the library can be linked to
#                    include(Check${name})

                    # If the test code compiled...
                    if(NOT ${name}_TEST_FAIL)
                        set(USE_INTERNAL_${name} 0)
                    else()
                        message(STATUS "===== ${name} test failed, using internal version =====")
                        set(USE_INTERNAL_${name} 1)
                    endif()

                else()
                    message(STATUS "===== ${name} not found, using internal version =====")
                    set(USE_INTERNAL_${name} 1)
                endif()

            endif()

      	    if(USE_INTERNAL_${name})
      		    include(Build${name})
      	    endif()

      	    # store include directory information
      	    include_directories(${${name}_INCLUDE_DIRS})
      	    LIST(APPEND MUQ_EXTERNAL_INCLUDES ${${name}_INCLUDE_DIRS})

      	    # store library information
      	    LIST(APPEND MUQ_LINK_LIBS ${${name}_LIBRARIES})
            message(STATUS "MUQ_LINK_LIBS: ${MUQ_LINK_LIBS}")
      	    LIST(APPEND MUQ_LINK_LIBS_STATIC ${${name}_LIBRARIES_STATIC})
            message(STATUS "MUQ_LINK_LIBS_STATIC: ${MUQ_LINK_LIBS_STATIC}")


            set(MUQ_HAS_${name} 1)
        else()
            set(MUQ_NEEDS_${name} OFF)
	          set(MUQ_HAS_${name} 0)
        endif()

endmacro(GetDependency)

file(MAKE_DIRECTORY ${CMAKE_INSTALL_PREFIX}/muq_external/)
file(MAKE_DIRECTORY ${CMAKE_INSTALL_PREFIX}/muq_external/include)
file(MAKE_DIRECTORY ${CMAKE_INSTALL_PREFIX}/muq_external/lib)


include_directories(${CMAKE_CURRENT_SOURCE_DIR}/external/include)

########################################################
##### LOOK FOR AND/OR BUILD REQUIRED DEPENDENCIES ######
########################################################
# Using the default find_package() for dependencies is necessary if linking to libraries will be needed;
# GetDependency is used otherwise
GetDependency(EIGEN3) # header only,
GetDependency(STANMATH) # header only
GetDependency(SUNDIALS) # this dep doesn't have a cmake config file so find_package() doesn't work
                        # GetDependency() calls the custom findSUNDIALS.cmake which sets the libs and includes manually

# these 3 are required if mpi is on, but the mpi build doesn't work for now so not modifying these yet
GetDependency(PARCER)
GetDependency(SPDLOG)
GetDependency(OTF2)

############################
##### LOOK FOR NLopt  ######
############################

find_package(NLopt REQUIRED)
if (NLOPT_FOUND)
    set(MUQ_HAS_NLOPT 1)
endif ()

########################################
##### LOOK FOR AND/OR BUILD HDF5  ######
########################################

set(HAVE_HDF5 1)
if (MUQ_FORCE_INTERNAL_HDF5)
    message(STATUS "Forcing internal HDF5 build.")
    set(BUILD_INTERNAL_HDF5 TRUE)
else()
    find_package(HDF5 REQUIRED COMPONENTS C CXX HL)
    if (HDF5_FOUND)
        set(MUQ_HAS_HDF5 1)
        set(BUILD_INTERNAL_HDF5 FALSE)
    else()
        message(STATUS "HDF5 not found, performing internal build.")
        set(BUILD_INTERNAL_HDF5 TRUE)
    endif()
endif()

if (BUILD_INTERNAL_HDF5)
    if(NOT DEFINED MUQ_INTERNAL_HDF5_VERSION)
        set(MUQ_INTERNAL_HDF5_VERSION "1.8.19")
    endif()
    include(FetchContent)
    FetchContent_Declare(HDF5
            GIT_REPOSITORY https://github.com/HDFGroup/hdf5.git
            GIT_TAG hdf5-1_8_19
    )
    FetchContent_MakeAvailable(HDF5)
    install(
            TARGETS hdf5
            EXPORT ${PROJECT_NAME}Targets
    )
endif()

#GetDependency(HDF5)

if(MUQ_USE_OPENMPI)
	find_package(ZLIB)
	include_directories(${ZLIB_INCLUDE_DIRS})
	LIST(APPEND MUQ_LINK_LIBS ${ZLIB_LIBRARIES})
	LIST(APPEND MUQ_LINK_LIBS_STATIC ${ZLIB_LIBRARIES_STATIC})
	LIST(APPEND MUQ_EXTERNAL_INCLUDES ${ZLIB_INCLUDE_DIRS})
	message("ZLIB_LIBRARIES" ${ZLIB_LIBRARIES})

endif()

############################################
##### LOOK FOR AND/OR BUILD NANOFLANN ######
############################################

GetDependency(NANOFLANN)
set(MUQ_NANOFLAN_PARAMS_COMPILES 0)
if(USE_INTERNAL_NANOFLANN)
    set(MUQ_NANOFLAN_PARAMS_COMPILES 0)
endif()
###############################################
##### LOOK FOR BOOST                     ######
###############################################

find_package(Boost REQUIRED COMPONENTS system filesystem graph)
if (Boost_FOUND)
    set(MUQ_HAS_BOOST 1)
endif ()

#list (FIND MUQ_REQUIRES BOOST dindex)
#if (${dindex} GREATER -1)
#    set(MUQ_NEEDS_BOOST ON)
#
#    find_package(BOOSTMUQ)
#    if(NOT DEFINED Boost_FOUND)
#	set(Boost_FOUND ${BOOST_FOUND})
#    endif()
#
#    if(Boost_FOUND)
#	# check to make sure the library can be linked to
#	include(CheckBoost)
#
#	if(NOT BOOST_TEST_FAIL)
#		set(USE_INTERNAL_BOOST 0)
#	else()
#		set(USE_INTERNAL_BOOST 1)
#	endif()
#
#    else()
#	set(USE_INTERNAL_BOOST 1)
#    endif()
#
#    if(USE_INTERNAL_BOOST)
#	include(BuildBoost)
#    endif()
#
#    # do we want to compile the python interface?
#    set(MUQ_PYTHON 0)
#    if(MUQ_USE_PYTHON)
#        set(MUQ_PYTHON 1)
#    endif()
#
#    # store include directory information
#    if(NOT DEFINED Boost_INCLUDE_DIRS)
#        set(Boost_INCLUDE_DIRS ${BOOST_INCLUDE_DIR})
#    endif()
#
#    include_directories(${Boost_INCLUDE_DIRS})
#    LIST(APPEND ${CMAKE_PROJECT_NAME}_EXTERNAL_INCLUDES ${Boost_INCLUDE_DIRS})
#
#    if(NOT DEFINED Boost_LIBRARIES)
#        set(Boost_LIBRARIES ${BOOST_LIBRARY})
#        set(Boost_LIBRARIES_STATIC ${BOOST_LIBRARIES_STATIC})
#    endif()
#
#    # store library information
#    LIST(APPEND ${CMAKE_PROJECT_NAME}_LINK_LIBS ${Boost_LIBRARIES})
#    LIST(APPEND ${CMAKE_PROJECT_NAME}_LINK_LIBS_STATIC ${Boost_LIBRARIES_STATIC})
#
#else()
#    set(MUQ_NEEDS_BOOST OFF)
#
#endif()

########################################
##### REMOVE DUPLICATE INCLUDES   ######
########################################

list( REMOVE_DUPLICATES MUQ_EXTERNAL_INCLUDES)
set(MUQ_EXTERNAL_INCLUDE_DIRS ${MUQ_EXTERNAL_INCLUDES} CACHE INTERNAL "List of external include directories for MUQ.")
