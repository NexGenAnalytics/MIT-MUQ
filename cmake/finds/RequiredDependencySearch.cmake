# define a macro to look for a package and install a local copy if we can't find it
macro (GetDependency name)
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
                if(NOT ${name}_FOUND)
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
      	    LIST(APPEND MUQ_LINK_LIBS_STATIC ${${name}_LIBRARIES_STATIC})

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
GetDependency(STANMATH)
GetDependency(PARCER)
GetDependency(SPDLOG)
GetDependency(OTF2)


###########################
##### LOOK FOR HDF5  ######
###########################

find_package(HDF5 REQUIRED COMPONENTS C CXX HL)

###########################
##### LOOK FOR ZLIB  ######
###########################

if(MUQ_USE_OPENMPI)
	find_package(ZLIB)
	include_directories(${ZLIB_INCLUDE_DIRS})
	LIST(APPEND MUQ_LINK_LIBS ${ZLIB_LIBRARIES})
	LIST(APPEND MUQ_LINK_LIBS_STATIC ${ZLIB_LIBRARIES_STATIC})
	LIST(APPEND MUQ_EXTERNAL_INCLUDES ${ZLIB_INCLUDE_DIRS})
	message("ZLIB_LIBRARIES" ${ZLIB_LIBRARIES})

endif()

############################
##### LOOK FOR NLOPT  ######
############################

find_package(NLopt REQUIRED)

###############################
##### LOOK FOR SUNDIALS  ######
###############################

find_package(SUNDIALS 5.5.0...<6.0.0 REQUIRED)
set(MUQ_HAS_SUNDIALS 1) # this is needed for preprocessor directives in the MUQ source code

#############################
##### LOOK FOR EIGEN3  ######
#############################

find_package(Eigen3 3.3 REQUIRED NO_MODULE)

###############################
##### LOOK FOR NANOFLANN ######
###############################

find_package(nanoflann REQUIRED)

###############################################
##### LOOK FOR BOOST                     ######
###############################################

set(BOOST_MIN_VERSION "1.56.0")
find_package(Boost ${BOOST_MIN_VERSION} COMPONENTS system filesystem graph regex)

########################################
##### REMOVE DUPLICATE INCLUDES   ######
########################################

list( REMOVE_DUPLICATES MUQ_EXTERNAL_INCLUDES)
set(MUQ_EXTERNAL_INCLUDE_DIRS ${MUQ_EXTERNAL_INCLUDES} CACHE INTERNAL "List of external include directories for MUQ.")
