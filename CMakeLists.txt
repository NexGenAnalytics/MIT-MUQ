cmake_minimum_required (VERSION 3.10)

LIST(APPEND CMAKE_MODULE_PATH ${CMAKE_CURRENT_SOURCE_DIR}/cmake/)

INCLUDE(CheckCXXSourceCompiles)
INCLUDE(CheckCXXSourceRuns)
INCLUDE(CheckFunctionExists)


# define the project name
project(MUQ)
set(${CMAKE_PROJECT_NAME}_MAJOR_VERSION 0)
set(${CMAKE_PROJECT_NAME}_MINOR_VERSION 4)
set(${CMAKE_PROJECT_NAME}_PATCH_VERSION 3)
set(${CMAKE_PROJECT_NAME}_VERSION ${${CMAKE_PROJECT_NAME}_MAJOR_VERSION}.${${CMAKE_PROJECT_NAME}_MINOR_VERSION}.${${CMAKE_PROJECT_NAME}_PATCH_VERSION})

message(STATUS "\n\nThis is MUQ, the MIT UQ Library!\n\n")

###############################################
# User options
###############################################
INCLUDE(SetupOptions)
INCLUDE(SetInstallPaths)

###############################################
# Compiler configuration
###############################################
INCLUDE(SetupCompiler)
INCLUDE(CompilerCheck)
include_directories(.)

###############################################
# Find source files and calculate required dependencies
###############################################
INCLUDE(CompileGroup)

# The ProcesCompileGroups function will add modules/CMakeLists.txt and then process any computed Compile/Test groups
INCLUDE(ProcessCompileGroups)

###############################################
# Manage dependencies
###############################################
INCLUDE(DependencySearch)

include_directories(. ./external/include)

# write a configure file containing information about what libraries are available
CONFIGURE_FILE(${CMAKE_CURRENT_SOURCE_DIR}/MUQ/config.h.in ${CMAKE_CURRENT_SOURCE_DIR}/MUQ/config.h)

###############################################
# set module dependencies
###############################################
INCLUDE(DependencyCheck)

###############################################
# Create BUILD targets for cmake
###############################################

# FRIZZI: we need to revise the logic for this
# Remove any compile groups with unsatisfied dependencies
INCLUDE(UpdateCompileGroups)

# Set up libraries and test executables based on the compile groups
INCLUDE(SetupCompileTargets)
INCLUDE(SetupTestTargets)

###############################################################
# INSTALLATION INFORMATION
###############################################################
get_directory_property(MUQ_COMPILE_DEFINITIONS COMPILE_DEFINITIONS)
message(STATUS "Storing MUQ_COMPILE_DEFINITIONS = ${MUQ_COMPILE_DEFINITIONS}")
message(STATUS "MUQ_LIBRARIES = ${MUQ_LIBRARIES}")

set(MUQ_LIBRARIES_CACHE ${MUQ_LIBRARIES} CACHE INTERNAL "List of MUQ libraries for linking in parent projects.")
set(MUQ_LINK_LIBRARIES_CACHE ${MUQ_LINK_LIBRARIES} CACHE INTERNAL "List of external libraries needed by MUQ for use in parent projects.")

set(MUQ_EXPORT_LIBRARIES "")
FOREACH(BUILT_LIB ${MUQ_LIBRARIES})
	list(APPEND MUQ_EXPORT_LIBRARIES ${CMAKE_INSTALL_PREFIX}/lib/${CMAKE_SHARED_LIBRARY_PREFIX}${BUILT_LIB}${CMAKE_SHARED_LIBRARY_SUFFIX})
ENDFOREACH()
message(STATUS "MUQ_EXPORT_LIBRARIES = ${MUQ_EXPORT_LIBRARIES}")

#export(PACKAGE ${CMAKE_PROJECT_NAME})

# Create the MUQConfig.cmake file
configure_file(cmake/MUQConfig.cmake.in "${PROJECT_BINARY_DIR}/${CMAKE_PROJECT_NAME}Config.cmake" @ONLY)
configure_file(cmake/MUQConfig.cmake.in "${PROJECT_BINARY_DIR}${CMAKE_FILES_DIRECTORY}/${CMAKE_PROJECT_NAME}Config.cmake" @ONLY)

# Install the MUQConfig.cmake file
install(FILES "${PROJECT_BINARY_DIR}${CMAKE_FILES_DIRECTORY}/${CMAKE_PROJECT_NAME}Config.cmake"
        DESTINATION "${CMAKE_INSTALL_PREFIX}/lib/cmake/${CMAKE_PROJECT_NAME}"
	    COMPONENT dev)

# Install the pymuq libraries
install(DIRECTORY "${CMAKE_CURRENT_SOURCE_DIR}/modules/Modeling/python/pymuqModeling/"
        DESTINATION "${PYTHON_INSTALL_PREFIX}/muq/pymuqModeling"
        COMPONENT dev)
install(DIRECTORY "${CMAKE_CURRENT_SOURCE_DIR}/modules/PythonPackage/"
        DESTINATION "${PYTHON_INSTALL_PREFIX}/muq/"
        COMPONENT dev)

# Install the export set for use with the install-tree
install(EXPORT ${CMAKE_PROJECT_NAME}Depends
        DESTINATION "${CMAKE_INSTALL_PREFIX}/lib/cmake/${CMAKE_PROJECT_NAME}"
	    COMPONENT dev)

# install all the MUQ headers -- may be better off installing for each package individually
install(DIRECTORY MUQ
        DESTINATION "${CMAKE_INSTALL_PREFIX}/include"
	      FILES_MATCHING PATTERN "*.h")

# install the python libraries -- Is there a way to avoid copying this for each library?
#install(DIRECTORY "${CMAKE_CURRENT_SOURCE_DIR}/modules/Approximation/python/pymuqApproximation/"
#        DESTINATION "${CMAKE_INSTALL_PREFIX}/lib/pymuqApproximation"
#        PATTERN "wrappers/" EXCLUDE)

install(DIRECTORY spdlog
        DESTINATION "${CMAKE_INSTALL_PREFIX}/include"
	    FILES_MATCHING PATTERN "*.h")

set(CMAKE_LIBRARY_OUTPUT_DIRECTORY ${INSTALL_LIB_DIR})
set(CMAKE_ARCHIVE_OUTPUT_DIRECTORY ${INSTALL_LIB_DIR})

###############################################
# Write summary to terminal
###############################################
# write a summary of the build configuration
include(MakeSummaryFile)
FILE(READ ${CMAKE_BINARY_DIR}/summary.log MUQ_LOG_SUMMARY)
MESSAGE("${MUQ_LOG_SUMMARY}")

###############################################
# Add a target to build documentation
###############################################
# add a target to generate API documentation with Doxygen
find_package(Doxygen)
if(DOXYGEN_FOUND)

	if(PYTHON_EXECUTABLE)
		set(Python3_EXECUTABLE ${PYTHON_EXECUTABLE})
		set(Python3_Interpreter_FOUND ON)
	else()
		find_package(Python3 COMPONENTS Interpreter)
	endif()

	if(Python3_Interpreter_FOUND)
		message(STATUS "PYTHON3 EXE: ${Python3_EXECUTABLE}")
		message(STATUS "To build documentation, run 'make doc'")
		option(MUQDOC_INCLUDE_ANALYTICS OFF)
		if(MUQDOC_INCLUDE_ANALYTICS)
			message(STATUS "Will include google analytics header information in doxygen output.")
			set(ANALYTICS_HEADER ${CMAKE_CURRENT_SOURCE_DIR}/documentation/doxFiles/analytics_header.html)
		else()
			set(ANALYTICS_HEADER ${CMAKE_CURRENT_SOURCE_DIR}/documentation/doxFiles/noanalytics_header.html)
		endif()

	    configure_file(${CMAKE_CURRENT_SOURCE_DIR}/documentation/muq.doxyfile.in ${CMAKE_CURRENT_BINARY_DIR}/muq.doxyfile @ONLY)
	    add_custom_target(doc
	        COMMAND ${Python3_EXECUTABLE} ${CMAKE_CURRENT_SOURCE_DIR}/documentation/scripts/doxygen_preprocessing.py ${CMAKE_CURRENT_SOURCE_DIR} ${CMAKE_CURRENT_BINARY_DIR}
	        COMMAND ${DOXYGEN_EXECUTABLE} ${CMAKE_CURRENT_BINARY_DIR}/muq.doxyfile
	        WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}
	        COMMENT "Generating API documentation with Doxygen" VERBATIM
	    )

	endif(Python3_Interpreter_FOUND)
endif(DOXYGEN_FOUND)
