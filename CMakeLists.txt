cmake_minimum_required (VERSION 3.10)

LIST(APPEND CMAKE_MODULE_PATH ${CMAKE_CURRENT_SOURCE_DIR}/cmake/)

INCLUDE(CheckCXXSourceCompiles)
INCLUDE(CheckCXXSourceRuns)
INCLUDE(CheckFunctionExists)

project(MUQ)
set(${CMAKE_PROJECT_NAME}_MAJOR_VERSION 0)
set(${CMAKE_PROJECT_NAME}_MINOR_VERSION 4)
set(${CMAKE_PROJECT_NAME}_PATCH_VERSION 3)
set(${CMAKE_PROJECT_NAME}_VERSION ${${CMAKE_PROJECT_NAME}_MAJOR_VERSION}.${${CMAKE_PROJECT_NAME}_MINOR_VERSION}.${${CMAKE_PROJECT_NAME}_PATCH_VERSION})

message(STATUS "\n\nThis is MUQ, the MIT UQ Library!\n\n")

###############################################
INCLUDE(SetupOptions)

if(PYTHON_INSTALL_PREFIX)
  message(STATUS
    "PYTHON_INSTALL_PREFIX was set by user to be ${PYTHON_INSTALL_PREFIX}.")
else()
  message(STATUS
    "PYTHON_INSTALL_PREFIX was not set by user, defaulting to CMAKE_INSTALL_PREFIX/python.")
  set(PYTHON_INSTALL_PREFIX ${CMAKE_INSTALL_PREFIX}/python)
endif()

message(STATUS "Python packages will be installed to ${PYTHON_INSTALL_PREFIX}.")


###############################################
# Compiler configuration
###############################################
INCLUDE(SetupCompiler)
INCLUDE(CompilerCheck)

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

# write a configure file containing information about what libraries are available
CONFIGURE_FILE(
  ${CMAKE_CURRENT_SOURCE_DIR}/MUQ/config.h.in
  ${CMAKE_CURRENT_SOURCE_DIR}/MUQ/config.h)

###############################################
# targets
###############################################
INCLUDE(SetupCompileTargets)
INCLUDE(SetupTestTargets)

###############################################################
# INSTALLATION INFORMATION
###############################################################

include(CMakePackageConfigHelpers)

configure_package_config_file(
  ${CMAKE_CURRENT_SOURCE_DIR}/cmake/MUQConfig.cmake.in
  "${CMAKE_CURRENT_BINARY_DIR}/MUQConfig.cmake"
  INSTALL_DESTINATION ${CMAKE_INSTALL_LIBDIR}/cmake/MUQ
)

install(FILES
  "${CMAKE_CURRENT_BINARY_DIR}/MUQConfig.cmake"
  DESTINATION ${CMAKE_INSTALL_LIBDIR}/cmake/MUQ
)

# Install the pymuq libraries
install(
  DIRECTORY "${CMAKE_CURRENT_SOURCE_DIR}/modules/Modeling/python/pymuqModeling/"
  DESTINATION "${PYTHON_INSTALL_PREFIX}/muq/pymuqModeling"
  COMPONENT dev
)

install(
  DIRECTORY "${CMAKE_CURRENT_SOURCE_DIR}/modules/PythonPackage/"
  DESTINATION "${PYTHON_INSTALL_PREFIX}/muq/"
  COMPONENT dev
)

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
