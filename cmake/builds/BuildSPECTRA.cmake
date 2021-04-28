include(ExternalProject)

set(SPECTRA_DEPENDS )
if( USE_INTERNAL_EIGEN3 )
  list(APPEND SPECTRA_DEPENDS EIGEN3)
endif()

set(SPECTRA_INSTALL_DIR ${CMAKE_INSTALL_PREFIX}/muq_external/)
ExternalProject_Add(
  SPECTRA
    PREFIX ${CMAKE_CURRENT_BINARY_DIR}/external/spectra
    GIT_REPOSITORY https://github.com/yixuan/spectra.git
    DEPENDS ${SPECTRA_DEPENDS}
    LOG_DOWNLOAD OFF
    LOG_UPDATE OFF
    LOG_PATCH OFF
    LOG_CONFIGURE OFF
    LOG_BUILD OFF
    LOG_INSTALL OFF
    LOG_TEST OFF
    BUILD_COMMAND ""
  	CONFIGURE_COMMAND ""
    INSTALL_COMMAND ${CMAKE_COMMAND} -E make_directory ${SPECTRA_INSTALL_DIR}/spectra/include && cp -r ${CMAKE_CURRENT_BINARY_DIR}/external/spectra/src/SPECTRA/include/Spectra ${SPECTRA_INSTALL_DIR}/spectra/include
  )

set_property( TARGET SPECTRA PROPERTY FOLDER "Externals")

set(SPECTRA_INCLUDE_DIRS ${SPECTRA_INSTALL_DIR}/spectra/include/)
message(STATUS "Adding ${SPECTRA_INSTALL_DIR}/spectra/include/ for a Spectra include directory.")
