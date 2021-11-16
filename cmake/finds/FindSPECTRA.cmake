
find_package(PkgConfig)

if(NOT DEFINED MUQ_SPECTRA_DIR)
	find_path(SPECTRA_INCLUDE_DIR Spectra)
else()
	find_path(SPECTRA_INCLUDE_DIR Spectra
	          HINTS ${MUQ_SPECTRA_DIR})
endif()
set(SPECTRA_INCLUDE_DIRS ${SPECTRA_INCLUDE_DIR} )

set(SPECTRA_FOUND 1)
if( NOT SPECTRA_INCLUDE_DIR )
	set(SPECTRA_FOUND 0)
endif()
