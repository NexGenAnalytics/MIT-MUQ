# PURPOSE:
# This file sets up the MUQ build targets (e.g., libMuqModeling).  Information is
# used from the compile groups that were processed in the ProcessCompileGroups.cmake
# file.
#

set(MUQ_LIBRARIES )
set(MUQ_PYTHON_LIBRARIES )

message("MUQ_LINK_LIBS = ${MUQ_LINK_LIBS}")
# Build all the targets
foreach(libName ${MUQ_TARGETS})

    list(LENGTH ${libName}_SOURCES strLength)
    if(${strLength} GREATER 0)

        message(STATUS "Creating ${libName} target.")
        string(REGEX MATCH "^pymuq" IsPythonWrapper ${libName})

        if(IsPythonWrapper)
            pybind11_add_module(${libName} SHARED NO_EXTRAS ${${libName}_SOURCES})
            list(APPEND MUQ_PYTHON_LIBRARIES ${libName})
        else()
            ADD_LIBRARY(${libName} SHARED ${${libName}_SOURCES})
            list(APPEND MUQ_LIBRARIES ${libName})
        endif()

        target_link_libraries(${libName} PUBLIC ${MUQ_LINK_LIBS})

        target_include_directories(
            ${libName}
            PUBLIC
            $<BUILD_INTERFACE:${CMAKE_SOURCE_DIR}>
            $<INSTALL_INTERFACE:${CMAKE_INSTALL_INCLUDEDIR}>
        )

        if(IsPythonWrapper)
            install(TARGETS ${libName}
                    EXPORT ${CMAKE_PROJECT_NAME}Depends
                    LIBRARY DESTINATION "${PYTHON_INSTALL_PREFIX}/muq"
                    ARCHIVE DESTINATION "${PYTHON_INSTALL_PREFIX}/muq")
        else()
            install(TARGETS ${libName}
                    EXPORT ${CMAKE_PROJECT_NAME}Depends
                    LIBRARY DESTINATION "${CMAKE_INSTALL_LIBDIR}"
                    ARCHIVE DESTINATION "${CMAKE_INSTALL_LIBDIR}")

            install(
                    EXPORT ${CMAKE_PROJECT_NAME}Depends
                    DESTINATION ${CMAKE_INSTALL_LIBDIR}/cmake/MUQ
                    NAMESPACE muq::
                    FILE "muqTargets.cmake"
            )            
        endif()
    endif()

endforeach()

INSTALL (
    DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}/MUQ
    DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}
    FILES_MATCHING PATTERN "*.h")

# If a group depends on an external library that is going to be built by MUQ, 
# then make sure we account for that dependency
foreach(group ${MUQ_GROUPS})
    if(${group}_IS_COMPILED)
        list(LENGTH ${group}_SOURCES strLength)

        # Add dependencies between different MUQ libraries
        foreach(depend ${${group}_REQUIRES_GROUPS})
            message(STATUS "Thinking about connection between ${${group}_LIBRARY} and ${${depend}_LIBRARY}")
            string(COMPARE EQUAL "${${depend}_LIBRARY}" "" result)
            if(NOT result)
                if(NOT ${${group}_LIBRARY} STREQUAL ${${depend}_LIBRARY})
                    IF( ${depend}_IS_COMPILED )
                        message(STATUS "Trying to add connection between ${${group}_LIBRARY} and ${${depend}_LIBRARY}")
                        target_link_libraries(${${group}_LIBRARY} PUBLIC ${${depend}_LIBRARY})
                        add_dependencies(${${group}_LIBRARY} ${${depend}_LIBRARY})
                    endif()
                endif()
              endif()
        endforeach()
    endif()
endforeach()

