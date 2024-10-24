SET(_log_summary  "${CMAKE_BINARY_DIR}/summary.log")
FILE(REMOVE ${_log_summary})

FILE(APPEND ${_log_summary}
"#############################################
#
#  MUQ configuration:
#        CMAKE_BUILD_TYPE:         ${CMAKE_BUILD_TYPE}
#        CMAKE_INSTALL_PREFIX:     ${CMAKE_INSTALL_PREFIX}
#        CMAKE_SOURCE_DIR:         ${CMAKE_SOURCE_DIR}
#        CMAKE_BINARY_DIR:         ${CMAKE_BINARY_DIR}
#        CMAKE_CXX_COMPILER name:  ${CMAKE_CXX_COMPILER_ID} ${CMAKE_CXX_COMPILER_VERSION} on platform ${CMAKE_SYSTEM_NAME} ${CMAKE_SYSTEM_PROCESSOR}
#        CMAKE_CXX_COMPILER path:  ${CMAKE_CXX_COMPILER}
#
"
)

IF(${CMAKE_BUILD_TYPE} MATCHES "Release")
FILE(APPEND ${_log_summary}
"#  Compiler flags used for this build:
#        CMAKE_CXX_FLAGS:        ${CMAKE_CXX_FLAGS} ${CMAKE_CXX_FLAGS_RELEASE}
"
  )
elseif(${CMAKE_BUILD_TYPE} MATCHES "Debug")
FILE(APPEND ${_log_summary}
"#  Compiler flags used for this build:
#        CMAKE_CXX_FLAGS:        ${CMAKE_CXX_FLAGS} ${CMAKE_CXX_FLAGS_DEBUG}
"
)
endif()

FILE(APPEND ${_log_summary}
"#
#  Compiler definitions used for this build:
#        COMPILE_DEFINITIONS:   ${MUQ_COMPILE_DEFINITIONS}
#
"
)

FILE(APPEND ${_log_summary} "#  MUQ links these dependencies: \n")
foreach(target ${MUQ_LINK_LIBS})
    # string(REPLACE "muq" "" moduleName ${target})
    FILE(APPEND ${_log_summary} "#        ${target}\n")
endforeach()
FILE(APPEND ${_log_summary}
"#
"
)

FILE(APPEND ${_log_summary} "#  MUQ Modules: \n")
foreach(target ${MUQ_TARGETS})
    string(REPLACE "muq" "" moduleName ${target})
    FILE(APPEND ${_log_summary} "#        ${moduleName}\n")
endforeach()

FILE(APPEND ${_log_summary} "#\n#  MUQ Libraries: \n")
FILE(APPEND ${_log_summary} "#        ${MUQ_LIBRARIES}\n")

FILE(APPEND ${_log_summary} "#############################################\n\n")
