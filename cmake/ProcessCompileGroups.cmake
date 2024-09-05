# PURPOSE:
# The code in this file loops through all of the compile groups, 
# extracts all the required and optional dependencies, and 
# collects all of the source files needed for the compile
# targets (e.g., libmuqModeling, etc...)
#

function(ForciblyEnable group)
  set(CAN_ENABLE_${group} ON)
  foreach(depend ${${group}_REQUIRES})
    if(NOT MUQ_USE_${depend})
      message(STATUS "    Cannot forcibly enable ${group} because of a disabled library dependency.")
      set(CAN_ENABLE_${group} OFF)
    endif()
  endforeach()

  if(CAN_ENABLE_${group})
    message(STATUS "  Forcibly enabling ${group}")
    set(MUQ_ENABLEGROUP_${group} ON CACHE INTERNAL "MUQ_ENABLEGROUP_${group}")

    foreach(depend ${${group}_REQUIRES_GROUPS})
        if(NOT MUQ_ENABLEGROUP_${depend})
            message(STATUS "    ${group} depends on the ${depend} group, but the ${depend} group was not enabled.")
            message(STATUS "    Turning the ${depend} group on.")
            set(MUQ_ENABLEGROUP_${depend} ON CACHE INTERNAL "MUQ_ENABLEGROUP_${depend}")
            ForciblyEnable(${depend})
        endif()
    endforeach()
  endif()
endfunction(ForciblyEnable)


# Initially, we have no targets to build
set(MUQ_TARGETS "" CACHE INTERNAL "List of MUQ libraries to build.")
set(MUQ_GROUPS "" CACHE INTERNAL "List of MUQ compile groups.")

# get info about every group
add_subdirectory(modules)

message("\n")
message("==========================================")
message("  PROCESSING COMPILE GROUPS               ")
message("==========================================")
set(MUQ_REQUIRES )
set(MUQ_DESIRES )
foreach(group ${MUQ_GROUPS})
      message(STATUS "${group} = ${MUQ_ENABLEGROUP_${group}}")

      set(CAN_ENABLE_${group} ON)
      foreach(depend ${${group}_REQUIRES})
        if(NOT MUQ_USE_${depend})
          message(STATUS "    Cannot enable ${group} because MUQ_USE_${depend}=OFF.")
          set(CAN_ENABLE_${group} OFF)
        endif()
      endforeach()

      if(CAN_ENABLE_${group})
        # Make sure all upstream dependency groups are enabled
        foreach(depend ${${group}_REQUIRES_GROUPS})
            if(MUQ_ENABLEGROUP_${group} AND NOT MUQ_ENABLEGROUP_${depend})
                message(STATUS "    ${group} depends on the ${depend} group, but the ${depend} group was not enabled.")
                ForciblyEnable(${depend})
            endif()
        endforeach()
      else()
        set(MUQ_ENABLEGROUP_${group} OFF)
      endif()

endforeach()


message("\n")
message("==========================================")
message("  PROCESSING TPLS                         ")
message("==========================================")
foreach(group ${MUQ_GROUPS})
  # if enabled, deal with this
  if(MUQ_ENABLEGROUP_${group})

      # loop over the requirements of this group
      foreach(depend ${${group}_REQUIRES})
        # if the TPL does NOT already exist, add it
        list (FIND MUQ_REQUIRES ${depend} dindex)
        if (${dindex} EQUAL -1)
          message(STATUS "Adding ${depend} to MUQ_REQUIRES because ${group} asked for it.")
          list(APPEND MUQ_REQUIRES ${depend})
        endif()
      endforeach()

      # Add to the list of desired (i.e., optional) external libraries
      foreach(depend ${${group}_DESIRES})
        list(APPEND MUQ_DESIRES ${depend})
      endforeach()

  endif()
endforeach()

# Remove duplicate requirements
list(REMOVE_DUPLICATES MUQ_REQUIRES)
if(MUQ_DESIRES)
    list(REMOVE_DUPLICATES MUQ_DESIRES)
endif()
MESSAGE(STATUS "List of TPLs needed = ${MUQ_REQUIRES}")

# Create a list of all MUQ libraries to build
set(MUQ_TARGETS )
foreach(group ${MUQ_GROUPS})
    if(MUQ_ENABLEGROUP_${group})
        message(STATUS "Adding target ${${group}_LIBRARY} for compile group ${group}")
        list(APPEND MUQ_TARGETS ${${group}_LIBRARY})
    endif()
endforeach()
list(REMOVE_DUPLICATES MUQ_TARGETS)
message("MUQ_TARGETS = ${MUQ_TARGETS}")

# Set up the source for each target library
foreach(target ${MUQ_TARGETS})
    set(${target}_SOURCES )

    foreach(group ${MUQ_GROUPS})
        if(MUQ_ENABLEGROUP_${group})
	          if(${${group}_LIBRARY} MATCHES ${target})

                # Check to see if a group has any source (e.g., *.cpp) files.  
                # Flag it as something that will be built if it does.
                list(LENGTH ${group}_SOURCES sources_length)
    	        if(sources_length GREATER 0)
                    set(${group}_IS_COMPILED ON CACHE INTERNAL 
                        "Whether or not the group ${group} is used in any library.")
    	        endif()

                list(APPEND ${target}_SOURCES ${${group}_SOURCES})
	          endif()
	      endif()
    endforeach()

    if(${target}_SOURCES)
        list(REMOVE_DUPLICATES ${target}_SOURCES)
    endif()
endforeach()
