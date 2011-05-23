message(STATUS "Entering package ${CMAKE_CURRENT_SOURCE_DIR}")

include_directories(${CMAKE_CURRENT_SOURCE_DIR}/include)

file(GLOB DIRECTORIES ${CMAKE_CURRENT_SOURCE_DIR}/*)
find_path(FLINK_PATH link_DchCalib.mk)
find_path(FBIN_PATH bin_DchCalib.mk PATHS ${DIRECTORIES})

#
# Recursively read link_<somepack>.mk and bin_<somepack>.mk to find the list of libraries a given
# bin or lib has to be linked to
#
foreach(type lib bin)
   string(TOUPPER ${type} PREFIX)
   execute_process(COMMAND ${SCRIPTS_PATH}/depresolve.py
                   ${RELEASE_ROOT} ${CMAKE_CURRENT_SOURCE_DIR} DchCalib ${type}
                   OUTPUT_VARIABLE ${PREFIX}_LINK_LIST ERROR_VARIABLE err
                   OUTPUT_STRIP_TRAILING_WHITESPACE)
   if(${err} MATCHES "ERROR")
       if(${type} MATCHES "link")
           message(FATAL_ERROR "Problems with lib .mk file: ${err}")
       elseif(${type} MATCHES "bin")
           message(WARNING "${err}, bin .mk file not found. Setting linker list to default lib one")

       else()
          message(FATAL_ERROR "${err}, Check dependency files!")
       endif()
   endif()
endforeach()

#
# Di nuovo una cosa piuttosto orrida. Rimuovo dalla lista delle librerie da
# linkare i segnaposto della lista di cui al precedente commento in italiano
#
foreach(ext LIB_LINK_LIST BIN_LINK_LIST)
   string(STRIP "${LIB_LINK_LIST}" LIB_LINK_LIST)
   string(STRIP "${BIN_LINK_LIST}" BIN_LINK_LIST)
   list(LENGTH ${ext} ll)
   if(${ll} EQUAL "1")
       if(NOT ${${ext}} MATCHES "^ERROR")
           SetSpecialLib(${ext})
       endif()
   else()
       SetSpecialLib(${ext})
   endif()
endforeach()
list(REVERSE LIB_LINK_LIST)
list(REVERSE BIN_LINK_LIST)

#
# To avoid listing a huge number of src files we "glob" it along with their
# headers
#
set(RELATIVE_SRC_PATH ${CMAKE_CURRENT_SOURCE_DIR}/src)
set(RELATIVE_HH_PATH ${CMAKE_CURRENT_SOURCE_DIR}/include)

file(GLOB DchCalib_SRCS src/*.c*) 
file(GLOB DchCalib_HH include/*.h*) 

set(DCHCALIB_CXX_FLAGS "-Wall")
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${DCHCALIB_CXX_FLAGS}")

#
# Generate library
#
if(VERBOSE)
   message(STATUS "LIB_LINK_LIST ${LIB_LINK_LIST}")
   message(STATUS "BIN_LINK_LIST ${BIN_LINK_LIST}")
endif()

add_library(DchCalib ${DchCalib_SRCS})

#
# If the bin_list of libraries to be linked against the exe is empty, we switch to
# the lib_list
#
if("${BIN_LINK_LIST}" STREQUAL "")
   set(BIN_LINK_LIST ${LIB_LINK_LIST})
   if(VERBOSE)
       message(STATUS "New BIN_LINK_LIST ${BIN_LINK_LIST}")
   endif()
endif()

target_link_libraries(DchCalib ${LIB_LINK_LIST})
