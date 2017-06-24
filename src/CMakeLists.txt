
#
# Source for the SIMPLE library.
# Better config and build alternative to perl compile process
#


#ADD_SUBDIRECTORY(${SRC}/defs ${LIB})
#ADD_SUBDIRECTORY(${SRC}/utils ${LIB})


file(GLOB DEFS_src "${CMAKE_SOURCE_DIR}/src/defs/*.f[90][08]" PARENT_SCOPE)
# message(STATUS "DEFS ${DEFS_src}")


# file(STRINGS "${CMAKE_SOURCE_DIR}/simple_utils/Makefile_target" UTILS_src REGEX "\(.*\)\\.o")
# message(STATUS "--${UTILS_src}--")
# string(REGEX REPLACE "\(.*:\)" "" UTILS_src ${UTILS_src})
# string(REGEX REPLACE "\\\\[;]?" "" UTILS_src ${UTILS_src})
# message(STATUS "UTILS2  --${UTILS_src}--")
# string(REGEX REPLACE "\(.*\).o" "${CMAKE_SOURCE_DIR}/simple_utils/\\1.f90" UTILS_src ${UTILS_src})
# message(STATUS "UTILS3 --${UTILS_src}--")

file(GLOB UTILS_src "${CMAKE_SOURCE_DIR}/src/simple_utils/*.f90" PARENT_SCOPE)
list(REMOVE_ITEM UTILS_src "${CMAKE_SOURCE_DIR}/src/simple_utils/simple_timer_cuda.f90")
# list(APPEND UTILS_src ${UTILS_src} "${CMAKE_SOURCE_DIR}/src/simple_utils/cuda/precision_m.f90")
# message(STATUS "UTILS3 --${UTILS_src}--")


##########################################################################
# Simple_main folder
##########################################################################
file(GLOB SIMPLE_src "${CMAKE_SOURCE_DIR}/src/simple_main/*.f[09][0-9]" PARENT_SCOPE)

list(REMOVE_ITEM SIMPLE_src
  "${CMAKE_SOURCE_DIR}/src/simple_main/simple_polarft_corrcalc_backup.f90"
  "${CMAKE_SOURCE_DIR}/src/simple_main/simple_polarft_corrcalc_dble.f90"
  "${CMAKE_SOURCE_DIR}/src/simple_main/simple_symsrcher.f90"
  "${CMAKE_SOURCE_DIR}/src/simple_main/simple_cftcc_scatsrch.f90")


# string(REGEX MATCHALL "(backup|dble|txt)"  ${SIMPLE_src})
# file(STRINGS "simple_main/Makefile_target" SIMPLE_src REGEX "^[^#]* \(.*\)\\.o")
# message(STATUS "${SIMPLE_src}")
# string(REGEX REPLACE "\([^:;]*\)\\.o" "simple_main/\\1.f90;" SIMPLE_src ${SIMPLE_src})
# string(REGEX REPLACE "\(.*:\)" "" SIMPLE_src ${SIMPLE_src})
# message(STATUS "SIMPLE_MAIN ${SIMPLE_src}")

set(CORELIB_src  ${DEFS_src} ${UTILS_src} ${SIMPLE_src})

# foreach(SOURCE ${CORELIB_src})
#   set(HASTIMER "")
#   file(STRINGS ${SOURCE} HASTIMER REGEX "^#include[ ]*[\"]simple_timer.h[\"].*$")
#   message(STATUS " ${SOURCE} HASTIMER --${HASTIMER}--")
#   if (NOT "${HASTIMER}x" MATCHES "x")
#     string(REGEX REPLACE ".*/\([^/]*\).(f90|f08)" "\\1.F08" SOURCE_F08 ${SOURCE})
#     message(STATUS " Replacing ${SOURCE} with ${SOURCE_F08}")
#     add_custom_target(${SOURCE_F08}
#       COMMAND "${CMAKE_CPP}" "${CMAKE_CPP_FLAGS}" "-E ${SOURCE} -o ${LIB}/${SOURCE_F08}"
#       DEPENDS ${SOURCE} "${CMAKE_SOURCE_DIR}/simple_utils/simple_timer.h"
#       WORKING_DIRECTORY ${LIB}
#       COMMENT "Preprocessing ${SOURCE} in CPP for variadic macros.")
#     string(REGEX REPLACE "${SOURCE}" "${LIB}/${SOURCE_F08}" CORELIB_src "${CORELIB_src}")
#   endif()
# endforeach(SOURCES)


if (USE_CUDA)
  list(APPEND UTILS_src ${UTILS_src} "${CMAKE_SOURCE_DIR}/src/simple_utils/simple_timer_cuda.f90")
endif()


list(APPEND SIMPLECORELIB  ${DEFS_src} ${UTILS_src} ${SIMPLE_src})
# message(STATUS "SIMPLE LIBRARY SOURCES: ${SIMPLECORELIB}")

if (${CMAKE_Fortran_COMPILER_ID} STREQUAL "Intel")
  add_library(${SIMPLELIB} SHARED ${SIMPLECORELIB} )
  SET_TARGET_PROPERTIES(${SIMPLELIB}
    PROPERTIES
    LINKER_LANGUAGE Fortran
    SOVERSION ${${PROJECT_NAME}_SOVERSION}
    )

else()
  add_library(${SIMPLELIB} SHARED ${SIMPLECORELIB})
  SET_TARGET_PROPERTIES(${SIMPLELIB}
    PROPERTIES
    LINKER_LANGUAGE Fortran
    SOVERSION ${${PROJECT_NAME}_SOVERSION}
    )
endif()

set_target_properties(${SIMPLELIB}
  PROPERTIES
  VERSION ${${PROJECT_NAME}_VERSION}
  COMPILE_FLAGS "${OpenMP_Fortran_FLAGS}"
  LINK_FLAGS "${OpenMP_Fortran_FLAGS}"
  LINK_LIBRARIES "${EXTRA_LIBS}"
  JOB_POOL_LINK NUM_JOBS
  )

# Make sure the compiler can find include files for our SIMPLE library
# when other libraries or executables link to Hello
target_include_directories (${SIMPLELIB} PUBLIC ${CMAKE_CURRENT_SOURCE_DIR})
