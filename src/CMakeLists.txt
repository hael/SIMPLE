
##########################################################################
# Source for the SIMPLE library.
# Better config and build alternative to old perl compile process
##########################################################################


file(GLOB DEFS_src "${CMAKE_CURRENT_SOURCE_DIR}/defs/*.f[90][08]" PARENT_SCOPE)
file(GLOB UTILS_src "${CMAKE_CURRENT_SOURCE_DIR}/simple_utils/*.f90" PARENT_SCOPE)
if ( USE_CUDA )
  file(GLOB CUDA_src "${CMAKE_CURRENT_SOURCE_DIR}/cuda/*.f[90][08]" PARENT_SCOPE)
##  list(REMOVE_ITEM UTILS_src "${CMAKE_CURRENT_SOURCE_DIR}/simple_utils/simple_timer_cuda.f90")
## else()
##   list(APPEND UTILS_src ${UTILS_src} "${CMAKE_SOURCE_DIR}/simple_utils/simple_timer_cuda.f90")
endif()

##########################################################################
# Simple_main folder
##########################################################################
file(GLOB SIMPLE_src "${CMAKE_CURRENT_SOURCE_DIR}/simple_main/*.f[09][0-9]" PARENT_SCOPE)


# Remove unused files
list(REMOVE_ITEM SIMPLE_src
  "${CMAKE_CURRENT_SOURCE_DIR}/simple_main/simple_polarft_corrcalc_backup.f90"
  "${CMAKE_CURRENT_SOURCE_DIR}/simple_main/simple_polarft_corrcalc_dble.f90"
  "${CMAKE_CURRENT_SOURCE_DIR}/simple_main/simple_cftcc_scatsrch.f90"
  "${CMAKE_CURRENT_SOURCE_DIR}/simple_main/simple_kninterpol_elemental.f90"
) #  "${CMAKE_CURRENT_SOURCE_DIR}/simple_main/simple_symsrcher.f90"


set(CORELIB_src  ${DEFS_src} ${UTILS_src} ${SIMPLE_src})

list(APPEND SIMPLECORELIB  ${DEFS_src} ${UTILS_src} ${SIMPLE_src})
# message(STATUS "SIMPLE LIBRARY SOURCES: ${SIMPLECORELIB}")

if (${CMAKE_Fortran_COMPILER_ID} STREQUAL "Intel")
  add_library(${SIMPLELIB} STATIC ${SIMPLECORELIB} )
  SET_TARGET_PROPERTIES(${SIMPLELIB}
    PROPERTIES
    LINKER_LANGUAGE Fortran
    SOVERSION ${${PROJECT_NAME}_VERSION}
    )
else()
  add_library(${SIMPLELIB} SHARED ${SIMPLECORELIB})
  SET_TARGET_PROPERTIES(${SIMPLELIB}
    PROPERTIES
    LINKER_LANGUAGE Fortran
    SOVERSION ${${PROJECT_NAME}_VERSION}
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

install(TARGETS ${SIMPLELIB}
  ARCHIVE DESTINATION lib
  LIBRARY DESTINATION lib
  )

# Make sure the compiler can find include files for our SIMPLE library
# when other libraries or executables link to it
target_include_directories (${SIMPLELIB} PUBLIC ${CMAKE_CURRENT_SOURCE_DIR})
