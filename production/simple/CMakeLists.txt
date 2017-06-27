########################################
# Set up how to compile the source files
########################################

# Add the source files
set(SIMPLE_EXE_src "${CMAKE_CURRENT_SOURCE_DIR}/simple_exec/simple_exec.f90")
set(SIMPLE_DISTR_src "${CMAKE_CURRENT_SOURCE_DIR}/simple_distr_exec/simple_distr_exec.f90")
#  link_directories(${LIB})

# Define the executable in terms of the source files


  set(HASTIMER "")
  file(STRINGS ${SIMPLE_EXE_src} HASTIMER REGEX "^#include[ ]*[\"]simple_timer.h[\"].*$")
  if (NOT "${HASTIMER}x" MATCHES "x")
    string(REGEX REPLACE "(f90|f08)" "F08" SOURCE_F08 ${SIMPLE_EXE_src})
    message(STATUS " Replacing ${SIMPLE_EXE_src} with ${SOURCE_F08}")
    add_custom_target(${SOURCE_F08} 
      COMMAND "${CMAKE_FortranCPP}" "${CMAKE_FortranCPP_FLAGS}" "-E ${SIMPLE_EXE_src} -o ${SOURCE_F08}"
      DEPENDS ${SIMPLE_EXE_src} "${CMAKE_SOURCE_DIR}/simple_utils/simple_timer.h"
      WORKING_DIRECTORY ${LIB}
      COMMENT "Preprocessing ${SIMPLE_EXE_src} in CPP for variadic macros."
      )
    set(SIMPLE_EXE_src  "${SOURCE_F08}")
  endif()
  file(STRINGS ${SIMPLE_DISTR_src} HASTIMER REGEX "^#include[ ]*[\"]simple_timer.h[\"].*$")
  if (NOT "${HASTIMER}x" MATCHES "x")
    string(REGEX REPLACE ".*/\([^/]*\).(f90|f08)" "\\1.F08" SOURCE_F08 ${SIMPLE_DISTR_src})

    message(STATUS " Replacing ${SIMPLE_DISTR_src} with ${SOURCE_F08}")
    add_custom_target(${SOURCE_F08} 
      COMMAND "${CMAKE_FortranCPP}" "${CMAKE_FortranCPP_FLAGS}" "-E ${SIMPLE_DISTR_src} -o ${SOURCE_F08}"
      DEPENDS ${SIMPLE_DISTR_src} "${CMAKE_SOURCE_DIR}/simple_utils/simple_timer.h"
      WORKING_DIRECTORY ${LIB}
      COMMENT "Preprocessing ${SIMPLE_DISTR_src} in CPP for variadic macros."
      )
    set(SIMPLE_EXE_src  "${SOURCE_F08}")

  endif()




add_executable(${SIMPLE_EXE} ${SIMPLE_EXE_src})
add_executable(${SIMPLE_DISTR_EXE} ${SIMPLE_DISTR_src})

#####################################################
# Add the needed libraries and special compiler flags
#####################################################

# This links foo to the bar library
target_link_libraries(${SIMPLE_EXE} LINK_PUBLIC  ${EXTRA_LIBS})
target_link_libraries(${SIMPLE_DISTR_EXE} LINK_PUBLIC ${EXTRA_LIBS})

IF(USE_OPENMP)
    set_target_properties(${SIMPLE_EXE} PROPERTIES
      COMPILE_FLAGS "${OpenMP_Fortran_FLAGS}"
      LINK_FLAGS "${OpenMP_Fortran_FLAGS}"
      JOB_POOL_COMPILE NUM_JOBS)
    set_target_properties(${SIMPLE_DISTR_EXE} PROPERTIES
      COMPILE_FLAGS "${OpenMP_Fortran_FLAGS}"
      LINK_FLAGS "${OpenMP_Fortran_FLAGS}"
      JOB_POOL_COMPILE NUM_JOBS)
ELSEIF(USE_MPI)
#    SET_TARGET_PROPERTIES(${SIMPLE_EXECEXE} PROPERTIES
#                          COMPILE_FLAGS "${MPI_Fortran_COMPILE_FLAGS}"
#                          LINK_FLAGS "${MPI_Fortran_LINK_FLAGS}")
#    INCLUDE_DIRECTORIES(${MPI_Fortran_INCLUDE_PATH})
#    TARGET_LINK_LIBRARIES(${SIMPLE_EXECEXE} ${MPI_Fortran_LIBRARIES})
ENDIF(USE_OPENMP)



# Export HEADER
# Enable portability to Windows
#

#include (GenerateExportHeader)
#GENERATE_EXPORT_HEADER(shared           # generates the export header shared_EXPORTS.h automatically
#  BASE_NAME shared
#  EXPORT_MACRO_NAME shared_EXPORTS
#  EXPORT_FILE_NAME shared_EXPORTS.h
#  STATIC_DEFINE SHARED_EXPORTS_BUILT_AS_STATIC)
#target_link_libraries(main shared)
# link our previously created shared library "shared" to the exectable "main"












#####################################
# Tell how to install this executable
#####################################


install(TARGETS ${SIMPLE_EXE} RUNTIME DESTINATION bin)
install(TARGETS ${SIMPLE_DISTR_EXE} RUNTIME DESTINATION bin)
