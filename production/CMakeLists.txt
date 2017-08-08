

# Add the source files
set(SIMPLE_EXEC_src "${CMAKE_CURRENT_SOURCE_DIR}/simple/simple_exec/simple_exec.f90")
set(SIMPLE_DISTR_src "${CMAKE_CURRENT_SOURCE_DIR}/simple/simple_distr_exec/simple_distr_exec.f90")


add_executable(${SIMPLE_SERIAL_EXE} ${SIMPLE_EXEC_src})
add_executable(${SIMPLE_DISTR_EXE} ${SIMPLE_DISTR_src})

# This links foo to the bar library
target_link_libraries(${SIMPLE_SERIAL_EXE} LINK_PUBLIC ${EXTRA_LIBS})
target_link_libraries(${SIMPLE_DISTR_EXE} LINK_PUBLIC ${EXTRA_LIBS})

IF(USE_OPENMP)
  set_target_properties(${SIMPLE_SERIAL_EXE}
    PROPERTIES
      COMPILE_FLAGS "${OpenMP_Fortran_FLAGS}"
      LINK_FLAGS "${OpenMP_Fortran_FLAGS}"
      JOB_POOL_COMPILE NUM_JOBS)
  set_target_properties(${SIMPLE_DISTR_EXE}
    PROPERTIES
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


install(TARGETS "${SIMPLE_SERIAL_EXE}"
  ARCHIVE DESTINATION "${CMAKE_INSTALL_BINDIR}"
  RUNTIME DESTINATION "${CMAKE_INSTALL_BINDIR}"
         CONFIGURATIONS DEBUG|RELEASE
)
install(TARGETS "${SIMPLE_DISTR_EXE}"
  ARCHIVE DESTINATION "${CMAKE_INSTALL_BINDIR}"
  RUNTIME DESTINATION "${CMAKE_INSTALL_BINDIR}"
  CONFIGURATIONS DEBUG|RELEASE
)


#################################################################
#
# ENABLE TESTING
#
#################################################################
if (${PROJECT_NAME}_BUILD_TESTS)
  if(VERBOSE)
    include(CTest)
  endif()
endif()
#  add_subdirectory(tests)


########################################
# Set up how to compile the tests
########################################
file(GLOB TESTDIRS
  LIST_DIRECTORIES true
  RELATIVE ${CMAKE_CURRENT_SOURCE_DIR}
  "tests/test_*"
  )
message( STATUS "Test directories ${TESTDIRS}")


foreach(TESTNAME ${TESTDIRS})
  string(REPLACE "tests/" "" TESTNAME "${TESTNAME}")
  set(TESTEXE simple_${TESTNAME})
  file(GLOB_RECURSE TESTSOURCES LIST_DIRECTORIES true
    "tests/${TESTNAME}/*.[fF][90][0538]")
  if("" STREQUAL "${TESTSOURCES}")
    # message(STATUS " skipping ${TESTNAME} -- no files in dir")
  else()

    add_executable(${TESTEXE} ${TESTSOURCES})
    target_link_libraries(${TESTEXE} ${EXTRA_LIBS})
    if(USE_OPENMP)
      set_target_properties(${TESTEXE} PROPERTIES
        COMPILE_FLAGS "${OpenMP_Fortran_FLAGS}"
        LINK_FLAGS "${OpenMP_Fortran_FLAGS}"
        LANGUAGE Fortran)
    endif(USE_OPENMP)
    install(TARGETS ${TESTEXE}
      ARCHIVE DESTINATION ${CMAKE_INSTALL_BINDIR}
      RUNTIME DESTINATION ${CMAKE_INSTALL_BINDIR})

    file(GLOB_RECURSE TESTTXTSOURCES LIST_DIRECTORIES true "tests/${TESTNAME}/*.txt")
    if(NOT "${TESTTXTSOURCES}" STREQUAL "")
    foreach(TXTNAME "${TESTTXTSOURCES}") 
       install(FILES ${TXTNAME}
       DESTINATION ${CMAKE_INSTALL_PREFIX}/tests/${TESTNAME}
       CONFIGURATIONS DEBUG|RELEASE
       PERMISSIONS  OWNER_READ OWNER_WRITE GROUP_READ GROUP_WRITE WORLD_READ WORLD_WRITE
       )
endforeach()
    endif()
  endif()
endforeach(TESTNAME)

#
# Default CTest tests
#
list(APPEND CTESTDIRS
  test_imgfile
  test_install
  test_parse
  test_speed
  test_srch
  test_timer
  test_units )

add_test(NAME test_parse
  CONFIGURATIONS Debug|Release
  WORKING_DIRECTORY ${CMAKE_INSTALL_BINDIR}
  COMMAND simple_test_parse nthr=4 verbose=yes) 

add_test(NAME test_install
  CONFIGURATIONS Debug|Release
  WORKING_DIRECTORY ${CMAKE_INSTALL_BINDIR}
  COMMAND simple_test_install nthr=4 verbose=yes) 


add_test(NAME test_units
  CONFIGURATIONS Debug|Release
  WORKING_DIRECTORY ${CMAKE_INSTALL_BINDIR}
  COMMAND simple_test_units nthr=4) 

add_test(NAME test_timer
  CONFIGURATIONS Debug|Release
  WORKING_DIRECTORY ${CMAKE_INSTALL_BINDIR}
  COMMAND simple_test_timer verbose=yes debug=yes) 


