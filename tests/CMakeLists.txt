
########################################
# Set up how to compile the tests
########################################
file(GLOB TESTDIRS
  LIST_DIRECTORIES true
  RELATIVE ${CMAKE_CURRENT_SOURCE_DIR}
  "test_*"
  )
message( STATUS "Test directories ${TESTDIRS}")


foreach(TESTNAME ${TESTDIRS})
  set(TESTEXE simple_${TESTNAME})
  file(GLOB_RECURSE TESTSOURCES LIST_DIRECTORIES true
    "${TESTNAME}/*.[fF][90][0538]")
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

    file(GLOB_RECURSE TESTTXTSOURCES LIST_DIRECTORIES true "${TESTNAME}/*.txt")
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
