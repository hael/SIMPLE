# Add the source files
add_executable(${SIMPLE_EXE} "${CMAKE_CURRENT_SOURCE_DIR}/simple_exec/simple_exec.f90")
add_executable(${SINGLE_EXE} "${CMAKE_CURRENT_SOURCE_DIR}/single_exec/single_exec.f90")
add_executable(${SIMPLE_PRIVATE_EXE}
"${CMAKE_CURRENT_SOURCE_DIR}/simple_private_exec/simple_private_exec.f90")
set(SIMPLE_MAIN_TARGETS ${SIMPLE_EXE} ${SINGLE_EXE} ${SIMPLE_PRIVATE_EXE})
# target_link_libraries(${SIMPLE_MAIN_TARGETS} LINK_PUBLIC ${EXTRA_LIBS})
target_link_libraries(${SIMPLE_EXE} LINK_PUBLIC ${EXTRA_LIBS})
target_link_libraries(${SINGLE_EXE} LINK_PUBLIC ${EXTRA_LIBS})
target_link_libraries(${SIMPLE_PRIVATE_EXE} LINK_PUBLIC ${EXTRA_LIBS})
if(USE_OPENMPI AND USE_MPI)
  set_target_properties(${SIMPLE_MAIN_TARGETS}
    PROPERTIES
    COMPILE_FLAGS "${OpenMP_Fortran_FLAGS} ${MPI_Fortran_FLAGS}"
    LINK_FLAGS "${OpenMP_Fortran_FLAGS} ${MPI_Fortran_LINK_FLAGS}"
    JOB_POOL_COMPILE NUM_JOBS)
  target_link_libraries(${SIMPLE_MAIN_TARGETS} LINK_PUBLIC  ${MPI_Fortran_LIBRARIES})
elseif(USE_OPENMP)
  set_target_properties(${SIMPLE_MAIN_TARGETS}
    PROPERTIES
    COMPILE_FLAGS "${OpenMP_Fortran_FLAGS}"
    LINK_FLAGS "${OpenMP_Fortran_FLAGS}"
    JOB_POOL_COMPILE NUM_JOBS)
elseif(USE_MPI)
  set_target_properties(${SIMPLE_MAIN_TARGETS}
    PROPERTIES
    COMPILE_FLAGS "${MPI_Fortran_FLAGS}"
    LINK_FLAGS "${MPI_Fortran_LINK_FLAGS}"
    JOB_POOL_COMPILE NUM_JOBS)
  target_link_libraries(${SIMPLE_MAIN_TARGETS} LINK_PUBLIC  ${MPI_Fortran_LIBRARIES})
endif()
install(TARGETS ${SIMPLE_MAIN_TARGETS}
  ARCHIVE DESTINATION "${CMAKE_INSTALL_BINDIR}"
  RUNTIME DESTINATION "${CMAKE_INSTALL_BINDIR}"
  CONFIGURATIONS DEBUG|RELEASE
)

# ENABLE TESTING
if (${PROJECT_NAME}_BUILD_TESTS)
  if(VERBOSE)
    include(CTest)
  endif()
endif()

# Set up how to compile the tests
file(GLOB TESTDIRS
  LIST_DIRECTORIES true
  RELATIVE ${CMAKE_CURRENT_SOURCE_DIR}
  "tests/test_*"
  )
set(DONT_CHECK_RANGE "-fno-range-check")
if(NOT USE_MPI)
  list(REMOVE_ITEM TESTDIRS "tests/test_mpi")
endif()
 list(REMOVE_ITEM TESTDIRS "tests/test_mpimp_benchmark")
if(NOT USE_NLOPT)
  list(REMOVE_ITEM TESTDIRS "tests/test_nlopt")
  list(REMOVE_ITEM TESTDIRS "tests/test_nlopt_lp_spectrum")
endif()

if(USE_MPI AND USE_OPENMP)
    set(EXTRA_FFLAGS "${EXTRA_FFLAGS} ${OpenMP_Fortran_FLAGS} ${MPI_Fortran_FLAGS}")
elseif(USE_MPI)
    set(EXTRA_FFLAGS "${EXTRA_FFLAGS} ${MPI_Fortran_FLAGS}")
elseif(USE_OPENMP)
    set(EXTRA_FFLAGS "${EXTRA_FFLAGS} ${OpenMP_Fortran_FLAGS}")
    if(${CMAKE_Fortran_COMPILER_ID} STREQUAL "GNU")
        set(EXTRA_FFLAGS "${EXTRA_FFLAGS} ${DONT_CHECK_RANGE}")
    endif()
endif()

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
    if(USE_MPI AND USE_OPENMP)
      target_link_libraries(${TESTEXE}  ${MPI_Fortran_LIBRARIES})
      set_target_properties(${TESTEXE} PROPERTIES
        COMPILE_FLAGS  "${EXTRA_FFLAGS}"
        LINK_FLAGS "${OpenMP_Fortran_FLAGS} ${MPI_Fortran_LINK_FLAGS}"
        LANGUAGE Fortran)
    elseif(USE_MPI)
      target_link_libraries(${TESTEXE}  ${MPI_Fortran_LIBRARIES})
      set_target_properties(${TESTEXE} PROPERTIES
        COMPILE_FLAGS  "${EXTRA_FFLAGS}"
        LINK_FLAGS "${MPI_Fortran_LINK_FLAGS}"
        LANGUAGE Fortran)
    elseif(USE_OPENMP)
      set_target_properties(${TESTEXE} PROPERTIES
          COMPILE_FLAGS "${EXTRA_FFLAGS}"
          LINK_FLAGS "${OpenMP_Fortran_FLAGS}"
          LANGUAGE Fortran)
    endif()
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
