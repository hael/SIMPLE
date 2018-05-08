# Add the source files
add_executable(${SIMPLE_EXE} "${CMAKE_CURRENT_SOURCE_DIR}/simple/simple_exec/simple_exec.f90")
add_executable(${SIMPLE_DISTR_EXE} "${CMAKE_CURRENT_SOURCE_DIR}/simple/simple_distr_exec/simple_distr_exec.f90")
add_executable(${SIMPLE_PRIVATE_EXE} "${CMAKE_CURRENT_SOURCE_DIR}/simple/simple_private_exec/simple_private_exec.f90")

# This links foo to the bar library
target_link_libraries(${SIMPLE_EXE} LINK_PUBLIC ${EXTRA_LIBS})
target_link_libraries(${SIMPLE_DISTR_EXE} LINK_PUBLIC ${EXTRA_LIBS})
target_link_libraries(${SIMPLE_PRIVATE_EXE} LINK_PUBLIC ${EXTRA_LIBS})

IF(USE_OPENMPI AND USE_MPI)
  set_target_properties(${SIMPLE_EXE}
    PROPERTIES
      COMPILE_FLAGS "${OpenMP_Fortran_FLAGS} ${MPI_Fortran_FLAGS}"
      LINK_FLAGS "${OpenMP_Fortran_FLAGS} ${MPI_Fortran_LINK_FLAGS}"
      JOB_POOL_COMPILE NUM_JOBS)
  set_target_properties(${SIMPLE_DISTR_EXE}
    PROPERTIES
      COMPILE_FLAGS "${OpenMP_Fortran_FLAGS} ${MPI_Fortran_FLAGS}"
      LINK_FLAGS "${OpenMP_Fortran_FLAGS} ${MPI_Fortran_LINK_FLAGS}"
      JOB_POOL_COMPILE NUM_JOBS)
  set_target_properties(${SIMPLE_PRIVATE_EXE}
    PROPERTIES
      COMPILE_FLAGS "${OpenMP_Fortran_FLAGS} ${MPI_Fortran_FLAGS}"
      LINK_FLAGS "${OpenMP_Fortran_FLAGS} ${MPI_Fortran_LINK_FLAGS}"
      JOB_POOL_COMPILE NUM_JOBS)

target_link_libraries(${SIMPLE_EXE} LINK_PUBLIC  ${MPI_Fortran_LIBRARIES})
target_link_libraries(${SIMPLE_DISTR_EXE} LINK_PUBLIC  ${MPI_Fortran_LIBRARIES})
target_link_libraries(${SIMPLE_PRIVATE_EXE} LINK_PUBLIC  ${MPI_Fortran_LIBRARIES})
ELSEIF(USE_OPENMP)
  set_target_properties(${SIMPLE_EXE}
    PROPERTIES
      COMPILE_FLAGS "${OpenMP_Fortran_FLAGS}"
      LINK_FLAGS "${OpenMP_Fortran_FLAGS}"
      JOB_POOL_COMPILE NUM_JOBS)
  set_target_properties(${SIMPLE_DISTR_EXE}
    PROPERTIES
      COMPILE_FLAGS "${OpenMP_Fortran_FLAGS}"
      LINK_FLAGS "${OpenMP_Fortran_FLAGS}"
      JOB_POOL_COMPILE NUM_JOBS)
  set_target_properties(${SIMPLE_PRIVATE_EXE}
  PROPERTIES
    COMPILE_FLAGS "${OpenMP_Fortran_FLAGS}"
    LINK_FLAGS "${OpenMP_Fortran_FLAGS}"
    JOB_POOL_COMPILE NUM_JOBS)
ELSEIF(USE_MPI)
     set_target_properties(${SIMPLE_EXE}
    PROPERTIES
      COMPILE_FLAGS "${MPI_Fortran_FLAGS}"
      LINK_FLAGS "${MPI_Fortran_LINK_FLAGS}"
      JOB_POOL_COMPILE NUM_JOBS)
  set_target_properties(${SIMPLE_DISTR_EXE}
    PROPERTIES
      COMPILE_FLAGS "${MPI_Fortran_FLAGS}"
      LINK_FLAGS "${MPI_Fortran_LINK_FLAGS}"
      JOB_POOL_COMPILE NUM_JOBS)
  set_target_properties(${SIMPLE_PRIVATE_EXE}
  PROPERTIES
    COMPILE_FLAGS "${MPI_Fortran_FLAGS}"
    LINK_FLAGS "${MPI_Fortran_LINK_FLAGS}"
    JOB_POOL_COMPILE NUM_JOBS)

target_link_libraries(${SIMPLE_EXE} LINK_PUBLIC  ${MPI_Fortran_LIBRARIES})
target_link_libraries(${SIMPLE_DISTR_EXE} LINK_PUBLIC  ${MPI_Fortran_LIBRARIES})
target_link_libraries(${SIMPLE_PRIVATE_EXE} LINK_PUBLIC  ${MPI_Fortran_LIBRARIES})
ENDIF()

install(TARGETS "${SIMPLE_EXE}"
  ARCHIVE DESTINATION "${CMAKE_INSTALL_BINDIR}"
  RUNTIME DESTINATION "${CMAKE_INSTALL_BINDIR}"
         CONFIGURATIONS DEBUG|RELEASE
)
install(TARGETS "${SIMPLE_DISTR_EXE}"
  ARCHIVE DESTINATION "${CMAKE_INSTALL_BINDIR}"
  RUNTIME DESTINATION "${CMAKE_INSTALL_BINDIR}"
  CONFIGURATIONS DEBUG|RELEASE
)

install(TARGETS "${SIMPLE_PRIVATE_EXE}"
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
if (CMAKE_Fortran_COMPILER_ID STREQUAL "PGI")
  list(REMOVE_ITEM TESTDIRS "tests/test_pgi")
else()
  list(REMOVE_ITEM TESTDIRS "tests/test_pgi"
    "tests/test_cufft")
endif()

if (CMAKE_Fortran_COMPILER_ID STREQUAL "Intel")
  list(REMOVE_ITEM TESTDIRS "tests/test_openacc"
   "tests/test_laplace_solver"
   "tests/test_accuracy")
else()
  list(REMOVE_ITEM TESTDIRS "tests/test_intel")
endif()


if(NOT USE_CUDA)
 list(REMOVE_ITEM TESTDIRS "tests/test_cuda")
endif()
if(NOT USE_OPENMP)
 list(REMOVE_ITEM TESTDIRS "tests/test_omp")
endif()
if(NOT USE_MPI)
 list(REMOVE_ITEM TESTDIRS "tests/test_mpi")
endif()

# if(USE_OPENACC)
#   message(STATUS "Searching for BLAS to test OpenACC")
#   set(BLA_F95 ON)
#   set(BLA_VENDOR "ATLAS")
#   find_package(BLAS QUIET)
#   if (BLAS_FOUND AND USE_OPENACC)
#     message(STATUS "Adding BLAS library to tests: L${BLAS_LINKER_FLAGS} ${BLAS_blas_LIBRARY}")
#     endif()
# endif()


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
        COMPILE_FLAGS  "${OpenMP_Fortran_FLAGS} ${MPI_Fortran_FLAGS} -fno-range-check"
        LINK_FLAGS "${OpenMP_Fortran_FLAGS} ${MPI_Fortran_LINK_FLAGS}"
        LANGUAGE Fortran)
    elseif(USE_MPI)
      target_link_libraries(${TESTEXE}  ${MPI_Fortran_LIBRARIES})
      set_target_properties(${TESTEXE} PROPERTIES
        COMPILE_FLAGS  "${MPI_Fortran_FLAGS} -fno-range-check"
        LINK_FLAGS "${MPI_Fortran_LINK_FLAGS}"
        LANGUAGE Fortran)
    elseif(USE_OPENMP)
      if(${CMAKE_Fortran_COMPILER_ID} STREQUAL "PGI")
	      set_target_properties(${TESTEXE} PROPERTIES
          COMPILE_FLAGS "${OpenMP_Fortran_FLAGS}"
          LINK_FLAGS "${OpenMP_Fortran_FLAGS}"
          LANGUAGE Fortran)
      else()
	      set_target_properties(${TESTEXE} PROPERTIES
          COMPILE_FLAGS "${OpenMP_Fortran_FLAGS} -fno-range-check"
          LINK_FLAGS "${OpenMP_Fortran_FLAGS}"
          LANGUAGE Fortran)
      endif()
    endif()
    if(USE_CUDA)
      file(GLOB_RECURSE TESTCUDASOURCES LIST_DIRECTORIES true
          "tests/${TESTNAME}/*.cu")
       if("" STREQUAL "${TESTCUDASOURCES}")
         # message(STATUS " skipping ${TESTNAME} -- no files in dir")
       else()
         CUDA_COMPILE( ${TESTEXE}_cuda-kernels ${TESTCUDASOURCES}
           OPTIONS  OPTIONS -DFLAG=2 "-DFLAG_OTHER=space in flag"
         DEBUG -g
         RELEASE --use_fast_math
         RELWITHDEBINFO --use_fast_math;-g
         MINSIZEREL --use_fast_math)
        target_link_libraries(${TESTEXE} ${TESTEXE}_generated_files)
        endif()
      target_link_libraries(${TESTEXE} cuda_kernels)
    endif()
    # if(BLAS_FOUND AND USE_OPENACC)
    #   message(STATUS "Adding  ${BLAS_blas_LIBRARY} library to ${TESTEXE}")
    #   target_link_libraries(${TESTEXE} ${BLAS_blas_LIBRARY})
    #   set_target_properties(${TESTEXE} PROPERTIES
    #     COMPILE_FLAGS "${OpenMP_Fortran_FLAGS} -fno-range-check"
    #     LINK_FLAGS "${OpenMP_Fortran_FLAGS} -L${BLAS_LINKER_FLAGS}"
    #     LANGUAGE Fortran)
    # else()
   # endif()
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
  WORKING_DIRECTORY ${CMAKE_INSTALL_BINDIR}
  COMMAND source add2.bashrc &&  simple_test_parse nthr=4 verbose=yes)

add_test(NAME test_install
  WORKING_DIRECTORY ${CMAKE_INSTALL_BINDIR}
  COMMAND source add2.bashrc  && simple_test_install nthr=4 verbose=yes)


add_test(NAME test_units
  WORKING_DIRECTORY ${CMAKE_INSTALL_BINDIR}
  COMMAND source add2.bashrc &&  simple_test_units nthr=4)

add_test(NAME test_timer
  WORKING_DIRECTORY ${CMAKE_INSTALL_BINDIR}
  COMMAND source add2.bashrc && ${CMAKE_COMMAND} -E simple_test_timer verbose=yes debug=yes)


add_test(NAME YeahhhTest COMMAND ${CMAKE_COMMAND} -E echo 'Yeahh')
