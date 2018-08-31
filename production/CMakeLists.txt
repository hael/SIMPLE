# Add the source files
add_executable(${SIMPLE_EXE} "${CMAKE_CURRENT_SOURCE_DIR}/simple/simple_exec/simple_exec.f90")
add_executable(${SIMPLE_DISTR_EXE} "${CMAKE_CURRENT_SOURCE_DIR}/simple/simple_distr_exec/simple_distr_exec.f90")
add_executable(${SIMPLE_PRIVATE_EXE} "${CMAKE_CURRENT_SOURCE_DIR}/simple/simple_private_exec/simple_private_exec.f90")
set(SIMPLE_MAIN_TARGETS ${SIMPLE_EXE} ${SIMPLE_DISTR_EXE} ${SIMPLE_PRIVATE_EXE})
# This links foo to the bar library

#target_link_libraries(${SIMPLE_MAIN_TARGETS} LINK_PUBLIC ${EXTRA_LIBS})
target_link_libraries(${SIMPLE_EXE} LINK_PUBLIC ${EXTRA_LIBS})
target_link_libraries(${SIMPLE_DISTR_EXE} LINK_PUBLIC ${EXTRA_LIBS})
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
set(EXTRA_FFLAGS "")
set(DONT_CHECK_RANGE "-fno-range-check")
if(NOT USE_CUDA)
 list(REMOVE_ITEM TESTDIRS "tests/test_cuda")
 list(REMOVE_ITEM TESTDIRS "tests/test_cublas")
else()
 list(REMOVE_ITEM TESTDIRS "tests/test_img_export")
endif()
if(NOT USE_OPENMP)
 list(REMOVE_ITEM TESTDIRS "tests/test_omp")
endif()
if(NOT USE_MPI)
  list(REMOVE_ITEM TESTDIRS "tests/test_mpi"
    "tests/test_mpimp_benchmark")
endif()
 list(REMOVE_ITEM TESTDIRS "tests/test_mpimp_benchmark")

# if(USE_OPENACC)
#   message(STATUS "Searching for BLAS to test OpenACC")
#   set(BLA_F95 ON)
#   set(BLA_VENDOR "ATLAS")
#   find_package(BLAS QUIET)
#   if (BLAS_FOUND AND USE_OPENACC)
#     message(STATUS "Adding BLAS library to tests: L${BLAS_LINKER_FLAGS} ${BLAS_blas_LIBRARY}")
#     endif()
# endif()


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
    if(USE_CUDA)
      file(GLOB_RECURSE TESTCUDASOURCES LIST_DIRECTORIES true
          "tests/${TESTNAME}/*.cu")
       if("" STREQUAL "${TESTCUDASOURCES}")
         # message(STATUS " skipping ${TESTNAME} -- no files in dir")
       else()
         CUDA_COMPILE( ${TESTEXE}_cuda_kernels ${TESTCUDASOURCES}
           OPTIONS
         DEBUG -g
         RELEASE --use_fast_math
         RELWITHDEBINFO --use_fast_math;-g
         MINSIZEREL --use_fast_math)
       target_link_libraries(${TESTEXE} ${TESTEXE}_cuda_kernels)
     endif()
     CUDA_ADD_CUBLAS_TO_TARGET(${TESTEXE})
     target_link_libraries(${TESTEXE} cuda_kernels stdc++)
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
  simple_test_ansi_colors
  simple_test_binoris
  simple_test_cartcorr_sanity
  simple_test_corrs2weights
  simple_test_extr_frac
  simple_test_fileio
  simple_test_inside_write
  simple_test_imgfile
  simple_test_list_files
  simple_test_maxnloc
  simple_test_multinomal
  simple_test_omp
  simple_test_opeaksio
  simple_test_order_corr
  simple_test_parse
  simple_test_phaseplate_correct_fsc
  simple_test_simd
  simple_test_sp_project
  simple_test_sym
  simple_test_syslib
  simple_test_timer
  simple_test_units )


foreach(CHECKUNIT ${CTESTDIRS})
add_test(NAME ${CHECKUNIT}
  WORKING_DIRECTORY ${CMAKE_BINARY_DIR}/Testing
  COMMAND ../bin/${CHECKUNIT}  nthr=4 verbose=yes)
add_dependencies( check ${CHECKUNIT} )
endforeach()
