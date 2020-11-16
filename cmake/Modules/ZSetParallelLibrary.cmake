# Turns on either OpenMP or MPI
set (OMP_NUM_PROCS 1 CACHE
  STRING "Number of processors OpenMP may use" FORCE)

unset (OpenMP_Fortran_FLAGS CACHE)
unset (GOMP_Fortran_LINK_FLAGS CACHE)
unset (MPI_FOUND CACHE)
unset (MPI_COMPILER CACHE)
unset (MPI_LIBRARY CACHE)
if (USE_OPENMP)
  # Find OpenMP
  include(ProcessorCount)
  ProcessorCount(N)
  if(NOT N EQUAL 0)
    set (OMP_NUM_PROCS ${N})
    set (OMP_NUM_PROCS_INTERNAL ${N})
  endif()
  message(STATUS "In ZSetParallelLibrary -- OMP_NUM_PROCS ${OMP_NUM_PROCS}")
  if (NOT OpenMP_Fortran_FLAGS)
    find_package (OpenMP_Fortran)
    if (NOT OpenMP_Fortran_FLAGS)
      message (FATAL_ERROR "Fortran compiler does not support OpenMP")
    endif (NOT OpenMP_Fortran_FLAGS)

    if(OpenMP_Fortran_VERSION)
       # set(OpenMP_Fortran_FLAGS "${OpenMP_Fortran_FLAGS} -DOPENMP_VERSION=${OpenMP_Fortran_VERSION}")
    endif(OpenMP_Fortran_VERSION)
  endif (NOT OpenMP_Fortran_FLAGS)
endif(USE_OPENMP)

if (USE_MPI)
  # Find MPI when FC != mpif90
  if (NOT MPI_Fortran_FOUND)
    find_package (MPI REQUIRED)
  endif (NOT MPI_Fortran_FOUND)
endif (USE_MPI)
