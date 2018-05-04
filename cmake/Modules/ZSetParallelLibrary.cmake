# Turns on either OpenMP or MPI
# If both are requested, the other is disabled
# When one is turned on, the other is turned off
# If both are off, we explicitly disable them just in case
# Turn off both OpenMP and MPI
set (OMP_NUM_PROCS 0 CACHE
  STRING "Number of processors OpenMP may use" FORCE)
unset (OpenMP_Fortran_FLAGS CACHE)
unset (GOMP_Fortran_LINK_FLAGS CACHE)
unset (MPI_FOUND CACHE)
unset (MPI_COMPILER CACHE)
unset (MPI_LIBRARY CACHE)
if (USE_OPENMP)
  # Find OpenMP
  if (NOT OpenMP_Fortran_FLAGS)
    find_package (OpenMP_Fortran)
    if (NOT OpenMP_Fortran_FLAGS)
      message (FATAL_ERROR "Fortran compiler does not support OpenMP")
    endif (NOT OpenMP_Fortran_FLAGS)
  endif (NOT OpenMP_Fortran_FLAGS)
endif()

if (USE_MPI)
  # Find MPI
  if (NOT MPI_Fortran_FOUND)
    find_package (MPI REQUIRED)
  endif (NOT MPI_Fortran_FOUND)

endif (USE_MPI)
