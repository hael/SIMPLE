# Find GNU/OpenCoarrays support for -fcoarray=lib builds.
#
# Provides:
#   Coarray_FOUND
#   Coarray::Coarray

find_package(MPI QUIET COMPONENTS Fortran)

find_program(Coarray_CAF_EXECUTABLE NAMES caf)
find_program(Coarray_CAFRUN_EXECUTABLE NAMES cafrun)

set(_Coarray_LIBRARY_HINTS)
set(_Coarray_CAF_MPI_LIBRARY_CANDIDATE)
if(Coarray_CAF_EXECUTABLE)
    foreach(_caf_show_arg IN ITEMS --show -show)
        execute_process(
            COMMAND "${Coarray_CAF_EXECUTABLE}" "${_caf_show_arg}"
            OUTPUT_VARIABLE _caf_show
            ERROR_VARIABLE _caf_show_error
            RESULT_VARIABLE _caf_show_result
            OUTPUT_STRIP_TRAILING_WHITESPACE
            ERROR_STRIP_TRAILING_WHITESPACE
        )
        string(STRIP "${_caf_show} ${_caf_show_error}" _caf_show)
        if(_caf_show_result EQUAL 0 AND _caf_show)
            break()
        endif()
    endforeach()
    if(_caf_show)
        separate_arguments(_caf_show_args UNIX_COMMAND "${_caf_show}")
        set(_Coarray_NEXT_IS_LIBRARY_DIR OFF)
        foreach(_caf_arg IN LISTS _caf_show_args)
            if(_Coarray_NEXT_IS_LIBRARY_DIR)
                list(APPEND _Coarray_LIBRARY_HINTS "${_caf_arg}")
                set(_Coarray_NEXT_IS_LIBRARY_DIR OFF)
            elseif(_caf_arg STREQUAL "-L")
                set(_Coarray_NEXT_IS_LIBRARY_DIR ON)
            elseif(_caf_arg MATCHES "^-L(.+)")
                list(APPEND _Coarray_LIBRARY_HINTS "${CMAKE_MATCH_1}")
            elseif(_caf_arg MATCHES "libcaf_mpi\\.(a|so|dylib)(\\.[0-9.]+)?$" AND EXISTS "${_caf_arg}")
                set(_Coarray_CAF_MPI_LIBRARY_CANDIDATE "${_caf_arg}")
            endif()
        endforeach()
    endif()
endif()

if(_Coarray_CAF_MPI_LIBRARY_CANDIDATE AND NOT Coarray_CAF_MPI_LIBRARY)
    set(Coarray_CAF_MPI_LIBRARY "${_Coarray_CAF_MPI_LIBRARY_CANDIDATE}" CACHE FILEPATH "OpenCoarrays caf_mpi library" FORCE)
endif()

find_library(Coarray_CAF_MPI_LIBRARY
    NAMES caf_mpi libcaf_mpi
    HINTS ${_Coarray_LIBRARY_HINTS}
)

set(Coarray_FOUND FALSE)
if(Coarray_CAF_MPI_LIBRARY AND MPI_Fortran_FOUND)
    set(Coarray_FOUND TRUE)
elseif(NOT Coarray_FIND_QUIETLY)
    if(NOT MPI_Fortran_FOUND)
        message(STATUS "Could NOT find Coarray: missing MPI Fortran")
    elseif(NOT Coarray_CAF_MPI_LIBRARY)
        message(STATUS "Could NOT find Coarray: missing caf_mpi/libcaf_mpi")
    endif()
endif()

if(Coarray_FIND_REQUIRED AND NOT Coarray_FOUND)
    message(FATAL_ERROR "Could NOT find Coarray (requires caf_mpi/libcaf_mpi and MPI Fortran)")
endif()

if(Coarray_FOUND AND NOT TARGET Coarray::Coarray)
    add_library(Coarray::Coarray INTERFACE IMPORTED)
    target_link_libraries(Coarray::Coarray INTERFACE
        "${Coarray_CAF_MPI_LIBRARY}"
        MPI::MPI_Fortran
    )
    target_compile_options(Coarray::Coarray INTERFACE
        $<$<COMPILE_LANGUAGE:Fortran>:-fcoarray=lib>
    )
    target_link_options(Coarray::Coarray INTERFACE
        $<$<LINK_LANGUAGE:Fortran>:-fcoarray=lib>
    )
endif()

mark_as_advanced(
    Coarray_CAF_EXECUTABLE
    Coarray_CAFRUN_EXECUTABLE
    Coarray_CAF_MPI_LIBRARY
)
