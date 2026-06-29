# Find GNU/OpenCoarrays support for -fcoarray=lib builds.
#
# Provides:
#   Coarray_FOUND
#   Coarray::Coarray

find_package(MPI QUIET COMPONENTS Fortran)

find_program(Coarray_CAF_EXECUTABLE NAMES caf)

set(_Coarray_LIBRARY_HINTS)
set(_Coarray_LIBRARY_FROM_CAF)

if(Coarray_CAF_EXECUTABLE)
    foreach(_caf_show_arg IN ITEMS --show -show)
        execute_process(
            COMMAND "${Coarray_CAF_EXECUTABLE}" "${_caf_show_arg}"
            OUTPUT_VARIABLE _caf_show
            ERROR_QUIET
            RESULT_VARIABLE _caf_show_result
            OUTPUT_STRIP_TRAILING_WHITESPACE
        )
        if(_caf_show_result EQUAL 0 AND _caf_show)
            break()
        endif()
        unset(_caf_show)
    endforeach()

    if(_caf_show)
        separate_arguments(_caf_show_args UNIX_COMMAND "${_caf_show}")
        foreach(_caf_arg IN LISTS _caf_show_args)
            if(_caf_arg MATCHES "^-L(.+)")
                list(APPEND _Coarray_LIBRARY_HINTS "${CMAKE_MATCH_1}")
            elseif(_caf_arg MATCHES "^/.*/libcaf_[^/]+\\.(a|so|dylib)(\\.[0-9.]+)?$" AND EXISTS "${_caf_arg}")
                get_filename_component(_caf_lib_dir "${_caf_arg}" DIRECTORY)
                list(APPEND _Coarray_LIBRARY_HINTS "${_caf_lib_dir}")
                set(_Coarray_LIBRARY_FROM_CAF "${_caf_arg}")
            endif()
        endforeach()
    endif()
endif()

if(_Coarray_LIBRARY_FROM_CAF AND NOT Coarray_CAF_LIBRARY)
    set(Coarray_CAF_LIBRARY "${_Coarray_LIBRARY_FROM_CAF}" CACHE FILEPATH "OpenCoarrays runtime library" FORCE)
endif()

find_library(Coarray_CAF_LIBRARY
    NAMES caf_mpi caf_openmpi caf_mpich libcaf_mpi libcaf_openmpi libcaf_mpich
    HINTS ${_Coarray_LIBRARY_HINTS}
)

set(Coarray_FOUND FALSE)
if(Coarray_CAF_LIBRARY AND MPI_Fortran_FOUND)
    set(Coarray_FOUND TRUE)
endif()

if(Coarray_FOUND AND NOT TARGET Coarray::Coarray)
    add_library(Coarray::Coarray INTERFACE IMPORTED)
    target_link_libraries(Coarray::Coarray INTERFACE
        "${Coarray_CAF_LIBRARY}"
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
    Coarray_CAF_LIBRARY
)
