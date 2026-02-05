# cmake/AutoDetectGCC.cmake
#
# Auto-detect a suitable GCC/GFortran toolchain if user has not set
# CMAKE_*_COMPILER explicitly. Must be included BEFORE project(...).
# If a compiler is already set (via cache or -D on command line), don't touch it
if(DEFINED CMAKE_Fortran_COMPILER AND NOT CMAKE_Fortran_COMPILER STREQUAL "")
    message(STATUS "AutoDetectGCC: CMAKE_Fortran_COMPILER already set to ${CMAKE_Fortran_COMPILER}")
    return()
endif()

# ---------------------------------------------------------------------------
# 1. Scan for candidate gfortran binaries
# ---------------------------------------------------------------------------
set(_gfortran_candidates
    gfortran-16
    gfortran-15
    gfortran-14
    gfortran
)
set(_best_gfortran "")
set(_best_version "")
foreach(_cand IN LISTS _gfortran_candidates)
    find_program(_fc_prog NAMES ${_cand})
    if(NOT _fc_prog)
        continue()
    endif()
    execute_process(
        COMMAND "${_fc_prog}" -dumpversion
        OUTPUT_VARIABLE _fc_ver
        ERROR_QUIET
        OUTPUT_STRIP_TRAILING_WHITESPACE
    )
    if(NOT _fc_ver)
        continue()
    endif()
    message(STATUS "AutoDetectGCC: candidate ${_fc_prog} (version ${_fc_ver})")
    # Require at least 14.2 here
    if(_fc_ver VERSION_LESS "14.2")
        message(STATUS "AutoDetectGCC: skipping ${_fc_prog} since version ${_fc_ver} is too old")
        continue()
    endif()
    if(NOT _best_gfortran OR _fc_ver VERSION_GREATER _best_version)
        set(_best_gfortran "${_fc_prog}")
        set(_best_version "${_fc_ver}")
    endif()
endforeach()

if(NOT _best_gfortran)
    message(STATUS "AutoDetectGCC: no suitable gfortran >= 8.5 found; CMake will use its default.")
    return()
endif()

message(STATUS "AutoDetectGCC: selected gfortran = ${_best_gfortran} (version ${_best_version})")

# ---------------------------------------------------------------------------
# 2. Derive gcc/g++ from same directory
# ---------------------------------------------------------------------------
get_filename_component(_fc_dir "${_best_gfortran}" DIRECTORY)
set(_gcc_candidates gcc-16 gcc-15 gcc-14 gcc)
set(_gxx_candidates g++-16 g++-15 g++-14 g++)
set(_best_gcc "")
foreach(_cand IN LISTS _gcc_candidates)
    find_program(_gcc_prog NAMES ${_cand} HINTS "${_fc_dir}")
    if(_gcc_prog)
        set(_best_gcc "${_gcc_prog}")
        break()
    endif()
endforeach()
set(_best_gxx "")
foreach(_cand IN LISTS _gxx_candidates)
    find_program(_gxx_prog NAMES ${_cand} HINTS "${_fc_dir}")
    if(_gxx_prog)
        set(_best_gxx "${_gxx_prog}")
        break()
    endif()
endforeach()

# ---------------------------------------------------------------------------
# 3. Set CMAKE_*_COMPILER cache entries BEFORE project()
# ---------------------------------------------------------------------------
set(CMAKE_Fortran_COMPILER "${_best_gfortran}" CACHE FILEPATH "Auto-detected gfortran" FORCE)
if(_best_gcc)
    set(CMAKE_C_COMPILER "${_best_gcc}" CACHE FILEPATH "Auto-detected gcc" FORCE)
endif()
if(_best_gxx)
    set(CMAKE_CXX_COMPILER "${_best_gxx}" CACHE FILEPATH "Auto-detected g++" FORCE)
endif()
message(STATUS "AutoDetectGCC: final choice:")
message(STATUS "  CMAKE_Fortran_COMPILER = ${CMAKE_Fortran_COMPILER}")
message(STATUS "  CMAKE_C_COMPILER       = ${CMAKE_C_COMPILER}")
message(STATUS "  CMAKE_CXX_COMPILER     = ${CMAKE_CXX_COMPILER}")
