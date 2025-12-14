# Compiler options for GCC/GFortran only
# Modern CMake 3.12+ approach

# Verify we're using GCC/GFortran
if(NOT CMAKE_Fortran_COMPILER_ID STREQUAL "GNU")
    message(FATAL_ERROR "Only GFortran is supported. Found: ${CMAKE_Fortran_COMPILER_ID}")
endif()

# On macOS, allow AppleClang for C/C++ if GFortran is being used
# GFortran needs compatible C/C++ compilers and AppleClang works fine
if(APPLE)
    if(CMAKE_C_COMPILER_ID MATCHES "Clang|AppleClang")
        message(STATUS "Using ${CMAKE_C_COMPILER_ID} for C (compatible with GFortran on macOS)")
    elseif(NOT CMAKE_C_COMPILER_ID STREQUAL "GNU")
        message(WARNING "C compiler is ${CMAKE_C_COMPILER_ID}. GCC or AppleClang recommended.")
    endif()
    
    if(CMAKE_CXX_COMPILER_ID MATCHES "Clang|AppleClang")
        message(STATUS "Using ${CMAKE_CXX_COMPILER_ID} for C++ (compatible with GFortran on macOS)")
    elseif(NOT CMAKE_CXX_COMPILER_ID STREQUAL "GNU")
        message(WARNING "C++ compiler is ${CMAKE_CXX_COMPILER_ID}. G++ or AppleClang recommended.")
    endif()
else()
    # On Linux, enforce GCC for C/C++
    if(NOT CMAKE_C_COMPILER_ID STREQUAL "GNU")
        message(FATAL_ERROR "Only GCC is supported for C. Found: ${CMAKE_C_COMPILER_ID}")
    endif()
    
    if(NOT CMAKE_CXX_COMPILER_ID STREQUAL "GNU")
        message(FATAL_ERROR "Only G++ is supported for C++. Found: ${CMAKE_CXX_COMPILER_ID}")
    endif()
endif()

# Require minimum GFortran version
if(CMAKE_Fortran_COMPILER_VERSION VERSION_LESS "4.9")
    message(FATAL_ERROR "GFortran 4.9 or higher required. Found: ${CMAKE_Fortran_COMPILER_VERSION}")
endif()

# Platform-specific definitions
if(APPLE)
    add_definitions(-D__FreeBSD__ -DMACOSX)
else()
    add_definitions(-D__linux__ -DLINUX)
endif()

# Fortran language features
add_definitions(-DUSE_F08=1)
add_definitions(-D_DEBUG=$<CONFIG:Debug>)
add_definitions(-DGNU)

# C/C++ standards
set(CMAKE_C_STANDARD 11)
set(CMAKE_C_STANDARD_REQUIRED ON)
set(CMAKE_CXX_STANDARD 14)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

# Platform-specific optimization
if(APPLE)
    # Handle Apple Silicon vs Intel
    execute_process(
        COMMAND sysctl -n machdep.cpu.brand_string
        OUTPUT_VARIABLE APPLE_CPU
        OUTPUT_STRIP_TRAILING_WHITESPACE
    )
    if(APPLE_CPU MATCHES "Apple M")
        set(ARCH_FLAG "-mtune=generic")
    else()
        set(ARCH_FLAG "-mtune=native")
    endif()
else()
    set(ARCH_FLAG "-march=native")
endif()

# Common Fortran flags (used in both Debug and Release)
# Note: -D__FILENAME__ uses Make variable expansion for source file name
set(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} -cpp -ffree-form -fimplicit-none -ffree-line-length-none -fno-second-underscore -Wall -Waliasing -Wampersand -Wsurprising -Wline-truncation -D__FILENAME__='\"\\$(notdir \\$<)\"'")

# Release flags for Fortran
set(CMAKE_Fortran_FLAGS_RELEASE "-O3 -funroll-loops ${ARCH_FLAG} -fPIC")

# Debug flags for Fortran
set(CMAKE_Fortran_FLAGS_DEBUG "-O0 -g -fbacktrace -fbounds-check -fcheck=all -Wuninitialized -Wunused -fPIC")

# C flags
set(CMAKE_C_FLAGS_RELEASE "-O3 ${ARCH_FLAG} -fPIC")
set(CMAKE_C_FLAGS_DEBUG "-O0 -g -fPIC")

# Suppress Clang deployment version warning on macOS
if(APPLE AND CMAKE_C_COMPILER_ID MATCHES "Clang")
    set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -Wno-overriding-deployment-version")
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wno-overriding-deployment-version")
endif()

# C++ flags
set(CMAKE_CXX_FLAGS_RELEASE "-O3 ${ARCH_FLAG} -fPIC")
set(CMAKE_CXX_FLAGS_DEBUG "-O0 -g -fPIC")

# Default to Release build if not specified
if(NOT CMAKE_BUILD_TYPE)
    set(CMAKE_BUILD_TYPE Release CACHE STRING "Build type" FORCE)
    set_property(CACHE CMAKE_BUILD_TYPE PROPERTY STRINGS Debug Release RelWithDebInfo)
endif()

# RPATH settings for shared libraries
if(APPLE)
    set(CMAKE_INSTALL_RPATH "@executable_path/../lib;${CMAKE_INSTALL_PREFIX}/${CMAKE_INSTALL_LIBDIR}")
    set(CMAKE_INSTALL_RPATH_USE_LINK_PATH TRUE)
else()
    set(CMAKE_INSTALL_RPATH "$ORIGIN/../lib:${CMAKE_INSTALL_PREFIX}/${CMAKE_INSTALL_LIBDIR}")
    set(CMAKE_INSTALL_RPATH_USE_LINK_PATH TRUE)
endif()