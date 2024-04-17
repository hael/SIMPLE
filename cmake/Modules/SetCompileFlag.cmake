#############################################################################
# Given a list of flags, this function will try each, one at a time,
# and choose the first flag that works.  If no flags work, then nothing
# will be set, unless the REQUIRED key is given, in which case an error
# will be given.
#
# Call is:
# SET_COMPILE_FLAG(FLAGVAR FLAGVAL (Fortran|C|CXX) <REQUIRED> flag1 flag2...)
#
# For example, if you have the flag CMAKE_C_FLAGS and you want to add
# warnings and want to fail if this is not possible, you might call this
# function in this manner:
# SET_COMPILE_FLAGS(CMAKE_C_FLAGS "${CMAKE_C_FLAGS}" C REQUIRED
#                   "-Wall"     # GNU
#                   "-warn all" # Intel
#                  )
# The optin "-Wall" will be checked first, and if it works, will be
# appended to the CMAKE_C_FLAGS variable.  If it doesn't work, then
# "-warn all" will be tried.  If this doesn't work then checking will
# terminate because REQUIRED was given.
#
# The reasong that the variable must be given twice (first as the name then
# as the value in quotes) is because of the way CMAKE handles the passing
# of variables in functions; it is difficult to extract a variable's
# contents and assign new values to it from within a function.
#############################################################################

include(${CMAKE_ROOT}/Modules/CheckCCompilerFlag.cmake)
include(${CMAKE_ROOT}/Modules/CheckCXXCompilerFlag.cmake)
# CMakeCheckCompilerFlagCommonPatterns.cmake
set(FC_flagcheck FALSE)
if(EXISTS ${CMAKE_ROOT}/Modules/CheckFortranCompilerFlag.cmake)
  include(${CMAKE_ROOT}/Modules/CheckFortranCompilerFlag.cmake)
  set(FC_flagcheck TRUE)
endif()

function(SET_COMPILE_FLAG FLAGVAR FLAGVAL LANG)
  # Do some up front setup if Fortran
  if(LANG STREQUAL "Fortran")
    # Create a list of error messages from compilers
    set(FAIL_REGEX
      "ignoring unknown option"             # Intel
      "invalid argument"                    # Intel
      "unrecognized .*option"               # GNU
      "[Uu]nknown switch"                   # Portland Group
      "ignoring unknown option"             # MSVC
      "warning D9002"                       # MSVC, any lang
      "[Uu]nknown option"                   # HP
      "[Ww]arning: [Oo]ption"               # SunPro
      "command option .* is not recognized" # XL
      )
  endif(LANG STREQUAL "Fortran")

  # Make a variable holding the flags.  Filter out REQUIRED if it is there
  set(FLAG_REQUIRED FALSE)
  set(FLAG_FOUND FALSE)
  unset(FLAGLIST)
  foreach(var ${ARGN})
    STRING(TOUPPER "${var}" UP)
    if(UP STREQUAL "REQUIRED")
      set(FLAG_REQUIRED TRUE)
    else()
      set(FLAGLIST ${FLAGLIST} "${var}")
    endif(UP STREQUAL "REQUIRED")
  endforeach(var ${ARGN})

  # Now, loop over each flag
  foreach(flag ${FLAGLIST})
    list(FIND ${FLAGVAR} "${flag}" FLAG_EXISTS_ALREADY)
    if(NOT FLAG_EXISTS_ALREADY)
      unset(FLAG_WORKS)
      # Check the flag for the given language
      if(LANG STREQUAL "C")
        CHECK_C_COMPILER_FLAG("${flag}" FLAG_WORKS)
      elseif(LANG STREQUAL "CXX")
        CHECK_CXX_COMPILER_FLAG("${flag}" FLAG_WORKS)
      elseif(LANG STREQUAL "Fortran")
	if(FC_flagcheck)
  	  CHECK_Fortran_COMPILER_FLAG("${flag}" FLAG_WORKS)
	else()
          # There is no nice function to do this for FORTRAN, so we must manually
          # create a test program and check if it compiles with a given flag.
          set(TESTFILE "${CMAKE_BINARY_DIR}${CMAKE_FILES_DIRECTORY}")
          set(TESTFILE "${TESTFILE}/CMakeTmp/testFortranFlags.f90")
          file(WRITE "${TESTFILE}"
            "
            program dummyprog
              i = 5
            end program dummyprog
            ")
          try_compile(FLAG_WORKS ${CMAKE_BINARY_DIR} ${TESTFILE}
            COMPILE_DEFINITIONS "${flag}" OUTPUT_VARIABLE OUTPUT)
          # Check that the output message doesn't match any errors
          foreach(rx ${FAIL_REGEX})
            if("${OUTPUT}" MATCHES "${rx}")
              set(FLAG_WORKS FALSE)
            endif("${OUTPUT}" MATCHES "${rx}")
          endforeach(rx ${FAIL_REGEX})
        endif(FC_flagcheck)
      else()
        message(FATAL_ERROR "Unknown language in SET_COMPILE_FLAGS: ${LANG}")
      endif(LANG STREQUAL "C")

      # If this worked, use these flags, otherwise use other flags
      if(FLAG_WORKS)
        # Append this flag to the end of the list that already exists
        set(${FLAGVAR} "${FLAGVAL} ${flag}" CACHE STRING
          "Set the ${FLAGVAR} flags" FORCE)
        set(FLAG_FOUND TRUE)
        break() # We found something that works, so exit
      endif(FLAG_WORKS)

    else()
      set(FLAG_FOUND TRUE)
      break()
    endif()
  endforeach(flag ${FLAGLIST})

  # Raise an error if no flag was found
  if(FLAG_REQUIRED AND NOT FLAG_FOUND)
    message(FATAL_ERROR "No compile flags were found")
  endif(FLAG_REQUIRED AND NOT FLAG_FOUND)

endfunction()

function(SET_PREPROCESSOR_FLAG FLAGVAR FLAGVAL LANG)

  # Do some up front setup if Fortran
  if(LANG STREQUAL "Fortran")
    # Create a list of error messages from compilers
    set(FAIL_REGEX
      "ignoring unknown option"             # Intel
      "invalid argument"                    # Intel
      "unrecognized .*option"               # GNU
      "[Uu]nknown switch"                   # Portland Group
      "ignoring unknown option"             # MSVC
      "warning D9002"                       # MSVC, any lang
      "[Uu]nknown option"                   # HP
      "[Ww]arning: [Oo]ption"               # SunPro
      "command option .* is not recognized" # XL
      )
  endif(LANG STREQUAL "Fortran")

  # Make a variable holding the flags.  Filter out REQUIRED if it is there
  set(FLAG_REQUIRED FALSE)
  set(FLAG_FOUND FALSE)
  unset(FLAGLIST)
  foreach(var ${ARGN})
    STRING(TOUPPER "${var}" UP)
    if(UP STREQUAL "REQUIRED")
      set(FLAG_REQUIRED TRUE)
    else()
      set(FLAGLIST ${FLAGLIST} "${var}")
    endif(UP STREQUAL "REQUIRED")
  endforeach(var ${ARGN})

  # Now, loop over each flag
  foreach(flag ${FLAGLIST})

    unset(FLAG_WORKS)

    # Check the flag for the given language
    if(LANG STREQUAL "C")
      CHECK_C_COMPILER_FLAG("${flag}" FLAG_WORKS)
    elseif(LANG STREQUAL "CXX")
      CHECK_CXX_COMPILER_FLAG("${flag}" FLAG_WORKS)
    elseif(LANG STREQUAL "Fortran")
      # There is no nice function to do this for FORTRAN, so we must manually
      # create a test program and check if it compiles with a given flag.
      set(TESTFILE "${CMAKE_BINARY_DIR}${CMAKE_FILES_DIRECTORY}")
      set(TESTFILE "${TESTFILE}/CMakeTmp/testFortranCPPFlags.F08")
      file(WRITE "${TESTFILE}"
        "
        #define c99_count(...)    _c99_count1 ( , ##__VA_ARGS__)/* */
        #define _c99_count1(...)  _c99_count2 (__VA_ARGS__,10,9,8,7,6,5,4,3,2,1,0)
        #define _c99_count2(_,x0,x1,x2,x3,x4,x5,x6,x7,x8,x9,n,...) n
        program dummyprog
         integer i
         integer,parameter :: nv=c99_count (__VA_ARGS__);
         character(255)::p_tokens= #__VA_ARGS__ ;
         i = 5
        end program dummyprog
        ")
      try_compile(FLAG_WORKS ${CMAKE_BINARY_DIR} ${TESTFILE}
        COMPILE_DEFINITIONS "${flag}" OUTPUT_VARIABLE OUTPUT)
      # Check that the output message doesn't match any errors
      foreach(rx ${FAIL_REGEX})
        if("${OUTPUT}" MATCHES "${rx}")
          set(FLAG_WORKS FALSE)
        endif("${OUTPUT}" MATCHES "${rx}")
      endforeach(rx ${FAIL_REGEX})
    else()
      message(FATAL_ERROR "Unknown language in SET_PREPROCESSOR_FLAG: ${LANG}")
    endif(LANG STREQUAL "C")

    # If this worked, use these flags, otherwise use other flags
    if(FLAG_WORKS)
      # Append this flag to the end of the list that already exists
      set(${FLAGVAR} "${FLAGVAL} ${flag}" CACHE STRING
        "Set the ${FLAGVAR} flags" FORCE)
      set(FLAG_FOUND TRUE)
      break() # We found something that works, so exit
    endif(FLAG_WORKS)

  endforeach(flag ${FLAGLIST})

  # Raise an error if no flag was found
  if(FLAG_REQUIRED AND NOT FLAG_FOUND)
    message(FATAL_ERROR "No compile flags were found")
  endif(FLAG_REQUIRED AND NOT FLAG_FOUND)

endfunction()
