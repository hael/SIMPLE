
########################################
# Set up how to compile the tests
########################################
file(GLOB TESTDIRS
  LIST_DIRECTORIES true
  RELATIVE ${CMAKE_CURRENT_SOURCE_DIR}
  test_*
  )
message( STATUS "Test directories ${TESTDIRS}")



foreach(TESTNAME ${TESTDIRS})
#  message(STATUS " jumping to test ${TESTNAME} ")
  set(TESTEXE simple_${TESTNAME})
#  message(STATUS "  test exec ${TESTEXE} ")
  file(GLOB_RECURSE TESTSOURCES LIST_DIRECTORIES true
    "${TESTNAME}/*.[fF][90][0538]")
  if("" STREQUAL "${TESTSOURCES}")
    # message(STATUS " skipping ${TESTNAME} -- no files in dir")
  else()

#       list(LENGTH ${TESTSOURCES} NUMSOURCES)
#     foreach(SOURCE ${TESTSOURCES})
#       set(HASTIMER "")
#       file(STRINGS ${SOURCE} HASTIMER REGEX "^#include[ ]*[\"]simple_timer.h[\"].*$")
#       message(STATUS " ${SOURCE} HASTIMER --${HASTIMER}--")
#       if ( "" STREQUAL "${HASTIMER}" )
#         message(STATUS " plain compile ${SOURCE} ")
#       else()
#         get_filename_component(base_name ${SOURCE} NAME_WE)
# #        string(REGEX REPLACE "${CMAKE_CURRENT_SOURCE_DIR}/${TESTNAME}/\([^/]*\).\(f90|f08\)"
#           set(SOURCE_F08 "${CMAKE_CURRENT_BINARY_DIR}/${base_name}.f08" ) #SOURCE_F08 ${SOURCE})
#         #list(APPEND F2008FILES ${F2008FILES} ${SOURCE_F08})
#         message(STATUS " Replacing ${SOURCE} with ${SOURCE_F08}")


#          add_custom_command(
#            OUTPUT ${SOURCE_F08}
#            COMMAND ${CMAKE_FCPP_COMPILER} -C -P -I${CMAKE_SOURCE_DIR}/simple_utils -I${LIB} -J${LIB} -E ${SOURCE} -o ${SOURCE_F08}
           
#            MAIN_DEPENDENCY ${SOURCE}
#            WORKING_DIRECTORY ${LIB}
#         #  DEPENDS "${CMAKE_SOURCE_DIR}/simple_utils/simple_timer.h"
#            COMMENT "Preprocessing ${SOURCE} in CPP for variadic macros."
#            VERBATIM
#            )
#         add_custom_target(${base_name} ALL
#           DEPENDS ${SOURCE_F08} #  "${CMAKE_SOURCE_DIR}/simple_utils/simple_timer.h"
#           COMMENT "Preprocessed ${SOURCE_F08}."
#           VERBATIM
#           )
#         list(REMOVE_ITEM TESTSOURCES ${SOURCE})
#         list(APPEND TESTSOURCES ${SOURCE_F08})

#         set_source_files_properties(${SOURCE_F08} PROPERTIES
#           GENERATED TRUE
#           Fortran_FORMAT FREE)

#       endif()
#     endforeach(SOURCE)


    add_executable(${TESTEXE} ${TESTSOURCES})
    target_link_libraries(${TESTEXE} ${EXTRA_LIBS})
    if(USE_OPENMP)
      set_target_properties(${TESTEXE} PROPERTIES
        COMPILE_FLAGS "${OpenMP_Fortran_FLAGS}"
        LINK_FLAGS "${OpenMP_Fortran_FLAGS}"
        LANGUAGE Fortran)
    endif(USE_OPENMP)

    install(TARGETS ${TESTEXE} RUNTIME DESTINATION bin)

endif()
endforeach(TESTNAME)


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
  CONFIGURATIONS Debug|Release
  WORKING_DIRECTORY bin
  COMMAND simple_test_parse nthr=4) 

add_test(NAME test_units
  CONFIGURATIONS Debug|Release
  WORKING_DIRECTORY bin
  COMMAND simple_test_units nthr=4) 

add_test(NAME test_timer
  CONFIGURATIONS Debug|Release
  WORKING_DIRECTORY bin
  COMMAND simple_test_timer ) 
