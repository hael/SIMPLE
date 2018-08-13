# CMake Configuration and build added by Izaak Beekman -- May, 2014

# Copy right (c) 2014, Izaak Beekman
# All rights reserved.

# This file is contributed to the JSON-Fortran project, and
# is licensed under the terms of JSON-Fortran license. The JSON-Fortran
# license is located in the LICENSE file which must be distributed with
# this software. The contributing author, Izaak Beekman, retains all
# rights permitted by the terms of the JSON-Fortran license.


 file ( GLOB JF_TEST_SRCS "${CMAKE_SOURCE_DIR}/production/tests/json_fortran/jf_test_*.F90" )
 file ( COPY "${CMAKE_SOURCE_DIR}/production/tests/json-fortran/files"
    DESTINATION "${CMAKE_BINARY_DIR}/Testing/" )
  set ( DATA_DIR "${CMAKE_BINARY_DIR}/Testing/files" )
  set_directory_properties ( PROPERTIES ADDITIONAL_MAKE_CLEAN_FILES
    "${DATA_DIR}/test2.json;${DATA_DIR}/test4.json;" )
  find_program ( JSONLINT jsonlint )
   find_program ( DIFF     diff )
    # Validate input
  if ( JSONLINT )
    file ( GLOB JSON_INPUTS "${DATA_DIR}/inputs/*.json" )
    file ( GLOB INVALID_JSON "${DATA_DIR}/inputs/*invalid*.json" "${DATA_DIR}/inputs/comments.json")

    list ( REMOVE_ITEM JSON_INPUTS ${INVALID_JSON} )
    list ( REMOVE_ITEM JSON_INPUTS "${DATA_DIR}/inputs/big.json" ) # This takes too long and is valid
                                                                   # JSON from a trusted source

    foreach ( VALID_JSON ${JSON_INPUTS} )
      get_filename_component ( TESTNAME "${VALID_JSON}" NAME )
      add_test ( NAME validate-${TESTNAME}
	WORKING_DIRECTORY "${DATA_DIR}/inputs"
	COMMAND ${JSONLINT} "--allow=nonescape-characters" "${VALID_JSON}" )
    endforeach ()

    foreach ( INVALID ${INVALID_JSON} )
      get_filename_component ( TESTNAME "${INVALID}" NAME )
      add_test ( NAME validate-${TESTNAME}
	WORKING_DIRECTORY "${DATA_DIR}/inputs"
	COMMAND ${JSONLINT} "${INVALID}" )
      set_property ( TEST validate-${TESTNAME}
	PROPERTY
	WILL_FAIL TRUE)
    endforeach ()
  endif ()
 set ( UNIT_TESTS '' )
  foreach ( UNIT_TEST ${JF_TEST_SRCS} )
    get_filename_component ( TEST ${UNIT_TEST} NAME_WE )
    if(MSVC_IDE)
    link_directories(${CMAKE_BINARY_DIR}/lib)
    endif()
    add_executable ( ${TEST} EXCLUDE_FROM_ALL ${UNIT_TEST} )
    target_link_libraries ( ${TEST} ${SIMPLELIB} )
    add_dependencies ( check_json ${TEST} )
    set_target_properties ( ${TEST}
      PROPERTIES
      RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/bin )
    add_test( NAME ${TEST}
      WORKING_DIRECTORY ${CMAKE_BINARY_DIR}/bin
      COMMAND ./${TEST})
    list ( APPEND UNIT_TESTS ${TEST} )
    if ( JSONLINT )
      set_property ( TEST ${TEST}
	APPEND
	PROPERTY DEPENDS validate-input1 validate-input2 )
    endif()
  endforeach ( UNIT_TEST )

  set_property ( TEST jf_test_03
    APPEND
    PROPERTY DEPENDS jf_test_02 )

  # Validate output
  if ( JSONLINT )
    file ( GLOB JSON_FILES "${DATA_DIR}/*.json" )
    foreach ( JSON_FILE ${JSON_FILES} )
      get_filename_component ( TESTNAME ${JSON_FILE} NAME )
      add_test ( NAME validate-output-${TESTNAME}
	WORKING_DIRECTORY "${DATA_DIR}"
	COMMAND ${JSONLINT} "--allow=nonescape-characters" ${TESTNAME} )
      set_property ( TEST validate-output-${TESTNAME}
	APPEND
	PROPERTY
	DEPENDS ${UNIT_TESTS}
	REQUIRED_FILES ${JSON_FILES} )
    endforeach ( JSON_FILE )
  endif ()

  # Check output for differences
  if ( DIFF )
    file ( GLOB JSON_FILES "${DATA_DIR}/*.json" )
    foreach ( JSON_FILE ${JSON_FILES} )
      get_filename_component ( JSON_STEM ${JSON_FILE} NAME_WE )
      add_test ( NAME regression-${JSON_STEM}.json
	WORKING_DIRECTORY "${DATA_DIR}"
	COMMAND ${DIFF} -q ${JSON_STEM}.json expected-outputs/${JSON_STEM}.json )
      set_property ( TEST regression-${JSON_STEM}.json
	APPEND
	PROPERTY
	DEPENDS ${UNIT_TESTS}
	REQUIRED_FILES ${JSON_FILES} )
    endforeach ( JSON_FILE )
  else ()
    message ( WARNING
      "For full test coverage diff, or a similar tool must be present on your system" )
  endif ()
