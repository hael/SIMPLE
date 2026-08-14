# production/CMakeLists.txt
#
# Build production executables and (optionally) tests.

# ------------------------------------------------------------------------------
# Main executables
# ------------------------------------------------------------------------------
set(SIMPLE_EXE          simple_exec)
set(SIMPLE_TEST_EXE     simple_test_exec)
set(SINGLE_EXE          single_exec)
set(STREAM_EXE          simple_stream)
set(SIMPLE_PRIVATE_EXE  simple_private_exec)
set(SIMPLE_WORKER_EXE   simple_persistent_worker)

add_executable(${SIMPLE_EXE}
    "${CMAKE_CURRENT_SOURCE_DIR}/simple_exec.f90")
add_executable(${SIMPLE_TEST_EXE}
    "${CMAKE_CURRENT_SOURCE_DIR}/simple_test_exec.f90")
add_executable(${SINGLE_EXE}
    "${CMAKE_CURRENT_SOURCE_DIR}/single_exec.f90")
add_executable(${STREAM_EXE}
    "${CMAKE_CURRENT_SOURCE_DIR}/simple_stream.f90")
add_executable(${SIMPLE_PRIVATE_EXE}
    "${CMAKE_CURRENT_SOURCE_DIR}/simple_private_exec.f90"
    "${CMAKE_CURRENT_SOURCE_DIR}/simple_private_exec_driver.f90")
add_executable(${SIMPLE_WORKER_EXE}
    "${CMAKE_CURRENT_SOURCE_DIR}/simple_persistent_worker.f90")

set(SIMPLE_EXECUTABLES
    ${SIMPLE_EXE}
    ${SIMPLE_TEST_EXE}
    ${SINGLE_EXE}
    ${STREAM_EXE}
    ${SIMPLE_PRIVATE_EXE}
    ${SIMPLE_WORKER_EXE}
)

foreach(exe IN LISTS SIMPLE_EXECUTABLES)
    set_target_properties(${exe} PROPERTIES
        LINKER_LANGUAGE Fortran
    )
    if(USE_COARRAYS AND exe STREQUAL SIMPLE_PRIVATE_EXE)
        set_target_properties(${exe} PROPERTIES
            Fortran_MODULE_DIRECTORY "${SIMPLE_COARRAY_MODULE_DIRECTORY}"
        )
        target_link_libraries(${exe}
            PRIVATE
                ${SIMPLE_COARRAYLIB}
                ${SIMPLE_LIBRARIES}
        )
    else()
        target_link_libraries(${exe}
            PRIVATE
                ${SIMPLELIB}
                ${SIMPLE_LIBRARIES}
        )
    endif()
    if( USE_OPENMP )
      target_link_options(${exe} PUBLIC -fopenmp)
    endif()
endforeach()

install(TARGETS ${SIMPLE_EXECUTABLES}
    RUNTIME DESTINATION ${CMAKE_INSTALL_BINDIR}
)

# ------------------------------------------------------------------------------
# Complete UI JSON validation
# ------------------------------------------------------------------------------
# The UI descriptor is a client contract.  Regenerate and validate the complete
# emitted registry whenever simple_exec is rebuilt, independently of the
# optional CTest suite. Use the serial front end so metadata validation does not
# initialize an optional coarray/MPI runtime.
find_package(Python3 REQUIRED COMPONENTS Interpreter)
set(UI_JSON_VALIDATION_OUTPUT "${CMAKE_CURRENT_BINARY_DIR}/simple_ui.json")
set(UI_JSON_VALIDATION_STAMP "${CMAKE_CURRENT_BINARY_DIR}/simple_ui.json.valid")
add_custom_command(
  OUTPUT "${UI_JSON_VALIDATION_STAMP}"
  COMMAND "${Python3_EXECUTABLE}"
    "${CMAKE_CURRENT_SOURCE_DIR}/validate_complete_ui_json.py"
    --executable "$<TARGET_FILE:${SIMPLE_EXE}>"
    --output "${UI_JSON_VALIDATION_OUTPUT}"
    --stamp "${UI_JSON_VALIDATION_STAMP}"
  DEPENDS
    ${SIMPLE_EXE}
    "${CMAKE_CURRENT_SOURCE_DIR}/validate_complete_ui_json.py"
  BYPRODUCTS "${UI_JSON_VALIDATION_OUTPUT}"
  COMMENT "Generating and validating the complete UI JSON registry"
  VERBATIM
)
add_custom_target(validate_complete_ui_json ALL
  DEPENDS "${UI_JSON_VALIDATION_STAMP}"
)

# The stream interface is a reviewed acquisition-facing screen. Its exact
# snapshot is deliberately opt-in: ordinary builds verify valid JSON without
# blocking routine wording or layout work on a separate approval step.
set(STREAM_UI_CONTRACT_OUTPUT "${CMAKE_CURRENT_BINARY_DIR}/simple_stream_ui_contract.json")
set(STREAM_UI_CONTRACT_STAMP "${CMAKE_CURRENT_BINARY_DIR}/simple_stream_ui_contract.valid")
set(STREAM_UI_CONTRACT_UPDATE_STAMP "${CMAKE_CURRENT_BINARY_DIR}/simple_stream_ui_contract.updated")
set(STREAM_UI_CONTRACT_REFERENCE "${CMAKE_CURRENT_SOURCE_DIR}/stream_ui_contract.json")
add_custom_command(
  OUTPUT "${STREAM_UI_CONTRACT_STAMP}"
  COMMAND "${Python3_EXECUTABLE}"
    "${CMAKE_CURRENT_SOURCE_DIR}/validate_complete_ui_json.py"
    --executable "$<TARGET_FILE:${SIMPLE_PRIVATE_EXE}>"
    --output "${STREAM_UI_CONTRACT_OUTPUT}"
    --stamp "${STREAM_UI_CONTRACT_STAMP}"
    --check-stream-contract
  DEPENDS
    ${SIMPLE_PRIVATE_EXE}
    "${CMAKE_CURRENT_SOURCE_DIR}/validate_complete_ui_json.py"
    "${STREAM_UI_CONTRACT_REFERENCE}"
  BYPRODUCTS "${STREAM_UI_CONTRACT_OUTPUT}"
  COMMENT "Reviewing the stream UI snapshot"
  VERBATIM
)
add_custom_target(validate_stream_ui_contract
  DEPENDS "${STREAM_UI_CONTRACT_STAMP}"
)
add_custom_target(update_stream_ui_contract
  COMMAND "${Python3_EXECUTABLE}"
    "${CMAKE_CURRENT_SOURCE_DIR}/validate_complete_ui_json.py"
    --executable "$<TARGET_FILE:${SIMPLE_PRIVATE_EXE}>"
    --output "${STREAM_UI_CONTRACT_OUTPUT}"
    --stamp "${STREAM_UI_CONTRACT_UPDATE_STAMP}"
    --update-stream-contract
  DEPENDS
    ${SIMPLE_PRIVATE_EXE}
    "${CMAKE_CURRENT_SOURCE_DIR}/validate_complete_ui_json.py"
  COMMENT "Updating the reviewed stream UI JSON reference"
  VERBATIM
)

# ------------------------------------------------------------------------------
# UI display-default audit
# ------------------------------------------------------------------------------
# This is intentionally opt-in: it is review evidence, not a runtime artifact,
# and it reports ambiguities for human resolution without failing normal builds.
set(UI_DEFAULT_AUDIT_JSON "${CMAKE_CURRENT_BINARY_DIR}/simple_ui_default_audit.json")
set(UI_DEFAULT_AUDIT_MARKDOWN "${CMAKE_CURRENT_BINARY_DIR}/simple_ui_default_audit.md")
add_custom_command(
  OUTPUT "${UI_DEFAULT_AUDIT_JSON}" "${UI_DEFAULT_AUDIT_MARKDOWN}"
  COMMAND "${Python3_EXECUTABLE}"
    "${CMAKE_SOURCE_DIR}/scripts/default_audit.py"
    --source-root "${CMAKE_SOURCE_DIR}"
    --ui-executable "$<TARGET_FILE:${SIMPLE_PRIVATE_EXE}>"
    --output-json "${UI_DEFAULT_AUDIT_JSON}"
    --output-markdown "${UI_DEFAULT_AUDIT_MARKDOWN}"
  DEPENDS
    ${SIMPLE_PRIVATE_EXE}
    "${CMAKE_SOURCE_DIR}/scripts/default_audit.py"
  COMMENT "Generating UI display-default audit report"
  VERBATIM
)
add_custom_target(audit_ui_defaults
  DEPENDS "${UI_DEFAULT_AUDIT_JSON}" "${UI_DEFAULT_AUDIT_MARKDOWN}"
)

# ------------------------------------------------------------------------------
# Tests
# ------------------------------------------------------------------------------
if(BUILD_TESTS)
  enable_testing()
  # Find all test "main" programs
  file(GLOB TEST_MAINS
    RELATIVE "${CMAKE_CURRENT_SOURCE_DIR}"
    "${CMAKE_CURRENT_SOURCE_DIR}/tests/simple_test_*.[fF]90"
  )
  # Optional filtering (adjust patterns to your actual filenames)
  if(NOT USE_MPI)
    list(FILTER TEST_MAINS EXCLUDE REGEX "mpi")
  endif()
  if(NOT USE_COARRAYS)
    list(FILTER TEST_MAINS EXCLUDE REGEX "coarray|coarrays")
  endif()
  # Tests whose fixtures require an external reference volume (vol1=)
  set(VOL1_FIXTURE_TESTS
    continuous_inplane_rotation2D_stage1_validation
    continuous_inplane_cc_grad
    continuous_inplane_refine3D
  )
  foreach(main_src IN LISTS TEST_MAINS)
    get_filename_component(main_we "${main_src}" NAME_WE)         # simple_test_ft_expanded
    string(REGEX REPLACE "^simple_test_" "" test_id "${main_we}") # ft_expanded
    set(test_name "${test_id}")              # ctest name
    set(test_exe  "simple_test_${test_id}")  # executable name (choose what you like)
    # Supporting sources for this test (e.g., simple_ft_expanded_tester.f90)
    file(GLOB test_deps
      "${CMAKE_CURRENT_SOURCE_DIR}/tests/simple_${test_id}_*.[fF]90"
    )
    # Build executable from main + deps
    set(test_sources
      "${CMAKE_CURRENT_SOURCE_DIR}/${main_src}"
      ${test_deps}
    )
    list(REMOVE_DUPLICATES test_sources)
    add_executable(${test_exe} ${test_sources})
    set_target_properties(${test_exe} PROPERTIES
      LINKER_LANGUAGE Fortran
    )
    if(USE_COARRAYS AND test_id STREQUAL "coarrays")
      set_target_properties(${test_exe} PROPERTIES
        Fortran_MODULE_DIRECTORY "${SIMPLE_COARRAY_MODULE_DIRECTORY}"
      )
      target_link_libraries(${test_exe}
        PRIVATE
          ${SIMPLE_COARRAYLIB}
          ${SIMPLE_LIBRARIES}
      )
    else()
      target_link_libraries(${test_exe}
        PRIVATE
          ${SIMPLELIB}
          ${SIMPLE_LIBRARIES}
      )
    endif()
    if( USE_OPENMP )
      target_link_options(${test_exe} PUBLIC -fopenmp)
    endif()
    install(TARGETS ${test_exe}
      RUNTIME DESTINATION ${CMAKE_INSTALL_BINDIR}
    )
    # Optional per-test data: tests/<prefix>*.txt
    file(GLOB test_data
      "${CMAKE_CURRENT_SOURCE_DIR}/tests/${main_we}*.txt"
      "${CMAKE_CURRENT_SOURCE_DIR}/tests/simple_${test_id}_*.txt"
    )
    if(test_data)
      install(FILES ${test_data}
        DESTINATION "${CMAKE_INSTALL_PREFIX}/tests/${test_id}"
      )
    endif()
    if(USE_COARRAYS AND test_id STREQUAL "coarrays")
      add_test(
        NAME ${test_name}
        COMMAND "${Coarray_CAFRUN_EXECUTABLE}" -np 2 $<TARGET_FILE:${test_exe}>
      )
    elseif(test_id IN_LIST VOL1_FIXTURE_TESTS)
      # These fixtures error stop without a reference volume; register them
      # only when one is supplied via -DSIMPLE_TEST_VOL1=<path.mrc>
      if(SIMPLE_TEST_VOL1)
        add_test(NAME ${test_name} COMMAND ${test_exe} vol1=${SIMPLE_TEST_VOL1})
      endif()
    else()
      add_test(NAME ${test_name} COMMAND ${test_exe})
    endif()
  endforeach()

endif()
