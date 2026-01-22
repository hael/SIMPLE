# production/CMakeLists.txt
#
# Build production executables and (optionally) tests.

# ------------------------------------------------------------------------------
# Main executables
# ------------------------------------------------------------------------------
set(SIMPLE_EXE          simple_exec)
set(SINGLE_EXE          single_exec)
set(STREAM_EXE          simple_stream)
set(SIMPLE_PRIVATE_EXE  simple_private_exec)

add_executable(${SIMPLE_EXE}
    "${CMAKE_CURRENT_SOURCE_DIR}/simple_exec.f90")
add_executable(${SINGLE_EXE}
    "${CMAKE_CURRENT_SOURCE_DIR}/single_exec.f90")
add_executable(${STREAM_EXE}
    "${CMAKE_CURRENT_SOURCE_DIR}/simple_stream.f90")
add_executable(${SIMPLE_PRIVATE_EXE}
    "${CMAKE_CURRENT_SOURCE_DIR}/simple_private_exec.f90")

set(SIMPLE_EXECUTABLES
    ${SIMPLE_EXE}
    ${SINGLE_EXE}
    ${STREAM_EXE}
    ${SIMPLE_PRIVATE_EXE}
)

foreach(exe IN LISTS SIMPLE_EXECUTABLES)
    target_link_libraries(${exe}
        PRIVATE
            ${SIMPLELIB}
            ${SIMPLE_LIBRARIES}
    )
endforeach()

install(TARGETS ${SIMPLE_EXECUTABLES}
    RUNTIME DESTINATION ${CMAKE_INSTALL_BINDIR}
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
    target_link_libraries(${test_exe}
      PRIVATE
        ${SIMPLELIB}
        ${SIMPLE_LIBRARIES}
    )
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
    add_test(NAME ${test_name} COMMAND ${test_exe})
  endforeach()
endif()
