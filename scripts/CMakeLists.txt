
find_package(Perl REQUIRED)
find_package(Python3 REQUIRED COMPONENTS Interpreter)

# Two generated records, rebuilt whenever their inputs change and not committed (.gitignore):
# the code map (scripts/generate_codeoverview.pl) and the test inventory
# (scripts/test_review_dossier.py, from the sources and the hand-written, committed review record
# doc/code_overview/test_review_record.md).
#  - Each is written in the build tree and copied into doc/ only when its content differs, so an
#    unchanged record keeps its timestamp and nothing downstream reruns; a stamp in the build tree
#    marks the run.
#  - Each is also made at configure time when it is missing (a fresh clone): the inventory's inputs
#    are a CONFIGURE_DEPENDS glob over doc/, and a file that first appears there during the build
#    makes the next make (the install after the fast gate) report a GLOB mismatch, reconfigure and
#    rescan the whole tree.
#  - The inventory reads the documentation, the code map among it, so it waits for the map.
set(SIMPLE_CODE_MAP       "${CMAKE_SOURCE_DIR}/doc/code_overview/code_base_map.md")
set(SIMPLE_TEST_INVENTORY "${CMAKE_SOURCE_DIR}/doc/code_overview/test_inventory.md")
set(SIMPLE_TEST_RECORD    "${CMAKE_SOURCE_DIR}/doc/code_overview/test_review_record.md")
set(SIMPLE_CODE_MAP_CMD "${PERL_EXECUTABLE}" "${CMAKE_SOURCE_DIR}/scripts/generate_codeoverview.pl"
  --quiet --root "${CMAKE_SOURCE_DIR}" --out)
set(SIMPLE_TEST_INVENTORY_CMD "${Python3_EXECUTABLE}" "${CMAKE_SOURCE_DIR}/scripts/test_review_dossier.py"
  --inventory-only --quiet --record "${SIMPLE_TEST_RECORD}" --inventory)
if(NOT EXISTS "${SIMPLE_CODE_MAP}")
  execute_process(COMMAND ${SIMPLE_CODE_MAP_CMD} "${SIMPLE_CODE_MAP}" OUTPUT_QUIET)
endif()
if(NOT EXISTS "${SIMPLE_TEST_INVENTORY}")
  execute_process(COMMAND ${SIMPLE_TEST_INVENTORY_CMD} "${SIMPLE_TEST_INVENTORY}" OUTPUT_QUIET)
endif()

file(GLOB_RECURSE SIMPLE_CODE_MAP_INPUTS CONFIGURE_DEPENDS
  LIST_DIRECTORIES false
  "${CMAKE_SOURCE_DIR}/src/*.f90"
  "${CMAKE_SOURCE_DIR}/src/*.F90"
  "${CMAKE_SOURCE_DIR}/src/*.inf"
  "${CMAKE_SOURCE_DIR}/production/*.f90"
  "${CMAKE_SOURCE_DIR}/production/*.F90"
  "${CMAKE_SOURCE_DIR}/production/*.inf"
  "${CMAKE_SOURCE_DIR}/scripts/*.f90"
  "${CMAKE_SOURCE_DIR}/scripts/*.F90"
  "${CMAKE_SOURCE_DIR}/scripts/*.inf"
)
set(SIMPLE_CODE_MAP_BUILD "${CMAKE_CURRENT_BINARY_DIR}/code_base_map.md")
set(SIMPLE_CODE_MAP_STAMP "${CMAKE_CURRENT_BINARY_DIR}/code_base_map.stamp")
add_custom_command(
  OUTPUT "${SIMPLE_CODE_MAP_STAMP}"
  COMMAND ${SIMPLE_CODE_MAP_CMD} "${SIMPLE_CODE_MAP_BUILD}"
  COMMAND "${CMAKE_COMMAND}" -E copy_if_different "${SIMPLE_CODE_MAP_BUILD}" "${SIMPLE_CODE_MAP}"
  COMMAND "${CMAKE_COMMAND}" -E touch "${SIMPLE_CODE_MAP_STAMP}"
  DEPENDS
    "${CMAKE_SOURCE_DIR}/scripts/generate_codeoverview.pl"
    ${SIMPLE_CODE_MAP_INPUTS}
  COMMENT "Regenerating the SIMPLE codebase map"
  VERBATIM
)
add_custom_target(generate_codebase_map ALL DEPENDS "${SIMPLE_CODE_MAP_STAMP}")

file(GLOB_RECURSE SIMPLE_TEST_INVENTORY_INPUTS CONFIGURE_DEPENDS
  LIST_DIRECTORIES false
  "${CMAKE_SOURCE_DIR}/src/*.f90"
  "${CMAKE_SOURCE_DIR}/production/*.f90"
  "${CMAKE_SOURCE_DIR}/scripts/*.py"
  "${CMAKE_SOURCE_DIR}/scripts/*.pl"
  "${CMAKE_SOURCE_DIR}/scripts/*.sh"
  "${CMAKE_SOURCE_DIR}/doc/*.md"
  "${CMAKE_SOURCE_DIR}/doc/*.txt"
  "${CMAKE_SOURCE_DIR}/.github/*.yml"
  "${CMAKE_SOURCE_DIR}/.github/*.yaml"
  "${CMAKE_SOURCE_DIR}/.github/*.md"
  "${CMAKE_SOURCE_DIR}/nice/*.md"
  "${CMAKE_SOURCE_DIR}/nice/*.txt"
  "${CMAKE_SOURCE_DIR}/nice/*.yml"
  "${CMAKE_SOURCE_DIR}/nice/*.yaml"
)
list(REMOVE_ITEM SIMPLE_TEST_INVENTORY_INPUTS "${SIMPLE_TEST_INVENTORY}")
set(SIMPLE_TEST_INVENTORY_BUILD "${CMAKE_CURRENT_BINARY_DIR}/test_inventory.md")
set(SIMPLE_TEST_INVENTORY_STAMP "${CMAKE_CURRENT_BINARY_DIR}/test_inventory.stamp")
add_custom_command(
  OUTPUT "${SIMPLE_TEST_INVENTORY_STAMP}"
  COMMAND ${SIMPLE_TEST_INVENTORY_CMD} "${SIMPLE_TEST_INVENTORY_BUILD}"
  COMMAND "${CMAKE_COMMAND}" -E copy_if_different "${SIMPLE_TEST_INVENTORY_BUILD}" "${SIMPLE_TEST_INVENTORY}"
  COMMAND "${CMAKE_COMMAND}" -E touch "${SIMPLE_TEST_INVENTORY_STAMP}"
  DEPENDS
    "${CMAKE_SOURCE_DIR}/scripts/test_review_dossier.py"
    "${CMAKE_SOURCE_DIR}/production/CMakeLists.txt"
    "${SIMPLE_TEST_RECORD}"
    ${SIMPLE_TEST_INVENTORY_INPUTS}
  COMMENT "Regenerating the SIMPLE test inventory"
  VERBATIM
)
add_custom_target(generate_test_inventory ALL DEPENDS "${SIMPLE_TEST_INVENTORY_STAMP}")
add_dependencies(generate_test_inventory generate_codebase_map)

install(DIRECTORY  ${CMAKE_CURRENT_SOURCE_DIR}
  DESTINATION  ${CMAKE_INSTALL_PREFIX}
  CONFIGURATIONS DEBUG|RELEASE
  PATTERN "backups/*" EXCLUDE
  PATTERN "*.txt"     EXCLUDE
  PATTERN "*.env"     EXCLUDE
  PATTERN "*.pl"
    PERMISSIONS OWNER_EXECUTE OWNER_WRITE OWNER_READ GROUP_EXECUTE GROUP_READ WORLD_EXECUTE WORLD_READ
  PATTERN "*.py"
    PERMISSIONS OWNER_EXECUTE OWNER_WRITE OWNER_READ GROUP_EXECUTE GROUP_READ WORLD_EXECUTE WORLD_READ
  )
