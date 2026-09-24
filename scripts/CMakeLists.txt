
find_package(Perl REQUIRED)
find_package(Python3 REQUIRED COMPONENTS Interpreter)

# Keep the two generated repository navigation records synchronized with the
# sources. CONFIGURE_DEPENDS makes a newly added source or documentation file
# reconfigure the glob; the custom commands otherwise run only when an input is
# newer than their output, so an incremental build does not repeatedly rescan
# the tree.
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
set(SIMPLE_CODE_MAP "${CMAKE_SOURCE_DIR}/doc/code_overview/code_base_map.md")
add_custom_command(
  OUTPUT "${SIMPLE_CODE_MAP}"
  COMMAND "${PERL_EXECUTABLE}" "${CMAKE_SOURCE_DIR}/scripts/generate_codeoverview.pl"
    --root "${CMAKE_SOURCE_DIR}"
    --out "${SIMPLE_CODE_MAP}"
  DEPENDS
    "${CMAKE_SOURCE_DIR}/scripts/generate_codeoverview.pl"
    ${SIMPLE_CODE_MAP_INPUTS}
  COMMENT "Regenerating the SIMPLE codebase map"
  VERBATIM
)
add_custom_target(generate_codebase_map ALL DEPENDS "${SIMPLE_CODE_MAP}")

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
set(SIMPLE_TEST_INVENTORY "${CMAKE_SOURCE_DIR}/doc/refactoring_notes/test_inventory.md")
list(REMOVE_ITEM SIMPLE_TEST_INVENTORY_INPUTS "${SIMPLE_TEST_INVENTORY}")
add_custom_command(
  OUTPUT "${SIMPLE_TEST_INVENTORY}"
  COMMAND "${Python3_EXECUTABLE}" "${CMAKE_SOURCE_DIR}/scripts/test_review_dossier.py"
    --inventory "${SIMPLE_TEST_INVENTORY}"
    --inventory-only
  DEPENDS
    "${CMAKE_SOURCE_DIR}/scripts/test_review_dossier.py"
    "${CMAKE_SOURCE_DIR}/production/CMakeLists.txt"
    ${SIMPLE_TEST_INVENTORY_INPUTS}
  COMMENT "Regenerating the SIMPLE test inventory"
  VERBATIM
)
add_custom_target(generate_test_inventory ALL DEPENDS "${SIMPLE_TEST_INVENTORY}")

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
