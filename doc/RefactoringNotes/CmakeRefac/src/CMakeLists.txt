# src/CMakeLists.txt

# Canonical target name (do NOT encode version)
set(SIMPLELIB simple)

# ----------------------------------------------------------------------
# Generated sources
# ----------------------------------------------------------------------
set(GENERATED_DIR "${CMAKE_CURRENT_BINARY_DIR}/generated")
file(MAKE_DIRECTORY "${GENERATED_DIR}")

set(SIMPLE_ARGS_F90 "${GENERATED_DIR}/simple_args.f90")
set(SIMPLE_VARLIST  "${GENERATED_DIR}/simple_varlist.txt")

add_custom_command(
  OUTPUT  ${SIMPLE_ARGS_F90} ${SIMPLE_VARLIST}
  COMMAND ${CMAKE_COMMAND} -E env SIMPLE_PATH=${CMAKE_BINARY_DIR}
          ${CMAKE_SOURCE_DIR}/scripts/simple_args_generator.pl
  DEPENDS
          ${CMAKE_SOURCE_DIR}/scripts/simple_args_generator.pl
          ${CMAKE_CURRENT_SOURCE_DIR}/main/simple_parameters.f90
  WORKING_DIRECTORY "${CMAKE_CURRENT_SOURCE_DIR}/main"
  COMMENT "Generating simple_args.f90"
)

set(GIT_HASH_F90 "${GENERATED_DIR}/git_commit.f90")

add_custom_command(
  OUTPUT ${GIT_HASH_F90}
  COMMAND ${CMAKE_COMMAND} -E env SIMPLE_PATH=${CMAKE_BINARY_DIR}
          ${CMAKE_SOURCE_DIR}/scripts/insert_git_commit_hash.pl
  DEPENDS ${CMAKE_SOURCE_DIR}/scripts/insert_git_commit_hash.pl
  WORKING_DIRECTORY "${CMAKE_CURRENT_SOURCE_DIR}/main"
  COMMENT "Embedding git commit hash"
)

# ----------------------------------------------------------------------
# Library
# ----------------------------------------------------------------------
add_library(${SIMPLELIB})

# --- Sources (use CONFIGURE_DEPENDS as interim step)
file(GLOB CONFIGURE_DEPENDS
  DEFS_src        "${CMAKE_CURRENT_SOURCE_DIR}/defs/*.[Ff]90"
  GNUFOR_src      "${CMAKE_CURRENT_SOURCE_DIR}/extlibs/gnufor/*.[Ff]90"
  JSON_src        "${CMAKE_CURRENT_SOURCE_DIR}/extlibs/json/*.[Ff]90"
  NETLIB_src      "${CMAKE_CURRENT_SOURCE_DIR}/extlibs/netlib/*.[Ff]90"
  UNIX_src        "${CMAKE_CURRENT_SOURCE_DIR}/extlibs/unix/src/*.[Ff]90"
  XML_src         "${CMAKE_CURRENT_SOURCE_DIR}/extlibs/xml/*/*.[Ff]90"
  UTILS_src       "${CMAKE_CURRENT_SOURCE_DIR}/utils/*.[Ff]90"
  FILEIO_src      "${CMAKE_CURRENT_SOURCE_DIR}/fileio/*.[Ff]90"
  SIMPLE_src      "${CMAKE_CURRENT_SOURCE_DIR}/main/*.[Ff]90"
)

set(JPG_src
  extlibs/jpg/stb_image_write.c
)

if(USE_LIBTIFF)
  list(APPEND JPG_src extlibs/tiff/libtiff_funs.c)
endif()

target_sources(${SIMPLELIB}
  PRIVATE
    ${DEFS_src}
    ${GNUFOR_src}
    ${JSON_src}
    ${NETLIB_src}
    ${UNIX_src}
    ${XML_src}
    ${UTILS_src}
    ${FILEIO_src}
    ${SIMPLE_src}
    ${JPG_src}
    main/interp/simple_kbinterpol_memo.c
    ${SIMPLE_ARGS_F90}
    ${GIT_HASH_F90}
)

# ----------------------------------------------------------------------
# Fortran modules (per-target!)
# ----------------------------------------------------------------------
set_target_properties(${SIMPLELIB} PROPERTIES
  VERSION   ${PROJECT_VERSION}
  SOVERSION ${PROJECT_VERSION_MAJOR}
  Fortran_MODULE_DIRECTORY
    ${PROJECT_BINARY_DIR}/modules/${SIMPLELIB}
)

target_include_directories(${SIMPLELIB}
  PUBLIC
    $<BUILD_INTERFACE:${PROJECT_BINARY_DIR}/modules/${SIMPLELIB}>
    $<INSTALL_INTERFACE:${CMAKE_INSTALL_INCLUDEDIR}/simple>
)

# ----------------------------------------------------------------------
# Dependencies
# ----------------------------------------------------------------------
target_link_libraries(${SIMPLELIB}
  PUBLIC
    ${SIMPLE_LIBRARIES}   # keep until Dependencies.cmake is modernized
)

# ----------------------------------------------------------------------
# Installation
# ----------------------------------------------------------------------
install(TARGETS ${SIMPLELIB}
  EXPORT SIMPLETargets
  LIBRARY DESTINATION ${CMAKE_INSTALL_LIBDIR}
  ARCHIVE DESTINATION ${CMAKE_INSTALL_LIBDIR}
)

install(DIRECTORY
  ${PROJECT_BINARY_DIR}/modules/${SIMPLELIB}/
  DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}/simple
  FILES_MATCHING PATTERN "*.mod"
)

install(FILES
  ${SIMPLE_VARLIST}
  DESTINATION ${CMAKE_INSTALL_LIBDIR}/simple
)
