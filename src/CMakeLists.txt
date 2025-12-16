# src/CMakeLists.txt
#
# Builds the SIMPLE shared library (${SIMPLELIB})

# ------------------------------------------------------------------------------
# Generated sources (simple_args + varlist)
# ------------------------------------------------------------------------------
set(SIMPLE_GEN_DIR "${SIMPLE_BINARY_LIB_DIR}")  # from root

add_custom_command(
    OUTPUT
        "${SIMPLE_GEN_DIR}/simple_args.f90"
        "${SIMPLE_GEN_DIR}/simple_varlist.txt"
    COMMAND SIMPLE_PATH=${CMAKE_BINARY_DIR}
            "${CMAKE_SOURCE_DIR}/scripts/simple_args_generator.pl"
    DEPENDS
        "${CMAKE_SOURCE_DIR}/scripts/simple_args_generator.pl"
        "${CMAKE_CURRENT_SOURCE_DIR}/main/simple_parameters.f90"
    WORKING_DIRECTORY "${CMAKE_CURRENT_SOURCE_DIR}/main"
    COMMENT "Generating simple_args.f90 and simple_varlist.txt"
    VERBATIM
)

add_custom_target(simple_args
    DEPENDS
        "${SIMPLE_GEN_DIR}/simple_args.f90"
        "${SIMPLE_GEN_DIR}/simple_varlist.txt"
)

# Git commit hash insertion (side-effect script)
add_custom_target(simple_git
    COMMAND SIMPLE_PATH=${CMAKE_BINARY_DIR}
            "${CMAKE_SOURCE_DIR}/scripts/insert_git_commit_hash.pl"
    DEPENDS "${CMAKE_SOURCE_DIR}/scripts/insert_git_commit_hash.pl"
    WORKING_DIRECTORY "${CMAKE_CURRENT_SOURCE_DIR}/main"
    COMMENT "Inserting git commit hash"
    VERBATIM
)

# ------------------------------------------------------------------------------
# Sources
# ------------------------------------------------------------------------------
file(GLOB DEFS_src               "defs/*.[Ff][90][08]")
file(GLOB GNUFOR_src             "extlibs/gnufor/*.[Ff][90][08]")
file(GLOB JSON_src               "extlibs/json/*.[Ff][90][08]")
file(GLOB NETLIB_src             "extlibs/netlib/*.[Ff][90][08]")
file(GLOB UNIX_src               "extlibs/unix/src/*.[Ff][90][08]")
file(GLOB XML_src                "extlibs/xml/*/*.[Ff][90][08]")
file(GLOB CPLOT2D_src            "extlibs/CPlot2D/*.cpp" "extlibs/CPlot2D/*.[Ff][90][08]")
file(GLOB_RECURSE UTILS_src      "utils/*.[Ff][90][08]")
file(GLOB FILEIO_src             "fileio/*.[Ff][90][08]")
file(GLOB STARFILE_src           "fileio/starfile/*.cpp")
file(GLOB_RECURSE SIMPLE_src     "main/*.[Ff][90][08]")

# JPG sources including C
file(GLOB JPG_src "extlibs/jpg/*.[Ff][90][08]")
list(APPEND JPG_src "extlibs/jpg/stb_image_write.c")
if(USE_LIBTIFF)
    list(APPEND JPG_src "extlibs/tiff/libtiff_funs.c")
endif()

# FILEIO C sources
list(APPEND FILEIO_src "fileio/simple_posix.c")

# OPERATORS C sources
set(OPERATORS_src "main/interp/simple_kbinterpol_memo.c")

# Generated simple_args.f90
list(APPEND SIMPLE_src "${SIMPLE_GEN_DIR}/simple_args.f90")

set(SIMPLE_ALL_SOURCES
    ${DEFS_src}
    ${GNUFOR_src}
    ${JSON_src}
    ${UTILS_src}
    ${FILEIO_src}
    ${STARFILE_src}
    ${CPLOT2D_src}
    ${NETLIB_src}
    ${XML_src}
    ${UNIX_src}
    ${JPG_src}
    ${OPERATORS_src}
    ${SIMPLE_src}
)

# ------------------------------------------------------------------------------
# Library target
# ------------------------------------------------------------------------------

add_library(${SIMPLELIB} STATIC ${SIMPLE_ALL_SOURCES})

add_dependencies(${SIMPLELIB} simple_args simple_git)

set_target_properties(${SIMPLELIB} PROPERTIES
    VERSION                 ${PROJECT_VERSION}
    SOVERSION               ${PROJECT_VERSION_MAJOR}
    Fortran_MODULE_DIRECTORY "${CMAKE_Fortran_MODULE_DIRECTORY}"
    POSITION_INDEPENDENT_CODE ON
)

target_link_libraries(${SIMPLELIB}
    PUBLIC
        ${SIMPLE_LIBRARIES}
)

target_include_directories(${SIMPLELIB}
    PUBLIC
        $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/inc>
        $<BUILD_INTERFACE:${CMAKE_Fortran_MODULE_DIRECTORY}>
        $<INSTALL_INTERFACE:include>
        $<BUILD_INTERFACE:${SIMPLE_GEN_DIR}>
)

# ------------------------------------------------------------------------------
# Installation
# ------------------------------------------------------------------------------
install(TARGETS ${SIMPLELIB}
    LIBRARY DESTINATION ${CMAKE_INSTALL_LIBDIR}
    ARCHIVE DESTINATION ${CMAKE_INSTALL_LIBDIR}
    RUNTIME DESTINATION ${CMAKE_INSTALL_BINDIR}
)
# Install Fortran modules
install(DIRECTORY "${CMAKE_Fortran_MODULE_DIRECTORY}/"
    DESTINATION "${CMAKE_INSTALL_INCLUDEDIR}/simple"
    FILES_MATCHING PATTERN "*.mod"
)
# Install varlist for command line argument error checking
install(FILES "${SIMPLE_GEN_DIR}/simple_varlist.txt"
    DESTINATION "${CMAKE_INSTALL_LIBDIR}/simple"
)
