# Build the SIMPLE library

# Generate simple_args.f90 from Perl script
add_custom_target(simple_args
    COMMAND SIMPLE_PATH=${CMAKE_BINARY_DIR} ${CMAKE_SOURCE_DIR}/scripts/simple_args_generator.pl
    DEPENDS "${CMAKE_SOURCE_DIR}/scripts/simple_args_generator.pl"
            "${CMAKE_CURRENT_SOURCE_DIR}/main/simple_parameters.f90"
    BYPRODUCTS "${CMAKE_BINARY_DIR}/lib/simple/simple_args.f90"
               "${CMAKE_BINARY_DIR}/lib/simple/simple_varlist.txt"
    WORKING_DIRECTORY "${CMAKE_CURRENT_SOURCE_DIR}/main"
    COMMENT "Generating simple_args.f90"
    USES_TERMINAL
)

# Insert git commit hash
add_custom_target(simple_git
    COMMAND SIMPLE_PATH=${CMAKE_BINARY_DIR} ${CMAKE_SOURCE_DIR}/scripts/insert_git_commit_hash.pl
    DEPENDS "${CMAKE_SOURCE_DIR}/scripts/insert_git_commit_hash.pl"
    WORKING_DIRECTORY "${CMAKE_CURRENT_SOURCE_DIR}/main"
    COMMENT "Inserting git commit hash"
    USES_TERMINAL
)

# Collect source files
file(GLOB DEFS_src       "defs/*.[Ff][90][08]")
file(GLOB GNUFOR_src     "extlibs/gnufor/*.[Ff][90][08]")
file(GLOB JSON_src       "extlibs/json/*.[Ff][90][08]")
file(GLOB NETLIB_src     "extlibs/netlib/*.[Ff][90][08]")
file(GLOB UNIX_src       "extlibs/unix/src/*.[Ff][90][08]")
file(GLOB XML_src        "extlibs/xml/*/*.[Ff][90][08]")
file(GLOB CPLOT2D_src    "extlibs/CPlot2D/*.cpp" "extlibs/CPlot2D/*.[Ff][90][08]")
file(GLOB_RECURSE UTILS_src "utils/*.[Ff][90][08]")
file(GLOB FILEIO_src     "fileio/*.[Ff][90][08]")
file(GLOB STARFILE_src   "fileio/starfile/*.cpp")
file(GLOB_RECURSE SIMPLE_src "main/*.[Ff][90][08]")

# JPG sources including C files
file(GLOB JPG_src "extlibs/jpg/*.[Ff][90][08]")
list(APPEND JPG_src "extlibs/jpg/stb_image_write.c")
if(USE_LIBTIFF)
    list(APPEND JPG_src "extlibs/tiff/libtiff_funs.c")
endif()

# FILEIO C sources
list(APPEND FILEIO_src "fileio/simple_posix.c")

# OPERATORS C sources
set(OPERATORS_src "main/interp/simple_kbinterpol_memo.c")

# Add generated simple_args
list(APPEND SIMPLE_src "${CMAKE_BINARY_DIR}/lib/simple/simple_args.f90")

# Combine all sources
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

# Create shared library
add_library(${SIMPLELIB} SHARED ${SIMPLE_ALL_SOURCES})

# Add dependencies on generated files
add_dependencies(${SIMPLELIB} simple_args simple_git)

# Set library properties
set_target_properties(${SIMPLELIB} PROPERTIES
    VERSION ${PROJECT_VERSION}
    SOVERSION ${PROJECT_VERSION_MAJOR}
    PUBLIC_HEADER ""
    Fortran_MODULE_DIRECTORY "${CMAKE_Fortran_MODULE_DIRECTORY}"
)

# Link libraries
target_link_libraries(${SIMPLELIB} PUBLIC ${SIMPLE_LIBRARIES})

# Include directories
target_include_directories(${SIMPLELIB} PUBLIC
    $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/inc>
    $<BUILD_INTERFACE:${CMAKE_Fortran_MODULE_DIRECTORY}>
    $<INSTALL_INTERFACE:include>
)

# Installation
install(TARGETS ${SIMPLELIB}
    LIBRARY DESTINATION ${CMAKE_INSTALL_LIBDIR}
    ARCHIVE DESTINATION ${CMAKE_INSTALL_LIBDIR}
    RUNTIME DESTINATION ${CMAKE_INSTALL_BINDIR}
)

# Install Fortran modules
install(DIRECTORY ${CMAKE_Fortran_MODULE_DIRECTORY}/
    DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}/simple
    FILES_MATCHING PATTERN "*.mod"
)

# Install varlist for testing
install(FILES "${CMAKE_BINARY_DIR}/lib/simple/simple_varlist.txt"
    DESTINATION ${CMAKE_INSTALL_LIBDIR}/simple
)