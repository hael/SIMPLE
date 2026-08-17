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
        "${CMAKE_CURRENT_SOURCE_DIR}/main/params/simple_parameters.f90"
    WORKING_DIRECTORY "${CMAKE_CURRENT_SOURCE_DIR}/main/params"
    COMMENT "Generating simple_args.f90 and simple_varlist.txt"
    VERBATIM
)

add_custom_target(simple_args
    DEPENDS
        "${SIMPLE_GEN_DIR}/simple_args.f90"
        "${SIMPLE_GEN_DIR}/simple_varlist.txt"
)

# Generated UI display-default lookup.  The generator reads the parameter and
# commander sources but emits a Fortran module consumed only by the UI layer;
# it never changes command execution.
find_package(Python3 REQUIRED COMPONENTS Interpreter)
set(SIMPLE_UI_DEFAULTS_FILE "${SIMPLE_GEN_DIR}/simple_ui_default_values.f90")
file(GLOB_RECURSE SIMPLE_UI_DEFAULTS_SOURCES CONFIGURE_DEPENDS
    "${CMAKE_SOURCE_DIR}/production/simple_exec.f90"
    "${CMAKE_SOURCE_DIR}/production/single_exec.f90"
    "${CMAKE_SOURCE_DIR}/production/simple_stream.f90"
    "${CMAKE_CURRENT_SOURCE_DIR}/defs/*.f90"
    "${CMAKE_CURRENT_SOURCE_DIR}/main/params/*.f90"
    "${CMAKE_CURRENT_SOURCE_DIR}/main/exec/*.f90"
    "${CMAKE_CURRENT_SOURCE_DIR}/main/commanders/*.f90"
    "${CMAKE_CURRENT_SOURCE_DIR}/main/stream/*.f90")
add_custom_command(
    OUTPUT "${SIMPLE_UI_DEFAULTS_FILE}"
    COMMAND "${Python3_EXECUTABLE}" "${CMAKE_SOURCE_DIR}/scripts/default_audit.py"
            --source-root "${CMAKE_SOURCE_DIR}"
            --output-fortran "${SIMPLE_UI_DEFAULTS_FILE}"
    DEPENDS
        "${CMAKE_SOURCE_DIR}/scripts/default_audit.py"
        ${SIMPLE_UI_DEFAULTS_SOURCES}
    COMMENT "Generating UI display-default lookup"
    VERBATIM
)
add_custom_target(simple_ui_defaults DEPENDS "${SIMPLE_UI_DEFAULTS_FILE}")

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
file(GLOB STARFILE_src           "extlibs/starfile/*.cpp")
file(GLOB UNIX_src               "extlibs/unix/src/*.[Ff][90][08]")
file(GLOB CURL_src               "extlibs/curl/src/*.[Ff][90][08]" "extlibs/curl/src/*.c")
file(GLOB XML_src                "extlibs/xml/*/*.[Ff][90][08]")
file(GLOB CPLOT2D_src            "extlibs/CPlot2D/*.cpp" "extlibs/CPlot2D/*.[Ff][90][08]")
file(GLOB_RECURSE UTILS_src      "utils/*.[Ff][90][08]")
file(GLOB FILEIO_src             "fileio/*.[Ff][90][08]")
file(GLOB_RECURSE SIMPLE_src     "main/*.[Ff][90][08]")

# JPG sources including C
file(GLOB JPG_src "extlibs/jpg/*.[Ff][90][08]")
list(APPEND JPG_src "extlibs/jpg/stb_image_write.c")
list(APPEND JPG_src "extlibs/tiff/libtiff_funs.c")

# FILEIO C sources


list(APPEND FILEIO_src "fileio/simple_posix.c")

# Keep warnings from bundled third-party sources out of SIMPLE's warning
# budget. Project-owned wrappers (for example fileio/simple_tifflib.f90) retain
# the normal warning policy so interoperability defects remain visible.
set(SIMPLE_VENDOR_SOURCES
    ${GNUFOR_src}
    ${JSON_src}
    ${NETLIB_src}
    ${STARFILE_src}
    ${UNIX_src}
    ${CURL_src}
    ${XML_src}
    ${CPLOT2D_src}
    ${JPG_src}
)
set_source_files_properties(${SIMPLE_VENDOR_SOURCES} PROPERTIES COMPILE_OPTIONS "-w")

# gfortran diagnoses private state in a parent module as unused even when its
# submodules consume it. The offload reconstructor similarly keeps debug state
# for a configuration selected by preprocessing. Suppress this warning only on
# those translation units so genuinely unused module state elsewhere remains
# visible to fix_warnings.py.
set(SIMPLE_CROSS_UNIT_STATE_SOURCES
    "${CMAKE_CURRENT_SOURCE_DIR}/main/class/simple_classaverager.f90"
    "${CMAKE_CURRENT_SOURCE_DIR}/main/image/simple_image.f90"
    "${CMAKE_CURRENT_SOURCE_DIR}/main/nu_filt/simple_nu_filter.f90"
    "${CMAKE_CURRENT_SOURCE_DIR}/main/volume/simple_reconstructor_openmpoffload.f90"
)
set_source_files_properties(${SIMPLE_CROSS_UNIT_STATE_SOURCES}
    PROPERTIES COMPILE_OPTIONS "-Wno-unused-value")

# These sources contain declarations shared by mutually exclusive CPU and GPU
# preprocessor branches. An ordinary debug build necessarily sees one branch's
# declarations as unused; retain them for the alternate build configuration.
set(SIMPLE_CONDITIONAL_BRANCH_SOURCES
    "${CMAKE_CURRENT_SOURCE_DIR}/main/pftc/simple_polarft_corr.f90"
    "${CMAKE_CURRENT_SOURCE_DIR}/main/volume/simple_reconstructor_openmpoffload.f90"
    "${CMAKE_CURRENT_SOURCE_DIR}/utils/simple_gpu_utils.f90"
)
set_property(SOURCE ${SIMPLE_CONDITIONAL_BRANCH_SOURCES} APPEND
    PROPERTY COMPILE_OPTIONS "-Wno-unused-variable")

# Private helpers in these modules are retained as test entry points, callback
# implementations, or optional workflow hooks even when the current executable
# set does not reference them directly.
set(SIMPLE_INTENTIONAL_HELPER_SOURCES
    "${CMAKE_CURRENT_SOURCE_DIR}/fileio/simple_class_sample_io.f90"
    "${CMAKE_CURRENT_SOURCE_DIR}/main/commanders/test/simple_commanders_test_fft.f90"
    "${CMAKE_CURRENT_SOURCE_DIR}/main/flex/simple_flex_projected_latent_model.f90"
    "${CMAKE_CURRENT_SOURCE_DIR}/main/image/simple_image.f90"
    "${CMAKE_CURRENT_SOURCE_DIR}/main/motion/simple_motion_align_hybrid.f90"
    "${CMAKE_CURRENT_SOURCE_DIR}/main/motion/simple_opt_image_weights.f90"
    "${CMAKE_CURRENT_SOURCE_DIR}/main/ori/simple_ori.f90"
    "${CMAKE_CURRENT_SOURCE_DIR}/main/ori/simple_ori_utils.f90"
    "${CMAKE_CURRENT_SOURCE_DIR}/main/stream/simple_stream_cluster2D_utils.f90"
    "${CMAKE_CURRENT_SOURCE_DIR}/main/stream/simple_stream_p03_initial_analysis.f90"
    "${CMAKE_CURRENT_SOURCE_DIR}/main/stream/simple_stream_p04_refpick_extract_new.f90"
    "${CMAKE_CURRENT_SOURCE_DIR}/main/stream/simple_stream_pool2D_utils.f90"
    "${CMAKE_CURRENT_SOURCE_DIR}/utils/comm/simple_ipc_tcp_socket_tester.f90"
    "${CMAKE_CURRENT_SOURCE_DIR}/utils/persistent_worker/simple_persistent_worker_server_tester.f90"
    "${CMAKE_CURRENT_SOURCE_DIR}/utils/simple_private_prgs.f90"
    "${CMAKE_CURRENT_SOURCE_DIR}/utils/simple_segmentation.f90"
)
set_property(SOURCE ${SIMPLE_INTENTIONAL_HELPER_SOURCES} APPEND
    PROPERTY COMPILE_OPTIONS "-Wno-unused-function")

# FoX exports a symbol named len with a declared type. Importing its umbrella
# module triggers a GNU diagnostic about the intrinsic even though SIMPLE does
# not redeclare or call that symbol.
set(SIMPLE_FOX_IMPORT_SOURCES
    "${CMAKE_CURRENT_SOURCE_DIR}/main/star/simple_relion.f90"
    "${CMAKE_CURRENT_SOURCE_DIR}/main/star/simple_starproject.f90"
    "${CMAKE_CURRENT_SOURCE_DIR}/main/star/simple_starproject_stream.f90"
    "${CMAKE_CURRENT_SOURCE_DIR}/main/strategies/parallelization/simple_preprocess_strategy.f90"
)
set_property(SOURCE ${SIMPLE_FOX_IMPORT_SOURCES} APPEND
    PROPERTY COMPILE_OPTIONS "-Wno-surprising")

# GNU 16 emits descriptor-flow false positives around allocatable assignment
# and move_alloc in these routines. The descriptors are initialized by the
# language runtime before the reported statements execute.
set(SIMPLE_ALLOCATABLE_DESCRIPTOR_SOURCES
    "${CMAKE_CURRENT_SOURCE_DIR}/fileio/simple_imghead.f90"
    "${CMAKE_CURRENT_SOURCE_DIR}/main/star/simple_starproject_stream.f90"
    "${CMAKE_CURRENT_SOURCE_DIR}/main/stream/simple_stream_watcher.f90"
)
set_property(SOURCE ${SIMPLE_ALLOCATABLE_DESCRIPTOR_SOURCES} APPEND
    PROPERTY COMPILE_OPTIONS "-Wno-maybe-uninitialized")

# Add Windows-only stubs
if(WIN32)
    list(APPEND FILEIO_src "fileio/symlink_stub.c")
    list(APPEND FILEIO_src "fileio/posix_stubs_win32.c")
endif()

# Generated simple_args.f90
list(APPEND SIMPLE_src "${SIMPLE_GEN_DIR}/simple_args.f90")
list(APPEND SIMPLE_src "${SIMPLE_UI_DEFAULTS_FILE}")

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
    ${CURL_src}
    ${JPG_src}
    ${SIMPLE_src}
)

# ------------------------------------------------------------------------------
# Library target
# ------------------------------------------------------------------------------

add_library(${SIMPLELIB} STATIC ${SIMPLE_ALL_SOURCES})

add_dependencies(${SIMPLELIB} simple_args simple_ui_defaults simple_git)
if(SIMPLE_DEPENDENCY_TARGETS)
    add_dependencies(${SIMPLELIB} ${SIMPLE_DEPENDENCY_TARGETS})
endif()

set_target_properties(${SIMPLELIB} PROPERTIES
    VERSION                 ${PROJECT_VERSION}
    SOVERSION               ${PROJECT_VERSION_MAJOR}
    Fortran_MODULE_DIRECTORY "${CMAKE_Fortran_MODULE_DIRECTORY}"
    POSITION_INDEPENDENT_CODE ON
)

target_link_libraries(${SIMPLELIB}
    PUBLIC
        ${SIMPLE_LIBRARIES}
        $<$<PLATFORM_ID:Windows>:ws2_32>
)

# GNU Fortran module files and derived-type initialization metadata are not safe
# to mix between normal and -fcoarray=lib compilation modes. Keep the primary
# library in normal mode for controllers and unit tests, and build a separate
# consistently coarray-enabled copy only for coarray image executables.
if(USE_COARRAYS)
    file(MAKE_DIRECTORY "${SIMPLE_COARRAY_MODULE_DIRECTORY}")

    add_library(${SIMPLE_COARRAYLIB} STATIC ${SIMPLE_ALL_SOURCES})

    add_dependencies(${SIMPLE_COARRAYLIB} simple_args simple_ui_defaults simple_git)
    if(SIMPLE_DEPENDENCY_TARGETS)
        add_dependencies(${SIMPLE_COARRAYLIB} ${SIMPLE_DEPENDENCY_TARGETS})
    endif()

    set_target_properties(${SIMPLE_COARRAYLIB} PROPERTIES
        VERSION                  ${PROJECT_VERSION}
        SOVERSION                ${PROJECT_VERSION_MAJOR}
        Fortran_MODULE_DIRECTORY "${SIMPLE_COARRAY_MODULE_DIRECTORY}"
        POSITION_INDEPENDENT_CODE ON
    )

    target_link_libraries(${SIMPLE_COARRAYLIB}
        PUBLIC
            ${SIMPLE_LIBRARIES}
            Coarray::Coarray
            $<$<PLATFORM_ID:Windows>:ws2_32>
    )

    target_include_directories(${SIMPLE_COARRAYLIB}
        PUBLIC
            $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/inc>
            $<BUILD_INTERFACE:${SIMPLE_COARRAY_MODULE_DIRECTORY}>
            $<INSTALL_INTERFACE:include>
            $<BUILD_INTERFACE:${SIMPLE_GEN_DIR}>
    )
endif()

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
