# CPack Generator
# common package information
set(CPACK_PACKAGE_NAME                     ${PROJECT_NAME})
set(CPACK_PACKAGE_NAME                     SIMPLE)
set(CPACK_PACKAGE_VENDOR                   "Elmlund Lab -- CSB, NCI, NIH -- 2024")
set(CPACK_PACKAGE_CONTACT                  "Hans Elmlund <Hans.Elmlund@nih.gov>")
set(CPACK_PACKAGE_HOMEPAGE_URL             "https://simplecryoem.com")
set(CPACK_PACKAGE_DESCRIPTION_SUMMARY      "SIMPLE is a SIngle-particle cryo electron Microscope Image Processing Engine")
set(CPACK_DEBIAN_PACKAGE_DESCRIPTION       ${CPACK_PACKAGE_DESCRIPTION_SUMMARY})
set(CPACK_PACKAGE_INSTALL_DIRECTORY        ${CPACK_PACKAGE_NAME})
set(CPACK_PACKAGE_VERSION                  "${${PROJECT_NAME}_VERSION}")
set(CPACK_PACKAGE_VERSION                  ${SIMPLE_VERSION})
set(CPACK_PACKAGE_VERSION_MAJOR            ${SIMPLE_VERSION_MAJOR})
set(CPACK_PACKAGE_VERSION_MINOR            ${SIMPLE_VERSION_MINOR})
set(CPACK_PACKAGE_VERSION_PATCH            ${SIMPLE_GIT_VERSION})
set(CPACK_PACKAGE_VERSION_MAJOR            ${PROJECT_VERSION_MAJOR})
set(CPACK_PACKAGE_VERSION_MINOR            ${PROJECT_VERSION_MINOR})
set(CPACK_PACKAGE_VERSION_PATCH            ${PROJECT_VERSION_PATCH})
set(CPACK_PACKAGE_FILE_NAME                ${CMAKE_PACKAGE_NAME}-${CPACK_PACKAGE_VERSION}-${CPACK_SYSTEM_NAME})
set(CPACK_PACKAGE_ICON                     "${CMAKE_SOURCE_DIR}/doc/SimpleManual/SIMPLE_logo")
set(CPACK_PACKAGE_CHECKSUM                 ${CMAKE_PACKAGE_FILE_NAME}-${CPACK_PACKAGE_CHECKSUM})
set(CPACK_RESOURCE_FILE_WELCOME            "${CMAKE_SOURCE_DIR}/WELCOME")
set(CPACK_RESOURCE_FILE_LICENSE            "${CMAKE_SOURCE_DIR}/LICENSE")
set(CPACK_PACKAGE_DESCRIPTION_FILE         "${CMAKE_SOURCE_DIR}/DESCRIPTION")
set(CPACK_RESOURCE_FILE_README             "${CMAKE_SOURCE_DIR}/README.txt")
set(CPACK_VERBATIM_VARIABLES               TRUE)
set(CPACK_SET_DESTDIR                      "/usr/local")
set(CPACK_INSTALL_PREFIX                   "${CMAKE_INSTALL_PREFIX}")
set(CPACK_OUTPUT_FILE_PREFIX               packages)
set(CPACK_PACKAGE_RELOCATABLE              TRUE)
set(CPACK_MONOLITHIC_INSTALL               TRUE)
set(CPACK_PACKAGE_INSTALL_DIRECTORY        "simple ${SIMPLE_VERSION}.${SIMPLE_GIT_VERSION}")
set(CPACK_SOURCE_PACKAGE_FILE_NAME         "simple-${CPACK_PACKAGE_VERSION_MAJOR}.${CPACK_PACKAGE_VERSION_MINOR}.${CPACK_PACKAGE_VERSION_PATCH}")
# Prefix Debug/Nightly release
set(CPACK_PACKAGE_FILE_NAME                "${CPACK_PACKAGE_NAME}-${CPACK_PACKAGE_VERSION_MAJOR}.${CPACK_PACKAGE_VERSION_MINOR}")
if(NIGHTLY)
    set(CPACK_PACKAGE_FILE_NAME "${CPACK_PACKAGE_FILE_NAME}-nightly")
    #execute_process(COMMAND "date +%Y.%m.%d" OUTPUT_VARIABLE NIGHTLY_DATE)
endif()
if(CMAKE_BUILD_TYPE MATCHES "Debug")
    set(CPACK_PACKAGE_FILE_NAME "${CPACK_PACKAGE_FILE_NAME}-dbg")
endif()
# Useful descriptions for components
set(CPACK_COMPONENT_LIBRARIES_DISPLAY_NAME "SIMPLE cryoEM library")
set(CPACK_COMPONENT_DOCUMENTATION_NAME     "Doxygen documentation")
set(CPACK_COMPONENT_HEADERS_NAME           "Developmental headers")
set(CPACK_COMPONENT_CMAKE_NAME             "CMake support")
# default settings
set(CPACK_GENERATOR_TGZ                    ON)
set(CPACK_SOURCE_TBZ2                      OFF)
set(CPACK_SOURCE_TZ                        OFF)
set(CPACK_SOURCE_TXZ                       OFF)
set(CPACK_SOURCE_ZIP                       ON)

# default package generators
if(APPLE)
    set(PACKAGE_GENERATOR "PackageMaker")
    set(PACKAGE_SOURCE_GENERATOR "TGZ;ZIP")
elseif(UNIX)
    set(PACKAGE_GENERATOR "DEB;RPM")
    set(PACKAGE_SOURCE_GENERATOR "TGZ;ZIP")
else()
    set(PACKAGE_GENERATOR "ZIP")
    set(PACKAGE_SOURCE_GENERATOR "ZIP")
endif()

if(WIN32)
    set(CPACK_GENERATOR ZIP WIX)
elseif(APPLE)
    set(CPACK_GENERATOR TGZ productbuild)
elseif(CMAKE_SYSTEM_NAME STREQUAL "Linux")
    set(CPACK_GENERATOR TGZ RPM)
else()
    set(CPACK_GENERATOR TGZ)
endif()

#--------------------------------
#get_filename_component(cpack_build_dir "${CMAKE_BINARY_DIR}" ABSOLUTE)
#get_filename_component(cpack_source_dir "${CMAKE_SOURCE_DIR}" ABSOLUTE)

# Mac OS X package
if(CPACK_GENERATOR MATCHES "PackageMaker|DragNDrop")
    set(CPACK_PACKAGE_FILE_NAME "${CPACK_PACKAGE_FILE_NAME}-${CPACK_PACKAGE_VERSION}")
    set(CPACK_PACKAGING_INSTALL_PREFIX /usr/local)
# Debian package
elseif(CPACK_GENERATOR MATCHES "DEB")
    set(CPACK_DEBIAN_PACKAGE_DEPENDS "gcc, build-essential, binutils, libfftw3-dev, gfortran, cmake, gnuplot")
    set(CPACK_DEBIAN_PACKAGE_OPTIONAL "nvidia-cuda-dev, libopenmpi-dev, libjpeg9-dev, libsqlite3-dev ")
    set(CPACK_DEBIAN_PACKAGE_SECTION "science")
    set(CPACK_DEBIAN_PACKAGE_PRIORITY "optional")
    set(CPACK_DEBIAN_PACKAGE_HOMEPAGE "https://simplecryoem.org")
    set(CPACK_DEBIAN_PACKAGE_ARCHITECTURE "amd64")
    set(CPACK_PACKAGE_FILE_NAME "${CPACK_PACKAGE_FILE_NAME}"
        "-${CPACK_PACKAGE_VERSION}"
        "-${CPACK_DEBIAN_PACKAGE_ARCHITECTURE}"
    )
# RPM package
elseif(CPACK_GENERATOR MATCHES "RPM")
    # https://github.com/pld-linux/hhvm
    set(CPACK_RPM_PACKAGE_REQUIRES "binutils-devel, cmake >= 2.8.7, "
        "gcc >= 6:4.6.0, libfftw3-devel, gnuplot   "
    )
    set(CPACK_RPM_PACKAGE_GROUP "Applications/Engineering")
    set(CPACK_RPM_PACKAGE_LICENSE "GPL3+ ")
    set(CPACK_RPM_PACKAGE_URL "https://simplecryoem.org")
    set(CPACK_RPM_PACKAGE_ARCHITECTURE "x86_64")
    set(CPACK_PACKAGE_FILE_NAME "${CPACK_PACKAGE_FILE_NAME}"
        "-${CPACK_PACKAGE_VERSION}"
        "-${CPACK_RPM_PACKAGE_ARCHITECTURE}"
    )
endif()
include(CPack)
message( STATUS "CPack setup complete")

#cpack_add_component_group(SIMPLE_SDK
#    DISPLAY_NAME SDK
#    DESCRIPTION "Developer tools, library, etc."
#    DEPENDS
#    REQUIRED
#    HIDDEN
#    GROUP
#    INSTALL_TYPES
#    DONWLOADED
#    ARCHIVE_FILE archivefilename
#    PARENT_GROUP MyProj_SDK
#)
#cpack_add_install_type(Full)
