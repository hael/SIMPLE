# This CMake script will delete build directories and files to bring the
# package back to it's distribution state

# We want to start from the top of the source dir, so if we are in build
# we want to start one directory up
get_filename_component(BASEDIR ${CMAKE_SOURCE_DIR} NAME)
#IF(${BASEDIR} STREQUAL "build")
#    SET(TOPDIR "${CMAKE_SOURCE_DIR}/..")
#ELSE()
set(TOPDIR "${CMAKE_SOURCE_DIR}")
#ENDIF()
message(STATUS "TOPDIR ${TOPDIR}")
macro(GET_PARENT_DIRECTORIES search_string return_list grandparents)
    file(GLOB_RECURSE new_list ${search_string})
    set(dir_list "")
    foreach(file_path ${new_list})
        get_filename_component(dir_path ${file_path} PATH)
        # Remove an extra directory component to return grandparent
        if(${grandparents})
            # Tack on a fake extension to trick CMake into removing a second
            # path component
            set(dir_path "${dir_path}.tmp")
            get_filename_component(dir_path ${dir_path} PATH)
        endif(${grandparents})
        set(dir_list ${dir_list} ${dir_path})
    endforeach()
    list(REMOVE_DUPLICATES dir_list)
    set(${return_list} ${dir_list})
endmacro()

# Find directories and files that we will want to remove
file(GLOB_RECURSE CMAKECACHE "${TOPDIR}/*CMakeCache.txt")
file(GLOB_RECURSE CMAKEINSTALL "${TOPDIR}/*cmake_install.cmake"
                               "${TOPDIR}/*install_manifest.txt"
)
file(GLOB_RECURSE MAKEFILE "${TOPDIR}/*Makefile")
file(GLOB_RECURSE MODFILES RELATIVE "${TOPDIR}" "*.mod")
file(GLOB_RECURSE CMAKETESTFILES "${TOPDIR}/*CTestTestfile.cmake")
set(TOPDIRECTORIES "${TOPDIR}/cmake" "${TOPDIR}/tests" "${TOPDIR}/bin/tests" "${TOPDIR}/bin/simple")

# CMake has trouble finding directories recursively, so locate these
# files and then save the parent directory of the files
get_parent_directories(Makefile.cmake CMAKEFILES 0)
get_parent_directories(CMakeDirectoryInformation.cmake SUBDIRECTORIES 0) 
get_parent_directories(LastTest.log CMAKETESTING 1)

# Place these files and directories into a list
set(DEL ${CMAKECACHE}
        ${CMAKEINSTALL}
        ${MAKEFILE}
        ${CMAKEFILES}
        ${CMAKETESTING}
        ${CMAKETESTFILES}
        ${SUBDIRECTORIES}
        ${TOPDIRECTORIES}
)

# If we are not in the build dir, delete that as well
if(NOT (${BASEDIR} STREQUAL "build"))
    file(GLOB BUILD "${TOPDIR}/build")
    set(DEL ${DEL} ${BUILD})
endif()

# Loop over the directories and delete each one
foreach(D ${DEL})
    if(EXISTS ${D})
        message(STATUS "Removing ${D}")
        file(REMOVE_RECURSE ${D})
    endif()
endforeach()
