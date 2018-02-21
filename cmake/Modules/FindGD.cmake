#.rst:
# FindGD
# -------
#
# Find libgd, the official reference library for the GD image format.
#
# Imported targets
# ^^^^^^^^^^^^^^^^
#
# This module defines the following :prop_tgt:`IMPORTED` target:
#
# ``GD::GD``
#   The libgd library, if found.
#
# Result variables
# ^^^^^^^^^^^^^^^^
#
# This module will set the following variables in your project:
#
# ``GD_INCLUDE_DIRS``
#   where to find gd.h, etc.
# ``GD_LIBRARIES``
#   the libraries to link against to use GD.
# ``GD_DEFINITIONS``
#   You should add_definitons(${GD_DEFINITIONS}) before compiling code
#   that includes gd library files.
# ``GD_FOUND``
#   If false, do not try to use GD.
# ``GD_VERSION_STRING``
#   the version of the GD library found (since CMake 2.8.8)
#
# Obsolete variables
# ^^^^^^^^^^^^^^^^^^
#
# The following variables may also be set, for backwards compatibility:
#
# ``GD_LIBRARY``
#   where to find the GD library.
# ``GD_INCLUDE_DIR``
#   where to find the GD headers (same as GD_INCLUDE_DIRS)
#
# Since GD depends on the PNG and ZLib compression library, none of the above
# will be defined unless PNG and ZLib can be found.

#=============================================================================
# Copyright 2002-2009 Kitware, Inc.
# Copyright 2016 Raumfeld
#
# Distributed under the OSI-approved BSD License (the "License");
# see accompanying file Copyright.txt for details.
#
# This software is distributed WITHOUT ANY WARRANTY; without even the
# implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
# See the License for more information.
#=============================================================================
# (To distribute this file outside of CMake, substitute the full
#  License text for the above reference.)

if(GD_FIND_QUIETLY)
  set(_FIND_PNG_ARG QUIET)
endif()
find_package(PNG ${_FIND_PNG_ARG})

if(PNG_FOUND)
  find_path(GD_GD_INCLUDE_DIR gd.h
    /usr/include
    /usr/local/include             # OpenBSD
    /opt/sw/include

  )
  message(STATUS "LibGD version path ${GD_GD_INCLUDE_DIR}")
  list(APPEND GD_NAMES gd libgd)
  # unset(GD_NAMES_DEBUG)
  # set(_GD_VERSION_SUFFIXES 31 21 20)
  # if (GD_FIND_VERSION MATCHES "^([0-9]+)\\.([0-9]+)(\\..*)?$")
  #   set(_GD_VERSION_SUFFIX_MIN "${CMAKE_MATCH_1}${CMAKE_MATCH_2}")
  #   if (GD_FIND_VERSION_EXACT)
  #     set(_GD_VERSION_SUFFIXES ${_GD_VERSION_SUFFIX_MIN})
  #   else ()
  #     string(REGEX REPLACE
  #         "${_GD_VERSION_SUFFIX_MIN}.*" "${_GD_VERSION_SUFFIX_MIN}"
  #         _GD_VERSION_SUFFIXES "${_GD_VERSION_SUFFIXES}")
  #   endif ()
  #   unset(_GD_VERSION_SUFFIX_MIN)
  # endif ()
  # foreach(v IN LISTS _GD_VERSION_SUFFIXES)
  #   list(APPEND GD_NAMES gd${v} libgd${v})
  #   list(APPEND GD_NAMES_DEBUG gd${v}d libgd${v}d)
  # endforeach()
  # unset(_GD_VERSION_SUFFIXES)
  # # For compatiblity with versions prior to this multi-config search, honor
  # # any GD_LIBRARY that is already specified and skip the search.

  find_library(GD_LIBRARY NAMES ${GD_NAMES}
    PATHS
        /usr/lib
        /usr/lib/x86_64-linux-gnu
        /usr/local/lib               # Homebrew
        /opt/local/lib               # Macports
        /usr/opt/local/lib
        /sw/lib                      # Fink
        )

   #  include(${CMAKE_ROOT}/Modules/SelectLibraryConfigurations.cmake)
   #  select_library_configurations(GD)
     mark_as_advanced(GD_LIBRARY)

   unset(GD_NAMES)

  # Set by select_library_configurations(), but we want the one from
  # find_package_handle_standard_args() below.
  unset(GD_FOUND)

  if (GD_LIBRARY AND GD_GD_INCLUDE_DIR)
      # gd.h includes png.h. Sigh.
      set(GD_INCLUDE_DIRS ${GD_GD_INCLUDE_DIR} ${PNG_INCLUDE_DIR} )
      set(GD_INCLUDE_DIR ${GD_INCLUDE_DIRS} ) # for backward compatiblity
      set(GD_LIBRARIES ${GD_LIBRARY} ${PNG_LIBRARY})

      # if (CYGWIN)
      #   if(BUILD_SHARED_LIBS)
      #      # No need to define GD_USE_DLL here, because it's default for Cygwin.
      #   else()
      #     set (GD_DEFINITIONS -DGD_STATIC)
      #   endif()
      # endif ()

      if(NOT TARGET GD::GD)
        add_library(GD::GD UNKNOWN IMPORTED)
        set_target_properties(GD::GD PROPERTIES
          INTERFACE_COMPILE_DEFINITIONS "${GD_DEFINITIONS}"
          INTERFACE_INCLUDE_DIRECTORIES "${GD_INCLUDE_DIRS}"
          INTERFACE_LINK_LIBRARIES PNG::PNG)
        if(EXISTS "${GD_LIBRARY}")
          set_target_properties(GD::GD PROPERTIES
            IMPORTED_LINK_INTERFACE_LANGUAGES "C"
            IMPORTED_LOCATION "${GD_LIBRARY}")
        endif()
        # if(EXISTS "${GD_LIBRARY_DEBUG}")
        #   set_property(TARGET GD::GD APPEND PROPERTY
        #     IMPORTED_CONFIGURATIONS DEBUG)
        #   set_target_properties(GD::GD PROPERTIES
        #     IMPORTED_LINK_INTERFACE_LANGUAGES_DEBUG "C"
        #     IMPORTED_LOCATION_DEBUG "${GD_LIBRARY_DEBUG}")
        # endif()
        # if(EXISTS "${GD_LIBRARY_RELEASE}")
        #   set_property(TARGET GD::GD APPEND PROPERTY
        #     IMPORTED_CONFIGURATIONS RELEASE)
        #   set_target_properties(GD::GD PROPERTIES
        #     IMPORTED_LINK_INTERFACE_LANGUAGES_RELEASE "C"
        #     IMPORTED_LOCATION_RELEASE "${GD_LIBRARY_RELEASE}")
        # endif()
      endif()
  endif ()

  # if (GD_GD_INCLUDE_DIR AND EXISTS "${GD_GD_INCLUDE_DIR}/gd.h")
  #   file(STRINGS "${GD_GD_INCLUDE_DIR}/gd.h" gd_version_str REGEX "^#define[ \t]+GD_.*_VERSION[ \t]+[0-9][ \t]+.*")
  #   message(STATUS "LibGD version ${gd_version_str}")
  #     string(REGEX REPLACE "^#define[ \t]+GD_.*_VERSION[ \t]+([0-9]).*" "\\1" GD_VERSION_STRING "${gd_version_str}")
  #     message(STATUS "LibGD version ${GD_VERSION_STRING}")
  #     unset(gd_version_str)
  # endif ()
endif()

# handle the QUIETLY and REQUIRED arguments and set GD_FOUND to TRUE if
# all listed variables are TRUE
include(${CMAKE_ROOT}/Modules/FindPackageHandleStandardArgs.cmake)
find_package_handle_standard_args(GD
                                  REQUIRED_VARS GD_LIBRARY GD_GD_INCLUDE_DIR )
#                                  VERSION_VAR GD_VERSION_STRING)

mark_as_advanced(GD_GD_INCLUDE_DIR GD_LIBRARY )
