# FindJPEG
# --------
#
# Find JPEG9
#
# Find the native JPEG includes and library This module defines
#
# ::
#
#   JPEG9_INCLUDE_DIR, where to find jpeglib.h, etc.
#   JPEG9_LIBRARIES, the libraries needed to use JPEG.
#   JPEG9_FOUND, If false, do not try to use JPEG.
#
# also defined, but not for general use are
#
# ::
#
#   JPEG9_LIBRARY, where to find the JPEG library.

#=============================================================================
# Copyright 2001-2009 Kitware, Inc.
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

find_path(JPEG9_INCLUDE_DIR jpeglib.h)
find_path(JPEG9_CONFIG_DIR jconfig.h)
set(JPEG9_NAMES ${JPEG9_NAMES} jpeg9 libjpeg.so.9)
find_library(JPEG9_LIBRARY NAMES ${JPEG9_NAMES} )

# handle the QUIETLY and REQUIRED arguments and set JPEG_FOUND to TRUE if
# all listed variables are TRUE
include(${CMAKE_ROOT}/Modules/FindPackageHandleStandardArgs.cmake)
find_package_handle_standard_args(JPEG9 DEFAULT_MSG JPEG9_LIBRARY JPEG9_INCLUDE_DIR)

file(STRINGS ${JPEG9_CONFIG_DIR}/jconfig.h JPEG9LIB_VERSION REGEX "^#define[ \t]+JPEG_LIB_VERSION")
message(STATUS "JPEG9_LIB_VERSION definition: ${JPEG9LIB_VERSION}")
string(REGEX REPLACE "^#define[ ]*JPEG9_LIB_VERSION[ \t]+\([0-9]\)\([0-9]\)" "\\1.\\2" JPEG9LIB_VERSION "${JPEG9LIB_VERSION}")
message(STATUS "JPEG9lib version ${JPEG9LIB_VERSION} FOUND")
 if ("${JPEG9LIB_VERSION}" VERSION_GREATER "8.0")
   set(JPEG9_FOUND OFF)
   message(STATUS "Jpeglib.h version not greater than 8.0 " )
 endif()

if(JPEG9_FOUND)
  set(JPEG9_LIBRARIES ${JPEG9_LIBRARY})
endif()

# Deprecated declarations.
set (NATIVE_JPEG9_INCLUDE_PATH ${JPEG9_INCLUDE_DIR} )
if(JPEG9_LIBRARY)
  get_filename_component (NATIVE_JPEG9_LIB_PATH ${JPEG9_LIBRARY} PATH)
endif()

mark_as_advanced(JPEG9_LIBRARY JPEG9_INCLUDE_DIR )
