# FindJPEG
# --------
#
# Find JPEG
#
# Find the native JPEG includes and library This module defines
#
# ::
#
#   JPEG_INCLUDE_DIR, where to find jpeglib.h, etc.
#   JPEG_LIBRARIES, the libraries needed to use JPEG.
#   JPEG_FOUND, If false, do not try to use JPEG.
#
# also defined, but not for general use are
#
# ::
#
#   JPEG_LIBRARY, where to find the JPEG library.

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

find_path(JPEG_INCLUDE_DIR jpeglib.h)
find_path(JPEG_CONFIG_DIR jconfig.h)
set(JPEG_NAMES ${JPEG_NAMES} jpeg9 libjpeg.so.9)
find_library(JPEG_LIBRARY NAMES ${JPEG_NAMES} )

# handle the QUIETLY and REQUIRED arguments and set JPEG_FOUND to TRUE if
# all listed variables are TRUE
include(${CMAKE_ROOT}/Modules/FindPackageHandleStandardArgs.cmake)
find_package_handle_standard_args(JPEG DEFAULT_MSG JPEG_LIBRARY JPEG_INCLUDE_DIR)

file(STRINGS ${JPEG_CONFIG_DIR}/jconfig.h JPEGLIB_VERSION REGEX "^#define[ \t]+JPEG_LIB_VERSION")
message(STATUS "JPEG_LIB_VERSION definition: ${JPEGLIB_VERSION}")
string(REGEX REPLACE "^#define[ ]*JPEG_LIB_VERSION[ \t]+\([0-9]\)\([0-9]\)" "\\1.\\2" JPEGLIB_VERSION "${JPEGLIB_VERSION}")
message(STATUS "JPEGlib version ${JPEGLIB_VERSION} FOUND")
 if ("${JPEGLIB_VERSION}" VERSION_GREATER "8.0")
   set(JPEG_FOUND OFF)
   message(STATUS "Jpeglib.h version not equal to or greater than 8.0 " )
 endif()

if(JPEG_FOUND)
  set(JPEG_LIBRARIES ${JPEG_LIBRARY})
endif()

# Deprecated declarations.
set (NATIVE_JPEG_INCLUDE_PATH ${JPEG_INCLUDE_DIR} )
if(JPEG_LIBRARY)
  get_filename_component (NATIVE_JPEG_LIB_PATH ${JPEG_LIBRARY} PATH)
endif()

mark_as_advanced(JPEG_LIBRARY JPEG_INCLUDE_DIR )
