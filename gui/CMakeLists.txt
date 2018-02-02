# cmake_minimum_required(VERSION 3.2 FATAL_ERROR)

# project (simple_gui)

# include(ExternalProject)

# ExternalProject_Add(
#     libjpeg
#     URL ${CMAKE_CURRENT_SOURCE_DIR}/ext/jpegsrc.v9b.tar.gz
#     SOURCE_DIR ${CMAKE_BINARY_DIR}/ext/src/jpeg
#     CONFIGURE_COMMAND ${CMAKE_BINARY_DIR}/ext/src/jpeg/configure --prefix=${CMAKE_BINARY_DIR}/ext --enable-shared=no
#     LOG_DOWNLOAD 1
#     LOG_UPDATE 1
#     LOG_CONFIGURE 1
#     LOG_BUILD 1
#     LOG_INSTALL 1
# )

# ExternalProject_Add(
#     libsqlite3
#     URL ${CMAKE_CURRENT_SOURCE_DIR}/ext/sqlite-autoconf-3210000.tar.gz
#     SOURCE_DIR ${CMAKE_BINARY_DIR}/ext/src/sqlite
#     CONFIGURE_COMMAND ${CMAKE_BINARY_DIR}/ext/src/sqlite/configure --prefix=${CMAKE_BINARY_DIR}/ext --enable-shared=no
#     LOG_DOWNLOAD 1
#     LOG_UPDATE 1
#     LOG_CONFIGURE 1
#     LOG_BUILD 1
#     LOG_INSTALL 1
# )

# add_subdirectory (ext)
# add_subdirectory (src)

##########################################################################
# Source for the SIMPLE GUI.
##########################################################################
INCLUDE(CMakeForceCompiler)
CMAKE_FORCE_C_COMPILER(gcc GNU)
CMAKE_FORCE_CXX_COMPILER(g++ GNU)

project(simple-gui CXX)

find_package (SQLite3 REQUIRED)
if (SQLITE3_FOUND)
  include_directories(${SQLITE3_INCLUDE_DIRS})
else()
    message(FATAL_ERROR "
FATAL: SQLite3 libraries and headers not found.
       Install libsqlite3-dev (Ubuntu/Debian) or relevant package.
")
endif (SQLITE3_FOUND)

find_package (JPEG REQUIRED)
if (JPEG_FOUND)
  include_directories(${JPEG_INCLUDE_DIR})
else()
    message(FATAL_ERROR "
FATAL: JPEG libraries and headers not found.
       Install libjpeg-dev (Ubuntu/Debian) or relevant package.
")
endif (JPEG_FOUND)

set (SIMPLE_GUI_SOURCES
  ${CMAKE_CURRENT_SOURCE_DIR}/src/simple_gui.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/src/mongoose.cpp  )
set_source_files_properties(${SIMPLE_GUI_SOURCES} PROPERTIES LANGUAGE CXX )

include_directories(${CMAKE_CURRENT_SOURCE_DIR}/src/)

add_executable(simple_gui ${SIMPLE_GUI_SOURCES})
target_link_libraries(simple_gui  LINK_PUBLIC ${SQLITE3_LIBRARIES})
target_link_libraries(simple_gui  LINK_PUBLIC ${JPEG_LIBRARIES})
set_target_properties(simple_gui
   PROPERTIES
   OUTPUT_NAME simple
   LANGUAGE CXX
   COMPILE_FLAGS "-std=c++11 -DSIMPLE_DIR=${CMAKE_INSTALL_PREFIX} ${OpenMP_CXX_FLAGS} -O3 "
   LINK_FLAGS "${OpenMP_CXX_FLAGS} -pthread -lm"
   JOB_POOL_COMPILE NUM_JOBS
   )


install(TARGETS simple_gui
  ARCHIVE DESTINATION "${CMAKE_INSTALL_BINDIR}"
  RUNTIME DESTINATION "${CMAKE_INSTALL_BINDIR}"
  CONFIGURATIONS DEBUG|RELEASE
  )

install(DIRECTORY ${CMAKE_SOURCE_DIR}/gui/www
  DESTINATION ${CMAKE_INSTALL_PREFIX}/
  FILE_PERMISSIONS  OWNER_EXECUTE OWNER_WRITE OWNER_READ
         GROUP_EXECUTE GROUP_READ WORLD_READ WORLD_EXECUTE
  DIRECTORY_PERMISSIONS OWNER_EXECUTE OWNER_WRITE OWNER_READ
         GROUP_EXECUTE GROUP_READ GROUP_WRITE
         WORLD_READ WORLD_EXECUTE
  CONFIGURATIONS DEBUG|RELEASE
  )
