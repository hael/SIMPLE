cmake_minimum_required(VERSION 3.2 FATAL_ERROR)

project (simple_gui)

include(ExternalProject)

ExternalProject_Add(
    libjpeg
    URL ${CMAKE_CURRENT_SOURCE_DIR}/ext/jpegsrc.v9b.tar.gz
    SOURCE_DIR ${CMAKE_BINARY_DIR}/ext/src/jpeg
    CONFIGURE_COMMAND ${CMAKE_BINARY_DIR}/ext/src/jpeg/configure --prefix=${CMAKE_BINARY_DIR}/ext --enable-shared=no
 #   BUILD_IN_SOURCE 1
    LOG_DOWNLOAD 1
    LOG_UPDATE 1
    LOG_CONFIGURE 1
    LOG_BUILD 1
    LOG_INSTALL 1
)

ExternalProject_Add(
    libsqlite3
    URL ${CMAKE_CURRENT_SOURCE_DIR}/ext/sqlite-autoconf-3210000.tar.gz
    SOURCE_DIR ${CMAKE_BINARY_DIR}/ext/src/sqlite
    CONFIGURE_COMMAND ${CMAKE_BINARY_DIR}/ext/src/sqlite/configure --prefix=${CMAKE_BINARY_DIR}/ext --enable-shared=no
 #   BUILD_IN_SOURCE 1
    LOG_DOWNLOAD 1
    LOG_UPDATE 1
    LOG_CONFIGURE 1
    LOG_BUILD 1
    LOG_INSTALL 1
)

ExternalProject_Add(
    web_root
    URL ${CMAKE_CURRENT_SOURCE_DIR}/www
    SOURCE_DIR ${CMAKE_BINARY_DIR}/www
    CONFIGURE_COMMAND ""
    BUILD_COMMAND ""
    INSTALL_COMMAND ""
 #   BUILD_IN_SOURCE 1
    LOG_DOWNLOAD 1
    LOG_UPDATE 1
    LOG_CONFIGURE 1
    LOG_BUILD 1
    LOG_INSTALL 1
)

add_subdirectory (ext)
add_subdirectory (src)
