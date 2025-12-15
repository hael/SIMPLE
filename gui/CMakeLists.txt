project(simple_gui LANGUAGES NONE)

include(ExternalProject)

# ------------------------------------------------------------------------------
# Options
# ------------------------------------------------------------------------------
option(GUI "Build and install GUI components" OFF)
option(BUILD_GUI "Build GUI from source (Linux)" OFF)
if(NOT GUI)
  return()
endif()

# ------------------------------------------------------------------------------
# Paths
# ------------------------------------------------------------------------------
set(GUI_DATA_DIR    "${CMAKE_BINARY_DIR}/gui_data")
set(GUI_BIN_DIR     "${CMAKE_BINARY_DIR}/bin")
set(NODEJS_DIR      "${GUI_BIN_DIR}/nodejs")
set(TUTORIALS_DIR   "${CMAKE_BINARY_DIR}/tutorials")

# ------------------------------------------------------------------------------
# Core GUI data (always installed when GUI=ON)
# ------------------------------------------------------------------------------
add_custom_target(gui_core ALL
  COMMENT "Installing GUI core assets"
  COMMAND ${CMAKE_COMMAND} -E copy_directory
          "${CMAKE_CURRENT_SOURCE_DIR}/client"
          "${GUI_DATA_DIR}/client"
  COMMAND ${CMAKE_COMMAND} -E copy_directory
          "${CMAKE_CURRENT_SOURCE_DIR}/server"
          "${GUI_DATA_DIR}/server"
  COMMAND ${CMAKE_COMMAND} -E copy_directory
          "${CMAKE_CURRENT_SOURCE_DIR}/client/public"
          "${GUI_DATA_DIR}/public"
  COMMAND ${CMAKE_COMMAND} -E tar xf
          "${CMAKE_CURRENT_SOURCE_DIR}/tutorials.tar.gz"
  COMMAND ${CMAKE_COMMAND} -E copy_directory
          "${GUI_BIN_DIR}/gui/tutorials"
          "${TUTORIALS_DIR}"
  COMMAND ${CMAKE_COMMAND} -E remove_directory
          "${GUI_BIN_DIR}/gui/tutorials"
)
install(DIRECTORY "${GUI_DATA_DIR}" DESTINATION .)
install(DIRECTORY "${TUTORIALS_DIR}" DESTINATION .)

# ------------------------------------------------------------------------------
# Node.js (platform-specific)
# ------------------------------------------------------------------------------
if(APPLE)
  set(NODEJS_URL
    "https://nodejs.org/dist/v21.7.1/node-v21.7.1-darwin-arm64.tar.gz"
  )
  set(NODEJS_MD5 "5e85b7df5071ad757f773f94dc886d68")
elseif(UNIX)
  set(NODEJS_URL
    "https://nodejs.org/dist/v16.13.0/node-v16.13.0-linux-x64.tar.gz"
  )
  set(NODEJS_MD5 "bb20c3f845896272f18bd8c91d1c9b07")
else()
  message(FATAL_ERROR "GUI is not supported on this platform")
endif()
ExternalProject_Add(nodejs
  URL             ${NODEJS_URL}
  URL_MD5         ${NODEJS_MD5}
  DOWNLOAD_DIR    "${CMAKE_BINARY_DIR}/external/tarballs"
  SOURCE_DIR      "${NODEJS_DIR}"
  CONFIGURE_COMMAND ""
  BUILD_COMMAND     ""
  INSTALL_COMMAND   ""
  TEST_COMMAND      ""
)
install(DIRECTORY "${NODEJS_DIR}"
        DESTINATION bin
        USE_SOURCE_PERMISSIONS)

# ------------------------------------------------------------------------------
# Multi-user GUI
# ------------------------------------------------------------------------------
add_custom_target(gui_multi ALL
  COMMENT "Installing GUI multi-user support"
  DEPENDS gui_core nodejs
  COMMAND ${CMAKE_COMMAND} -E env
          "PATH=${NODEJS_DIR}/bin:$ENV{PATH}"
          npm install -g --unsafe-perm binding
  COMMAND ${CMAKE_COMMAND} -E copy
          "${CMAKE_CURRENT_SOURCE_DIR}/server/multi_user/simple_multiuser.sh"
          "${GUI_BIN_DIR}/simple_multiuser"
  WORKING_DIRECTORY
          "${GUI_DATA_DIR}/server/multi_user"
)
install(PROGRAMS
  "${GUI_BIN_DIR}/simple_multiuser"
  DESTINATION bin
)

# ------------------------------------------------------------------------------
# Single-user GUI
# ------------------------------------------------------------------------------
if(APPLE)
  set(PREBUILT_GUI
      "${CMAKE_CURRENT_SOURCE_DIR}/precompiled/simple_gui.mac.tar.gz")
  if(NOT EXISTS "${PREBUILT_GUI}")
    message(FATAL_ERROR
      "Prebuilt macOS GUI not found. Re-run cmake with -DGUI=OFF")
  endif()
  add_custom_target(gui_single ALL
    COMMENT "Installing precompiled macOS GUI"
    DEPENDS gui_core
    COMMAND ${CMAKE_COMMAND} -E tar xf "${PREBUILT_GUI}"
    COMMAND ${CMAKE_COMMAND} -E copy
            "${CMAKE_CURRENT_SOURCE_DIR}/precompiled/simple_mac.sh"
            "${GUI_BIN_DIR}/simple"
  )
  install(PROGRAMS "${GUI_BIN_DIR}/simple" DESTINATION bin)
elseif(UNIX)
  if(BUILD_GUI)
    add_custom_target(gui_single ALL
      COMMENT "Building single-user GUI from source"
      DEPENDS gui_core nodejs
      COMMAND ${CMAKE_COMMAND} -E env
              "PATH=${NODEJS_DIR}/bin:$ENV{PATH}"
              npm install
      COMMAND ${CMAKE_COMMAND} -E env
              "PATH=${NODEJS_DIR}/bin:$ENV{PATH}"
              npm run build
      COMMAND ${CMAKE_COMMAND} -E copy
              "${GUI_DATA_DIR}/server/single_user/dist/simple_gui.linux"
              "${GUI_BIN_DIR}/simple"
      WORKING_DIRECTORY
              "${GUI_DATA_DIR}/server/single_user"
    )
  else()
    set(PREBUILT_GUI
        "${CMAKE_CURRENT_SOURCE_DIR}/precompiled/simple_gui.linux")
    if(NOT EXISTS "${PREBUILT_GUI}")
      message(FATAL_ERROR
        "Prebuilt Linux GUI not found. Re-run cmake with -DGUI=OFF")
    endif()
    add_custom_target(gui_single ALL
      COMMENT "Installing precompiled Linux GUI"
      DEPENDS gui_core
      COMMAND ${CMAKE_COMMAND} -E copy
              "${PREBUILT_GUI}"
              "${GUI_BIN_DIR}/simple"
      COMMAND ${CMAKE_COMMAND} -E tar xf
              "${CMAKE_CURRENT_SOURCE_DIR}/precompiled/node_modules.linux.tar.gz"
      COMMAND ${CMAKE_COMMAND} -E copy_directory
              "${GUI_BIN_DIR}/gui/node_modules"
              "${GUI_BIN_DIR}/node_modules"
    )
    install(DIRECTORY
      "${GUI_BIN_DIR}/node_modules"
      DESTINATION bin
    )
  endif()
  install(PROGRAMS "${GUI_BIN_DIR}/simple" DESTINATION bin)
endif()
