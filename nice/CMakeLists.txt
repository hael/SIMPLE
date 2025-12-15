# ------------------------------------------------------------------------------
# NICE (optional Django-based service)
# ------------------------------------------------------------------------------
option(NICE "Enable NICE service" OFF)
if(NOT NICE)
  return()
endif()

# ------------------------------------------------------------------------------
# Python
# ------------------------------------------------------------------------------
find_package(Python3 3.10 REQUIRED COMPONENTS Interpreter)

# ------------------------------------------------------------------------------
# Paths
# ------------------------------------------------------------------------------
set(NICE_BUILD_DIR   "${CMAKE_CURRENT_BINARY_DIR}")
set(NICE_VENV_DIR    "${NICE_BUILD_DIR}/venv")

# ------------------------------------------------------------------------------
# Copy sources and substitute configuration
# ------------------------------------------------------------------------------
add_custom_target(nice_src ALL
  COMMENT "Preparing NICE source tree"
  COMMAND ${CMAKE_COMMAND} -E copy_directory_if_different
          "${CMAKE_CURRENT_SOURCE_DIR}/example_templates"
          "${NICE_BUILD_DIR}/example_templates"
  COMMAND ${CMAKE_COMMAND} -E copy_directory_if_different
          "${CMAKE_CURRENT_SOURCE_DIR}/example_systemd"
          "${NICE_BUILD_DIR}/example_systemd"
  COMMAND ${CMAKE_COMMAND} -E copy_directory_if_different
          "${CMAKE_CURRENT_SOURCE_DIR}/example_proxy"
          "${NICE_BUILD_DIR}/example_proxy"
  COMMAND ${CMAKE_COMMAND} -E copy_directory_if_different
          "${CMAKE_CURRENT_SOURCE_DIR}/nice_project"
          "${NICE_BUILD_DIR}/nice_project"
  COMMAND ${CMAKE_COMMAND} -E copy_directory_if_different
          "${CMAKE_CURRENT_SOURCE_DIR}/nice_lite"
          "${NICE_BUILD_DIR}/nice_lite"
  COMMAND ${CMAKE_COMMAND} -E copy_directory_if_different
          "${CMAKE_CURRENT_SOURCE_DIR}/test_display"
          "${NICE_BUILD_DIR}/test_display"
  COMMAND ${CMAKE_COMMAND} -E copy
          "${CMAKE_CURRENT_SOURCE_DIR}/manage.py"
          "${NICE_BUILD_DIR}/manage.py"
  COMMAND ${CMAKE_COMMAND}
          -D NICE_FQDN=${NICE_FQDN}
          -D SIMPLE_PATH=${CMAKE_BINARY_DIR}
          -D SRC_PATH=${CMAKE_CURRENT_SOURCE_DIR}
          -P "${CMAKE_CURRENT_LIST_DIR}/substitute_nice_src.cmake"
)

# ------------------------------------------------------------------------------
# Python virtual environment
# ------------------------------------------------------------------------------
set(NICE_PY_PACKAGES
  django==5.2.5
  pytailwindcss==0.2.0
  dirsync==2.2.6
)

if(NOT APPLE)
  list(APPEND NICE_PY_PACKAGES
    pysqlite3-binary==0.5.4
    gunicorn==23.0.0
  )
endif()

add_custom_target(nice_venv ALL
  COMMENT "Creating NICE Python virtual environment"
  COMMAND ${Python3_EXECUTABLE} -m venv "${NICE_VENV_DIR}"
  COMMAND "${NICE_VENV_DIR}/bin/python3" -m pip install --upgrade pip
  COMMAND "${NICE_VENV_DIR}/bin/pip3" install
          --isolated
          --require-virtualenv
          --no-cache-dir
          --no-input
          ${NICE_PY_PACKAGES}
)

# ------------------------------------------------------------------------------
# Build NICE assets and database
# ------------------------------------------------------------------------------
add_custom_target(nice_build ALL
  COMMENT "Building NICE web assets"
  DEPENDS nice_src nice_venv
  COMMAND "${NICE_VENV_DIR}/bin/tailwindcss"
          -i nice_lite/templates/tailwind.css
          -o nice_lite/static/nice_lite/nice_lite.css
          --minify
  COMMAND "${NICE_VENV_DIR}/bin/python3"
          "${NICE_BUILD_DIR}/manage.py" makemigrations
  COMMAND "${NICE_VENV_DIR}/bin/python3"
          "${NICE_BUILD_DIR}/manage.py" migrate
  COMMAND "${NICE_VENV_DIR}/bin/python3"
          "${NICE_BUILD_DIR}/manage.py" collectstatic --no-input
)
