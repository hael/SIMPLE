# Run at build time by the fortls target (scripts/CMakeLists.txt): installs the pinned fortls in a
# virtual environment of the build tree, then registers the build tree's Claude Code plugin for this
# checkout (.claude/settings.local.json, never committed). Each step is skipped when already done;
# offline or without the claude CLI the build goes on with a warning.
# Inputs: PYTHON, VENV, MARKETPLACE (the build tree's marketplace directory), SOURCE_DIR.
set(FORTLS_VERSION 3.2.2)
set(fortls_exe "${VENV}/bin/fortls")
set(installed_version "")
if(EXISTS "${fortls_exe}")
  execute_process(COMMAND "${fortls_exe}" --version OUTPUT_VARIABLE installed_version
    OUTPUT_STRIP_TRAILING_WHITESPACE ERROR_QUIET)
endif()
if(NOT installed_version STREQUAL FORTLS_VERSION)
  execute_process(COMMAND "${PYTHON}" -m venv "${VENV}" RESULT_VARIABLE rc OUTPUT_QUIET ERROR_VARIABLE err)
  if(rc EQUAL 0)
    execute_process(COMMAND "${VENV}/bin/python" -m pip install --quiet --disable-pip-version-check
      "fortls==${FORTLS_VERSION}" RESULT_VARIABLE rc OUTPUT_QUIET ERROR_VARIABLE err)
  endif()
  if(NOT rc EQUAL 0)
    message(WARNING "fortls ${FORTLS_VERSION} was not installed, so Claude Code has no Fortran code "
      "intelligence (offline? configure with -DFORTLS=OFF to skip): ${err}")
    return()
  endif()
endif()

find_program(CLAUDE_CLI claude)
if(NOT CLAUDE_CLI)
  return()
endif()
# registered already when the checkout's local settings name this build tree's marketplace
set(local_settings "${SOURCE_DIR}/.claude/settings.local.json")
if(EXISTS "${local_settings}")
  file(READ "${local_settings}" settings_text)
  string(FIND "${settings_text}" "${MARKETPLACE}" at)
  if(at GREATER -1)
    return()
  endif()
endif()
# a marketplace of the same name from another build tree is replaced
execute_process(COMMAND "${CLAUDE_CLI}" plugin marketplace remove simple-build
  WORKING_DIRECTORY "${SOURCE_DIR}" OUTPUT_QUIET ERROR_QUIET)
execute_process(COMMAND "${CLAUDE_CLI}" plugin marketplace add "${MARKETPLACE}" --scope local
  WORKING_DIRECTORY "${SOURCE_DIR}" RESULT_VARIABLE rc OUTPUT_QUIET ERROR_VARIABLE err)
if(rc EQUAL 0)
  execute_process(COMMAND "${CLAUDE_CLI}" plugin install fortls-lsp@simple-build --scope local
    WORKING_DIRECTORY "${SOURCE_DIR}" RESULT_VARIABLE rc OUTPUT_QUIET ERROR_VARIABLE err)
endif()
if(NOT rc EQUAL 0)
  message(WARNING "The Claude Code plugin fortls-lsp was not registered: ${err}")
else()
  message(STATUS "Claude Code plugin fortls-lsp registered for this checkout (start a new session to load it)")
endif()
