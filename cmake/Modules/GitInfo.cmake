set(SIMPLE_GIT_VERSION "not known")
if(EXISTS ${CMAKE_CURRENT_SOURCE_DIR}/.git)
  find_package(Git)
  
  if(GIT_FOUND)
    EXECUTE_PROCESS(
      COMMAND ${GIT_EXECUTABLE} describe --abbrev=8  --always
      WORKING_DIRECTORY "${CMAKE_CURRENT_SOURCE_DIR}"
      OUTPUT_VARIABLE "SIMPLE_GIT_VERSION"
      ERROR_QUIET
      OUTPUT_STRIP_TRAILING_WHITESPACE)
  else(GIT_FOUND)
    #MESSAGE( STATUS "Git not found: ${SIMPLE_GIT_VERSION}" )
  endif(GIT_FOUND)
  set(SIMPLE_GIT_VERSION "git${SIMPLE_GIT_VERSION}")
else()
  set(SIMPLE_GIT_VERSION "-${SIMPLE_VERSION}")
endif(EXISTS ${CMAKE_CURRENT_SOURCE_DIR}/.git)

message( STATUS "Git version: ${SIMPLE_GIT_VERSION}" )
message( STATUS "Git version template : ${CMAKE_SOURCE_DIR}/cmake/Modules/GitVersion.h.in" )
message( STATUS "Git version template : ${CMAKE_INSTALL_PREFIX}/lib/simple/SimpleGitVersion.h")
configure_file(${CMAKE_SOURCE_DIR}/cmake/Modules/GitVersion.h.in  ${CMAKE_BINARY_DIR}/lib/simple/SimpleGitVersion.h @ONLY)
