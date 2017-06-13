project(simple_gui)
enable_language(CXX)

add_executable(simple_ws src/simple_ws.cpp
  src/lodepng.cpp
  src/base64.cpp)

set_target_properties(simple_ws PROPERTIES
  COMPILE_FLAGS "${OpenMP_Fortran_FLAGS} -O3 -DSIMPLE_DIR=${CMAKE_BINARY_DIR}"
  LINK_FLAGS "${OpenMP_Fortran_FLAGS} "
  LANGUAGE CXX)

# POSIX threads

# SET_COMPILE_FLAG(CMAKE_CXX_FLAGS_RELEASE "${CMAKE_CXX_FLAGS_RELEASE}"
#   CXX
#   "-mp"                 # PGI
#   "-pthreads"           # GNU
#   "-threads"            # Intel
#   )


add_executable(simple src/simple.cpp
  src/lodepng.cpp
  src/base64.cpp
  )
set_target_properties(simple PROPERTIES
  COMPILE_FLAGS "${OpenMP_Fortran_FLAGS} -O3 -DSIMPLE_DIR=${CMAKE_BINARY_DIR}"
  LINK_FLAGS "${OpenMP_Fortran_FLAGS}"
  LANGUAGE CXX)

install(DIRECTORY DESTINATION  ${CMAKE_INSTALL_PREFIX}/gui/bin)
install(TARGETS simple_ws simple RUNTIME DESTINATION ${CMAKE_INSTALL_PREFIX}/gui/bin)
if(MACOSX)
  install(FILES src/ext/websocketd-mac DESTINATION ${CMAKE_INSTALL_PREFIX}/gui/bin RENAME websocketd)
else()
  install(FILES src/ext/websocketd-linux DESTINATION ${CMAKE_INSTALL_PREFIX}/gui/bin RENAME websocketd)
endif()
install(DIRECTORY etc www DESTINATION ${CMAKE_INSTALL_PREFIX}/gui/)
#  PERMISSIONS OWNER_EXECUTE OWNER_WRITE OWNER_READ GROUP_EXECUTE GROUP_READ)

