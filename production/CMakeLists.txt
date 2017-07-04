

add_subdirectory(simple)

#################################################################
#
# ENABLE TESTING
#
#################################################################
if (${PROJECT_NAME}_BUILD_TESTS)
  if(VERBOSE)
    include(CTest)
  endif()
  add_subdirectory(tests)
endif()
