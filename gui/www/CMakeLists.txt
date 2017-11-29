add_custom_command(OUTPUT ${CMAKE_BINARY_DIR}/www
		   COMMAND mkdir -p ${CMAKE_BINARY_DIR}/www
		   COMMAND cp -r * ${CMAKE_BINARY_DIR}/www/
		   WORKING_DIRECTORY "${CMAKE_CURRENT_SOURCE_DIR}"
		   COMMENT "Setting up web directory ...")
add_custom_target(simple_web ALL DEPENDS ${CMAKE_BINARY_DIR}/www)
