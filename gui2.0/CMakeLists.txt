project(simple_gui)

cmake_minimum_required(VERSION 3.9)
set(BUILD_GUI TRUE)

if(${BUILD_GUI})
	find_program(NODE_PATH node)
	find_program(NPM_PATH npm)

	if(${NODE_PATH} STREQUAL "NODE_PATH-NOTFOUND" )
		message(FATAL_ERROR "Could not locate an executable for node. Make sure node is installed or build without gui using the flag")
	else()
		message("-- Located node executable at ${NODE_PATH}")
	endif()

	if(${NPM_PATH} STREQUAL "NPM_PATH-NOTFOUND" )
                message(FATAL_ERROR "Could not locate an executable for npm. Make sure npm is installed or build without gui using the flag")
	 else()
                message("-- Located npm executable at ${NPM_PATH}")
        endif()

	message("-- gui build environment is complete")
	
	add_custom_target(
		GUI ALL
		COMMAND ${CMAKE_COMMAND} -E copy_directory ${CMAKE_SOURCE_DIR}/src ${CMAKE_CURRENT_BINARY_DIR}/src
		COMMAND ${CMAKE_COMMAND} -E copy_directory ${CMAKE_SOURCE_DIR}/external ${CMAKE_CURRENT_BINARY_DIR}/external
		COMMAND ${CMAKE_COMMAND} -E copy ${CMAKE_SOURCE_DIR}/tsconfig.json ${CMAKE_CURRENT_BINARY_DIR}
		COMMAND ${CMAKE_COMMAND} -E copy ${CMAKE_SOURCE_DIR}/package.json ${CMAKE_CURRENT_BINARY_DIR}
		COMMAND ${CMAKE_COMMAND} -E copy ${CMAKE_SOURCE_DIR}/users.htpasswd ${CMAKE_CURRENT_BINARY_DIR}
		COMMAND ${CMAKE_COMMAND} -E copy ${CMAKE_SOURCE_DIR}/simple.sqlite ${CMAKE_CURRENT_BINARY_DIR}
		COMMAND ${CMAKE_COMMAND} -E copy ${CMAKE_SOURCE_DIR}/install ${CMAKE_CURRENT_BINARY_DIR}
		COMMAND ${NPM_PATH} install

	)
	
endif()



