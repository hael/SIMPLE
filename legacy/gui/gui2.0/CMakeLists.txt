project(simple_gui)

cmake_minimum_required(VERSION 3.2)
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
		gui ALL
		COMMAND ${CMAKE_COMMAND} -E remove_directory ${CMAKE_CURRENT_BINARY_DIR}/src
		COMMAND ${CMAKE_COMMAND} -E remove_directory ${CMAKE_CURRENT_BINARY_DIR}/dist
		COMMAND ${CMAKE_COMMAND} -E make_directory ${CMAKE_CURRENT_BINARY_DIR}/dist
		COMMAND ${CMAKE_COMMAND} -E make_directory ${CMAKE_CURRENT_BINARY_DIR}/dist/views
		COMMAND ${CMAKE_COMMAND} -E make_directory ${CMAKE_CURRENT_BINARY_DIR}/dist/public
		COMMAND ${CMAKE_COMMAND} -E make_directory ${CMAKE_CURRENT_BINARY_DIR}/dist/public/css
		COMMAND ${CMAKE_COMMAND} -E make_directory ${CMAKE_CURRENT_BINARY_DIR}/dist/public/img
		COMMAND ${CMAKE_COMMAND} -E make_directory ${CMAKE_CURRENT_BINARY_DIR}/dist/public/js
		COMMAND ${CMAKE_COMMAND} -E copy_directory ${CMAKE_SOURCE_DIR}/gui2.0/src ${CMAKE_CURRENT_BINARY_DIR}/src
		COMMAND ${CMAKE_COMMAND} -E copy_directory ${CMAKE_SOURCE_DIR}/gui2.0/external ${CMAKE_CURRENT_BINARY_DIR}/external
		COMMAND ${CMAKE_COMMAND} -E copy ${CMAKE_SOURCE_DIR}/gui2.0/tsconfig.json ${CMAKE_CURRENT_BINARY_DIR}
		COMMAND ${CMAKE_COMMAND} -E copy ${CMAKE_SOURCE_DIR}/gui2.0/package.json ${CMAKE_CURRENT_BINARY_DIR}
		COMMAND ${CMAKE_COMMAND} -E copy ${CMAKE_SOURCE_DIR}/gui2.0/binding.gyp ${CMAKE_CURRENT_BINARY_DIR}
		COMMAND ${CMAKE_COMMAND} -E copy ${CMAKE_SOURCE_DIR}/gui2.0/concat.sh ${CMAKE_CURRENT_BINARY_DIR}
		COMMAND ${NPM_PATH} install
		COMMAND ${CMAKE_COMMAND} -E copy_directory ${CMAKE_CURRENT_BINARY_DIR}/build/modules ${CMAKE_CURRENT_BINARY_DIR}/dist/modules
		COMMAND ${CMAKE_COMMAND} -E copy_directory ${CMAKE_CURRENT_BINARY_DIR}/build/server ${CMAKE_CURRENT_BINARY_DIR}/dist/server
		COMMAND ${CMAKE_COMMAND} -E copy ${CMAKE_CURRENT_BINARY_DIR}/build/inferno.js ${CMAKE_CURRENT_BINARY_DIR}/dist
		COMMAND ${NPM_PATH} run compile
		COMMAND ${CMAKE_COMMAND} -E copy_directory ${CMAKE_CURRENT_BINARY_DIR}/build/Release ${CMAKE_CURRENT_BINARY_DIR}/dist/Release
		COMMAND ${NPM_PATH} run concat
		COMMAND ${CMAKE_COMMAND} -E remove_directory ${CMAKE_CURRENT_BINARY_DIR}/src
		COMMAND ${CMAKE_COMMAND} -E remove_directory ${CMAKE_CURRENT_BINARY_DIR}/build
		COMMAND ${CMAKE_COMMAND} -E copy ${CMAKE_SOURCE_DIR}/gui2.0/users.htpasswd ${CMAKE_CURRENT_BINARY_DIR}/dist/server
		#COMMAND ${CMAKE_COMMAND} -E copy ${CMAKE_SOURCE_DIR}/gui2.0/simple.sqlite ${CMAKE_CURRENT_BINARY_DIR}/dist
		COMMAND ${CMAKE_COMMAND} -E copy ${CMAKE_SOURCE_DIR}/gui2.0/src/modules/simple/simple_user_interface.json ${CMAKE_CURRENT_BINARY_DIR}/dist/modules/simple
		COMMAND ${CMAKE_COMMAND} -E copy_directory ${CMAKE_SOURCE_DIR}/gui2.0/src/bin ${CMAKE_BINARY_DIR}/bin
		COMMAND ${CMAKE_COMMAND} -E remove ${CMAKE_CURRENT_BINARY_DIR}/concat.sh
		COMMAND ${CMAKE_COMMAND} -E remove ${CMAKE_CURRENT_BINARY_DIR}/tsconfig.json
		COMMAND ${CMAKE_COMMAND} -E remove ${CMAKE_CURRENT_BINARY_DIR}/package.json
		COMMAND ${CMAKE_COMMAND} -E remove ${CMAKE_CURRENT_BINARY_DIR}/binding.gyp
		COMMAND ${CMAKE_COMMAND} -E remove ${CMAKE_CURRENT_BINARY_DIR}/package-lock.json
	)
	
endif()


