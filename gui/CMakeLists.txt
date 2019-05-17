project(simple_gui)

cmake_minimum_required(VERSION 3.2)
set(BUILD_GUI TRUE)
set(BUILD_GUI_MULTIUSER TRUE)
set(BUILD_GUI_APPLE FALSE)
set(BUILD_GUI_UNIX FALSE)

if(${BUILD_GUI})
	find_program(NODE_PATH node)
	find_program(NPM_PATH npm)
	
	if(${APPLE})
		if(EXISTS ${CMAKE_CURRENT_SOURCE_DIR}/dist/simple_gui.mac)
			message("-- Located pre-built gui binary for mac")
			set(BUILD_GUI_APPLE TRUE)
		else()
               		 message(WARNING "Could not locate a prebuilt gui binary for your architecture")
        	endif()
	elseif(${UNIX})
		if(EXISTS ${CMAKE_CURRENT_SOURCE_DIR}/dist/simple_gui.linux)
                	message("-- Located pre-built gui binary for linux")
			set(BUILD_GUI_UNIX TRUE)
		else()
                        message(WARNING "Could not locate a prebuilt gui binary for your architecture")
                endif()
	else()
		message(WARNING "Could not locate a prebuilt gui binary for your architecture")
        endif()

	if(${NODE_PATH} STREQUAL "NODE_PATH-NOTFOUND" )
		message(WARNING "Could not locate an executable for node. GUI will be installed without multiuser support")
		set(BUILD_GUI_MULTIUSER FALSE)
	else()
		message("-- Located node executable at ${NODE_PATH}")
	endif()

	if(${NPM_PATH} STREQUAL "NPM_PATH-NOTFOUND" )
		message(WARNING "Could not locate an executable for npm. GUI will be installed without multiuser support")
		set(BUILD_GUI_MULTIUSER FALSE)
	else()
                message("-- Located npm executable at ${NPM_PATH}")
        endif()

	if(${BUILD_GUI_APPLE})
		add_custom_target(
			gui_apple ALL
			COMMENT "Installing precompiled GUI ... "
			COMMAND ${CMAKE_COMMAND} -E copy ${CMAKE_CURRENT_SOURCE_DIR}/dist/simple_gui.mac ${CMAKE_BINARY_DIR}/bin/simple
		)
	endif()

	if(${BUILD_GUI_UNIX})
                add_custom_target(
                        gui_unix ALL
                        COMMENT "Installing precompiled GUI ... "
                        COMMAND ${CMAKE_COMMAND} -E copy ${CMAKE_CURRENT_SOURCE_DIR}/dist/simple_gui.linux ${CMAKE_BINARY_DIR}/bin/simple
                )
    endif()
    
    if(${BUILD_GUI_MULTIUSER})
                add_custom_target(
                        gui_multiuser ALL
                        COMMENT "Installing multiuser GUI ... "
                        COMMAND ${CMAKE_COMMAND} -E copy_directory ${CMAKE_CURRENT_SOURCE_DIR}/core ${CMAKE_BINARY_DIR}/bin/gui/core
                        COMMAND ${CMAKE_COMMAND} -E copy_directory ${CMAKE_CURRENT_SOURCE_DIR}/external ${CMAKE_BINARY_DIR}/bin/gui/external
                        COMMAND ${CMAKE_COMMAND} -E copy_directory ${CMAKE_CURRENT_SOURCE_DIR}/icons ${CMAKE_BINARY_DIR}/bin/gui/icons
                        COMMAND ${CMAKE_COMMAND} -E copy_directory ${CMAKE_CURRENT_SOURCE_DIR}/public ${CMAKE_BINARY_DIR}/bin/gui/public
                        COMMAND ${CMAKE_COMMAND} -E copy ${CMAKE_CURRENT_SOURCE_DIR}/package.json ${CMAKE_BINARY_DIR}/bin/gui/package.json
                        COMMAND ${NPM_PATH} install
                        COMMAND ${NPM_PATH} rebuild
                )
    endif()
		#endif()
		#if(${BUILD_GUI_APPLE})
		#        COMMAND ${CMAKE_COMMAND} -E copy ${CMAKE_CURRENT_SOURCE_DIR}/dist/simple_gui.mac ${CMAKE_BINARY_DIR}/bin/simple
		#endif()
		#COMMAND ${NPM_PATH} install
		#COMMAND ${NPM_PATH} run compile
		#COMMAND ${NPM_PATH} run concat
		#COMMAND ${CMAKE_COMMAND} -E remove_directory ${CMAKE_CURRENT_BINARY_DIR}/src
		#COMMAND ${CMAKE_COMMAND} -E remove_directory ${CMAKE_CURRENT_BINARY_DIR}/build
		#COMMAND ${CMAKE_COMMAND} -E copy ${CMAKE_SOURCE_DIR}/gui2.0/users.htpasswd ${CMAKE_CURRENT_BINARY_DIR}/dist/server
		#COMMAND ${CMAKE_COMMAND} -E copy ${CMAKE_SOURCE_DIR}/gui2.0/simple.sqlite ${CMAKE_CURRENT_BINARY_DIR}/dist
		#COMMAND ${CMAKE_COMMAND} -E copy ${CMAKE_SOURCE_DIR}/gui2.0/src/modules/simple/simple_user_interface.json ${CMAKE_CURRENT_BINARY_DIR}/dist/modules/simple
		#COMMAND ${CMAKE_COMMAND} -E copy_directory ${CMAKE_SOURCE_DIR}/gui2.0/src/bin ${CMAKE_BINARY_DIR}/bin
		#COMMAND ${CMAKE_COMMAND} -E remove ${CMAKE_CURRENT_BINARY_DIR}/concat.sh
		#COMMAND ${CMAKE_COMMAND} -E remove ${CMAKE_CURRENT_BINARY_DIR}/tsconfig.json
		#COMMAND ${CMAKE_COMMAND} -E remove ${CMAKE_CURRENT_BINARY_DIR}/package.json
		#COMMAND ${CMAKE_COMMAND} -E remove ${CMAKE_CURRENT_BINARY_DIR}/binding.gyp
		#COMMAND ${CMAKE_COMMAND} -E remove ${CMAKE_CURRENT_BINARY_DIR}/package-lock.json
	
endif()


