project(simple_gui)

cmake_minimum_required(VERSION 3.2)

if(${GUI})
	if(${APPLE})
		if(EXISTS ${CMAKE_CURRENT_SOURCE_DIR}/precompiled/simple_gui.mac.tar.gz)
			message("-- Located pre-built gui binary for mac")
			add_custom_target(
				gui_apple ALL
				COMMENT "Installing precompiled GUI ... "
				COMMAND ${CMAKE_COMMAND} -E tar xf ${CMAKE_CURRENT_SOURCE_DIR}/precompiled/simple_gui.mac.tar.gz
				COMMAND ${CMAKE_COMMAND} -E copy_directory ${CMAKE_BINARY_DIR}/bin/gui/mac ${CMAKE_BINARY_DIR}/bin/mac
				COMMAND ${CMAKE_COMMAND} -E create_symlink ${CMAKE_BINARY_DIR}/bin/mac/simple_gui.app/Contents/MacOS/simple_gui ${CMAKE_BINARY_DIR}/bin/simple
				COMMAND ${CMAKE_COMMAND} -E copy ${CMAKE_CURRENT_SOURCE_DIR}/precompiled/simple_gui_multiuser.mac ${CMAKE_BINARY_DIR}/bin/simple_multiuser
				COMMAND ${CMAKE_COMMAND} -E copy_directory ${CMAKE_CURRENT_SOURCE_DIR}/client ${CMAKE_BINARY_DIR}/gui_data/client
				COMMAND ${CMAKE_COMMAND} -E copy_directory ${CMAKE_CURRENT_SOURCE_DIR}/server ${CMAKE_BINARY_DIR}/gui_data/server
				COMMAND ${CMAKE_COMMAND} -E copy_directory ${CMAKE_CURRENT_SOURCE_DIR}/client/public ${CMAKE_BINARY_DIR}/gui_data/public
				COMMAND ${CMAKE_COMMAND} -E tar xf ${CMAKE_CURRENT_SOURCE_DIR}/precompiled/node_modules.mac.tar.gz
                		COMMAND ${CMAKE_COMMAND} -E copy_directory ${CMAKE_BINARY_DIR}/bin/gui/node_modules ${CMAKE_BINARY_DIR}/bin/node_modules
			)
		else()
			message(FATAL_ERROR "Could not locate a prebuilt gui binary for your architecture. Please rerun cmake with -DGUI=OFF")
		endif()
	elseif(${UNIX})
		if(EXISTS ${CMAKE_CURRENT_SOURCE_DIR}/precompiled/simple_gui.linux)
            message("-- Located pre-built gui binary for linux")
			add_custom_target(
                gui_unix ALL
                COMMENT "Installing precompiled GUI ... "
                COMMAND ${CMAKE_COMMAND} -E copy ${CMAKE_CURRENT_SOURCE_DIR}/precompiled/simple_gui.linux ${CMAKE_BINARY_DIR}/bin/simple
                COMMAND ${CMAKE_COMMAND} -E copy ${CMAKE_CURRENT_SOURCE_DIR}/precompiled/simple_gui_multiuser.linux ${CMAKE_BINARY_DIR}/bin/simple_multiuser
		COMMAND ${CMAKE_COMMAND} -E copy_directory ${CMAKE_CURRENT_SOURCE_DIR}/client ${CMAKE_BINARY_DIR}/gui_data/client
		COMMAND ${CMAKE_COMMAND} -E copy_directory ${CMAKE_CURRENT_SOURCE_DIR}/server ${CMAKE_BINARY_DIR}/gui_data/server
		COMMAND ${CMAKE_COMMAND} -E copy_directory ${CMAKE_CURRENT_SOURCE_DIR}/client/public ${CMAKE_BINARY_DIR}/gui_data/public
		COMMAND ${CMAKE_COMMAND} -E tar xf ${CMAKE_CURRENT_SOURCE_DIR}/precompiled/node_modules.linux.tar.gz
		COMMAND ${CMAKE_COMMAND} -E copy_directory ${CMAKE_BINARY_DIR}/bin/gui/node_modules ${CMAKE_BINARY_DIR}/bin/node_modules
            )
		else()
            message(FATAL_ERROR "Could not locate a prebuilt gui binary for your architecture. Please rerun cmake with -DGUI=OFF")
        endif()
	else()
		message(FATAL_ERROR "Could not locate a prebuilt gui binary for your architecture. Please rerun cmake with -DGUI=OFF")
    endif()
    
endif()

