project(simple_gui)

cmake_minimum_required(VERSION 3.2)

if(NOT ${BUILD_GUI} AND ${GUI})
	if(${APPLE})
		if(EXISTS ${CMAKE_CURRENT_SOURCE_DIR}/dist/simple_gui.mac)
			message("-- Located pre-built gui binary for mac")
			add_custom_target(
				gui_apple ALL
				COMMENT "Installing precompiled GUI ... "
				COMMAND ${CMAKE_COMMAND} -E copy ${CMAKE_CURRENT_SOURCE_DIR}/dist/simple_gui.mac ${CMAKE_BINARY_DIR}/bin/simple
			)
		else()
			message(FATAL_ERROR "Could not locate a prebuilt gui binary for your architecture. Please rerun cmake with either -DGUI=OFF or -DBUILD_GUI=ON")
		endif()
	elseif(${UNIX})
		if(EXISTS ${CMAKE_CURRENT_SOURCE_DIR}/dist/simple_gui.linux)
            message("-- Located pre-built gui binary for linux")
			add_custom_target(
                gui_unix ALL
                COMMENT "Installing precompiled GUI ... "
                COMMAND ${CMAKE_COMMAND} -E copy ${CMAKE_CURRENT_SOURCE_DIR}/dist/simple_gui.linux ${CMAKE_BINARY_DIR}/bin/simple
            )
		else()
            message(FATAL_ERROR "Could not locate a prebuilt gui binary for your architecture. Please rerun cmake with either -DGUI=OFF or -DBUILD_GUI=ON")
        endif()
	else()
		message(FATAL_ERROR "Could not locate a prebuilt gui binary for your architecture. Please rerun cmake with either -DGUI=OFF or -DBUILD_GUI=ON")
    endif()
endif()

if(${MULTIUSER_GUI} AND ${GUI})
	find_program(NODE_PATH node)
	find_program(NPM_PATH npm)

	if(${NODE_PATH} STREQUAL "NODE_PATH-NOTFOUND" )
		message(FATAL_ERROR "Could not locate an executable for node. Cannot compile GUI. Please install node if you wish to compile the GUI")
	else()
		message("-- Located node executable at ${NODE_PATH}")
	endif()

	if(${NPM_PATH} STREQUAL "NPM_PATH-NOTFOUND" )
		message(WARNING "Could not locate an executable for npm. Please install npm if you wish to compile the GUI")
	else()
        message("-- Located npm executable at ${NPM_PATH}")
    endif()

    add_custom_target(
		gui_compile ALL
		COMMENT "Compiling GUI ... "
		COMMAND ${CMAKE_COMMAND} -E copy_directory ${CMAKE_CURRENT_SOURCE_DIR}/core ${CMAKE_BINARY_DIR}/bin/gui/core
        COMMAND ${CMAKE_COMMAND} -E copy_directory ${CMAKE_CURRENT_SOURCE_DIR}/external ${CMAKE_BINARY_DIR}/bin/gui/external
        COMMAND ${CMAKE_COMMAND} -E copy_directory ${CMAKE_CURRENT_SOURCE_DIR}/icons ${CMAKE_BINARY_DIR}/bin/gui/icons
        COMMAND ${CMAKE_COMMAND} -E copy_directory ${CMAKE_CURRENT_SOURCE_DIR}/public ${CMAKE_BINARY_DIR}/bin/gui/public
        COMMAND ${CMAKE_COMMAND} -E copy ${CMAKE_CURRENT_SOURCE_DIR}/package.json ${CMAKE_BINARY_DIR}/bin/gui/package.json
        COMMAND ${NPM_PATH} install > ${CMAKE_BINARY_DIR}/bin/gui/install.log
        COMMAND ${NPM_PATH} rebuild > ${CMAKE_BINARY_DIR}/bin/gui/install.log
    )
endif()


if(${BUILD_GUI} AND ${GUI})
	find_program(NODE_PATH node)
	find_program(NPM_PATH npm)

	if(${NODE_PATH} STREQUAL "NODE_PATH-NOTFOUND" )
		message(FATAL_ERROR "Could not locate an executable for node. Cannot compile GUI. Please install node if you wish to compile the GUI")
	else()
		message("-- Located node executable at ${NODE_PATH}")
	endif()

	if(${NPM_PATH} STREQUAL "NPM_PATH-NOTFOUND" )
		message(WARNING "Could not locate an executable for npm. Please install npm if you wish to compile the GUI")
	else()
        message("-- Located npm executable at ${NPM_PATH}")
    endif()
	if(${APPLE})
		add_custom_target(
			gui_compile ALL
			COMMENT "Compiling GUI ... "
			COMMAND ${CMAKE_COMMAND} -E copy_directory ${CMAKE_CURRENT_SOURCE_DIR}/core ${CMAKE_BINARY_DIR}/bin/gui/core
			COMMAND ${CMAKE_COMMAND} -E copy_directory ${CMAKE_CURRENT_SOURCE_DIR}/external ${CMAKE_BINARY_DIR}/bin/gui/external
			COMMAND ${CMAKE_COMMAND} -E copy_directory ${CMAKE_CURRENT_SOURCE_DIR}/icons ${CMAKE_BINARY_DIR}/bin/gui/icons
			COMMAND ${CMAKE_COMMAND} -E copy_directory ${CMAKE_CURRENT_SOURCE_DIR}/public ${CMAKE_BINARY_DIR}/bin/gui/public
			COMMAND ${CMAKE_COMMAND} -E copy ${CMAKE_CURRENT_SOURCE_DIR}/package.json ${CMAKE_BINARY_DIR}/bin/gui/package.json
			COMMAND ${NPM_PATH} install > ${CMAKE_BINARY_DIR}/bin/gui/install.log
			COMMAND ${NPM_PATH} rebuild > ${CMAKE_BINARY_DIR}/bin/gui/install.log
			COMMAND ${NPM_PATH} run build > ${CMAKE_BINARY_DIR}/bin/gui/install.log
			
			COMMAND ${CMAKE_COMMAND} -E create_symlink ${CMAKE_CURRENT_SOURCE_DIR}/dist/simple_gui.mac ${CMAKE_BINARY_DIR}/bin/simple
		)
	elseif(${UNIX})
		add_custom_target(
			gui_compile ALL
			COMMENT "Compiling GUI ... "
			COMMAND ${CMAKE_COMMAND} -E copy_directory ${CMAKE_CURRENT_SOURCE_DIR}/core ${CMAKE_BINARY_DIR}/bin/gui/core
			COMMAND ${CMAKE_COMMAND} -E copy_directory ${CMAKE_CURRENT_SOURCE_DIR}/external ${CMAKE_BINARY_DIR}/bin/gui/external
			COMMAND ${CMAKE_COMMAND} -E copy_directory ${CMAKE_CURRENT_SOURCE_DIR}/icons ${CMAKE_BINARY_DIR}/bin/gui/icons
			COMMAND ${CMAKE_COMMAND} -E copy_directory ${CMAKE_CURRENT_SOURCE_DIR}/public ${CMAKE_BINARY_DIR}/bin/gui/public
			COMMAND ${CMAKE_COMMAND} -E copy ${CMAKE_CURRENT_SOURCE_DIR}/package.json ${CMAKE_BINARY_DIR}/bin/gui/package.json
			COMMAND ${NPM_PATH} install > ${CMAKE_BINARY_DIR}/bin/gui/install.log
			COMMAND ${NPM_PATH} rebuild > ${CMAKE_BINARY_DIR}/bin/gui/install.log
			COMMAND ${NPM_PATH} run build > ${CMAKE_BINARY_DIR}/bin/gui/install.log
			COMMAND ${CMAKE_COMMAND} -E copy ${CMAKE_BINARY_DIR}/bin/gui/dist/simple_gui.linux ${CMAKE_BINARY_DIR}/bin/simple
		)
	endif()
endif()

