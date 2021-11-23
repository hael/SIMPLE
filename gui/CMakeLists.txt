project(simple_gui)

cmake_minimum_required(VERSION 3.2)

include(ExternalProject)

if(${GUI})

	add_custom_target(
                        gui_core ALL
                        COMMENT "Installing GUI core ... "
                        COMMAND ${CMAKE_COMMAND} -E copy_directory ${CMAKE_CURRENT_SOURCE_DIR}/client ${CMAKE_BINARY_DIR}/gui_data/client
                        COMMAND ${CMAKE_COMMAND} -E copy_directory ${CMAKE_CURRENT_SOURCE_DIR}/server ${CMAKE_BINARY_DIR}/gui_data/server
                        COMMAND ${CMAKE_COMMAND} -E copy_directory ${CMAKE_CURRENT_SOURCE_DIR}/client/public ${CMAKE_BINARY_DIR}/gui_data/public
                        COMMAND ${CMAKE_COMMAND} -E tar xf ${CMAKE_CURRENT_SOURCE_DIR}/tutorials.tar.gz
                        COMMAND ${CMAKE_COMMAND} -E copy_directory ${CMAKE_BINARY_DIR}/bin/gui/tutorials ${CMAKE_BINARY_DIR}/tutorials
                        COMMAND ${CMAKE_COMMAND} -E remove_directory ${CMAKE_BINARY_DIR}/bin/gui/tutorials
        )
	
	install(DIRECTORY ${CMAKE_BINARY_DIR}/gui_data DESTINATION ${CMAKE_INSTALL_PREFIX})
	install(DIRECTORY ${CMAKE_BINARY_DIR}/tutorials DESTINATION ${CMAKE_INSTALL_PREFIX})

	if(${APPLE})

		ExternalProject_Add(
                        nodejs_apple
                        URL https://nodejs.org/dist/v16.13.0/node-v16.13.0-darwin-x64.tar.gz
                        URL_MD5 fa4e4210611a00d34a50dd2a6ce0b738
                        DOWNLOAD_DIR ${CMAKE_BINARY_DIR}/external/tarballs
                        SOURCE_DIR ${CMAKE_BINARY_DIR}/bin/nodejs
                        CONFIGURE_COMMAND ""
                        BUILD_COMMAND ""
                        INSTALL_COMMAND ""
                        TEST_COMMAND ""
                )

		install(DIRECTORY ${CMAKE_BINARY_DIR}/bin/nodejs DESTINATION ${CMAKE_INSTALL_PREFIX}/bin USE_SOURCE_PERMISSIONS)

		add_custom_target(
                        gui_multi_apple ALL
                        COMMENT "Installing GUI multiuser ... "
                        DEPENDS gui_core nodejs_apple
                        COMMAND ${CMAKE_COMMAND} -E env "PATH=${CMAKE_BINARY_DIR}/bin/nodejs/bin:$ENV{PATH}" npm install
			COMMAND ${CMAKE_COMMAND} -E copy ${CMAKE_CURRENT_SOURCE_DIR}/server/multi_user/simple_multiuser.sh ${CMAKE_BINARY_DIR}/bin/simple_multiuser
                        WORKING_DIRECTORY ${CMAKE_BINARY_DIR}/gui_data/server/multi_user
                )

	install(PROGRAMS ${CMAKE_BINARY_DIR}/bin/simple_multiuser DESTINATION ${CMAKE_INSTALL_PREFIX}/bin)

		if(EXISTS ${CMAKE_CURRENT_SOURCE_DIR}/precompiled/simple_gui.mac.tar.gz)
			message("-- Located pre-built gui binary for mac")
			add_custom_target(
				gui_single_apple ALL
				DEPENDS gui_core
				COMMENT "Installing precompiled single-user GUI ... "
				COMMAND ${CMAKE_COMMAND} -E tar xf ${CMAKE_CURRENT_SOURCE_DIR}/precompiled/simple_gui.mac.tar.gz
				COMMAND ${CMAKE_COMMAND} -E copy ${CMAKE_CURRENT_SOURCE_DIR}/precompiled/simple_mac.sh ${CMAKE_BINARY_DIR}/bin/simple
			)

		install(PROGRAMS ${CMAKE_BINARY_DIR}/bin/simple DESTINATION ${CMAKE_INSTALL_PREFIX}/bin)
		
		else()
			message(FATAL_ERROR "Could not locate a prebuilt gui binary for your architecture. Please rerun cmake with -DGUI=OFF")
		endif()
	elseif(${UNIX})

		ExternalProject_Add(
			nodejs_unix
    			URL https://nodejs.org/dist/v16.13.0/node-v16.13.0-linux-x64.tar.gz
			URL_MD5 bb20c3f845896272f18bd8c91d1c9b07
			DOWNLOAD_DIR ${CMAKE_BINARY_DIR}/external/tarballs
			SOURCE_DIR ${CMAKE_BINARY_DIR}/bin/nodejs
			CONFIGURE_COMMAND ""
			BUILD_COMMAND ""
			INSTALL_COMMAND ""
    			TEST_COMMAND ""
		)

		install(DIRECTORY ${CMAKE_BINARY_DIR}/bin/nodejs DESTINATION ${CMAKE_INSTALL_PREFIX}/bin USE_SOURCE_PERMISSIONS)

		add_custom_target(
                        gui_multi_unix ALL
			COMMENT "Installing GUI multiuser ... "
			DEPENDS gui_core nodejs_unix
			COMMAND ${CMAKE_COMMAND} -E env "PATH=${CMAKE_BINARY_DIR}/bin/nodejs/bin:$ENV{PATH}" npm install
			COMMAND ${CMAKE_COMMAND} -E copy ${CMAKE_CURRENT_SOURCE_DIR}/server/multi_user/simple_multiuser.sh ${CMAKE_BINARY_DIR}/bin/simple_multiuser
			WORKING_DIRECTORY ${CMAKE_BINARY_DIR}/gui_data/server/multi_user
		)

		install(PROGRAMS ${CMAKE_BINARY_DIR}/bin/simple_multiuser DESTINATION ${CMAKE_INSTALL_PREFIX}/bin)
		
		if(${BUILD_GUI})
			add_custom_target(
                 	       	gui_single_unix ALL
                	       	COMMENT "Installing GUI singleuser ... "
                        	DEPENDS gui_core nodejs_unix
                        	COMMAND ${CMAKE_COMMAND} -E env "PATH=${CMAKE_BINARY_DIR}/bin/nodejs/bin:$ENV{PATH}" npm install
				COMMAND ${CMAKE_COMMAND} -E env "PATH=${CMAKE_BINARY_DIR}/bin/nodejs/bin:$ENV{PATH}" npm run build
				COMMAND ${CMAKE_COMMAND} -E copy ${CMAKE_BINARY_DIR}/gui_data/server/single_user/dist/simple_gui.linux ${CMAKE_BINARY_DIR}/bin/simple
                        	WORKING_DIRECTORY ${CMAKE_BINARY_DIR}/gui_data/server/single_user
                	)
			
			install(PROGRAMS ${CMAKE_BINARY_DIR}/bin/simple DESTINATION ${CMAKE_INSTALL_PREFIX}/bin)

		else()
			if(EXISTS ${CMAKE_CURRENT_SOURCE_DIR}/precompiled/simple_gui.linux)
                		message("-- Located pre-built single-user gui binary for linux")
		
				add_custom_target(
                			gui_single_unix ALL
					DEPENDS gui_core
                			COMMENT "Installing precompiled single-user GUI ... "
                			COMMAND ${CMAKE_COMMAND} -E copy ${CMAKE_CURRENT_SOURCE_DIR}/precompiled/simple_gui.linux ${CMAKE_BINARY_DIR}/bin/simple
					COMMAND ${CMAKE_COMMAND} -E tar xf ${CMAKE_CURRENT_SOURCE_DIR}/precompiled/node_modules.linux.tar.gz
					COMMAND ${CMAKE_COMMAND} -E copy_directory ${CMAKE_BINARY_DIR}/bin/gui/node_modules ${CMAKE_BINARY_DIR}/bin/node_modules
            			)
				
				install(PROGRAMS ${CMAKE_BINARY_DIR}/bin/simple DESTINATION ${CMAKE_INSTALL_PREFIX}/bin)
				install(DIRECTORY ${CMAKE_BINARY_DIR}/bin/gui/node_modules DESTINATION ${CMAKE_INSTALL_PREFIX}/bin)
				
			else()
            			message(FATAL_ERROR "Could not locate a prebuilt single-user gui binary for your architecture. Please rerun cmake with -DGUI=OFF")
        		endif()
		endif()
	else()
		message(FATAL_ERROR "Could not locate a prebuilt single-user gui binary for your architecture. Please rerun cmake with -DGUI=OFF")
    endif()

endif()

