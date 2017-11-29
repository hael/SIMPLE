add_custom_command(OUTPUT Makefile
		   COMMAND ./configure --prefix=${CMAKE_CURRENT_BINARY_DIR} > /tmp/libjpeg_config.log  
		   WORKING_DIRECTORY "${CMAKE_CURRENT_SOURCE_DIR}"
		   COMMENT "Configuring libjpeg ...")

add_custom_target(simple_libjpeg_configure ALL DEPENDS Makefile)

add_custom_command(OUTPUT ${CMAKE_CURRENT_BINARY_DIR}/../ext/jpeg-9b/lib/libjpeg.a
                   COMMAND make install > /tmp/libjpeg_build.log
                   WORKING_DIRECTORY "${CMAKE_CURRENT_SOURCE_DIR}"
                   COMMENT "Compiling libjpeg ...")

add_custom_target(simple_libjpeg ALL DEPENDS ${CMAKE_CURRENT_BINARY_DIR}/../ext/jpeg-9b/lib/libjpeg.a)
