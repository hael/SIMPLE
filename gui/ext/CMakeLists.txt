add_custom_command(OUTPUT ${CMAKE_CURRENT_BINARY_DIR}/../ext/jpeg-9b/lib/libjpeg.a
		   COMMAND tar -xf ${CMAKE_CURRENT_SOURCE_DIR}/jpegsrc.v9b.tar.gz
		   COMMAND cd ${CMAKE_CURRENT_BINARY_DIR}/jpeg-9b/ && ./configure --prefix=${CMAKE_CURRENT_BINARY_DIR}/jpeg-9b > /tmp/libjpeg_config.log 
		   COMMAND cd ${CMAKE_CURRENT_BINARY_DIR}/jpeg-9b/ && make install > /tmp/libjpeg_build.log
		   WORKING_DIRECTORY "${CMAKE_CURRENT_BINARY_DIR}"
		   COMMENT "Configuring libjpeg ...")

add_custom_target(simple_libjpeg ALL DEPENDS ${CMAKE_CURRENT_BINARY_DIR}/../ext/jpeg-9b/lib/libjpeg.a)

