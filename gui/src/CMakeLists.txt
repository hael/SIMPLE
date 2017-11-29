add_custom_command(OUTPUT ${CMAKE_BINARY_DIR}/bin/simple
		   COMMAND g++ simple_gui.cpp ../ext/mongoose/mongoose.c ${CMAKE_CURRENT_BINARY_DIR}/../ext/jpeg-9b/lib/libjpeg.a -O3 -pthread -lsqlite3 -o ${CMAKE_BINARY_DIR}/bin/simple -I ${CMAKE_CURRENT_BINARY_DIR}/../ext/jpeg-9b/include/ -std=c++11 > /tmp/gui_build.log
		   WORKING_DIRECTORY "${CMAKE_CURRENT_SOURCE_DIR}"
		   DEPENDS simple_libjpeg
		   COMMENT "Compiling GUI ...")
add_custom_target(simple_GUI ALL DEPENDS ${CMAKE_BINARY_DIR}/bin/simple)
