#add_custom_command(OUTPUT ${CMAKE_BINARY_DIR}/bin/simple
#		   COMMAND g++ simple_gui.cpp ../ext/mongoose/mongoose.c ${CMAKE_CURRENT_BINARY_DIR}/../ext/jpeg-9b/lib/libjpeg.a -O3 -pthread -lsqlite3 -o ${CMAKE_BINARY_DIR}/bin/simple -I ${CMAKE_CURRENT_BINARY_DIR}/../ext/jpeg-9b/include/ -std=c++11 > /tmp/gui_build.log
#		   WORKING_DIRECTORY "${CMAKE_CURRENT_SOURCE_DIR}"
#		   DEPENDS ${CMAKE_CURRENT_BINARY_DIR}/../ext/jpeg-9b/lib/libjpeg.a
#		   COMMENT "Compiling GUI ...")
#add_custom_target(simple_GUI ALL DEPENDS ${CMAKE_BINARY_DIR}/bin/simple)
#set(CMAKE_CXX_COMPILER "g++-4.9")
#set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -O3 -pthread")
#
#set(CMAKE_CXX_STANDARD 11)
#enable_language( CXX C)
#add_executable (simple simple_gui.cpp)
#target_include_directories (simple PUBLIC ${CMAKE_BINARY_DIR}/ext/include)

set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -O3 -pthread -std=c++11")
set(LIBRARY_OUTPUT_PATH "${CMAKE_BINARY_DIR}/bin")
add_executable (simple simple_gui.cpp)
add_dependencies (simple libjpeg libsqlite3 web_root mongoose)
target_include_directories (simple PUBLIC ${CMAKE_BINARY_DIR}/ext/include)
target_link_libraries (simple ${CMAKE_BINARY_DIR}/ext/lib/libsqlite3.a ${CMAKE_BINARY_DIR}/ext/lib/libjpeg.a ${CMAKE_BINARY_DIR}/ext/lib/libmongoose.a dl)
install (TARGETS simple DESTINATION ${CMAKE_BINARY_DIR}/bin)

