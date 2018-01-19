set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -O3 -pthread -std=c++11")
set(LIBRARY_OUTPUT_PATH "${CMAKE_BINARY_DIR}/bin")
add_executable (simple simple_gui.cpp)
add_dependencies (simple libjpeg libsqlite3 mongoose)
target_include_directories (simple PUBLIC ${CMAKE_BINARY_DIR}/ext/include)
target_link_libraries (simple ${CMAKE_BINARY_DIR}/ext/lib/libsqlite3.a ${CMAKE_BINARY_DIR}/ext/lib/libjpeg.a ${CMAKE_BINARY_DIR}/ext/lib/libmongoose.a dl)
install (TARGETS simple DESTINATION ${CMAKE_BINARY_DIR}/bin)

