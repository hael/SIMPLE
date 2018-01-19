set(LIBRARY_OUTPUT_PATH "${CMAKE_BINARY_DIR}/ext/lib")
add_library (mongoose mongoose.c)
file(COPY mongoose.h DESTINATION ${CMAKE_BINARY_DIR}/ext/include)
