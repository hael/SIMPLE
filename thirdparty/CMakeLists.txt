# JSON-Fortran
include(ExternalProject)

set(JSON_LIB ${CMAKE_CURRENT_BINARY_DIR}/src/libjson-build/lib) 

ExternalProject_Add( libjson
  URL ${CMAKE_SOURCE_DIR}/thirdparty/json-fortran-6.1.0.tar.gz
  PREFIX ${CMAKE_CURRENT_BINARY_DIR}/../
  CMAKE_ARGS -DCMAKE_MODULE_PATH=${CMAKE_MODULE_PATH}
  CMAKE_ARGS -Wno-dev
  CMAKE_ARGS -DCMAKE_INSTALL_PREFIX=${CMAKE_CURRENT_BINARY_DIR}/..
  CMAKE_ARGS -DUSE_GNU_INSTALL_CONVENTION:BOOL=TRUE
  CMAKE_ARGS -DSKIP_DOC_GEN=on
  CMAKE_ARGS -DCMAKE_Fortran_FLAGS=${CMAKE_Fortran_FLAGS}
  CMAKE_ARGS -DCMAKE_Fortran_FLAGS_DEBUG=${CMAKE_Fortran_FLAGS_DEBUG}
  CMAKE_ARGS -DCMAKE_Fortran_FLAGS_RELEASE=${CMAKE_Fortran_FLAGS_RELEASE}
  CMAKE_ARGS -DCMAKE_INSTALL_LIBDIR=${CMAKE_BINARY_DIR}/lib/
  CMAKE_ARGS -DCMAKE_INSTALL_INCLUDEDIR=${CMAKE_BINARY_DIR}/lib/
  BUILD_COMMAND make -j 
 # BUILD_COMMAND cp ${JSON_LIB}/libjsonfortran.* ${CMAKE_BINARY_DIR}/lib/ &&
#  BUILD_COMMAND cp ${JSON_LIB}/../include/*.mod ${CMAKE_BINARY_DIR}/lib/
  INSTALL_DIR ${CMAKE_CURRENT_BINARY_DIR}/../
  #INSTALL_COMMAND ""
  )	

#  find_package(jsonfortran-gnu REQUIRED PATH ${CMAKE_BINARY_DIR}/lib/cmake/jsonfortran-gnu-6.1.0)
#if (jsonfortran_INCLUDE_DIRS)

#include ( ${jsonfortran-gnu_CONFIG_INSTALL_DIR}/jsonfortran-gnu-config.cmake )


#find_library(JSONFORTRAN_LIBRARY 
##  NAME jsonfortran libjsonfortran
#  HINTS ../lib
#  PATH ${CMAKE_BINARY_DIR}/lib
#  )


#cp ${JSON_LIB}/libjsonfortran.a ${CMAKE_BINARY_DIR}/lib/
set(JSON_INC ${CMAKE_CURRENT_BINARY_DIR}/src/libjson-build/include) 
set(JSONFORTRAN_LIBRARY  jsonfortran)
set(EXTRA_LIBS ${EXTRA_LIBS} jsonfortran)
include_directories(${JSON_INC})
#target_link_libraries(${SIMPLELIB}  LINK_PUBLIC ${JSON_LIBRARIES})
mark_as_advanced (JSONFORTRAN_LIBRARY EXTRA_LIBS)
message(STATUS " JSON_Fortran completed  ${JSONFORTRAN_LIBRARY}")
message(STATUS " JSON_Fortran completed  ${EXTRA_LIBS}")