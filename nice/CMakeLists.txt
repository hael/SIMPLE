cmake_minimum_required(VERSION 3.12)

project(NICE
        VERSION ${SIMPLE_VERSION}
        )

include(ExternalProject)

if(${NICE})

	ExternalProject_Add(
            nice
            URL ${CMAKE_CURRENT_SOURCE_DIR}/nice-2.0.3.tar.gz
            CMAKE_ARGS -DINSTALL_PREFIX=${CMAKE_BINARY_DIR}       
    	)
    
endif()

