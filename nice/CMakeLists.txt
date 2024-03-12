include(ExternalProject)

if(${NICE})

	ExternalProject_Add(
        nice
        URL ${CMAKE_CURRENT_SOURCE_DIR}/nice-2.0.1.tar.gz
        CMAKE_ARGS
			-DINSTALL_PREFIX=${CMAKE_BINARY_DIR}       
    )
    
endif()

