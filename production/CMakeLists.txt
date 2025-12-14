# Build production executables

# Define executable names
set(SIMPLE_EXE simple_exec)
set(SINGLE_EXE single_exec)
set(STREAM_EXE simple_stream)
set(SIMPLE_PRIVATE_EXE simple_private_exec)

# Create executables
add_executable(${SIMPLE_EXE} 
    "${CMAKE_CURRENT_SOURCE_DIR}/simple_exec/simple_exec.f90")
add_executable(${SINGLE_EXE} 
    "${CMAKE_CURRENT_SOURCE_DIR}/single_exec/single_exec.f90")
add_executable(${STREAM_EXE} 
    "${CMAKE_CURRENT_SOURCE_DIR}/simple_stream/simple_stream.f90")
add_executable(${SIMPLE_PRIVATE_EXE}
    "${CMAKE_CURRENT_SOURCE_DIR}/simple_private_exec/simple_private_exec.f90")

# List of all executables for easy iteration
set(SIMPLE_EXECUTABLES 
    ${SIMPLE_EXE} 
    ${SINGLE_EXE} 
    ${STREAM_EXE} 
    ${SIMPLE_PRIVATE_EXE}
)

# Link all executables to the SIMPLE library
foreach(exe ${SIMPLE_EXECUTABLES})
    target_link_libraries(${exe} PRIVATE ${SIMPLELIB} ${SIMPLE_LIBRARIES})
endforeach()

# Install executables
install(TARGETS ${SIMPLE_EXECUTABLES}
    RUNTIME DESTINATION ${CMAKE_INSTALL_BINDIR}
)

# ============================================================================
# Tests
# ============================================================================
if(BUILD_TESTS)
    # Find all test directories
    file(GLOB TEST_DIRS
        LIST_DIRECTORIES true
        RELATIVE ${CMAKE_CURRENT_SOURCE_DIR}
        "tests/test_*"
    )
    
    # Remove MPI tests if MPI not enabled
    if(NOT USE_MPI)
        list(REMOVE_ITEM TEST_DIRS "tests/test_mpi")
    endif()
    
    foreach(test_dir ${TEST_DIRS})
        # Extract test name
        string(REPLACE "tests/" "" test_name "${test_dir}")
        set(test_exe "simple_${test_name}")
        
        # Find Fortran sources for this test
        file(GLOB_RECURSE test_sources
            "${CMAKE_CURRENT_SOURCE_DIR}/${test_dir}/*.[fF][90][0538]"
        )
        
        # Only build if sources exist
        if(test_sources)
            add_executable(${test_exe} ${test_sources})
            target_link_libraries(${test_exe} PRIVATE ${SIMPLELIB} ${SIMPLE_LIBRARIES})
            
            # Install test executable
            install(TARGETS ${test_exe}
                RUNTIME DESTINATION ${CMAKE_INSTALL_BINDIR}
            )
            
            # Install test data files if they exist
            file(GLOB test_data "${CMAKE_CURRENT_SOURCE_DIR}/${test_dir}/*.txt")
            if(test_data)
                install(FILES ${test_data}
                    DESTINATION ${CMAKE_INSTALL_PREFIX}/tests/${test_name}
                )
            endif()
            
            # Add as CTest
            add_test(NAME ${test_name} COMMAND ${test_exe})
        endif()
    endforeach()
endif()