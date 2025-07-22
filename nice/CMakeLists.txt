include(ExternalProject)

if(${NICE})

	ExternalProject_Add(
        nice
        URL ${CMAKE_CURRENT_SOURCE_DIR}/nice_pkg_7619eed.tar.gz
        CMAKE_ARGS -DINSTALL_PREFIX=${CMAKE_BINARY_DIR}/nice
        INSTALL_COMMAND ${CMAKE_COMMAND} -E copy_directory ${CMAKE_BINARY_DIR}/nice/nice-prefix/src/nice-build/lib ${CMAKE_BINARY_DIR}/nice/lib && ${CMAKE_COMMAND} -E copy_directory ${CMAKE_BINARY_DIR}/nice/nice-prefix/src/nice-build/var ${CMAKE_BINARY_DIR}/nice/var && ${CMAKE_COMMAND} -E copy_directory ${CMAKE_BINARY_DIR}/nice/nice-prefix/src/nice-build/etc ${CMAKE_BINARY_DIR}/nice/etc && ${CMAKE_COMMAND} -E copy_directory ${CMAKE_BINARY_DIR}/nice/nice-prefix/src/nice-build/bin ${CMAKE_BINARY_DIR}/nice/bin && ${CMAKE_COMMAND} -E copy_directory ${CMAKE_BINARY_DIR}/nice/nice-prefix/src/nice-build/share ${CMAKE_BINARY_DIR}/nice/share && make install 
    )

    write_file(${CMAKE_BINARY_DIR}/bin/nice "#!/usr/bin/env bash\nroot=`dirname -- \"$0\"`\n$root/../nice/bin/nice local")
    file(CHMOD ${CMAKE_BINARY_DIR}/bin/nice PERMISSIONS OWNER_READ OWNER_WRITE OWNER_EXECUTE GROUP_READ GROUP_EXECUTE WORLD_READ WORLD_EXECUTE)
    
endif()

