set(POST_INSTALL_SCRIPT ${CMAKE_SOURCE_DIR}/cmake/PostInstall/PostInstallMessage.cmake)
install(
    CODE 
    "
    execute_process(COMMAND \"${CMAKE_COMMAND}\"
                            -DCMAKE_SYSTEM_NAME=${CMAKE_SYSTEM_NAME}
                            -DCMAKE_INSTALL_PREFIX=${CMAKE_INSTALL_PREFIX}
                            -DBUILD_SHARED_LIBS=${BUILD_SHARED_LIBS}
                            -P \"${POST_INSTALL_SCRIPT}\")
    "
)
