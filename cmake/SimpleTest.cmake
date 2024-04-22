# CTest - Simple Test Sete
enable_testing()
# emulate GNU Autotools `make check`
add_custom_target(check COMMAND ${CMAKE_CTEST_COMMAND})

add_test(
    NAME somethingworks
    COMMAND /path/to/tool someArg more Arg
)

