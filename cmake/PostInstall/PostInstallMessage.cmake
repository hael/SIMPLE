##
# Print some post install messages for the user
##
message("\n\n")
message("Installation complete.")
message("==========================================================================")
message("Please ensure the following variables are set properly in add2.*rc file:  ")
message("    SIMPLE_EMAIL  SIMPLE_QSYS  SIMPLE_PATH  SIMPLE_SOURCE_PATH            ")
message("To use SIMPLE, append the relevant add2.* to your HOME shell rc file:     ")
message(" bash$ cat add2.bashrc >> ~/.bashrc                                       ")
message(" tcsh$ cat add2.tcshrc >> ~/.tcshrc                                       ")

if (${CMAKE_SYSTEM_NAME} MATCHES "Darwin")
	message("Fixing library names with install_name_tool")
#	execute_process(COMMAND install_name_tool -add_rpath ${CMAKE_INSTALL_PREFIX}/lib ${CMAKE_INSTALL_PREFIX}/bin/simple_exec)
#	execute_process(COMMAND install_name_tool -add_rpath ${CMAKE_INSTALL_PREFIX}/lib ${CMAKE_INSTALL_PREFIX}/bin/simple_distr_exec)
endif()
message("==========================================================================")
message("For minimal installation to work correctly add:\n${CMAKE_INSTALL_PREFIX}/bin and\n${CMAKE_INSTALL_PREFIX}/scripts\n to your PATH environment variable.")
if(BUILD_SHARED_LIBS)
message("For shared or dynamic builds: ")
if ("${CMAKE_SYSTEM_NAME}" MATCHES "Darwin")
	message("MacOS users please add ${CMAKE_INSTALL_PREFIX}/lib to one of: \n\tDYLD_LIBRARY_PATH or\n\tDYLD_FALLBACK_LIBRARY_PATH or\n\tRPATH ")
else()
	message("Please add ${CMAKE_INSTALL_PREFIX}/lib to LD_LIBRARY_PATH")
endif()
endif()
message("==========================================================================")
message("\n\n")


