##
# Print some post install messages for the user
##
message("\n\n")
message("Installation complete.")
message("==========================================================================")
message("Please ensure the following variables are set properly in add2.*rc file:  ")
message("    SIMPLE_EMAIL  SIMPLE_QSYS  SIMPLE_PATH                                ")
message("To use SIMPLE, append the relevant add2.* to your HOME shell rc file:     ")
message(" bash$ cat add2.bashrc >> ~/.bashrc                                       ")
message(" tcsh$ cat add2.tcshrc >> ~/.tcshrc                                       ")
message("==========================================================================")
message("For minimal installation to work correctly add:\n${CMAKE_INSTALL_PREFIX}/bin and\n${CMAKE_INSTALL_PREFIX}/scripts\n to your PATH environment variable.")
#if(BUILD_SHARED_LIBS)
# message("For shared or dynamic builds: ")
#message("SHARED LIB: ${CMAKE_INSTALL_PREFIX}/lib appended LD_LIBRARY_PATH")
#endif()
message("==========================================================================")
message("\n\n")
