# The informational executable requires GLib 2.64; EL8's interface library does not.
set(source "${SOURCE_DIR}/src/CMakeLists.txt")
file(READ "${source}" contents)
set(command_block [=[configure_file(gtk-fortran.f90 gtk-fortran.out.f90)
add_executable(${gtk_V_fortran} "gtk-fortran.out.f90")
target_link_libraries(${gtk_V_fortran} gtk-fortran_static ${GTK_LIBRARIES})
set_target_properties(${gtk_V_fortran} PROPERTIES
  Fortran_MODULE_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/modules/ )
add_test(${gtk_V_fortran} ./${gtk_V_fortran})]=])
set(install_block [=[install(TARGETS ${gtk_V_fortran}
  DESTINATION ${CMAKE_INSTALL_BINDIR}
  PERMISSIONS OWNER_READ OWNER_WRITE OWNER_EXECUTE GROUP_READ GROUP_EXECUTE WORLD_READ WORLD_EXECUTE
)]=])
foreach(block IN ITEMS command_block install_block)
    string(FIND "${contents}" "${${block}}" position)
    if(position EQUAL -1)
        message(FATAL_ERROR "Pinned gtk-fortran source changed: cannot find ${block} in ${source}")
    endif()
    string(REPLACE "${${block}}" "" contents "${contents}")
endforeach()
file(WRITE "${source}" "${contents}")
