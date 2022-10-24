# Source for the SIMPLE library.
add_custom_target(simple_args
  COMMAND SIMPLE_PATH=${CMAKE_BINARY_DIR} ${CMAKE_SOURCE_DIR}/scripts/simple_args_generator.pl
  DEPENDS "${CMAKE_SOURCE_DIR}/scripts/simple_args_generator.pl" "${CMAKE_CURRENT_SOURCE_DIR}/main/simple_parameters.f90"
  BYPRODUCTS "${CMAKE_BINARY_DIR}/lib/simple/simple_args.f90" "${CMAKE_BINARY_DIR}/lib/simple/simple_varlist.txt"
  WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}/main
  USES_TERMINAL
  )
add_custom_target(simple_git
  COMMAND SIMPLE_PATH=${CMAKE_BINARY_DIR} ${CMAKE_SOURCE_DIR}/scripts/insert_git_commit_hash.pl
  DEPENDS "${CMAKE_SOURCE_DIR}/scripts/insert_git_commit_hash.pl"
  WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}/main
  USES_TERMINAL
)
file(GLOB DEFS_src   "${CMAKE_CURRENT_SOURCE_DIR}/defs/*.f[90][08]"   PARENT_SCOPE)
file(GLOB FILEIO_src "${CMAKE_CURRENT_SOURCE_DIR}/fileio/*.f[90][08]" PARENT_SCOPE)
list(APPEND FILEIO_src
  "${CMAKE_CURRENT_SOURCE_DIR}/fileio/simple_posix.c"
  )
file(GLOB OPERATORS_src "${CMAKE_CURRENT_SOURCE_DIR}/ops/*.f90" PARENT_SCOPE)
list(APPEND OPERATORS_src
  "${CMAKE_CURRENT_SOURCE_DIR}/ops/simple_kbinterpol_memo.c"
  )
file(GLOB UTILS_src    "${CMAKE_CURRENT_SOURCE_DIR}/utils/*.f90" "${CMAKE_CURRENT_SOURCE_DIR}/utils/*/*.f90" PARENT_SCOPE)
file(GLOB STARFILE_src "${CMAKE_CURRENT_SOURCE_DIR}/starfile/*.cpp" PARENT_SCOPE)
file(GLOB CPLOT2D_src  "${CMAKE_CURRENT_SOURCE_DIR}/CPlot2D/*.cpp"  PARENT_SCOPE)
file(GLOB XML_src      "${CMAKE_CURRENT_SOURCE_DIR}/xml/*/*.[Ff]90" PARENT_SCOPE)

file(GLOB MAIN_src "${CMAKE_CURRENT_SOURCE_DIR}/main/*/*(^commander)*.f[90][08]" "${CMAKE_CURRENT_SOURCE_DIR}/main/*(^commander)*.f[90][08]" PARENT_SCOPE)
file(GLOB JPG_src "${CMAKE_CURRENT_SOURCE_DIR}/jpg/*.f[90][08]" PARENT_SCOPE)
list(APPEND JPG_src
  "${CMAKE_CURRENT_SOURCE_DIR}/jpg/stb_image_write.c")
if(USE_LIBTIFF)
    list(APPEND JPG_src "${CMAKE_CURRENT_SOURCE_DIR}/tiff/libtiff_funs.c")
endif()
file(GLOB SIMPLE_src "${CMAKE_CURRENT_SOURCE_DIR}/main/*/*.f[09][0-9]" "${CMAKE_CURRENT_SOURCE_DIR}/main/*.f[09][0-9]" PARENT_SCOPE)

# add the simple_args module to SIMPLE_src
list(APPEND SIMPLE_src
  "${CMAKE_BINARY_DIR}/lib/simple/simple_args.f90"
  "${CMAKE_BINARY_DIR}/lib/simple/simple_modules.f90"
  )
list(APPEND SIMPLECORELIB  ${DEFS_src} ${UTILS_src} ${FILEIO_src} ${STARFILE_src} ${CPLOT2D_src} ${XML_src} ${JPG_src} ${OPERATORS_src} ${SIMPLE_src})

add_custom_target(simple_modules
  COMMAND  echo Creating simple_module \; ( echo module simple_modules\; for file in ${DEFS_src} ${UTILS_src} ${FILEIO_src} ${JPG_src} ${OPERATORS_src} ${MAIN_src}\; do if \\[ \"\$\$\{file\#\#*.\}\" = \"f90\" \\] &&  test ! \"\$\${file}\" != \"\$\${file%simple_modules*}\" \; then      echo use \$\$\(basename \$\$file .f90\)\;fi\;done\;echo contains\;echo end module ) >> ${CMAKE_BINARY_DIR}/lib/simple/simple_modules_tmp.f90
  COMMAND if \\[ ! -f ${CMAKE_BINARY_DIR}/lib/simple/simple_modules.f90 \\]\;then mv ${CMAKE_BINARY_DIR}/lib/simple/simple_modules_tmp.f90 ${CMAKE_BINARY_DIR}/lib/simple/simple_modules.f90\;else cmp --silent ${CMAKE_BINARY_DIR}/lib/simple/simple_modules_tmp.f90 ${CMAKE_BINARY_DIR}/lib/simple/simple_modules.f90 || cp ${CMAKE_BINARY_DIR}/lib/simple/simple_modules_tmp.f90 ${CMAKE_BINARY_DIR}/lib/simple/simple_modules.f90\; rm -f  ${CMAKE_BINARY_DIR}/lib/simple/simple_modules_tmp.f90\;fi
  BYPRODUCTS "${CMAKE_BINARY_DIR}/lib/simple/simple_modules.f90"
  WORKING_DIRECTORY ${CMAKE_BINARY_DIR}
  USES_TERMINAL
  )

  if(BUILD_SHARED_LIBS)
    add_library(${SIMPLELIB}  ${SIMPLECORELIB})
  else()
    add_library(${SIMPLELIB} STATIC ${SIMPLECORELIB})
  endif()

add_dependencies(${SIMPLELIB} simple_args simple_git simple_modules)
set_target_properties(${SIMPLELIB}
  PROPERTIES
  LINKER_LANGUAGE Fortran
  SOVERSION ${${PROJECT_NAME}_VERSION}
  VERSION ${${PROJECT_NAME}_VERSION}
  COMPILE_FLAGS "${OpenMP_Fortran_FLAGS}"
  LINK_FLAGS "${OpenMP_Fortran_FLAGS}"
  LINK_LIBRARIES "${EXTRA_LIBS}"
  JOB_POOL_LINK NUM_JOBS
  )

install(TARGETS ${SIMPLELIB}
  ARCHIVE DESTINATION "${CMAKE_INSTALL_LIBDIR}"
  LIBRARY DESTINATION "${CMAKE_INSTALL_LIBDIR}"
  RUNTIME DESTINATION "${CMAKE_INSTALL_LIBDIR}"
  )
#install varlist for testing in production
install(FILES "${CMAKE_BINARY_DIR}/lib/simple/simple_varlist.txt"
  DESTINATION "${CMAKE_INSTALL_LIBDIR}/simple/")
# Make sure the compiler can find include files for our SIMPLE library
# when other libraries or executables link to it
target_include_directories (${SIMPLELIB} PUBLIC ${CMAKE_CURRENT_SOURCE_DIR})
