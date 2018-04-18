
##########################################################################
# Source for the SIMPLE library.
# Better config and build alternative to old perl compile process
##########################################################################

add_custom_target(simple_args
  COMMAND SIMPLE_PATH=${CMAKE_BINARY_DIR} ${CMAKE_SOURCE_DIR}/scripts/simple_args_generator.pl
  DEPENDS "${CMAKE_SOURCE_DIR}/scripts/simple_args_generator.pl" "${CMAKE_CURRENT_SOURCE_DIR}/main/simple_params.f90"
  BYPRODUCTS "${CMAKE_BINARY_DIR}/lib/simple/simple_args.f90" "${CMAKE_BINARY_DIR}/lib/simple/simple_varlist.txt"
  WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}/main
  USES_TERMINAL
  )

add_custom_target(simple_gen_doc_f90
  COMMAND [ -x "${CMAKE_INSTALL_PREFIX}/bin/simple_exec" ] && SIMPLE_PATH=${CMAKE_INSTALL_PREFIX} ${CMAKE_SOURCE_DIR}/scripts/gen_simple_doc_primer.pl f90 > ${CMAKE_BINARY_DIR}/lib/simple/simple_gen_doc.f90
  DEPENDS "${CMAKE_SOURCE_DIR}/scripts/gen_simple_doc_primer.pl"
  BYPRODUCTS "${CMAKE_BINARY_DIR}/lib/simple/simple_gen_doc.f90"
  WORKING_DIRECTORY ${CMAKE_SOURCE_DIR}
  COMMENT "Create new simple_gen_doc.f90 ... requires existing build of simple"
  USES_TERMINAL
  )

# if(NOT EXISTS "${CMAKE_BINARY_DIR}/lib/simple/simple_gen_doc.f90")
# add_custom_target(simple_gen_doc
#   COMMAND cp "${CMAKE_CURRENT_SOURCE_DIR}/main/simple_gen_doc.f90" "${CMAKE_BINARY_DIR}/lib/simple/simple_gen_doc.f90"
#   DEPENDS "${CMAKE_CURRENT_SOURCE_DIR}/main/simple_gen_doc.f90"
#   WORKING_DIRECTORY ${CMAKE_SOURCE_DIR}
#   USES_TERMINAL
#   )
# else()
# add_custom_target(simple_gen_doc
#   COMMAND echo "Leaving gen_doc alone"
#   DEPENDS "${CMAKE_SOURCE_DIR}/scripts/gen_simple_doc_primer.pl"
#   BYPRODUCTS "${CMAKE_BINARY_DIR}/lib/simple/simple_gen_doc.f90"
#   WORKING_DIRECTORY ${CMAKE_SOURCE_DIR}
#   USES_TERMINAL
#   )
# endif()

file(GLOB DEFS_src "${CMAKE_CURRENT_SOURCE_DIR}/defs/*.f[90][08]" PARENT_SCOPE)
file(GLOB FILEIO_src "${CMAKE_CURRENT_SOURCE_DIR}/fileio/*.f[90][08]" PARENT_SCOPE)
 list(APPEND FILEIO_src
    "${CMAKE_CURRENT_SOURCE_DIR}/fileio/simple_posix.c"
    )
file(GLOB OPERATORS_src "${CMAKE_CURRENT_SOURCE_DIR}/ops/*.f90" PARENT_SCOPE)
file(GLOB UTILS_src "${CMAKE_CURRENT_SOURCE_DIR}/utils/*.f90" PARENT_SCOPE)
if ( USE_CUDA )
  file(GLOB CUDA_src "${CMAKE_CURRENT_SOURCE_DIR}/cuda/*.f[90][08]" PARENT_SCOPE)
  list(APPEND UTILS_src ${UTILS_src} ${CUDA_src})
##  list(REMOVE_ITEM UTILS_src "${CMAKE_CURRENT_SOURCE_DIR}/utils/simple_timer_cuda.f90")
## else()
##   list(APPEND UTILS_src ${UTILS_src} "${CMAKE_SOURCE_DIR}/utils/simple_timer_cuda.f90")
endif()
file(GLOB MAIN_src "${CMAKE_CURRENT_SOURCE_DIR}/main/*(^commander)*.f[90][08]" PARENT_SCOPE)
#
  file(GLOB IMG_src "${CMAKE_CURRENT_SOURCE_DIR}/img/*.f[90][08]" PARENT_SCOPE)
  list(APPEND IMG_src
    "${CMAKE_CURRENT_SOURCE_DIR}/img/stb_image_write.c"
    )
  if ( BUILD_WITH_LIBGD )
    list(APPEND IMG_src
    "${CMAKE_CURRENT_SOURCE_DIR}/img/cgd.c"
    "${CMAKE_CURRENT_SOURCE_DIR}/img/cjpg.c"
    )
endif()


########################################################################
# main folder
##########################################################################
file(GLOB SIMPLE_src "${CMAKE_CURRENT_SOURCE_DIR}/main/*.f[09][0-9]" PARENT_SCOPE)

# add the simple_args module to SIMPLE_src
list(APPEND SIMPLE_src
  "${CMAKE_BINARY_DIR}/lib/simple/simple_args.f90"
   "${CMAKE_BINARY_DIR}/lib/simple/simple_modules.f90"
)

#set(CORELIB_src  ${DEFS_src} ${UTILS_src} ${SIMPLE_src})
list(APPEND SIMPLECORELIB  ${DEFS_src} ${UTILS_src} ${FILEIO_src} ${IMG_src} ${OPERATORS_src} ${SIMPLE_src})
# message(STATUS "SIMPLE LIBRARY SOURCES: ${SIMPLECORELIB}")

##########################################################################
# Remove unused files here
##########################################################################
#list(REMOVE_ITEM SIMPLECORELIB
#    "${CMAKE_CURRENT_SOURCE_DIR}/fileio/simple_png.f90"
#)

add_custom_target(simple_modules
  COMMAND  rm -f ${CMAKE_BINARY_DIR}/lib/simple/simple_modules.f90\;  (echo module simple_modules \; for file in  ${DEFS_src} ${UTILS_src} ${FILEIO_src} ${IMG_src} ${OPERATORS_src} ${MAIN_src}\; do    if  \\[ \"\$\$\{file\#\#*.\}\" = \"f90\" \\] &&  test ! \"\$\${file}\" != \"\$\${file%simple_modules*}\" \; then      echo use \$\$\(basename \$\$file .f90\)\;    fi\;  done \;  echo contains\;  echo end module ) >> ${CMAKE_BINARY_DIR}/lib/simple/simple_modules.f90
  BYPRODUCTS "${CMAKE_BINARY_DIR}/lib/simple/simple_modules.f90"
  WORKING_DIRECTORY ${CMAKE_BINARY_DIR}
  USES_TERMINAL
  )



if (${CMAKE_Fortran_COMPILER_ID} STREQUAL "Intel")
  add_library(${SIMPLELIB} STATIC ${SIMPLECORELIB} )
  SET_TARGET_PROPERTIES(${SIMPLELIB}
    PROPERTIES
    LINKER_LANGUAGE Fortran
    SOVERSION ${${PROJECT_NAME}_VERSION}
    )
else()
    add_library(${SIMPLELIB}  ${SIMPLECORELIB})
endif()

add_dependencies(${SIMPLELIB} simple_args simple_modules)
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
