
##########################################################################
# Source for the SIMPLE library.
# Better config and build alternative to old perl compile process
##########################################################################

add_custom_target(simple_args
  COMMAND SIMPLE_PATH=${CMAKE_BINARY_DIR} ${CMAKE_SOURCE_DIR}/scripts/simple_args_generator.pl
  DEPENDS "${CMAKE_SOURCE_DIR}/scripts/simple_args_generator.pl" "${CMAKE_CURRENT_SOURCE_DIR}/main/simple_parameters.f90"
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

file(GLOB DEFS_src "${CMAKE_CURRENT_SOURCE_DIR}/defs/*.f[90][08]" PARENT_SCOPE)
file(GLOB FILEIO_src "${CMAKE_CURRENT_SOURCE_DIR}/fileio/*.f[90][08]" PARENT_SCOPE)
list(APPEND FILEIO_src
  "${CMAKE_CURRENT_SOURCE_DIR}/fileio/simple_posix.c"
  )
file(GLOB OPERATORS_src "${CMAKE_CURRENT_SOURCE_DIR}/ops/*.f90" PARENT_SCOPE)
list(APPEND OPERATORS_src
  "${CMAKE_CURRENT_SOURCE_DIR}/ops/simple_kbinterpol_memo.c"
  )
file(GLOB UTILS_src "${CMAKE_CURRENT_SOURCE_DIR}/utils/*.f90" PARENT_SCOPE)
#list(REMOVE_ITEM  UTILS_src
  #"${CMAKE_CURRENT_SOURCE_DIR}/utils/simple_qsort_mt.f90"
#  )
if (NOT ${CMAKE_Fortran_COMPILER_ID} STREQUAL "Intel")
  list(REMOVE_ITEM  UTILS_src
    "${CMAKE_CURRENT_SOURCE_DIR}/utils/simple_intel_qsort.f90"
    )
endif()

file(GLOB STARFILE_src "${CMAKE_CURRENT_SOURCE_DIR}/starfile/*.cpp" PARENT_SCOPE)
file(GLOB CPLOT2D_src "${CMAKE_CURRENT_SOURCE_DIR}/CPlot2D/*.cpp" PARENT_SCOPE)
file(GLOB CUDA_src "${CMAKE_CURRENT_SOURCE_DIR}/cuda/*.[Ff][90][08]" PARENT_SCOPE)
message(STATUS "SIMPLE CUDA LIBRARY ${CUDA_src}")

file(GLOB MAIN_src "${CMAKE_CURRENT_SOURCE_DIR}/main/*(^commander)*.f[90][08]" PARENT_SCOPE)
file(GLOB JPG_src "${CMAKE_CURRENT_SOURCE_DIR}/jpg/*.f[90][08]" PARENT_SCOPE)
list(APPEND JPG_src
  "${CMAKE_CURRENT_SOURCE_DIR}/jpg/stb_image_write.c")
list(REMOVE_ITEM JPG_src
  "${CMAKE_CURRENT_SOURCE_DIR}/jpg/simple_openslide.f90"
  )
if(USE_LIBTIFF)
list(APPEND JPG_src
  "${CMAKE_CURRENT_SOURCE_DIR}/jpg/write-bigtiff.c"
  "${CMAKE_CURRENT_SOURCE_DIR}/jpg/write_tiff.cpp"
  "${CMAKE_CURRENT_SOURCE_DIR}/jpg/ctif.c"
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
list(APPEND SIMPLECORELIB  ${DEFS_src} ${UTILS_src} ${FILEIO_src} ${STARFILE_src} ${CPLOT2D_src} ${JPG_src} ${OPERATORS_src} ${SIMPLE_src})
if ( USE_CUDA )
  list(APPEND SIMPLECORELIB ${CUDA_src})

endif()
# message(STATUS "SIMPLE LIBRARY SOURCES: ${SIMPLECORELIB}")

##########################################################################
# Remove unused files here
##########################################################################
#list(REMOVE_ITEM SIMPLECORELIB
#  "${CMAKE_CURRENT_SOURCE_DIR}/main/simple_pftcc_shsrch.f90"
#  )

add_custom_target(simple_modules
  COMMAND  echo Creating simple_module \; ( echo module simple_modules\; for file in ${DEFS_src} ${UTILS_src} ${FILEIO_src} ${JPG_src} ${OPERATORS_src} ${MAIN_src}\; do if \\[ \"\$\$\{file\#\#*.\}\" = \"f90\" \\] &&  test ! \"\$\${file}\" != \"\$\${file%simple_modules*}\" \; then      echo use \$\$\(basename \$\$file .f90\)\;fi\;done\;echo contains\;echo end module ) >> ${CMAKE_BINARY_DIR}/lib/simple/simple_modules_tmp.f90
  COMMAND if \\[ ! -f ${CMAKE_BINARY_DIR}/lib/simple/simple_modules.f90 \\]\;then mv ${CMAKE_BINARY_DIR}/lib/simple/simple_modules_tmp.f90 ${CMAKE_BINARY_DIR}/lib/simple/simple_modules.f90\;else cmp --silent ${CMAKE_BINARY_DIR}/lib/simple/simple_modules_tmp.f90 ${CMAKE_BINARY_DIR}/lib/simple/simple_modules.f90 || cp ${CMAKE_BINARY_DIR}/lib/simple/simple_modules_tmp.f90 ${CMAKE_BINARY_DIR}/lib/simple/simple_modules.f90\; rm -f  ${CMAKE_BINARY_DIR}/lib/simple/simple_modules_tmp.f90\;fi
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
  if(BUILD_SHARED_LIBS)
    add_library(${SIMPLELIB}  ${SIMPLECORELIB})
  else()
    add_library(${SIMPLELIB} STATIC ${SIMPLECORELIB})
  endif()
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


if(USE_CUDA)
  file(GLOB CUKERNELS_src "${CMAKE_CURRENT_SOURCE_DIR}/cuda/kernels/*.cu" PARENT_SCOPE)
  if(NOT CUDA_TOOLKIT_SDK_DIR)
    list(REMOVE_ITEM CUKERNELS_src
      "${CMAKE_CURRENT_SOURCE_DIR}/cuda/kernels/blockreduction.cu"
      )
  endif()
  if(BUILD_SHARED_LIBS)
    CUDA_ADD_LIBRARY( cuda_kernels  ${CUKERNELS_src}
      OPTIONS DEBUG -g
      RELEASE --use_fast_math
      RELWITHDEBINFO --use_fast_math;-g
      MINSIZEREL --use_fast_math
      )
  else()
    CUDA_ADD_LIBRARY( cuda_kernels   ${CUKERNELS_src}
      STATIC
      OPTIONS DEBUG -g
      RELEASE --use_fast_math
      RELWITHDEBINFO --use_fast_math;-g
      MINSIZEREL --use_fast_math
      )
  endif()

endif()


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
