

add_latex_document(simple2-5_usermanual.tex
   BIBFILES Prime2bibfile.bib
#   IMAGE_DIRS ../SIMPLE_logo ../CPUtopo
#   "${CMAKE_CURRENT_SOURCE_DIR}/SIMPLE_logo"
#   "${CMAKE_CURRENT_SOURCE_DIR}/CPUtopo"
   IMAGES ../SIMPLE_logo/rawlogo.png ../CPUtopo/cputopo.ai
   TARGET_NAME Simple2.5_User_Manual
   FORCE_PDF
   EXCLUDE_FROM_ALL
   EXCLUDE_FROM_DEFAULTS)

 install(FILES ${CMAKE_CURRENT_BINARY_DIR}/simple2-5_usermanual.pdf
   DESTINATION "${CMAKE_INSTALL_PREFIX}/doc"
   CONFIGURATIONS RELEASE
   RENAME Simple2-5_User_Manual.pdf OPTIONAL)
