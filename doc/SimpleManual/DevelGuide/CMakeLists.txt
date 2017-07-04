

add_latex_document(simple3-0_devel.tex
   BIBFILES Prime2bibfile.bib
   IMAGES ../SIMPLE_logo/rawlogo.png ../Commander/Commander.png ../LayeredArch/layered_arch.png ../SuppleDesign/SuppleDesign.png
   TARGET_NAME DeveloperGuide-Simple3
   FORCE_PDF
   EXCLUDE_FROM_ALL
   EXCLUDE_FROM_DEFAULTS)

 install(FILES ${CMAKE_CURRENT_BINARY_DIR}/simple3-0_devel.pdf
   DESTINATION "${CMAKE_INSTALL_PREFIX}/doc"
   CONFIGURATIONS RELEASE)
