
# Do not use BIBFILES unless you plan to add citations
add_latex_document(simple3-0_devel.tex
   IMAGES ../SIMPLE_logo/rawlogo.png ./Commander/Commander.png ./LayeredArch/layered_arch.png ./SuppleDesign/SuppleDesign.png
   TARGET_NAME SimpleDeveloperGuide
   FORCE_PDF
   EXCLUDE_FROM_ALL
   EXCLUDE_FROM_DEFAULTS)

 install(FILES ${CMAKE_CURRENT_BINARY_DIR}/simple3-0_devel.pdf
   DESTINATION "${CMAKE_INSTALL_PREFIX}/doc"
   CONFIGURATIONS RELEASE
   RENAME SimpleDeveloperGuide.pdf OPTIONAL)
