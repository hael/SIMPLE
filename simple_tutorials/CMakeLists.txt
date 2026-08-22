
add_latex_document(Tutorials.tex
   BIBFILES Prime2bibfile.bib
   IMAGES ./SIMPLE_logo/rawlogo.png  ./CPUtopo/cputopo.ai    ./pspec/pspec.png           ./straightvsintg/straightvsintg.ai           ./spectra/pspecs.ai           ./astigctf/astigctf.ai          ./trpv1ptcls.png  ./trpv1iters/trpv1iters.ai   ./ribocavgs.png  ./Pfrib3D/Pfrib3D.ai  ./bgal/bgal_cavgs.png  ./bgal/bgal_cavgs_sel.png  ./bgal/bgal_maps.ai
   TARGET_NAME SimpleTutorial-March2017
   FORCE_PDF
   EXCLUDE_FROM_ALL
   EXCLUDE_FROM_DEFAULTS)

 install(FILES ${CMAKE_CURRENT_BINARY_DIR}/SimpleTutorial-March2017.pdf
   DESTINATION "${CMAKE_INSTALL_PREFIX}/doc"
   CONFIGURATIONS RELEASE
   OPTIONAL)
