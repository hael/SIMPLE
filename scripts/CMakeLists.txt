
include(${CMAKE_ROOT}/Modules/FindPerl.cmake)


list(APPEND SIMPLE_PERL_SCRIPTS
cleanSIMPLEout.pl
convert_frealign2simple.pl
convert_helical2box.pl
convert_relion2simple.pl
convert_simple2relion.py
cp_filetab.pl
dm42mrc.pl
emx2simple.pl
exec_batchbox.pl
exec_ctffind4.1.X.pl
exec_ctffind.pl
gen_simple_doc_primer.pl
get_PixSz_XML.pl
kill.pl
parse_ctffind3.pl
parse_ctffind4.pl
partition.pl
prime_cleanup.pl
puts_gitignore.pl
rawctf2simplectf.pl
relion2emanbox.pl
relion2simplectf.pl
relion2simple.pl
restart.pl
sam2emanbox.pl
simmovies.pl
simple_args_generator.pl
split_filetab.pl
tab2space.pl)

foreach(PERL_SCRIPT ${SIMPLE_PERL_SCRIPTS})
  install(FILES ${CMAKE_CURRENT_SOURCE_DIR}/${PERL_SCRIPT}
    DESTINATION ${SCRP_BIN}
    PERMISSIONS OWNER_READ GROUP_READ WORLD_READ OWNER_EXECUTE GROUP_EXECUTE WORLD_EXECUTE
    )
endforeach()
