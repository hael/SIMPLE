rem Make print_bib and select_bib programs
rem
gfortran -o print_bib xmlparse.f90 read_xml_prims.f90 write_xml_prims.f90 bibliography.f90 print_bib.f90
gfortran -o select_bib xmlparse.f90 select_bib.f90
