# Main program is simple_test_foofactory.f90 
# process_fsource called with arg simple_test_foofactory.f90 
# simple_test_foofactory.f90 Uses Module simple_foofactory
# simple_test_foofactory.f90 Uses Module simple_foo
# Full list of modules in simple_test_foofactory.f90: simple_foofactory simple_foo 
# Uses simple_foofactory which is in ../simple_foofactory.f90
# process_fsource called with arg ../simple_foofactory.f90 
# ../simple_foofactory.f90 Uses Module simple_foo
# ../simple_foofactory.f90 Uses Module simple_foor
# ../simple_foofactory.f90 Uses Module simple_fooi
# Full list of modules in ../simple_foofactory.f90: simple_foo simple_foor simple_fooi 
# Uses simple_foo which is in ../simple_foo.f90
# process_fsource called with arg ../simple_foo.f90 
# Full list of modules in ../simple_foo.f90:  
# Uses simple_foor which is in ../simple_foor.f90
# process_fsource called with arg ../simple_foor.f90 
# ../simple_foor.f90 Uses Module simple_foo
# Full list of modules in ../simple_foor.f90: simple_foo 
# Uses simple_foo which is in ../simple_foo.f90
# ../simple_foo.f90 already in list
# Uses simple_fooi which is in ../simple_fooi.f90
# process_fsource called with arg ../simple_fooi.f90 
# ../simple_fooi.f90 Uses Module simple_foo
# Full list of modules in ../simple_fooi.f90: simple_foo 
# Uses simple_foo which is in ../simple_foo.f90
# ../simple_foo.f90 already in list
# Uses simple_foo which is in ../simple_foo.f90
# ../simple_foo.f90 already in list

# ------------------Macro-Defs---------------------
F90=gfortran 

# -------------------End-macro-Defs---------------------------

# Here is the link step 
simple_test_foofactory:simple_foo.o simple_foor.o simple_fooi.o simple_foofactory.o simple_test_foofactory.o 
	 $(F90) -o simple_test_foofactory simple_foo.o simple_foor.o simple_fooi.o simple_foofactory.o simple_test_foofactory.o  -fimplicit-none -fall-intrinsics -ftree-vectorize -O3 -L/sw/lib/ -lfftw3f -lfftw3f_threads -fopenmp

# Here are the compile steps
 
simple_foo.o:../simple_foo.f90  
	 $(F90) -c  -fimplicit-none -fall-intrinsics -ftree-vectorize -O3 -L/sw/lib/ -lfftw3f -lfftw3f_threads -fopenmp ../simple_foo.f90

simple_foor.o:../simple_foor.f90 simple_foo.o 
	 $(F90) -c  -fimplicit-none -fall-intrinsics -ftree-vectorize -O3 -L/sw/lib/ -lfftw3f -lfftw3f_threads -fopenmp ../simple_foor.f90

simple_fooi.o:../simple_fooi.f90 simple_foo.o 
	 $(F90) -c  -fimplicit-none -fall-intrinsics -ftree-vectorize -O3 -L/sw/lib/ -lfftw3f -lfftw3f_threads -fopenmp ../simple_fooi.f90

simple_foofactory.o:../simple_foofactory.f90 simple_foo.o simple_foor.o simple_fooi.o 
	 $(F90) -c  -fimplicit-none -fall-intrinsics -ftree-vectorize -O3 -L/sw/lib/ -lfftw3f -lfftw3f_threads -fopenmp ../simple_foofactory.f90

simple_test_foofactory.o:simple_test_foofactory.f90 simple_foofactory.o simple_foo.o 
	 $(F90) -c  -fimplicit-none -fall-intrinsics -ftree-vectorize -O3 -L/sw/lib/ -lfftw3f -lfftw3f_threads -fopenmp simple_test_foofactory.f90
# This entry allows you to type " make clean " to get rid of
# all object and module files 
clean:
	rm -f -r f_{files,modd}* *~ *.g90 *.o *.mod *.M *.d V*.inc *.vo \
	V*.f *.dbg album F.err
  
