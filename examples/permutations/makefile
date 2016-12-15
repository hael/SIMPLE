# Using compile cmd gfortran
# Using openmp
# Using opt
# Using m64
# Main program is permutations.f90 
# process_fsource called with arg permutations.f90 
# Full list of modules in permutations.f90:  

# ------------------Macro-Defs---------------------
F90=gfortran 

# -------------------End-macro-Defs---------------------------

# Here is the link step 
permutations:permutations.o 
	 $(F90) -o permutations permutations.o   -fimplicit-none -fall-intrinsics -ftree-vectorize -O3 -L/sw/lib -lfftw3f -lfftw3f_threads

# Here are the compile steps
 
permutations.o:permutations.f90  
	 $(F90) -c -fimplicit-none -fall-intrinsics -ftree-vectorize -O3 -L/sw/lib -lfftw3f -lfftw3f_threads permutations.f90
# This entry allows you to type " make clean " to get rid of
# all object and module files 
clean:
	rm -f -r f_{files,modd}* *~ *.g90 *.o *.mod *.M *.d V*.inc *.vo \
	V*.f *.dbg album F.err
  
