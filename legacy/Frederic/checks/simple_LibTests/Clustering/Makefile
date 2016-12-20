##############################################
# settings
#
# COMPILER = gfortran-4.9
# CFLAGS   = -ffree-form -cpp -O3 -fno-second-underscore -fopenmp -DBENCH
# LINK     = -lblas -llapack -lstdc++ -lpthread -lm
#
# .SUFFIXES: .f90 .o
# .f90.o:
#        $(COMPILER) $(CFLAGS) -c $<
#
##############################################
# Massive settings
# module load fftw/3.3.4-gcc ;
# module load gcc/4.9.1 ;
# module load lapack/3.4.2 ;
# module load openmpi ;
#
COMPILER_GFORTRAN = gfortran-4.9
FFTW_LIB = /usr/lib/x86_64-linux-gnu
CFLAGS   = -ffree-form -cpp -O3 -fno-second-underscore -fopenmp -DLINUX -DBENCH
LINK     = -L $(FFTW_LIB) -lfftw3                   \
           -L $(FFTW_LIB) -lfftw3f                  \
           -lblas -llapack -lstdc++ -lpthread -lm

SIMPLE_SOURCE=../../../
XC_DIR=$(SIMPLE_SOURCE)/checks/simple_LibTests/Clustering
COMDIR = $(SIMPLE_SOURCE)/obj/GFORTgpu
OBJDIR = $(SIMPLE_SOURCE)/obj/GFORTgpu

.SUFFIXES: .f90 .o 
.f90.o:
	$(COMPILER_GFORTRAN) $(CFLAGS) -I $(OBJDIR) -c $<

COMMON =

#Prarralle IO
MAIN =
#Executables
XCS = 
#Compilations
all: $(XCS)

mc: $(COMMON) $(MAIN)
	$(COMPILER_GFORTRAN) $(CFLAGS) -o $(XC_DIR)/mc \
        $(COMMON) $(MAIN) $(LINK)

clean :
	rm -rf *.o *.d *~* *#* *.spi *.mrc prime2D_* *.txt *.log fort.2500
cleanall:




