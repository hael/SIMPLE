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
XC_DIR=$(SIMPLE_SOURCE)/checks/simple_LibTests/BenchMarking
COMDIR = $(SIMPLE_SOURCE)/obj/GFORTgpu
OBJDIR = $(SIMPLE_SOURCE)/obj/GFORTgpu

.SUFFIXES: .f90 .o 
.f90.o:
	$(COMPILER_GFORTRAN) $(CFLAGS) -I $(OBJDIR) -c $<

COMMON = $(OBJDIR)/simple_defs.o \
         $(OBJDIR)/matrixGetter.o \
         $(OBJDIR)/simple_testfunction.o \
         $(OBJDIR)/greeting_version.o \
         $(OBJDIR)/simple_timing.o \
         $(OBJDIR)/timming.o \
         $(OBJDIR)/timming_c.o \
         $(OBJDIR)/get_cpu_time_c.o \
         $(OBJDIR)/timestamp.o \
         $(OBJDIR)/initialu.o \
         $(OBJDIR)/su3random.o \
         $(OBJDIR)/su2random.o \
         $(OBJDIR)/simple_sorting.o \
         $(OBJDIR)/simple_textHandler.o \
         $(OBJDIR)/get_systemQuery_cpu.o \
         $(OBJDIR)/get_systemQuery_cpu_c.o \
         $(OBJDIR)/simple_systemQuery_cpu.o \
         $(OBJDIR)/simple_file_utils.o \
         $(OBJDIR)/simple_yaml_output.o \
         $(OBJDIR)/simple_yaml_strings.o \
         $(OBJDIR)/strlen.o \
         $(OBJDIR)/simple_eglossary.o \
         $(OBJDIR)/simple_err_defs.o \
         $(OBJDIR)/simple_error_handling.o \
         $(OBJDIR)/simple_eglossary_lowlev.o \
         $(OBJDIR)/simple_file_highlev.o \
         $(OBJDIR)/simple_dynamic_memory.o \
         $(OBJDIR)/simple_memory_profiling.o \
         $(OBJDIR)/simple_params.o \
         $(OBJDIR)/simple_cmdline.o \
         $(OBJDIR)/simple_jiffys.o \
         $(OBJDIR)/simple_math.o \
         $(OBJDIR)/simple_ori.o \
         $(OBJDIR)/simple_imghead.o \
         $(OBJDIR)/simple_imgheadrec.o \
         $(OBJDIR)/simple_args.o \
         $(OBJDIR)/simple_rnd.o \
         $(OBJDIR)/simple_hash.o \
         $(OBJDIR)/simple_softmax_weights.o \
         $(OBJDIR)/simple_ran_tabu.o \
         $(OBJDIR)/simple_stat.o \
         $(OBJDIR)/simple_AVratios.o \
         $(OBJDIR)/simple_strings.o \
         $(OBJDIR)/gnufor2.o \
         $(OBJDIR)/simple_oris.o \
         $(OBJDIR)/simple_cmd_dict.o \
         $(OBJDIR)/simple_convergence.o \
         $(OBJDIR)/simple_opt_spec.o \
         $(OBJDIR)/simple_simplex_opt.o \
         $(OBJDIR)/simple_opt_subs.o \
         $(OBJDIR)/simple_optimizer.o

#Data analyzer
MAIN = simple_BenchMrkDtAnlzr.o
#Executables
XCS = mcBenchMrkDtAnlzr

#Compilations
all: $(XCS)

mcBenchMrkDtAnlzr: $(COMMON) $(MAIN)
	$(COMPILER_GFORTRAN) $(CFLAGS) -o $(XC_DIR)/$(XCS) \
        $(COMMON) $(MAIN) $(LINK)

clean :
	rm -f f *.o *.pyc *.d *~* *#* $(XC_DIR)/$(XCS) \
 z_prime3D.asc z_prime3D_nrotnk*

cleanall:
#Not needed yet but if need be will be implemented
	rm -rf *.o *.pyc *.d *~* *#* *.asc *_lp*_bz*_smpd*_Bunch \
 *_lp*_bz*_smpd*_msk* *.mrc msk*_*\
 $(XC_DIR)/$(XCS) z_prime3D.asc z_prime3D_nrotnk*
