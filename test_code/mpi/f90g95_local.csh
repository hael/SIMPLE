FLNAME=`basename $1 .f90`

GFORTRAN=gfortran-4.9

SOURCE_DIR=/home/frederic/Monash/SourceCode/Simple/Restructured/HansVersion/Simple_Restruct.projet

MODDIR=$SOURCE_DIR/obj/GFORTgpu
OBJDIR=$SOURCE_DIR/obj/GFORTgpu

MAGMADIR=/opt/Devel_tools/magma-1.6.1
MAGMADIR_LIB=$MAGMADIR/lib

CUDADIR=/usr/local/cuda

SHARE_LIB=/opt/Devel_tools/NVIDIA_GPU_Computing_SDK/shared/lib/linux

OPTION="-ffree-form -cpp -O3 -fPIC -fno-second-underscore -fopenmp -DBENCH -DMAGMA"

FFTW_LIB=/usr/lib/x86_64-linux-gnu

$GFORTRAN $OPTION -I $OBJDIR -I $MODDIR -o $FLNAME \
                                           $FLNAME.f90 \
                                           $OBJDIR/simple_defs.o \
                                           $OBJDIR/matrixGetter.o \
                                           $OBJDIR/initialu.o \
                                           $OBJDIR/su3random.o \
                                           $OBJDIR/su2random.o \
                                           $OBJDIR/simple_SU3_tester.o \
                                           $OBJDIR/simple_timing.o \
                                           $OBJDIR/simple_random.o \
                                           $OBJDIR/simple_sorting.o \
                                           $OBJDIR/timestamp.o \
                                           $OBJDIR/timming.o \
                                           $OBJDIR/timming_c.o \
                                           $OBJDIR/get_cpu_time_c.o \
                                           $OBJDIR/simple_textHandler.o \
                                           $OBJDIR/greeting_version.o \
                                           $OBJDIR/simple_testfunction.o \
                                           $OBJDIR/simple_matrixHandler_mpi_cpu.o \
                                           $OBJDIR/simple_math_gpu_c.o \
                                           $OBJDIR/simple_math.o \
                                           -L $FFTW_LIB -lfftw3 \
                                           -L $FFTW_LIB -lfftw3f \
                                           -L $FFTW_LIB -lfftw3_threads \
                                           -lblas \
                                           -llapack \
                                           -lstdc++ \
                                           -lpthread \
                                           -lm


