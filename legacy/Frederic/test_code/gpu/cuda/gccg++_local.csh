FLNAME=`basename $1 .cpp`

GFORTRAN=gcc-4.9

SOURCE_DIR=/home/frederic/Monash/SourceCode/Simple/Restructured/HansVersion/Simple_Restruct.projet

MODDIR=$SOURCE_DIR/include/simple
OBJDIR=$SOURCE_DIR/obj/GFORTgpu

MAGMADIR=/opt/Devel_tools/magma-1.6.1
MAGMADIR_LIB=$MAGMADIR/lib

CUDADIR=/usr/local/cuda
CUDAINC=$CUDADIR/include

SHARE_LIB=/opt/Devel_tools/NVIDIA_GPU_Computing_SDK/shared/lib/linux

PLATFORM=`uname`
if [ "$PLATFORM" = "Linux" ]
then
 PLAT="-DLINUX"
 elif [ "$PLATFORM" = "Darwin" ]
then
 PLAT="-DMACOSX"
else
echo "You seem to be compiling/running in wonderland"
fi

OPTION="-cpp -O3 -fopenmp $PLAT -DOPENMP -DBENCH -DCUDA -DMAGMA"

$GFORTRAN $OPTION \
          -I $OBJDIR -I $MODDIR -I $CUDAINC -o $FLNAME \
                                               $FLNAME.cpp \
                                               $OBJDIR/timming.o \
                                               $OBJDIR/timming_c.o \
                                               -L $FFTW_LIB -lfftw3 \
                                               -L $FFTW_LIB -lfftw3f \
                                               -L $FFTW_LIB -lfftw3_threads \
                                               -lblas \
                                               -llapack \
                                               -lstdc++ \
                                               -lrt \
                                               -lpthread \
                                               -lcuda \
                                               -lcublas \
                                               -lcudart \
                                               -lcufft \
                                               -lcudnn \
                                               -L $CUDADIR/lib64 \
                                               -lmagma \
                                               -L $MAGMADIR_LIB \
                                               -lm
