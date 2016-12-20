FLNAME=`basename $1 .f90`

MPIF90=/usr/lib64/mpi/gcc/openmpi/bin/mpif90

MPIDIR=/usr/lib64/mpi/gcc/openmpi/include
MPEDIR=/opt/Devel_tools/mpich2-1.5/src/mpe2

MODDIR=/home/Frederic-PC.linux64/SourceCode/PhysSim/include/QCD
OBJDIR=/home/Frederic-PC.linux64/SourceCode/PhysSim/obj/GFORTgpu

OPTION="-ffree-form -cpp -O3 -fno-second-underscore -DBENCH -DMAGMA -DCUDA -DPHYSSYM_MPI"

$MPIF90 $OPTION -I $OBJDIR -I $MODDIR -I $MPIDIR -o $FLNAME               \
                                         $FLNAME.f90                      \
                                        $OBJDIR/physsym_mpi_defs.o       \
                                         $OBJDIR/matrixGetter.o           \
                                         $OBJDIR/matrixMultiplier_cpu.o   \
                                         $OBJDIR/initialu.o               \
                                         $OBJDIR/su3random.o              \
                                         $OBJDIR/su2random.o              \
                                         $OBJDIR/random.o                 \
                                         $OBJDIR/physsym_sorting.o        \
                                         $OBJDIR/timestamp.o              \
                                         $OBJDIR/physsym_textHandler.o    \
                                         $OBJDIR/timming_c.o              \
                                         $OBJDIR/physsym_timing.o         \
                                         $OBJDIR/greeting_version.o
