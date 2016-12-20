FLNAME=`basename $1 .f90`

GFORTRAN=gfortran

MODDIR=/home/frederic/SourceCode/PhysSim/include/QCD
#/home/Frederic-PC.linux64/SourceCode/PhysSim/include/QCD
OBJDIR=/home/frederic/SourceCode/PhysSim/obj/GFORTgpu
#/home/Frederic-PC.linux64/SourceCode/PhysSim/obj/GFORTgpu

OPTION="-ffree-form -cpp -O3 -fno-second-underscore -DBENCH -DMAGMA -DCUDA"

$GFORTRAN $OPTION -I $OBJDIR -I $MODDIR -o $FLNAME                           \
                                         $FLNAME.f90                         \
                                         $OBJDIR/matrixGetter.o              \
                                         $OBJDIR/physsym_matrixHandler_cpu.o \
                                         $OBJDIR/initialu.o                  \
                                         $OBJDIR/su3random.o                 \
                                         $OBJDIR/su2random.o                 \
                                         $OBJDIR/random.o                    \
                                         $OBJDIR/physsym_sorting.o           \
                                         $OBJDIR/timestamp.o                 \
                                         $OBJDIR/physsym_textHandler.o       \
                                         $OBJDIR/timming_c.o                 \
                                         $OBJDIR/physsym_timing.o            \
                                         $OBJDIR/greeting_version.o
