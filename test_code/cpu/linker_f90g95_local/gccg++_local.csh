FLNAME=`basename $1 .cpp`

GFORTRAN=gcc-4.9

SOURCE_DIR=/home/frederic/Monash/SourceCode/Simple/Restructured/HansVersion/Simple_Restruct.projet

MODDIR=$SOURCE_DIR/include/simple
OBJDIR=$SOURCE_DIR/obj/GFORTgpu

OPTION="-cpp -O3 -DBENCH"

$GFORTRAN $OPTION -I $OBJDIR -I $MODDIR -o $FLNAME                         \
                                           $FLNAME.cpp
