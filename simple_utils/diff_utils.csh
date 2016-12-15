#!/bin/bash
set MY_HOME_SIDE = /home/frederic/SimpleSideBrch 

set DIFF_PATH = $MY_HOME_SIDE/21Mar16/Simple_Restruct.projet/simple_utils

echo "::::::::::::::"
echo "$DIFF_PATH/diff_utils.csh ./diff_utils.csh"
echo "::::::::::::::"
diff diff_utils.csh $DIFF_PATH
echo "::::::::::::::"
echo "$DIFF_PATH/gnufor2.f90 ./gnufor2.f90"
echo "::::::::::::::"
diff gnufor2.f90 $DIFF_PATH
echo "::::::::::::::"
echo "$DIFF_PATH/Makefile_target ./Makefile_target"
echo "::::::::::::::"
diff Makefile_target $DIFF_PATH
echo "::::::::::::::"
echo ""
echo "::::::::::::::"
diff simple_jiffys.f90 $DIFF_PATH
echo "::::::::::::::"
echo ""
echo "::::::::::::::"
diff simple_math.f90 $DIFF_PATH
echo "::::::::::::::"
echo ""
echo "::::::::::::::"
diff simple_rnd.f90 $DIFF_PATH
echo "::::::::::::::"
echo ""
echo "::::::::::::::"
diff simple_stat.f90 $DIFF_PATH
echo "::::::::::::::"
echo ""
echo "::::::::::::::"
diff simple_strings.f90 $DIFF_PATH
echo "::::::::::::::"
echo ""
echo "::::::::::::::"
diff simple_syscalls.f90 $DIFF_PATH
