#!/bin/bash
set MY_HOME_SIDE = /home/frederic/SimpleSideBrch/Laptop/10Jul16/Simple_Restruct.projet

set DIFF_PATH = $MY_HOME_SIDE/defs

echo "::::::::::::::"
echo "diff_defs.csh"
echo "::::::::::::::"
diff diff_defs.csh $DIFF_PATH
echo "::::::::::::::"
echo "./Makefile_target $DIFF_PATH/Makefile_target "
echo "::::::::::::::"
diff Makefile_target $DIFF_PATH
echo "::::::::::::::"
echo "nr.f90"
echo "::::::::::::::"
diff nr.f90 $DIFF_PATH
echo "::::::::::::::"
echo "nrtype.f90"
echo "::::::::::::::"
diff nrtype.f90 $DIFF_PATH
echo "::::::::::::::"
echo "nrutil.f90"
echo "::::::::::::::"
diff nrutil.f90 $DIFF_PATH
echo "::::::::::::::"
echo "simple_cuda_defs.f90"
echo "::::::::::::::"
diff simple_cuda_defs.f90 $DIFF_PATH
echo "::::::::::::::"
echo "simple_deepLearning_defs.f90"
echo "::::::::::::::"
diff simple_deepLearning_defs.f90 $DIFF_PATH
echo "::::::::::::::"
echo "simple_defs.f90"
echo "::::::::::::::"
diff simple_defs.f90 $DIFF_PATH
echo "::::::::::::::"
echo "simple_err_defs.f90"
echo "::::::::::::::"
diff simple_err_defs.f90 $DIFF_PATH
echo "::::::::::::::"
echo "simple_fftw3.f90"
echo "::::::::::::::"
diff simple_fftw3.f90 $DIFF_PATH
echo "::::::::::::::"
echo "simple_file_defs.f90"
echo "::::::::::::::"
diff simple_file_defs.f90 $DIFF_PATH
echo "::::::::::::::"
echo "simple_lattice_defs.f90"
echo "::::::::::::::"
diff simple_lattice_defs.f90 $DIFF_PATH
echo "::::::::::::::"
echo "simple_mpi_defs.f90"
echo "::::::::::::::"
diff simple_mpi_defs.f90 $DIFF_PATH
echo "::::::::::::::"
echo "simple_OpenCL_defs.f90"
echo "::::::::::::::"
diff simple_OpenCL_defs.f90 $DIFF_PATH
echo "::::::::::::::"
echo "simple_resources_defs.f90"
echo "::::::::::::::"
diff simple_resources_defs.f90 $DIFF_PATH
echo "::::::::::::::"
echo "simple_units.f90"
echo "::::::::::::::"
diff simple_units.f90 $DIFF_PATH
