#!/bin/bash
set MY_HOME_SRC = /home/frederic/Monash/SourceCode/Simple/Restructured/HansVersion
set MY_HOME_SIDE = /home/frederic/SimpleSideBrch 

set DIFF_PATH = $MY_HOME_SIDE/21Mar16/Simple_Restruct.projet/test_code/cpu
set DIFF_PATH_TIMER = $MY_HOME_SRC/backups/Simple_Restruct.projet/test_code/cpu

echo "::::::::::::::"
echo "clean_All_Testing_cpu.csh"
echo "::::::::::::::"
diff clean_All_Testing_cpu.csh $DIFF_PATH
echo "::::::::::::::"
echo "compile_All_Testing_cpu.csh"
echo "::::::::::::::"
diff compile_All_Testing_cpu.csh $DIFF_PATH
echo "::::::::::::::"
echo "diff_test_TestCPU.csh"
echo "::::::::::::::"
diff diff_test_TestCPU.csh $DIFF_PATH
echo "::::::::::::::"
echo "diff_timer_tester.csh"
echo "::::::::::::::"
diff diff_timer_tester.csh $DIFF_PATH
echo "::::::::::::::"
echo "f90g95_local.csh"
echo "::::::::::::::"
diff f90g95_local.csh $DIFF_PATH
echo "::::::::::::::"
echo "functiontest_2d_cpu.m"
echo "::::::::::::::"
diff functiontest_2d_cpu.m $DIFF_PATH
echo "::::::::::::::"
echo "gccg++_local.csh"
echo "::::::::::::::"
diff gccg++_local.csh $DIFF_PATH
echo "::::::::::::::"
echo "Makefile_target"
echo "::::::::::::::"
diff Makefile_target $DIFF_PATH
echo "::::::::::::::"
echo "run_All_Testing_cpu.csh"
echo "::::::::::::::"
diff run_All_Testing_cpu.csh $DIFF_PATH
echo "::::::::::::::"
echo "simple_contopt_tester.f90"
echo "::::::::::::::"
diff simple_contopt_tester.f90 $DIFF_PATH
echo "::::::::::::::"
echo "simple_polarft_tester.f90"
echo "::::::::::::::"
diff simple_polarft_tester.f90 $DIFF_PATH
echo "::::::::::::::"
echo "testing_BFGSNum_2D_cpu.f90"
echo "::::::::::::::"
diff testing_BFGSNum_2D_cpu.f90 $DIFF_PATH
echo "::::::::::::::"
echo "testing_congruence.f90"
echo "::::::::::::::"
diff testing_congruence.f90 $DIFF_PATH
echo "::::::::::::::"
echo "testing_corr_cpu.f90"
echo "::::::::::::::"
diff testing_corr_cpu.f90 $DIFF_PATH
echo "::::::::::::::"
echo "testing_corr_gpu.f90"
echo "::::::::::::::"
diff testing_corr_gpu.f90 $DIFF_PATH
echo "::::::::::::::"
echo "testing_cshift.f90"
echo "::::::::::::::"
diff testing_cshift.f90 $DIFF_PATH
echo "::::::::::::::"
echo "testing_DataParallelIO.f90"
echo "::::::::::::::"
diff testing_DataParallelIO.f90 $DIFF_PATH
echo "::::::::::::::"
echo "testing_dgetri_cpu.f90"
echo "::::::::::::::"
diff testing_dgetri_cpu.f90 $DIFF_PATH
echo "::::::::::::::"
echo "testing_FFT_1D_Z2Z_cpu.f90"
echo "::::::::::::::"
diff testing_FFT_1D_Z2Z_cpu.f90 $DIFF_PATH
echo "::::::::::::::"
echo "testing_FFT_2D_C2C_cpu.f90"
echo "::::::::::::::"
diff testing_FFT_2D_C2C_cpu.f90 $DIFF_PATH
echo "::::::::::::::"
echo "testing_FFT_2D_C2S_cpu.f90"
echo "::::::::::::::"
diff testing_FFT_2D_C2S_cpu.f90 $DIFF_PATH
echo "::::::::::::::"
echo "testing_FFT_2D_D2Z_cpu.f90"
echo "::::::::::::::"
diff testing_FFT_2D_D2Z_cpu.f90 $DIFF_PATH
echo "::::::::::::::"
echo "testing_FFT_2D_S2C_cpu.f90"
echo "::::::::::::::"
diff testing_FFT_2D_S2C_cpu.f90 $DIFF_PATH
echo "::::::::::::::"
echo "testing_FFT_2D_Z2D_cpu.f90"
echo "::::::::::::::"
diff testing_FFT_2D_Z2D_cpu.f90 $DIFF_PATH
echo "::::::::::::::"
echo "testing_FFT_2D_Z2Z_cpu.f90"
echo "::::::::::::::"
diff testing_FFT_2D_Z2Z_cpu.f90 $DIFF_PATH
echo "::::::::::::::"
echo "testing_FFT_3D_C2C_cpu.f90"
echo "::::::::::::::"
diff testing_FFT_3D_C2C_cpu.f90 $DIFF_PATH
echo "::::::::::::::"
echo "testing_FFT_3D_C2S_cpu.f90"
echo "::::::::::::::"
diff testing_FFT_3D_C2S_cpu.f90 $DIFF_PATH
echo "::::::::::::::"
echo "testing_FFT_3D_D2Z_cpu.f90"
echo "::::::::::::::"
diff testing_FFT_3D_D2Z_cpu.f90 $DIFF_PATH
echo "::::::::::::::"
echo "testing_FFT_3D_S2C_cpu.f90"
echo "::::::::::::::"
diff testing_FFT_3D_S2C_cpu.f90 $DIFF_PATH
echo "::::::::::::::"
echo "testing_FFT_3D_Z2D_cpu.f90"
echo "::::::::::::::"
diff testing_FFT_3D_Z2D_cpu.f90 $DIFF_PATH
echo "::::::::::::::"
echo "testing_FFT_3D_Z2Z_cpu.f90"
echo "::::::::::::::"
diff testing_FFT_3D_Z2Z_cpu.f90 $DIFF_PATH
echo "::::::::::::::"
echo "testing_NumDiff_1D_cpu.f90"
echo "::::::::::::::"
diff testing_NumDiff_1D_cpu.f90 $DIFF_PATH
echo "::::::::::::::"
echo "testing_NumGrad_2D_cpu.f90"
echo "::::::::::::::"
diff testing_NumGrad_2D_cpu.f90 $DIFF_PATH
echo "::::::::::::::"
echo "testing_RowColMajCpp.cpp"
echo "::::::::::::::"
diff testing_RowColMajCpp.cpp $DIFF_PATH
echo "::::::::::::::"
echo "testing_RowColMajF90.f90"
echo "::::::::::::::"
diff testing_RowColMajF90.f90 $DIFF_PATH
echo "::::::::::::::"
echo "testing_sgetri_cpu.f90"
echo "::::::::::::::"
diff testing_sgetri_cpu.f90 $DIFF_PATH
echo "::::::::::::::"
echo "testing_SMatmult_gpu.f90"
echo "::::::::::::::"
diff testing_SMatmult_gpu.f90 $DIFF_PATH
echo "::::::::::::::"
echo "testing_summing.f90"
echo "::::::::::::::"
diff testing_summing.f90 $DIFF_PATH
echo "::::::::::::::"
echo "testing_systemQuery_cpu.f90"
echo "::::::::::::::"
diff testing_systemQuery_cpu.f90 $DIFF_PATH
echo "::::::::::::::"
echo "testing_yaml.f90"
echo "::::::::::::::"
diff testing_yaml.f90 $DIFF_PATH
echo "::::::::::::::"
echo "testing_zgetri_cpu.f90"
echo "::::::::::::::"
diff testing_zgetri_cpu.f90 $DIFF_PATH
echo "::::::::::::::"
echo "timer_tester_c.f90"
echo "::::::::::::::"
diff timer_tester_c.f90 $DIFF_PATH
echo "::::::::::::::"
echo "timer_tester.f90"
echo "::::::::::::::"
diff timer_tester.f90 $DIFF_PATH
echo "::::::::::::::"
echo "timer_tester.f90 with $DIFF_PATH_TIMER"
echo "::::::::::::::"
diff timer_tester_c.f90 $DIFF_PATH_TIMER
echo "::::::::::::::"
echo "timer_tester.f90 with $DIFF_PATH_TIMER"
echo "::::::::::::::"
diff timer_tester.f90   $DIFF_PATH_TIMER
