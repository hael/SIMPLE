#!/bin/bash
set MY_HOME_SIDE = /home/frederic/SimpleSideBrch 

set DIFF_PATH = $MY_HOME_SIDE/21Mar16/Simple_Restruct.projet/test_code/gpu/cuda

echo "::::::::::::::"
echo "clean_All_Testing_gpu.csh"
echo "::::::::::::::"
diff clean_All_Testing_gpu.csh $DIFF_PATH;
echo "::::::::::::::"
echo "compile_All_Testing_gpu.csh"
echo "::::::::::::::"
diff compile_All_Testing_gpu.csh $DIFF_PATH
echo "::::::::::::::"
echo "compile_cufft_test.csh"
echo "::::::::::::::"
diff compile_cufft_test.csh $DIFF_PATH
echo "::::::::::::::"
echo "diff_test_TestCuda.csh"
echo "::::::::::::::"
diff diff_test_TestCuda.csh $DIFF_PATH
echo "::::::::::::::"
echo "f90g95_local.csh"
echo "::::::::::::::"
diff f90g95_local.csh $DIFF_PATH
echo "::::::::::::::"
echo "f90g95_zgetri_local"
echo "::::::::::::::"
diff f90g95_zgetri_local $DIFF_PATH
echo "::::::::::::::"
echo "functiontest_2d.m"
echo "::::::::::::::"
diff functiontest_2d.m $DIFF_PATH
echo "::::::::::::::"
echo "Makefile_target"
echo "::::::::::::::"
diff Makefile_target $DIFF_PATH
echo "::::::::::::::"
echo "my_cufft_1D_C2C.cu"
echo "::::::::::::::"
diff my_cufft_1D_C2C.cu $DIFF_PATH
echo "::::::::::::::"
echo "my_cufft_2D_Z2Z.cu"
echo "::::::::::::::"
diff my_cufft_2D_Z2Z.cu $DIFF_PATH
echo "::::::::::::::"
echo "polarft_corr_Hadmr_tesla_F_N.cu"
echo "::::::::::::::"
diff polarft_corr_Hadmr_tesla_F_N.cu $DIFF_PATH
echo "::::::::::::::"
echo "polarft_corr_Hadmr_tesla_gpu.cu"
echo "::::::::::::::"
diff polarft_corr_Hadmr_tesla_gpu.cu $DIFF_PATH
echo "::::::::::::::"
echo "polarft_corr_Hadmr_tesla_N_N.cu"
echo "::::::::::::::"
diff polarft_corr_Hadmr_tesla_N_N.cu $DIFF_PATH
echo "::::::::::::::"
echo "polarft_corr_Hadmr_tesla_P_N.cu"
echo "::::::::::::::"
diff polarft_corr_Hadmr_tesla_P_N.cu $DIFF_PATH
echo "::::::::::::::"
echo "polarft_corr_Hadmr_tesla_X_N.cu"
echo "::::::::::::::"
diff polarft_corr_Hadmr_tesla_X_N.cu $DIFF_PATH
echo "::::::::::::::"
echo "polarft_corr_Helpr_tesla_gpu.cu"
echo "::::::::::::::"
diff polarft_corr_Helpr_tesla_gpu.cu $DIFF_PATH
echo "::::::::::::::"
echo "polarft_gencorrAll_tesla_gpu.cu"
echo "::::::::::::::"
diff polarft_gencorrAll_tesla_gpu.cu $DIFF_PATH
echo "::::::::::::::"
echo "polarft_gencorrAll_tesla_Z_N.cu"
echo "::::::::::::::"
diff polarft_gencorrAll_tesla_Z_N.cu $DIFF_PATH
echo "::::::::::::::"
echo "polarft_multi-GPUs_tesla_gpu.cu"
echo "::::::::::::::"
diff polarft_multi-GPUs_tesla_gpu.cu $DIFF_PATH
echo "::::::::::::::"
echo "polarft_multi-GPUs_tesla_Z_M.cu"
echo "::::::::::::::"
diff polarft_multi-GPUs_tesla_Z_M.cu $DIFF_PATH
echo "::::::::::::::"
echo "run_All_Testing_gpu.csh"
echo "::::::::::::::"
diff run_All_Testing_gpu.csh $DIFF_PATH
echo "::::::::::::::"
echo "simple_ErrorHandler.cpp"
echo "::::::::::::::"
diff simple_ErrorHandler.cpp $DIFF_PATH
echo "::::::::::::::"
echo "simple_polarft_corr_gpu_c.c"
echo "::::::::::::::"
diff simple_polarft_corr_gpu_c.c $DIFF_PATH
echo "::::::::::::::"
echo "simple_polarft_corr_gpu.cpp"
echo "::::::::::::::"
diff simple_polarft_corr_gpu.cpp $DIFF_PATH
echo "::::::::::::::"
echo "simple_polarft_gpu_c.cpp"
echo "::::::::::::::"
diff simple_polarft_gpu_c.cpp $DIFF_PATH
echo "::::::::::::::"
echo "simple_testfunction.f90"
echo "::::::::::::::"
diff simple_testfunction.f90 $DIFF_PATH
echo "::::::::::::::"
echo "testing_corr_gpu.f90"
echo "::::::::::::::"
diff testing_corr_gpu.f90 $DIFF_PATH
echo "::::::::::::::"
echo "testing_corr_gpu-v2.f90"
echo "::::::::::::::"
diff testing_corr_gpu-v2.f90 $DIFF_PATH
echo "::::::::::::::"
echo "testing_cuFFT_1D_Z2Z_gpu.f90"
echo "::::::::::::::"
diff testing_cuFFT_1D_Z2Z_gpu.f90 $DIFF_PATH
echo "::::::::::::::"
echo "testing_cuFFT_2D_C2C_gpu.f90"
echo "::::::::::::::"
diff testing_cuFFT_2D_C2C_gpu.f90 $DIFF_PATH
echo "::::::::::::::"
echo "testing_cuFFT_2D_C2S_gpu.f90"
echo "::::::::::::::"
diff testing_cuFFT_2D_C2S_gpu.f90 $DIFF_PATH
echo "::::::::::::::"
echo "testing_cuFFT_2D_D2Z_gpu.f90"
echo "::::::::::::::"
diff testing_cuFFT_2D_D2Z_gpu.f90 $DIFF_PATH
echo "::::::::::::::"
echo "testing_cuFFT_2D_S2C_gpu.f90"
echo "::::::::::::::"
diff testing_cuFFT_2D_S2C_gpu.f90 $DIFF_PATH
echo "::::::::::::::"
echo "testing_cuFFT_2D_Z2D_gpu.f90"
echo "::::::::::::::"
diff testing_cuFFT_2D_Z2D_gpu.f90 $DIFF_PATH
echo "::::::::::::::"
echo "testing_cuFFT_2D_Z2Z_gpu.f90"
echo "::::::::::::::"
diff testing_cuFFT_2D_Z2Z_gpu.f90 $DIFF_PATH
echo "::::::::::::::"
echo "testing_cuFFT_3D_C2C_gpu.f90"
echo "::::::::::::::"
diff testing_cuFFT_3D_C2C_gpu.f90 $DIFF_PATH
echo "::::::::::::::"
echo "testing_cuFFT_3D_C2S_gpu.f90"
echo "::::::::::::::"
diff testing_cuFFT_3D_C2S_gpu.f90 $DIFF_PATH
echo "::::::::::::::"
echo "testing_cuFFT_3D_D2Z_gpu.f9"
echo "::::::::::::::"
diff testing_cuFFT_3D_D2Z_gpu.f90 $DIFF_PATH
echo "::::::::::::::"
echo "testing_cuFFT_3D_S2C_gpu.f90"
echo "::::::::::::::"
diff testing_cuFFT_3D_S2C_gpu.f90 $DIFF_PATH
echo "::::::::::::::"
echo "testing_cuFFT_3D_Z2D_gpu.f90"
echo "::::::::::::::"
diff testing_cuFFT_3D_Z2D_gpu.f90 $DIFF_PATH
echo "::::::::::::::"
echo "testing_cuFFT_3D_Z2Z_gpu.f90"
echo "::::::::::::::"
diff testing_cuFFT_3D_Z2Z_gpu.f90 $DIFF_PATH
echo "::::::::::::::"
echo "testing_deviceQuery_gpu.f90"
echo "::::::::::::::"
diff testing_deviceQuery_gpu.f90 $DIFF_PATH
echo "::::::::::::::"
echo "testing_dgetri_gpu.f90"
echo "::::::::::::::"
diff testing_dgetri_gpu.f90 $DIFF_PATH
echo "::::::::::::::"
echo "testing_discrepancy.f90"
echo "::::::::::::::"
diff testing_discrepancy.f90 $DIFF_PATH
echo "::::::::::::::"
echo "testing_gen_polar_coords.f90"
echo "::::::::::::::"
diff testing_gen_polar_coords.f90 $DIFF_PATH
echo "::::::::::::::"
echo "testing_GnCrA_3D_gpu.f90"
echo "::::::::::::::"
diff testing_GnCrA_3D_gpu.f90 $DIFF_PATH
echo "::::::::::::::"
echo "testing_Hadmr_3D_gpu.f90"
echo "::::::::::::::"
diff testing_Hadmr_3D_gpu.f90 $DIFF_PATH
echo "::::::::::::::"
echo "testing_MultiGPU_gpu.f90"
echo "::::::::::::::"
diff testing_MultiGPU_gpu.f90 $DIFF_PATH
echo "::::::::::::::"
echo "testing_zgetri_gpu.f90"
echo "::::::::::::::"
diff testing_zgetri_gpu.f90 $DIFF_PATH
echo "::::::::::::::"
echo "Times_magma_dgetri_block_gpu_v2.log"
echo "::::::::::::::"
diff Times_magma_dgetri_block_gpu_v2.log $DIFF_PATH
echo "::::::::::::::"
echo "Times_magma_zgetri_block_gpu_v2.log"
echo "::::::::::::::"
diff Times_magma_zgetri_block_gpu_v2.log $DIFF_PATH
echo "::::::::::::::"
echo "zz2dgemm_ElmtWs_tesla_gpu.cu"
echo "::::::::::::::"
diff zz2dgemm_ElmtWs_tesla_gpu.cu $DIFF_PATH
echo "::::::::::::::"
echo "zz2dgemm_ElmtWs_tesla_N_N.cu"
echo "::::::::::::::"
diff zz2dgemm_ElmtWs_tesla_N_N.cu $DIFF_PATH
echo "::::::::::::::"
echo "zz2dgemm_ElmtWs_tesla_N_T.cu"
echo "::::::::::::::"
diff zz2dgemm_ElmtWs_tesla_N_T.cu $DIFF_PATH
echo "::::::::::::::"
echo "zz2dgemm_ElmtWs_tesla_sumsq_gpu.cu"
echo "::::::::::::::"
diff zz2dgemm_ElmtWs_tesla_sumsq_gpu.cu $DIFF_PATH
echo "::::::::::::::"
echo "zz2dgemm_ElmtWs_tesla_sumsq_N_N.cu"
echo "::::::::::::::"
diff zz2dgemm_ElmtWs_tesla_sumsq_N_N.cu $DIFF_PATH
echo "::::::::::::::"
echo "zz2dgemm_ElmtWs_tesla_sumsq_T_T.cu"
echo "::::::::::::::"
diff zz2dgemm_ElmtWs_tesla_sumsq_T_T.cu $DIFF_PATH
echo "::::::::::::::"
echo "zz2dgemm_ElmtWs_tesla_T_N.c"
echo "::::::::::::::"
diff zz2dgemm_ElmtWs_tesla_T_N.cu $DIFF_PATH
echo "::::::::::::::"
echo "zz2dgemm_ElmtWs_tesla_T_T.cu"
echo "::::::::::::::"
diff zz2dgemm_ElmtWs_tesla_T_T.cu $DIFF_PATH
