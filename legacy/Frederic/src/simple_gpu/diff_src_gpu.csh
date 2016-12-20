#!/bin/bash
set MY_HOME_SIDE = /home/frederic/SimpleSideBrch 

set DIFF_PATH = $MY_HOME_SIDE/11Apr16/Simple_Restruct.projet/src/simple_gpu

echo "::::::::::::::"
echo ""
echo "::::::::::::::"
diff diff_src_gpu.csh $DIFF_PATH
echo "::::::::::::::"
echo ""
echo "::::::::::::::"
diff gen_polar_coords_gpu.cu $DIFF_PATH
echo "::::::::::::::"
echo ""
echo "::::::::::::::"
diff invert_gpu_utils.f90 $DIFF_PATH
echo "::::::::::::::"
echo ""
echo "::::::::::::::"
diff magma_dgetri_blocked_gpu-v2.cpp $DIFF_PATH
echo "::::::::::::::"
echo ""
echo "::::::::::::::"
diff magma_get_getri_nb_gpu.cpp $DIFF_PATH
echo "::::::::::::::"
echo ""
echo "::::::::::::::"
diff magma_invert_gpu-v2.cpp $DIFF_PATH
echo "::::::::::::::"
echo ""
echo "::::::::::::::"
diff magma_zgetri_blocked_gpu-v2.cpp $DIFF_PATH
echo "::::::::::::::"
echo ""
echo "::::::::::::::"
diff Makefile_target $DIFF_PATH
echo "::::::::::::::"
echo ""
echo "::::::::::::::"
diff matmul_cuda_device.c $DIFF_PATH
echo "::::::::::::::"
echo ""
echo "::::::::::::::"
diff matmul_cuda_transfers.c $DIFF_PATH
echo "::::::::::::::"
echo ""
echo "::::::::::::::"
diff matmul_gpu_utils.f90 $DIFF_PATH
echo "::::::::::::::"
echo ""
echo "::::::::::::::"
diff simple_math_gpu_c.cpp $DIFF_PATH
echo "::::::::::::::"
echo ""
echo "::::::::::::::"
diff simple_math_gpu.f90 $DIFF_PATH
echo "::::::::::::::"
echo ""
echo "::::::::::::::"
diff simple_polarft_gpu.f90 $DIFF_PATH
echo "::::::::::::::"
echo ""
echo "::::::::::::::"
diff matvec_cuda_device.c $DIFF_PATH
echo "::::::::::::::"
echo ""
echo "::::::::::::::"
diff matvec_cuda_transfers.c $DIFF_PATH
echo "::::::::::::::"
echo ""
echo "::::::::::::::"
diff matvec_gpu_utils.f90 $DIFF_PATH
echo "::::::::::::::"
echo "simple_cuda.f90"
echo "::::::::::::::"
diff simple_cuda.f90 $DIFF_PATH
echo "::::::::::::::"
echo "simple_mesh1D.cpp"
echo "::::::::::::::"
diff simple_mesh1D.cpp $DIFF_PATH
echo "::::::::::::::"
echo "simple_mesh1DV.cpp"
echo "::::::::::::::"
diff simple_mesh1DV.cpp $DIFF_PATH
echo "::::::::::::::"
echo "simple_mesh3D.cpp"
echo "::::::::::::::"
diff simple_mesh3D.cpp $DIFF_PATH
echo "::::::::::::::"
echo "simple_pfts_Sizes.cpp"
echo "::::::::::::::"
diff simple_pfts_Sizes.cpp $DIFF_PATH
echo "::::::::::::::"
echo "Global_polarft.cpp"
echo "::::::::::::::"
diff Global_polarft.cpp $DIFF_PATH
echo "::::::::::::::"
echo "Global_polarft.cpp"
echo "::::::::::::::"
diff simple_img_2D_cart_Sizes.cpp $DIFF_PATH
echo "::::::::::::::"
echo "carte2D_ftExt-corr_gpu.cu"
echo "::::::::::::::"
carte2D_ftExt-corr_gpu.cu $DIFF_PATH
echo "::::::::::::::"
echo "carte2D_ftExt-corr_C_N.cu"
echo "::::::::::::::"
carte2D_ftExt-corr_C_N.cu $DIFF_PATH
echo "::::::::::::::"
echo "carte2D_ftExt-corr_C_F.cu"
echo "::::::::::::::"
carte2D_ftExt-corr_C_F.cu $DIFF_PATH

