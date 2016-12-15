######################## Fortran Code ########################
./f90g95_local.csh ./test_code/gpu/cuda/testing_corr_gpu-v2.f90;
./f90g95_local.csh ./test_code/gpu/cuda/testing_gen_polar_coords.f90;
./f90g95_local.csh ./test_code/gpu/cuda/testing_MultiGPU_gpu.f90;
./f90g95_local.csh ./test_code/gpu/cuda/testing_Hadmr_3D_gpu.f90;
./f90g95_local.csh ./test_code/gpu/cuda/testing_GnCrA_3D_gpu.f90;
./f90g95_local.csh ./test_code/gpu/cuda/testing_DeepLearning_f90.f90;

./f90g95_local.csh ./test_code/gpu/cuda/testing_2Dcarte_corr_gpu.f90;
./f90g95_local.csh ./test_code/gpu/cuda/testing_2DshCar_corr_gpu.f90;

./f90g95_local.csh ./test_code/gpu/cuda/testing_eglossary.f90;
./f90g95_local.csh ./test_code/gpu/cuda/testing_int2char_f90.f90;
./f90g95_local.csh ./test_code/gpu/cuda/testing_readWriteImg.f90;
./f90g95_local.csh ./test_code/gpu/cuda/testing_DataParallelIO.f90;
############################## C++ Code ########################
./gccg++_local.csh ./test_code/gpu/cuda/testing_DeepLearning_cpp.cpp;
./gccg++_local.csh ./test_code/gpu/cuda/testing_i2c_indx_cpp.cpp;
./gccg++_local.csh ./test_code/gpu/cuda/testing_i2c_posi_cpp.cpp;
###################### CUDA code with nvcc #####################
nvcc -I /usr/local/cuda-7.0/include/ my_cufft_1D_C2C.cu -L /usr/local/cuda-7.0/lib64/ -lcufft  -o mccufft_1D_C2C ;
nvcc -I /usr/local/cuda-7.0/include/ my_cufft_2D_Z2Z.cu -L /usr/local/cuda-7.0/lib64/ -lcufft  -o mccufft_2D_Z2Z;
