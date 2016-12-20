nvcc -I /usr/local/cuda-7.0/include/ my_cufft_1D_C2C.cu -L /usr/local/cuda-7.0/lib64/ -lcufft  -o mccufft_1D_C2C ;
nvcc -I /usr/local/cuda-7.0/include/ my_cufft_2D_Z2Z.cu -L /usr/local/cuda-7.0/lib64/ -lcufft  -o mccufft_2D_Z2Z 
