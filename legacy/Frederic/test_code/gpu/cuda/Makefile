##############################################
# settings
#
# COMPILER = gfortran-4.9
# CFLAGS   = -ffree-form -cpp -O3 -fno-second-underscore -fopenmp -DBENCH
# LINK     = -lblas -llapack -lstdc++ -lpthread -lm
#
# .SUFFIXES: .f90 .o
# .f90.o:
#        $(COMPILER) $(CFLAGS) -c $<
#
##############################################
# Massive settings
# module load fftw/3.3.4-gcc ;
# module load gcc/4.9.1 ;
# module load lapack/3.4.2 ;
# module load openmpi ;
#
COMPILER_GFORTRAN = gfortran-4.9
COMPILER_GCC = g++-4.9
FFTW_LIB = /usr/lib/x86_64-linux-gnu

MAGMADIR=/opt/Devel_tools/magma-1.6.1
MAGMADIR_LIB=$(MAGMADIR)/lib
CUDADIR=/usr/local/cuda
CUDAINC=$(CUDADIR)/include
SHARE_LIB=/opt/Devel_tools/NVIDIA_GPU_Computing_SDK/shared/lib/linux

CFLAGS   = -ffree-form -cpp -fPIC -fno-second-underscore -O3 -fopenmp -DLINUX -DOPENMP -DBENCH -DCUDA -DMAGMA
CCFLAGS  = -cpp -O3 -fopenmp -DLINUX -DOPENMP -DBENCH -DCUDA -DMAGMA
LINK     = -L $(FFTW_LIB) -lfftw3                   \
           -L $(FFTW_LIB) -lfftw3f                  \
           -L $(FFTW_LIB) -lfftw3_threads           \
           -lblas -llapack -lstdc++ -lpthread -lrt  \
           -L $(CUDADIR)/lib64                      \
           -lcuda -lcublas -lcudart -lcufft -lcudnn \
           -lmagma -L $(MAGMADIR_LIB) -lm
LINK_CPP =

SIMPLE_SOURCE=../../..
XC_DIR=$(SIMPLE_SOURCE)/test_code/gpu/cuda
COMDIR = $(SIMPLE_SOURCE)/obj/GFORTgpu
OBJDIR = $(SIMPLE_SOURCE)/obj/GFORTgpu
SIMPLE_INC = $(SIMPLE_SOURCE)/include/simple

.SUFFIXES: .f90 .cpp .o 
.f90.o:
	$(COMPILER_GFORTRAN) $(CFLAGS) -I $(OBJDIR) -c $<
.cpp.o:
	$(COMPILER_GCC) $(CCFLAGS) -I $(OBJDIR) -I $(CUDAINC) \
        -I $(SIMPLE_INC) -c $<

COMMON = $(OBJDIR)/simple_defs.o \
         $(OBJDIR)/matrixGetter.o \
         $(OBJDIR)/simple_testfunction.o \
         $(OBJDIR)/greeting_version.o \
         $(OBJDIR)/simple_timing.o \
         $(OBJDIR)/timming.o \
         $(OBJDIR)/timming_c.o \
         $(OBJDIR)/get_cpu_time_c.o \
         $(OBJDIR)/timestamp.o \
         $(OBJDIR)/initialu.o \
         $(OBJDIR)/su3random.o \
         $(OBJDIR)/su2random.o \
         $(OBJDIR)/simple_sorting.o \
         $(OBJDIR)/simple_textHandler.o \
         $(OBJDIR)/get_systemQuery_cpu.o \
         $(OBJDIR)/get_systemQuery_cpu_c.o \
         $(OBJDIR)/simple_systemQuery_cpu.o \
         $(OBJDIR)/simple_file_utils.o \
         $(OBJDIR)/simple_DataParallel_IO.o \
         $(OBJDIR)/simple_file_defs.o \
         $(OBJDIR)/simple_resources_defs.o \
         $(OBJDIR)/simple_fftw3.o \
         $(OBJDIR)/fft123D_cpu.o \
         $(OBJDIR)/get_fft123D_cpu.o \
         $(OBJDIR)/get_fft123D_cpu_c.o \
         $(OBJDIR)/fft123D_gpu.o \
         $(OBJDIR)/get_fft123D_gpu.o \
         $(OBJDIR)/get_fft123D_gpu_c.o \
         $(OBJDIR)/simple_deviceQuery_gpu.o \
         $(OBJDIR)/get_deviceQuery_gpu.o \
         $(OBJDIR)/get_deviceQuery_gpu_c.o \
         $(OBJDIR)/simple_cuda.o \
         $(OBJDIR)/simple_cuda_defs.o \
         $(OBJDIR)/fortran.o \
         $(OBJDIR)/simple_fortran.o \
         $(OBJDIR)/polarft_corr_Helpr_gpu.o \
         $(OBJDIR)/polarft_corr_Helpr_gpu_c.o \
         $(OBJDIR)/Global_polarft.o \
         $(OBJDIR)/simple_pfts_Sizes.o \
         $(OBJDIR)/simple_mesh3D.o \
         $(OBJDIR)/simple_mesh1D.o \
         $(OBJDIR)/simple_mesh1DV.o \
         $(OBJDIR)/simple_img_2D_cart_Sizes.o \
         $(OBJDIR)/invert_gpu_utils.o \
         $(OBJDIR)/invert_cpu_utils.o \
         $(OBJDIR)/magma_invert_gpu-v2.o \
         $(OBJDIR)/magma_zgetri_blocked_gpu-v2.o \
         $(OBJDIR)/magma_dgetri_blocked_gpu-v2.o \
         $(OBJDIR)/magma_get_getri_nb_gpu.o \
         $(OBJDIR)/simple_math_gpu.o \
         $(OBJDIR)/simple_math_gpu_c.o \
         $(OBJDIR)/matmul_gpu_utils.o \
         $(OBJDIR)/matvec_gpu_utils.o \
         $(OBJDIR)/simple_ErrorHandler.o \
         $(OBJDIR)/matmul_cuda_transfers.o \
         $(OBJDIR)/matmul_cuda_device.o \
         $(OBJDIR)/matvec_cuda_transfers.o \
         $(OBJDIR)/matvec_cuda_device.o \
         $(OBJDIR)/simple_testfuns.o \
         $(OBJDIR)/gen_polar_coords_gpu.o \
         $(OBJDIR)/simple_math.o \
         $(OBJDIR)/simple_jiffys.o \
         $(OBJDIR)/simple_rnd.o \
         $(OBJDIR)/zz2dgemm_ElmtWs_tesla_N_N.o \
         $(OBJDIR)/zz2dgemm_ElmtWs_tesla_T_N.o \
         $(OBJDIR)/zz2dgemm_ElmtWs_tesla_N_T.o \
         $(OBJDIR)/zz2dgemm_ElmtWs_tesla_T_T.o \
         $(OBJDIR)/zz2dgemm_ElmtWs_tesla_gpu.o \
         $(OBJDIR)/zz2dgemm_ElmtWs_tesla_sumsq_N_N.o \
         $(OBJDIR)/zz2dgemm_ElmtWs_tesla_sumsq_T_T.o \
         $(OBJDIR)/zz2dgemm_ElmtWs_tesla_sumsq_gpu.o \
         $(OBJDIR)/simple_polarft_gpu_c.o \
         $(OBJDIR)/simple_polarft_corr_gpu_c.o \
         $(OBJDIR)/simple_polarft_corr_gpu.o \
         $(OBJDIR)/polarft_corr_Hadmr_gpu.o \
         $(OBJDIR)/polarft_corr_Hadmr_N_N.o \
         $(OBJDIR)/polarft_corr_Hadmr_F_N.o \
         $(OBJDIR)/polarft_corr_Hadmr_P_N.o \
         $(OBJDIR)/polarft_corr_Hadmr_X_N.o \
         $(OBJDIR)/polarft_gencorrAll_gpu.o \
         $(OBJDIR)/polarft_gencorrAll_Z_N.o \
         $(OBJDIR)/polarft_multi-GPUs_gpu.o \
         $(OBJDIR)/polarft_multi-GPUs_Z_M.o \
         $(OBJDIR)/simple_ft_expanded.o \
         $(OBJDIR)/carte2D_ftExt-corr_gpu.o \
         $(OBJDIR)/carte2D_ftExt-corr_C_N.o \
         $(OBJDIR)/carte2D_ftExt-corr_C_F.o \
         $(OBJDIR)/simple_image.o \
         $(OBJDIR)/simple_ftiter.o \
         $(OBJDIR)/gnufor2.o \
         $(OBJDIR)/simple_stat.o \
         $(OBJDIR)/simple_winfuns.o \
         $(OBJDIR)/simple_online_var.o \
         $(OBJDIR)/simple_imgfile.o \
         $(OBJDIR)/simple_imghead.o \
         $(OBJDIR)/simple_imgheadrec.o \
         $(OBJDIR)/simple_deepLearning.o \
         $(OBJDIR)/simple_deepLearning_defs.o \
         $(OBJDIR)/simple_cudnn_fortran.o \
         $(OBJDIR)/simple_random.o \
         $(OBJDIR)/strlen.o \
         $(OBJDIR)/simple_eglossary.o \
         $(OBJDIR)/simple_err_defs.o \
         $(OBJDIR)/simple_error_handling.o \
         $(OBJDIR)/simple_eglossary_lowlev.o \
         $(OBJDIR)/simple_yaml_output.o \
         $(OBJDIR)/simple_yaml_strings.o \
         $(OBJDIR)/simple_file_highlev.o \
         $(OBJDIR)/simple_dynamic_memory.o \
         $(OBJDIR)/simple_cmdline.o \
         $(OBJDIR)/simple_ctf.o \
         $(OBJDIR)/simple_oris.o \
         $(OBJDIR)/simple_params.o \
         $(OBJDIR)/simple_build.o \
         $(OBJDIR)/simple_args.o \
         $(OBJDIR)/simple_ori.o \
         $(OBJDIR)/simple_opt_spec.o \
         $(OBJDIR)/simple_simplex_opt.o \
         $(OBJDIR)/simple_softmax_weights.o \
         $(OBJDIR)/simple_ran_tabu.o \
         $(OBJDIR)/simple_eo_reconstructor.o \
         $(OBJDIR)/simple_ppca.o \
         $(OBJDIR)/simple_centre_clust.o \
         $(OBJDIR)/simple_spatial_median.o \
         $(OBJDIR)/simple_sym.o \
         $(OBJDIR)/simple_projector.o \
         $(OBJDIR)/simple_reconstructor.o \
         $(OBJDIR)/simple_comlin.o \
         $(OBJDIR)/simple_hash.o \
         $(OBJDIR)/simple_bfgs_opt.o \
         $(OBJDIR)/simple_gridding.o \
         $(OBJDIR)/simple_polarft.o \
         $(OBJDIR)/simple_strings.o \
         $(OBJDIR)/simple_opt_subs.o \
         $(OBJDIR)/simple_optimizer.o \
         $(OBJDIR)/simple_memory_profiling.o \
         $(OBJDIR)/simple_AVratios.o \
         $(OBJDIR)/simple_estimate_ssnr.o \
         $(OBJDIR)/simple_cmd_dict.o \
         $(OBJDIR)/simple_convergence.o \
         $(OBJDIR)/simple_chash.o \
         $(OBJDIR)/simple_map_reduce.o \
         $(OBJDIR)/polarft_krnl-Opti_gpu.o \
         $(OBJDIR)/polarft_krnl-Opti_O_N.o


COMMON_C = $(OBJDIR)/timming.o \
           $(OBJDIR)/timming_c.o

#Prarralle IO
MAIN = testing_corr_gpu-v2.o
#Inverstion Lapacks
MAIN_GENPOLAR = testing_gen_polar_coords.o
#Polar FT testers
MAIN_MULTIGPU = testing_MultiGPU_gpu.o
MAIN_HADMR_3D = testing_Hadmr_3D_gpu.o
MAIN_GNCRA_3D = testing_GnCrA_3D_gpu.o
#motion correction testers
MAIN_2DCARTE_CORR = testing_2Dcarte_corr_gpu.o
MAIN_2DSHCAR_CORR = testing_2DshCar_corr_gpu.o
#DeepLearning
MAIN_DEEPLEARNING_F90 = testing_DeepLearning_f90.o
MAIN_DEEPLEARNING_CPP = testing_DeepLearning_cpp.o
#Prarralle IO
MAIN_DATAPARALLELIO = testing_DataParallelIO.o
#eglossary
MAIN_EGLOSSARY  = testing_eglossary.o
MAIN_READWRITEIMG  = testing_readWriteImg.o
#characters and integer converters
MAIN_INT2CHAR_F90 = testing_int2char_f90.o
MAIN_I2C_INDX_CPP = testing_i2c_indx_cpp.o
MAIN_I2C_POSI_CPP = testing_i2c_posi_cpp.o
#YAML
MAIN_YAML = testing_yaml.o

#Executables
XCS = mccorr_gpu-v2
XCS_GENPOLAR = mcgen_polar_coords
#Polar FT testers
XCS_MULTIGPU = mcMultiGPU_gpu
XCS_HADMR_3D = mcHadmr_3D_gpu
XCS_GNCRA_3D = mcGnCrA_3D_gpu
#motion correction testers
XCS_2DCARTE_CORR = mc2Dcarte_corr_gpu
XCS_2DSHCAR_CORR = mc2DshCar_corr_gpu
#DeepLearning
XCS_DEEPLEARNING_F90 = mcDeepLearning_f90
XCS_DEEPLEARNING_CPP = mcDeepLearning_cpp
#parallel IO
XCS_DATAPARALLELIO = mcDataParallelIO
#eglossary
XCS_EGLOSSARY  = mceglossary
XCS_READWRITEIMG  = mcreadWriteImg
#characters and integer converters
XCS_INT2CHAR_F90 = mcint2char_f90
XCS_I2C_INDX_CPP = mci2c_indx_cpp
XCS_I2C_POSI_CPP = mci2c_posi_cpp
#YAML
XCS_YAML = mcyaml

#Compilations
all: $(XCS) $(XCS_GENPOLAR) $(XCS_MULTIGPU) $(XCS_HADMR_3D) $(XCS_GNCRA_3D) \
     $(XCS_2DCARTE_CORR) $(XCS_2DSHCAR_CORR) $(XCS_DEEPLEARNING_F90) \
     $(XCS_DEEPLEARNING_CPP) $(XCS_DATAPARALLELIO) $(XCS_READWRITEIMG) \
     $(XCS_INT2CHAR_F90) $(XCS_EGLOSSARY) $(XCS_I2C_INDX_CPP) \
     $(XCS_I2C_POSI_CPP) $(XCS_YAML)

mccorr_gpu-v2: $(COMMON) $(MAIN)
	$(COMPILER_GFORTRAN) $(CFLAGS) -o $(XC_DIR)/$(XCS) \
        $(COMMON) $(MAIN) $(LINK)

mcgen_polar_coords: $(COMMON) $(MAIN_GENPOLAR)
	$(COMPILER_GFORTRAN) $(CFLAGS) -o $(XC_DIR)/$(XCS_GENPOLAR) \
        $(COMMON) $(MAIN_GENPOLAR) $(LINK)

mcMultiGPU_gpu: $(COMMON) $(MAIN_MULTIGPU)
	$(COMPILER_GFORTRAN) $(CFLAGS) -o $(XC_DIR)/$(XCS_MULTIGPU) \
        $(COMMON) $(MAIN_MULTIGPU) $(LINK)

mcHadmr_3D_gpu: $(COMMON) $(MAIN_HADMR_3D)
	$(COMPILER_GFORTRAN) $(CFLAGS) -o $(XC_DIR)/$(XCS_HADMR_3D) \
        $(COMMON) $(MAIN_HADMR_3D) $(LINK)

mcGnCrA_3D_gpu: $(COMMON) $(MAIN_GNCRA_3D)
	$(COMPILER_GFORTRAN) $(CFLAGS) -o $(XC_DIR)/$(XCS_GNCRA_3D) \
        $(COMMON) $(MAIN_GNCRA_3D) $(LINK)

mc2Dcarte_corr_gpu: $(COMMON) $(MAIN_2DCARTE_CORR)
	$(COMPILER_GFORTRAN) $(CFLAGS) -o $(XC_DIR)/$(XCS_2DCARTE_CORR) \
        $(COMMON) $(MAIN_2DCARTE_CORR) $(LINK)

mc2DshCar_corr_gpu: $(COMMON) $(MAIN_2DSHCAR_CORR)
	$(COMPILER_GFORTRAN) $(CFLAGS) -o $(XC_DIR)/$(XCS_2DSHCAR_CORR) \
        $(COMMON) $(MAIN_2DSHCAR_CORR) $(LINK)

mcDeepLearning_f90: $(COMMON) $(MAIN_DEEPLEARNING_F90)
	$(COMPILER_GFORTRAN) $(CFLAGS) -o $(XC_DIR)/$(XCS_DEEPLEARNING_F90) \
        $(COMMON) $(MAIN_DEEPLEARNING_F90) $(LINK)

mcDeepLearning_cpp: $(MAIN_DEEPLEARNING_CPP)
	$(COMPILER_GCC) $(CCFLAGS) -o $(XC_DIR)/$(XCS_DEEPLEARNING_CPP) \
        $(MAIN_DEEPLEARNING_CPP) $(LINK)

mcDataParallelIO : $(COMMON) $(MAIN_DATAPARALLELIO)
	$(COMPILER_GFORTRAN) $(CFLAGS) -o $(XC_DIR)/$(XCS_DATAPARALLELIO) \
        $(COMMON) $(MAIN_DATAPARALLELIO) $(LINK)

mceglossary: $(COMMON) $(MAIN_EGLOSSARY)
	$(COMPILER_GFORTRAN) $(CFLAGS) -o $(XC_DIR)/$(XCS_EGLOSSARY) \
        $(COMMON) $(MAIN_EGLOSSARY) $(LINK)

mcreadWriteImg: $(COMMON) $(MAIN_READWRITEIMG)
	$(COMPILER_GFORTRAN) $(CFLAGS) -o $(XC_DIR)/$(XCS_READWRITEIMG) \
        $(COMMON) $(MAIN_READWRITEIMG) $(LINK)

mcint2char_f90: $(COMMON) $(MAIN_INT2CHAR_F90)
	$(COMPILER_GFORTRAN) $(CFLAGS) -o $(XC_DIR)/$(XCS_INT2CHAR_F90) \
        $(COMMON) $(MAIN_INT2CHAR_F90) $(LINK)

mci2c_indx_cpp: $(COMMON_C) $(MAIN_I2C_INDX_CPP)
	$(COMPILER_GCC) $(CCFLAGS) -o $(XC_DIR)/$(XCS_I2C_INDX_CPP) \
        $(COMMON_C) $(MAIN_I2C_INDX_CPP) $(LINK)

mci2c_posi_cpp: $(COMMON_C) $(MAIN_I2C_POSI_CPP)
	$(COMPILER_GCC) $(CCFLAGS) -o $(XC_DIR)/$(XCS_I2C_POSI_CPP) \
        $(COMMON_C) $(MAIN_I2C_POSI_CPP) $(LINK)

mcyaml: $(COMMON) $(MAIN_YAML)
	$(COMPILER_GFORTRAN) $(CFLAGS) -o $(XC_DIR)/$(XCS_YAML) \
        $(COMMON) $(MAIN_YAML) $(LINK)

clean :
	rm -f *.o *.d *~* *#* *.pyc *.asc dwt.bmp orig.bmp recon.bmp *.yaml  \
 $(XC_DIR)/$(XCS) $(XC_DIR)/$(XCS_GENPOLAR) \
 $(XC_DIR)/$(XCS_MULTIGPU) $(XC_DIR)/$(XCS_HADMR_3D) \
 $(XC_DIR)/$(XCS_GNCRA_3D) $(XC_DIR)/$(XCS_2DCARTE_CORR) \
 $(XC_DIR)/$(XCS_2DSHCAR_CORR) $(XC_DIR)/$(XCS_DEEPLEARNING_F90) \
 $(XC_DIR)/$(XCS_DEEPLEARNING_CPP) $(XC_DIR)/$(XCS_DATAPARALLELIO) \
 $(XC_DIR)/$(XCS_EGLOSSARY) $(XC_DIR)/$(XCS_READWRITEIMG) \
 $(XC_DIR)/$(XCS_INT2CHAR_F90) $(XC_DIR)/$(XCS_I2C_INDX_CPP) \
 $(XC_DIR)/$(XCS_I2C_POSI_CPP) $(XC_DIR)/$(XCS_YAML)

cleanall:
######################## Fortran Code ########################
#./f90g95_local.csh ./test_code/cpu/testing_DataParallelIO.f90;

#./f90g95_local.csh ./test_code/cpu/testing_sgetri_cpu.f90;
#./f90g95_local.csh ./test_code/cpu/testing_dgetri_cpu.f90;
#./f90g95_local.csh ./test_code/cpu/testing_zgetri_cpu.f90;

#./f90g95_local.csh ./test_code/cpu/testing_FFT_3D_C2S_cpu.f90;
#./f90g95_local.csh ./test_code/cpu/testing_FFT_3D_S2C_cpu.f90;
#./f90g95_local.csh ./test_code/cpu/testing_FFT_3D_C2C_cpu.f90;
#./f90g95_local.csh ./test_code/cpu/testing_FFT_3D_D2Z_cpu.f90;
#./f90g95_local.csh ./test_code/cpu/testing_FFT_3D_Z2D_cpu.f90;
#./f90g95_local.csh ./test_code/cpu/testing_FFT_3D_Z2Z_cpu.f90;

#./f90g95_local.csh ./test_code/cpu/testing_FFT_2D_C2S_cpu.f90;
#./f90g95_local.csh ./test_code/cpu/testing_FFT_2D_S2C_cpu.f90;
#./f90g95_local.csh ./test_code/cpu/testing_FFT_2D_C2C_cpu.f90;
#./f90g95_local.csh ./test_code/cpu/testing_FFT_2D_D2Z_cpu.f90;
#./f90g95_local.csh ./test_code/cpu/testing_FFT_2D_Z2D_cpu.f90;
#./f90g95_local.csh ./test_code/cpu/testing_FFT_2D_Z2Z_cpu.f90;

#./f90g95_local.csh ./test_code/cpu/testing_FFT_1D_Z2Z_cpu.f90;
############################
#./f90g95_local.csh ./test_code/cpu/testing_systemQuery_cpu.f90;
#./f90g95_local.csh ./test_code/cpu/testing_yaml.f90;

#./f90g95_local.csh ./test_code/cpu/testing_summing.f90;
#./f90g95_local.csh ./test_code/cpu/testing_RowColMajF90.f90;
#./f90g95_local.csh ./test_code/cpu/testing_congruence.f90;
#./f90g95_local.csh ./test_code/cpu/testing_BFGSNum_2D_cpu.f90;
#./f90g95_local.csh ./test_code/cpu/testing_corr_cpu.f90;
#./f90g95_local.csh ./test_code/cpu/testing_corr_gpu.f90;
#./f90g95_local.csh ./test_code/cpu/testing_cshift.f90;
#./f90g95_local.csh ./test_code/cpu/testing_NumDiff_1D_cpu.f90;
#./f90g95_local.csh ./test_code/cpu/testing_NumGrad_2D_cpu.f90;
#./f90g95_local.csh ./test_code/cpu/timer_tester_c.f90;
#./f90g95_local.csh ./test_code/cpu/timer_tester.f90;

#./f90g95_local.csh ./test_code/gpu/cuda/testing_SMatmult_cpu.f90;
#./f90g95_local.csh ./test_code/gpu/cuda/testing_SMatvect_cpu.f90;

############################## C++ Code ########################
#./gccg++_local.csh ./test_code/cpu/testing_RowColMajCpp.cpp;




