FLNAME=/Users/hael/src/fortran/simple_wfrederic/Simple_Restruct.projet/production/simple_tests/simple_test_ncomps/simple_test_ncomps

GFORTRAN=/sw/bin/gfortran

SOURCE_DIR=/Users/hael/src/fortran/simple_wfrederic/Simple_Restruct.projet

MODDIR=$SOURCE_DIR/obj/GFORTgpu
OBJDIR=$SOURCE_DIR/obj/GFORTgpu

CUDADIR=

SHARE_LIB=/

OPTION="-fimplicit-none -fall-intrinsics -ffree-form -cpp -fpic -fno-second-underscore -g -Og -Wall -fbounds-check -DMACOSX -fopenmp      "

FFTW_LIB=/sw/lib

$GFORTRAN $OPTION $LDFLAGS -I $OBJDIR -I $MODDIR -o $FLNAME \
                                              $FLNAME.f90 \
                                              $OBJDIR/simple_params.o \
                                              $OBJDIR/simple_defs.o \
                                              $OBJDIR/simple_jiffys.o \
                                              $OBJDIR/simple_imghead.o \
                                              $OBJDIR/simple_imgheadrec.o \
                                              $OBJDIR/simple_ori.o \
                                              $OBJDIR/simple_hash.o \
                                              $OBJDIR/simple_strings.o \
                                              $OBJDIR/simple_rnd.o \
                                              $OBJDIR/simple_math.o \
                                              $OBJDIR/simple_simplex_opt.o \
                                              $OBJDIR/simple_optimizer.o \
                                              $OBJDIR/simple_opt_spec.o \
                                              $OBJDIR/simple_opt_subs.o \
                                              $OBJDIR/simple_ran_tabu.o \
                                              $OBJDIR/simple_stat.o \
                                              $OBJDIR/simple_AVratios.o \
                                              $OBJDIR/simple_cmdline.o \
                                              $OBJDIR/simple_cmd_dict.o \
                                              $OBJDIR/simple_chash.o \
                                              $OBJDIR/simple_args.o \
                                              $OBJDIR/simple_map_reduce.o \
                                              $OBJDIR/simple_image.o \
                                              $OBJDIR/simple_ftiter.o \
                                              $OBJDIR/simple_fftw3.o \
                                              $OBJDIR/gnufor2.o \
                                              $OBJDIR/simple_cuda_defs.o \
                                              $OBJDIR/simple_imgfile.o \
                                              $OBJDIR/simple_sll.o \
                                              $OBJDIR/simple_arr.o \
                                              $OBJDIR/simple_winfuns.o \
                                              $OBJDIR/simple_online_var.o \
                                              $OBJDIR/simple_timing.o \
                                              $OBJDIR/simple_sorting.o \
                                              $OBJDIR/simple_error_handling.o \
                                              $OBJDIR/simple_eglossary_lowlev.o \
                                              $OBJDIR/simple_yaml_strings.o \
                                              $OBJDIR/simple_file_defs.o \
                                              $OBJDIR/simple_err_defs.o \
                                              $OBJDIR/simple_eglossary.o \
                                              $OBJDIR/get_systemQuery_cpu.o \
                                              $OBJDIR/get_systemQuery_cpu_c.o \
                                              $OBJDIR/fortran.o \
                                              $OBJDIR/simple_cudnn_fortran.o \
                                              $OBJDIR/get_fft123D_cpu_c.o \
                                              $OBJDIR/get_fft123D_cpu.o \
                                              $OBJDIR/strlen.o \
                                              $OBJDIR/get_cpu_time_c.o \
                                              $OBJDIR/timestamp.o \
                                              $OBJDIR/timming.o \
                                              $OBJDIR/timming_c.o \
                                              $OBJDIR/polarft_corr_Helpr_gpu_c.o \
                                              -L $FFTW_LIB -lfftw3 \
                                              -L $FFTW_LIB -lfftw3f \
                                              -lblas \
                                              -llapack \
                                              -lstdc++ \
                                              -lpthread \
                                              -lm
mv $FLNAME $SOURCE_DIR/bin/bin_tests
