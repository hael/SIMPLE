FLNAME=`basename $1 .f90`

GFORTRAN=gfortran

SOURCE_DIR=/Users/hael/src/fortran/simple_wfrederic/Simple_Restruct.projet/

MODDIR=$SOURCE_DIR/obj/GFORTgpu
OBJDIR=$SOURCE_DIR/obj/GFORTgpu

MAGMADIR=/opt/Devel_tools/magma-1.6.1
MAGMADIR_LIB=$MAGMADIR/lib

CUDADIR=/usr/local/cuda

SHARE_LIB=/opt/Devel_tools/NVIDIA_GPU_Computing_SDK/shared/lib/linux

PLATFORM=`uname`
if [ "$PLATFORM" = "Linux" ]
then
 PLAT="-DLINUX"
 elif [ "$PLATFORM" = "Darwin" ]
then
 PLAT="-DMACOSX"
else
echo "You seem to be compiling/running in wonderland"
fi

OPTION="-ffree-form -cpp -O3 -fPIC -fno-second-underscore -fopenmp $PLAT -DBENCH"

FFTW_LIB=/sw/lib

$GFORTRAN $OPTION -I $OBJDIR -I $MODDIR -o $FLNAME \
                                           $FLNAME.f90 \
                                           $OBJDIR/simple_cmdline.o \
                                           $OBJDIR/simple_defs.o \
                                           $OBJDIR/simple_args.o \
                                           $OBJDIR/simple_jiffys.o \
                                           $OBJDIR/simple_comlin.o \
                                           $OBJDIR/simple_image.o \
                                           $OBJDIR/simple_ftiter.o \
                                           $OBJDIR/simple_math.o \
                                           $OBJDIR/simple_rnd.o \
                                           $OBJDIR/simple_fftw3.o \
                                           $OBJDIR/gnufor2.o \
                                           $OBJDIR/simple_stat.o \
                                           $OBJDIR/simple_sll.o \
                                           $OBJDIR/simple_arr.o \
                                           $OBJDIR/simple_ran_tabu.o \
                                           $OBJDIR/simple_winfuns.o \
                                           $OBJDIR/simple_online_var.o \
                                           $OBJDIR/simple_oris.o \
                                           $OBJDIR/simple_ori.o \
                                           $OBJDIR/simple_hash.o \
                                           $OBJDIR/simple_strings.o \
                                           $OBJDIR/simple_softmax_weights.o \
                                           $OBJDIR/simple_opt_spec.o \
                                           $OBJDIR/simple_simplex_opt.o \
                                           $OBJDIR/simple_optimizer.o \
                                           $OBJDIR/simple_opt_subs.o \
                                           $OBJDIR/simple_hac.o \
                                           $OBJDIR/simple_pair_dtab.o \
                                           $OBJDIR/simple_kmeans.o \
                                           $OBJDIR/simple_projector.o \
                                           $OBJDIR/simple_params.o \
                                           $OBJDIR/simple_units.o \
                                           $OBJDIR/simple_ctf.o \
                                           $OBJDIR/simple_gridding.o \
                                           $OBJDIR/simple_polarft.o \
                                           $OBJDIR/simple_reconstructor.o \
                                           $OBJDIR/simple_procimgfile.o \
                                           $OBJDIR/simple_sym.o \
                                           $OBJDIR/simple_eo_reconstructor.o \
                                           $OBJDIR/simple_masker.o \
                                           $OBJDIR/simple_opt_factory.o \
                                           $OBJDIR/simple_bfgs_opt.o \
                                           $OBJDIR/simple_powell_opt.o \
                                           $OBJDIR/simple_oasis_opt.o \
                                           $OBJDIR/simple_testfuns.o \
                                           $OBJDIR/simple_particle_swarm_opt.o \
                                           $OBJDIR/simple_de_opt.o \
                                           $OBJDIR/matrixGetter.o \
                                           $OBJDIR/initialu.o \
                                           $OBJDIR/su3random.o \
                                           $OBJDIR/su2random.o \
                                           $OBJDIR/simple_SU3_tester.o \
                                           $OBJDIR/simple_timing.o \
                                           $OBJDIR/simple_random.o \
                                           $OBJDIR/simple_sorting.o \
                                           $OBJDIR/timestamp.o \
                                           $OBJDIR/timming.o \
                                           $OBJDIR/timming_c.o \
                                           $OBJDIR/get_cpu_time_c.o \
                                           $OBJDIR/simple_textHandler.o \
                                           $OBJDIR/greeting_version.o \
                                           $OBJDIR/simple_imgfile.o \
                                           $OBJDIR/simple_syscalls.o \
                                           $OBJDIR/simple_bforce_opt.o \
                                           $OBJDIR/simple_imghead.o \
                                           $OBJDIR/simple_imgheadrec.o \
                                           $OBJDIR/get_systemQuery_cpu.o \
                                           $OBJDIR/simple_systemQuery_cpu.o \
                                           $OBJDIR/get_systemQuery_cpu_c.o \
                                           $OBJDIR/fft123D_cpu.o \
                                           $OBJDIR/get_fft123D_cpu.o \
                                           $OBJDIR/get_fft123D_cpu_c.o \
                                           $OBJDIR/simple_testfunction.o \
                                           $OBJDIR/simple_file_utils.o \
                                           $OBJDIR/simple_DataParallel_IO.o \
                                           $OBJDIR/simple_file_defs.o \
                                           $OBJDIR/simple_resources_defs.o \
                                           $OBJDIR/invert_cpu_utils.o \
                                           $OBJDIR/simple_cuda.o \
                                           $OBJDIR/simple_deviceQuery_gpu.o \
                                           $OBJDIR/get_deviceQuery_gpu.o \
                                           $OBJDIR/get_deviceQuery_gpu_c.o \
                                           $OBJDIR/simple_cuda_defs.o \
                                           $OBJDIR/strlen.o \
                                           $OBJDIR/simple_eglossary.o \
                                           $OBJDIR/simple_err_defs.o \
                                           $OBJDIR/simple_error_handling.o \
                                           $OBJDIR/simple_eglossary_lowlev.o \
                                           $OBJDIR/simple_yaml_output.o \
                                           $OBJDIR/simple_yaml_strings.o \
                                           $OBJDIR/simple_file_highlev.o \
                                           $OBJDIR/simple_dynamic_memory.o \
                                           $OBJDIR/simple_chash.o \
                                           $OBJDIR/simple_map_reduce.o \
                                           -L $FFTW_LIB -lfftw3 \
                                           -L $FFTW_LIB -lfftw3f \
                                           -L $FFTW_LIB -lfftw3_threads \
                                           -lblas \
                                           -llapack \
                                           -lstdc++ \
                                           -lpthread \
                                           -lm

#mv $FLNAME $SOURCE_DIR/bin
