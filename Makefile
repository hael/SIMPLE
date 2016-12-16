include Makefile_macros

# Source directories.

SRCDIRS = checks defs/gpu/OpenCL defs/gpu/cuda defs/gpu defs production simple_utils/common/mpe simple_utils/common simple_utils src/lapack/clapack src/lapack/lapack95 src/lapack src/simple_gpu/OpenCL src/simple_gpu/cuda src/simple_gpu src/simple_main/extra src/simple_main src test_code/cpu test_code/gpu/OpenCL test_code/gpu/cuda test_code/gpu test_code/mpi test_code

# Search directories.

vpath %.c   $(SRCDIRS)
vpath %.cpp $(SRCDIRS)
vpath %.cu  $(SRCDIRS)
vpath %.f   $(SRCDIRS)
vpath %.for $(SRCDIRS)
vpath %.t90 $(SRCDIRS)
vpath %.mpi.f90 $(SRCDIRS)
vpath %.f90 $(SRCDIRS)
vpath %.o   $(OBJDIR)

# Implicit rules.

%.o:%.c
	$(CCLIB) -o $(OBJDIR)/$@ $<

%.o:%.cpp
	$(CPPCLIB) -o $(OBJDIR)/$@ $<

%.o:%.cu
	$(NVCCCLIB) -o $(OBJDIR)/$@ $<

%.o:%.f
	$(F77CLIB) -o $(OBJDIR)/$@ $<

%.o:%.for
	$(F90CLIB77) -o $(OBJDIR)/$@ $<; $(F90POST)

%.o:%.t90
	./template $<
	$(F90CLIB) -o $(OBJDIR)/$@ $(basename $<).f90; $(F90POST)

%.o:%.mpi.f90
	$(MPIF90CLIB) -o $(OBJDIR)/$@ $<; $(MPIF90POST)

%.o:%.f90
	$(F90CLIB) -o $(OBJDIR)/$@ $<; $(F90POST)

# Targets.

.PHONY: checks_code OpenCL_defs_code cuda_defs_code gpu_code defs_code production_code utils_common_mpe_code utils_common_code utils_code clapack_code lapack95_code lapack_code src_gpu_kernel_OpenCL_code src_gpu_kernel_cuda_code src_gpu_code simple_extra_code Simple_code src_code cpu_test_code openCL_test_code cuda_test_code gpu_test_code mpi_test_code test_code default all makedir checkclean clean cleanall robodoc check_news finalCmpl wc tar detar;

default:
	@./makemake > /dev/null; \
	make --no-print-directory all

all: makedir checkclean checks_code OpenCL_defs_code cuda_defs_code gpu_code defs_code production_code utils_common_mpe_code utils_common_code utils_code clapack_code lapack95_code lapack_code src_gpu_kernel_OpenCL_code src_gpu_kernel_cuda_code src_gpu_code simple_extra_code Simple_code src_code cpu_test_code openCL_test_code cuda_test_code gpu_test_code mpi_test_code test_code finalCmpl;

checks_code: ;

OpenCL_defs_code: OpenCL_cpp_defs ;

OpenCL_cpp_defs: ;

cuda_defs_code: cuda_defs ;

cuda_defs: fortran.o              \
           simple_fortran.o       \
           simple_cudnn_fortran.o ;

gpu_code: ;

defs_code: numRecp                      \
           ps_defs                      \
           ps_mpi_defs                  \
           ps_cuda_defs                 \
           ps_dpl_defs                  \
           ps_openCL_defs               \
           ps_lattice_defs              \
           ps_file_defs                 \
           ps_reso_defs                 \
           ps_err_defs                  \
           s_defs                       \
           s_hdf5                       \
           s_h5_ncem                    \
           s_h5_rdwt                    \
           fftw3                        ;

numRecp: nrtype.o                       \
         nr.o                           \
         nrutil.o                       ;

ps_defs: simple_defs.o                  ;

ps_mpi_defs:  simple_mpi_defs.o         ;

ps_cuda_defs: simple_cuda_defs.o        ;

ps_dpl_defs: simple_deepLearning_defs.o ;

ps_openCL_defs: simple_OpenCL_defs.o    ;

ps_lattice_defs: simple_lattice_defs.o  ;

ps_file_defs: simple_file_defs.o        ;

ps_reso_defs: simple_resources_defs.o   ;

ps_err_defs: simple_err_defs.o          ;

s_defs:                                 ;

s_hdf5:                                 ;

s_h5_ncem:                              ;
s_h5_rdwt:                              ;

fftw3: simple_fftw3.o                   ;

production_code: ;

utils_common_mpe_code: mpe_code ;

mpe_code: mpe_log.o             \
          clog_commset.o        ;

utils_common_code: string_handling    \
                   s_testfunc         \
                   textHandler        \
                   greetings          \
                   matrix_getter      \
                   sortters           \
                   su3_matrix         \
                   yaml_strings       \
                   eglossary          \
                   malloc             \
                   profiling          \
                   file_handler       \
                   yaml_modules       \
                   memory             \
                   timers             \
                   file_high_lev      \
                   buffers            \
                   random_code        \
                   matrix_multiplier  \
                   matrix_testers     \
                   deviceQuery        \
                   systemQuery        \
                   cuda_fft           ;

string_handling: strlen.o             ;

s_testfunc: simple_testfuns.o         \
            simple_testfunction.o     ;

textHandler: simple_textHandler.o     ;

greetings: greeting_version.o         ;


yaml_strings: simple_yaml_strings.o   ;

eglossary: simple_eglossary_lowlev.o  \
           simple_error_handling.o    \
           simple_eglossary.o         ;

malloc: simple_module_file_malloc.o   ;

profiling: simple_memory_profiling.o  ;

file_handler: simple_file_utils.o     \
              simple_DataParallel_IO.o;

yaml_modules: simple_yaml_output.o    ;

memory: simple_dynamic_memory.o       ;

sortters: simple_sorting.o            ;

timers: timming.o                     \
        timming_c.o                   \
        simple_timing.o               \
        get_cpu_time_c.o              \
        timestamp.o                   ;

file_high_lev: simple_file_highlev.o  ;

su3_matrix: initialu.o                \
            su2random.o               \
            su3random.o               ;

buffers:                              ;
matrix_getter: matrixGetter.o         ;
random_code: simple_random.o          ;
matrix_multiplier:                    ;

matrix_testers: simple_SU3_tester.o   ;

deviceQuery: simple_deviceQuery_gpu.o \
             get_deviceQuery_gpu.o    \
             get_deviceQuery_gpu_c.o  ;

systemQuery: simple_systemQuery_cpu.o \
             get_systemQuery_cpu.o    \
             get_systemQuery_cpu_c.o  ;

cuda_fft: fft123D_cpu.o               \
          get_fft123D_cpu_c.o         \
          get_fft123D_cpu.o           \
          fft123D_gpu.o               \
          get_fft123D_gpu_c.o         \
          get_fft123D_gpu.o           ;

utils_code: s_utes  \
            s_utils ;

s_utes: simple_imgheadrec.o \
        simple_imghead.o    \
        simple_jiffys.o     \
        simple_rnd.o        \
        simple_math.o       \
        gnufor2.o           \
        simple_imgfile.o    \
        simple_stat.o       \
        simple_online_var.o \
        simple_ran_tabu.o   \
        simple_strings.o    ;


s_utils: simple_syscalls.o  ;

clapack_code: ;

lapack95_code: lapack95_lib ;

lapack95_lib: 

	@make --makefile=src/lapack/lapack95/Makefile_lapack95

lapack_code: ;

src_gpu_kernel_OpenCL_code: kernel_opencl_code ;

kernel_opencl_code:  ;

src_gpu_kernel_cuda_code: kernel_cuda ;

kernel_cuda: ;

src_gpu_code: cuda_gpu       \
              invert_gpu     \
              simple_gpu     \
              math_gpu       \
              polarft_gpu    \
              genericClasses \
              cartsian2D     \
              deepLearning   ;

cuda_gpu: simple_cuda.o ;

invert_gpu: invert_gpu_utils.o            \
            magma_zgetri_blocked_gpu-v2.o \
            magma_dgetri_blocked_gpu-v2.o \
            magma_get_getri_nb_gpu.o      \
            magma_invert_gpu-v2.o         ;

simple_gpu: ;

math_gpu: matmul_cpu_utils.o      \
          matmul_gpu_utils.o      \
          matvec_gpu_utils.o      \
          matvec_cpu_utils.o      \
          invert_cpu_utils.o      \
          matmul_cuda_transfers.o \
          matmul_cuda_device.o    \
          matvec_cuda_transfers.o \
          matvec_cuda_device.o    \
          simple_math_gpu.o       \
          simple_math_gpu_c.o     \
          gen_polar_coords_gpu.o  ;

polarft_gpu: ;

genericClasses: Global_polarft.o           \
                simple_img_2D_cart_Sizes.o \
                simple_mesh1D.o            \
                simple_mesh1DV.o           \
                simple_mesh3D.o            \
                simple_pfts_Sizes.o        ;

cartsian2D: carte2D_ftExt-corr_gpu.o \
            carte2D_ftExt-corr_C_N.o \
            carte2D_ftExt-corr_C_F.o ;

deepLearning: simple_deepLearning.o  ;

simple_extra_code: ;

Simple_code: opt            \
             general        \
             simple         \
             build          \
             comlin_corr    \
             comlin_srch    \
             s_other        \
             matcher        \
             comlin_sym     \
             cluster_cavg   \
             qsys           \
             masters        \
             highlev        ;

opt: simple_opt_spec.o               \
     simple_optimizer.o              \
     simple_opt_subs.o               \
     simple_simplex_opt.o            \
     simple_sa_opt.o                 \
     simple_powell_opt.o             \
     simple_oasis_opt.o              \
     simple_particle_swarm_opt.o     \
     simple_de_opt.o                 \
     simple_bforce_opt.o             \
     simple_bfgs_opt.o               \
     simple_opt_factory.o            ;

general: simple_hash.o            \
         simple_chash.o           \
         simple_softmax_weights.o \
         simple_ori.o             \
         simple_nrtxtfile.o       \
         simple_args.o            \
         simple_gen_doc.o         \
         simple_cmd_dict.o        \
         simple_cmdline.o         \
         simple_AVratios.o        \
         simple_map_reduce.o      \
         simple_params.o          \
         simple_oris.o            \
         simple_ftiter.o          \
         simple_arr.o             \
         simple_sll.o             \
         simple_winfuns.o         ;

simple: simple_image.o                  \
        simple_ft_expanded.o            \
        simple_ctf.o                    \
        simple_sym.o                    \
        simple_corrmat.o                \
        simple_gridding.o               \
        simple_polarft.o                \
        simple_polarft_corrcalc.o       \
        simple_projector.o              \
        simple_volpft_corrcalc.o        \
        simple_ft_shsrch.o              \
        simple_ftexp_shsrch.o           \
        simple_polarft_shsrch.o         \
        simple_volpft_srch.o            \
        simple_pftcc_shsrch.o           \
        simple_pftcc_inplsrch.o         \
        simple_nnimgs.o                 \
        simple_filterer.o               \
        simple_unblur.o                 \
        simple_convergence.o            \
        simple_prime_srch.o             \
        simple_prime2D_srch.o           \
        simple_prime3D_srch.o           \
        simple_procimgfile.o            ;

build: simple_pair_dtab.o               \
       simple_cluster_valid.o           \
       simple_hac.o                     \
       simple_comlin.o                  \
       simple_kmeans.o                  \
       simple_centre_clust.o            \
       simple_shc_cluster.o             \
       simple_denspeak_cluster.o        \
       simple_aff_prop.o                \
       simple_clusterer.o               \
       simple_spatial_median.o          \
       simple_ppca.o                    \
       simple_estimate_ssnr.o           \
       simple_reconstructor.o           \
       simple_simulator.o               \
       simple_eo_reconstructor.o        \
       simple_build.o                   ;
       

comlin_corr: simple_comlin_corr.o       ;

comlin_srch: simple_comlin_srch.o       ;

matcher: simple_picker.o                \
         simple_hadamard_common.o       \
         simple_cartft_corrcalc.o       \
         simple_hadamard2D_matcher.o    \
         simple_hadamard3D_matcher.o    \
         simple_crossover.o             \
         simple_scatter_orisrch.o       \
         simple_cftcc_srch.o            \
         simple_cont3D_matcher.o        ;

s_other: simple_masker.o                \
         simple_spatial_median.o        ;
         
comlin_sym: simple_comlin_symsrch.o     \
            simple_comlin_sym.o         ;

cluster_cavg: simple_matchpursuit.o     \
              simple_cavgppca.o         \
              simple_refinecluster.o    \
              simple_classrefine.o      \
              simple_clustercavg.o      ;

qsys: simple_qsys_base.o                \
      simple_qsys_slurm.o               \
      simple_qsys_local.o               \
      simple_qsys_factory.o             \
      simple_qsys_ctrl.o                \
      simple_qsys_funs.o                ;

masters: simple_rec_master.o             \
         simple_symsrcher.o              \
         simple_commander_base.o         \
         simple_commander_checks.o       \
         simple_commander_comlin.o       \
         simple_commander_distr.o        \
         simple_commander_imgproc.o      \
         simple_commander_mask.o         \
         simple_commander_misc.o         \
         simple_commander_oris.o         \
         simple_commander_preproc.o      \
         simple_commander_prime2D.o      \
         simple_commander_prime3D.o      \
         simple_commander_rec.o          \
         simple_commander_sim.o          \
         simple_commander_volops.o       \
         simple_commander_distr_wflows.o ;

highlev: simple_highlev.o               ;

src_code: ;

cpu_test_code: s_test        \
               cpu_test      ;

s_test: simple_optimiser_tester.o       \
        simple_scatter_orisrch_tester.o \
        simple_prime2D_srch_tester.o    \
        simple_wiener2D_tester.o        \
        simple_inpl_srch_tester.o       \
        simple_ft_expanded_tester.o     \
        simple_volpft_srch_tester.o     \
        simple_classrefine_tester.o     \

cpu_test: simple_prime3D_srch_tester.o ;

cpu_test:  ;

openCL_test_code: openCL_test      ;

openCL_test: getGPU_interface.o    \
             testing_openCL_defs.o ;

cuda_test_code: cuda_test      \
                testfunc       \
                classes        \
                corr_calc      \
                error_handlers \
                wavelets       \
                corr_calc_opt  ;

cuda_test: simple_polarft_gpu_c.o            \
           zz2dgemm_ElmtWs_tesla_N_N.o       \
           zz2dgemm_ElmtWs_tesla_T_N.o       \
           zz2dgemm_ElmtWs_tesla_N_T.o       \
           zz2dgemm_ElmtWs_tesla_T_T.o       \
           zz2dgemm_ElmtWs_tesla_gpu.o       \
           zz2dgemm_ElmtWs_tesla_sumsq_N_N.o \
           zz2dgemm_ElmtWs_tesla_sumsq_T_T.o \
           zz2dgemm_ElmtWs_tesla_sumsq_gpu.o ;

testfunc: ;

classes: ;

corr_calc: simple_polarft_corr_gpu_c.o \
           simple_polarft_corr_gpu.o   \
           polarft_corr_Helpr_gpu.o    \
           polarft_corr_Helpr_gpu_c.o  \
           polarft_corr_Hadmr_gpu.o    \
           polarft_corr_Hadmr_N_N.o    \
           polarft_corr_Hadmr_F_N.o    \
           polarft_corr_Hadmr_P_N.o    \
           polarft_corr_Hadmr_X_N.o    \
           polarft_gencorrAll_gpu.o    \
           polarft_gencorrAll_Z_N.o    \
           polarft_multi-GPUs_gpu.o    \
           polarft_multi-GPUs_Z_M.o    ;

error_handlers: simple_ErrorHandler.o  ;

wavelets: Global_wavelets.o      \
          waveletsD_c.o          \
          waveletsD.o            \
          wavelets_Coefficient.o \
          wavelets_FilterCoefs.o \
          wavelets_TransForm1D.o \
          wavelets_TransForm2D.o \
          wavelets_TransFormCm.o \
          wavelets_Convoluters.o ;

corr_calc_opt: polarft_krnl-Opti_gpu.o    \
               polarft_krnl-Opti_O_N.o    ;

gpu_test_code: ;

mpi_test_code: mpi_test      ;

mpi_test: simple_matrixHandler_mpi_cpu.o;

test_code: ;
 
makedir:
	@if [ ! -d $(OBJDIR) ]; \
	then \
		echo "Creating directory "$(OBJDIR)"."; \
		mkdir $(OBJDIR); \
	fi; \
	if [ ! -d $(MODDIR) ]; \
	then \
		echo "Creating directory "$(MODDIR)"."; \
		mkdir $(MODDIR); \
	fi

checkclean:
	@if [ -e CLEAN ]; \
	then \
		clean=1; \
		if [ -e $(OBJDIR)/LASTCLEAN ]; \
		then \
			if cmp -s CLEAN $(OBJDIR)/LASTCLEAN; \
			then \
				clean=0; \
			fi; \
		fi; \
		if [ $$clean -eq 1 ]; \
		then \
			make -s clean; \
			cp CLEAN $(OBJDIR)/LASTCLEAN; \
		fi \
	fi; \

clean:
	@echo "Cleaning .o and .mod files in "$(OBJDIR)"."; \
	rm -f $(OBJDIR)/*.o; \
	rm -f $(MODDIR)/*.mod; \
	rm -f $(MODDIR)/*.h

cleanall:
	@echo "Cleaning .o and .mod files in "$(OBJDIR)" and "$(OBJDIR)"/lapack/lapack95."; \
	rm -f $(OBJDIR)/*.o; \
	rm -f $(OBJDIR)/lapack/lapack95/*.o; \
	rm -f $(MODDIR)/*.mod; \
	rm -f $(MODDIR)/*.h; \
        cd ./checks/simple_LibTests/cpu_test_code; \
        pwd; \
        make clean; \
        pwd; \
        cd ./../gpu_test_code; \
        pwd; \
        make clean; \
        pwd; \
        cd ./../Clustering/; \
        pwd; \
        make clean; \
        pwd; \
        cd ./../1WCM/; \
        pwd; \
        make clean; \
        pwd; \
        cd ../../../bin/bin_tests; \
        pwd; \
        rm simple_test_*; \
        pwd; \
        cd ../; \
        pwd; \
        rm simple_*; \
        pwd
robodoc:
	@echo "Building RoBoDoc documentation."; \
	make check_news > doc/robodoc/.lastnews; \
	robodoc --rc doc/robodoc/robodoc.rc; \
	cd robodoc; \
	./makemdoc;
check:
	@echo "Running building in check scenarios for simple models."; \
        export PYTHONPATH=/Applications/EMAN//EMAN.1.9/lib:/Applications/EMAN2//lib:/Applications/EMAN2//bin:/Applications/EMAN2//extlib/site-packages:/Applications/EMAN2//extlib/site-packages/ipython-1.2.1-py2.7.egg:; \
        export LD_LIBRARY_PATH=; \
        python ./checks/simple_LibTests/report.py --use_gpu=$(use_gpu) --bench_gpu=$(bench_gpu) --help=$(help)
bench:
	@echo "Running building in check scenarios for simple models."; \
        export PYTHONPATH=/Applications/EMAN//EMAN.1.9/lib:/Applications/EMAN2//lib:/Applications/EMAN2//bin:/Applications/EMAN2//extlib/site-packages:/Applications/EMAN2//extlib/site-packages/ipython-1.2.1-py2.7.egg:; \
        export LD_LIBRARY_PATH=; \
        python ./checks/simple_LibTests/bench.py --use_gpu=$(use_gpu) --bench_gpu=$(bench_gpu) --fix_gpu=$(fix_gpu) --set_gpu=$(set_gpu) --help=$(help)
check_help:
	@echo "Running building in check scenarios for simple models."; \
        export PYTHONPATH=/Applications/EMAN//EMAN.1.9/lib:/Applications/EMAN2//lib:/Applications/EMAN2//bin:/Applications/EMAN2//extlib/site-packages:/Applications/EMAN2//extlib/site-packages/ipython-1.2.1-py2.7.egg:; \
        export LD_LIBRARY_PATH=; python ./checks/simple_LibTests/report.py --help
check_cpu:
	@echo "moving to directory: ./checks/simple_LibTests/cpu_test_code and compiling check codes..."; \
        cd ./checks/simple_LibTests/cpu_test_code; \
        pwd; \
        make; \
        ./run_Testing_cpu.csh; \
        cd ../../../
clean_check_cpu:
	@echo "moving to directory: ./checks/simple_LibTests/cpu_test_code and compiling check codes..."; \
        cd ./checks/simple_LibTests/cpu_test_code; \
        pwd; \
        make clean; \
        cd ../Clustering/; \
        make clean; \
        pwd; \
        cd ../1WCM/; \
        make clean; \
        pwd
check_gpu:
	@echo "moving to directory: ./checks/simple_LibTests/gpu_test_code and compiling check codes..."; \
        cd ./checks/simple_LibTests/gpu_test_code; \
        pwd; \
        make; \
        ./run_Testing_gpu.csh; \
        cd ../../../
clean_check_gpu:
	@echo "moving to directory: ./checks/simple_LibTests/gpu_test_code and compiling check codes..."; \
        cd ./checks/simple_LibTests/gpu_test_code; \
        pwd; \
        make clean
finalCmpl:
	@news=1; \
	if [ -e news/.lastnews ]; \
	then \
		if cmp -s news/news news/.lastnews; \
		then \
			news=0; \
		fi; \
	fi; \
	if [ $$news -eq 1 ]; \
	then \
		echo ""; \
		tput bold; \
		tput setaf 4; \
		echo "**************************************************************"; \
		echo "* The Simple-v3.0 news:                                      *"; \
		echo "* Date: 16 May 2016                                          *"; \
		echo "* Date: currently partially completed,                       *"; \
		echo "*       completion anticipated at latest January 2017        *"; \
		echo "* 1.) GPU accelerated motion correction                      *"; \
		echo "* 2.) GPU accelerated PRIME-2D/3D                            *"; \
		echo "* 3.) GPU accelerated multi-particle reconstruction          *"; \
		echo "*     (heterogeneity analysis)                               *"; \
		echo "**************************************************************"; \
		tput sgr0; \
		echo ""; \
	fi; \

check_news:
	@news=1; \
	if [ -e news/.lastnews ]; \
	then \
		if cmp -s news/news news/.lastnews; \
		then \
			news=0; \
		fi; \
	fi; \
	if [ $$news -eq 1 ]; \
	then \
		echo ""; \
		tput bold; \
		tput setaf 4; \
		echo "**************************************************************"; \
		echo "* The Simple-v2.1 news:                                      *"; \
		echo "* Date: 16 May 2016                                          *"; \
                echo "* New 2D alignment/clustering code (PRIME2D) that            *"; \
		echo "* 1.) generates improved 2D class averages from large        *"; \
		echo "*     single-particle cryo-EM datasets                       *"; \
		echo "* 2.) can be used in conjunction wiht PRIME3D to obtain a    *"; \
		echo "*     reliable 3D starting model in a rapid and unbiased     *"; \
		echo "*     fashion                                                *"; \
		echo "* 3.) overcomes inherent limitations in widely used          *"; \
		echo "*     clustering  approaches                                 *"; \
		echo "*                                                            *"; \
		echo "* Our next release: Simple-v3.0                              *"; \
		echo "* Date: currently partially completed,                       *"; \
		echo "*       completion anticipated at latest January 2017        *"; \
		echo "* 1.) GPU accelerated motion correction                      *"; \
		echo "* 2.) GPU accelerated PRIME-2D/3D                            *"; \
		echo "* 3.) GPU accelerated multi-particle reconstruction          *"; \
		echo "*     (heterogeneity analysis)                               *"; \
		echo "**************************************************************"; \
		tput sgr0; \
		echo ""; \
	fi; \

wc:
	@echo "Count the number of lines, words and characters in the modules."; \
	unset list; \
	DIRS="$(SRCDIRS)"; \
	for dir in $(SRCDIRS); \
	do \
		for ext in "*.c" "*.cpp" "*.cu" "*.h" "*.f90" "*.f" "*.for"; \
		do \
			list=$$list" "$$(find $$dir"/" -maxdepth 1 -name "$$ext"); \
		done; \
	done; \
	wc $$list

tar:
	@echo "Creating a tar archive 'modules.tar.gz' with the modules."; \
	unset list; \
	for dir in $(SRCDIRS); \
	do \
		for ext in "*.c" "*.cpp" "*.cu" "*.h" "*.f90" ".t90" "*.f" "*.for" "Makefile_*" "*.par*"; \
		do \
			list=$$list" "$$(find $$dir"/" -maxdepth 1 -name "$$ext"); \
		done; \
	done; \
	tar czvf modules.tar.gz makemake template macros/Makefile_macros* $$list

detar:
	@echo "Creating a tar backup 'modules.old.tar.gz' with the modules before unpacking the tar archive 'modules.tar.gz'."; \
	unset list; \
	for dir in $(SRCDIRS); \
	do \
		for ext in "*.c" "*.cpp" "*.cu" "*.h" "*.f90" ".t90" "*.f" "*.for" "Makefile_*" "*.par*"; \
		do \
			list=$$list" "$$(find $$dir"/" -maxdepth 1 -name "$$ext"); \
		done; \
	done; \
	tar czvf modules.old.tar.gz makemake template macros/Makefile_macros* $$list; \
	tar xzvf modules.tar.gz
