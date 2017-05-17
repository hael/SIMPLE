include Makefile_macros

# Source directories.

SRCDIRS = defs production simple_utils/cuda simple_utils src/simple_main src test_code/cpu test_code

# Search directories.

vpath %.c   $(SRCDIRS)
vpath %.cpp $(SRCDIRS)
vpath %.cu  $(SRCDIRS)
vpath %.f   $(SRCDIRS)
vpath %.for $(SRCDIRS)
vpath %.f08 $(SRCDIRS)
vpath %.t90 $(SRCDIRS)
vpath %.mpi.f90 $(SRCDIRS)
vpath %.f90 $(SRCDIRS)
vpath %.o   $(OBJDIR)
vpath %.i   $(OBJDIR)

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

%.o:%.f08
	$(FPP) $(FPPFLAGS)  -o $(OBJDIR)/$(basename $@).f08 $<
	$(F90CLIB) -o $(OBJDIR)/$@ $(OBJDIR)/$(basename $@).f08

%.o:%.t90
	./template $<
	$(F90CLIB) -o $(OBJDIR)/$@ $(basename $<).f90; $(F90POST)

%.o:%.mpi.f90
	$(MPIF90CLIB) -o $(OBJDIR)/$@ $<; $(MPIF90POST)

%.o:%.f90
	$(F90CLIB) -o $(OBJDIR)/$@ $<; $(F90POST)

# Targets.

.PHONY: defs_code production_code cuda_mods utils_code Simple_code src_code cpu_test_code test_code default all makedir checkclean clean cleanall wc;

default:
	@./makemake > /dev/null; \
	make --no-print-directory all $(MFLAGS)

all: makedir checkclean defs_code production_code cuda_mods utils_code Simple_code src_code cpu_test_code test_code;

defs_code: simple_defs.o      \
           simple_defs_conv.o \
           simple_fftw3.o     ;

production_code: ;

cuda_mods: precision_m.o;

utils_code: s_utes                \
            s_utils               ;

s_utes: simple_arr.o              \
        simple_sll.o              \
        simple_strings.o          \
        simple_imgheadrec.o       \
        simple_imghead.o          \
        simple_filehandling.o     \
        simple_jiffys.o           \
        simple_math.o             \
        simple_magic_boxes.o      \
        simple_rnd.o              \
        gnufor2.o                 \
        simple_imgfile.o          \
        simple_stat.o             \
        simple_online_var.o       \
        simple_ran_tabu.o         \
        simple_testfuns.o         ;

s_utils: simple_timer_omp.o    \
         simple_timer.o       \
	       simple_syscalls.o  ;

Simple_code: opt            \
             general        \
             simple         \
             build          \
             s_other        \
             qsys           \
             matcher        \
             comlin_srch    \
             cluster_cavg   \
             masters        ;

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
         simple_sauron.o          \
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
         simple_kbinterpol.o      \
         simple_winfuns.o         \
         simple_fsc_compare.o     \
         simple_restart.o         ;

simple: simple_image.o                  \
        simple_ft_expanded.o            \
        simple_ctf.o                    \
        simple_sym.o                    \
        simple_corrmat.o                \
        simple_gridding.o               \
        simple_polarft_corrcalc.o       \
        simple_projector.o              \
        simple_mask_projector.o         \
        simple_pftcc_opt.o              \
        simple_simplex_pftcc_opt.o      \
        simple_projector_hlev.o         \
        simple_volpft_corrcalc.o        \
        simple_ft_shsrch.o              \
        simple_ftexp_shsrch.o           \
        simple_volpft_srch.o            \
        simple_pftcc_shsrch.o           \
        simple_pftcc_inplsrch.o         \
        simple_pftcc_srch.o             \
        simple_filterer.o               \
        simple_unblur.o                 \
        simple_convergence.o            \
        simple_convergence_perptcl.o    \
        simple_prime_srch.o             \
        simple_prime2D_srch.o           \
        simple_prime3D_srch.o           \
        simple_cont3D_srch.o            \
        simple_cont3D_greedysrch.o      \
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
       simple_combinatorics.o           \
       simple_ppca.o                    \
       simple_estimate_ssnr.o           \
       simple_reconstructor.o           \
       simple_simulator.o               \
       simple_eo_reconstructor.o        \
       simple_build.o                   ;

comlin_srch: simple_comlin_srch.o       ;

matcher: simple_picker.o                  \
         simple_tseries_tracker.o         \
         simple_hadamard_common.o         \
         simple_cartft_corrcalc.o         \
         simple_scatter_orisrch.o         \
         simple_cftcc_shsrch.o            \
         simple_cont3D_matcher.o          \
         simple_hadamard2D_matcher.o      \
         simple_hadamard3D_matcher.o      ;

s_other: simple_masker.o                  ;

cluster_cavg: simple_matchpursuit.o       \
              simple_cavgppca.o           ;

qsys: simple_qsys_base.o                  \
      simple_qsys_slurm.o                 \
      simple_qsys_sge.o                   \
      simple_qsys_pbs.o                   \
      simple_qsys_local.o                 \
      simple_qsys_factory.o               \
      simple_qsys_ctrl.o                  \
      simple_qsys_funs.o                  \
      simple_qsys_env.o                   ;

masters: simple_rec_master.o              \
         simple_commander_base.o          \
         simple_commander_volops.o        \
         simple_commander_checks.o        \
         simple_commander_comlin.o        \
         simple_commander_distr.o         \
         simple_commander_imgproc.o       \
         simple_scaler.o                  \
         simple_commander_mask.o          \
         simple_commander_misc.o          \
         simple_commander_oris.o          \
         simple_unblur_iter.o             \
         simple_ctffind_iter.o            \
         simple_pick_iter.o               \
         simple_commander_preproc.o       \
         simple_commander_prime2D.o       \
         simple_commander_prime3D.o       \
         simple_commander_rec.o           \
         simple_commander_sim.o           \
         simple_commander_tseries.o       \
         simple_commander_stream_wflows.o \
         simple_commander_distr_wflows.o  \
         simple_commander_hlev_wflows.o   ;
         

src_code: ;

cpu_test_code: simple_optimiser_tester.o       \
               simple_scatter_orisrch_tester.o \
               simple_wiener2D_tester.o        \
               simple_ft_expanded_tester.o     \
               simple_speedtester.o            \
               simple_volpft_srch_tester.o     \
               simple_timer_basic_test.o       \
               simple_timer_profile_test.o     \
               simple_timer_omp_test.o         ;

#Variadic macros cannot be done with gfortran alone,
# cpp produces the *.i, gfortran can compile the rest
# simple_timer_profile_test.o: simple_timer_profile_test.f08

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
	fi;

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
	fi;

clean:
	@echo "Cleaning .o and .mod files in "$(OBJDIR)"."; \
	rm -f $(OBJDIR)/*.o; \
	rm -f $(MODDIR)/*.mod; \
	rm -f $(MODDIR)/*.h;

cleanall:
	@echo "Cleaning .o and .mod files in "$(OBJDIR)"."; \
    rm -f $(OBJDIR)/*.o; \
    rm -f $(MODDIR)/*.mod; \
    rm -f $(MODDIR)/*.h; \
    cd ./bin/bin_tests; \
    pwd; \
    rm simple_test_*; \
    pwd; \
    cd ../; \
    pwd; \
    rm simple_*; \
    pwd;

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
	wc $$list;

