#!/bin/bash
set MY_HOME_SIDE = ~/Monash/SourceCode/Simple/Restructured/HansVersion/backups/Simple_Restruct.projet

set DIFF_PATH = $MY_HOME_SIDE/src/simple_main

echo "::::::::::::::"
echo "$DIFF_PATH/diff_src_cpu.csh ./diff_src_cpu.csh "
echo "::::::::::::::"
diff diff_src_cpu.csh $DIFF_PATH
echo "::::::::::::::"
echo "$DIFF_PATH/.gitignore ./.gitignore"
echo "::::::::::::::"
diff .gitignore $DIFF_PATH
echo "::::::::::::::"
echo ""
echo "::::::::::::::"
diff Makefile_target $DIFF_PATH
echo "::::::::::::::"
echo ""
echo "::::::::::::::"
diff simple_aff_prop.f90 $DIFF_PATH
echo "::::::::::::::"
echo ""
echo "::::::::::::::"
diff simple_arr.f90 $DIFF_PATH
echo "::::::::::::::"
echo ""
echo "::::::::::::::"
diff simple_bfgs_opt.f90 $DIFF_PATH
echo "::::::::::::::"
echo ""
echo "::::::::::::::"
diff simple_bforce_opt.f90 $DIFF_PATH
echo "::::::::::::::"
echo ""
echo "::::::::::::::"
diff simple_build.f90 $DIFF_PATH
echo "::::::::::::::"
echo ""
echo "::::::::::::::"
diff simple_cavgppca.f90 $DIFF_PATH
echo "::::::::::::::"
echo ""
echo "::::::::::::::"
diff simple_centre_clust.f90 $DIFF_PATH
echo "::::::::::::::"
echo ""
echo "::::::::::::::"
diff simple_clustercavg.f90 $DIFF_PATH
echo "::::::::::::::"
echo ""
echo "::::::::::::::"
diff simple_cluster_valid.f90 $DIFF_PATH
echo "::::::::::::::"
echo ""
echo "::::::::::::::"
diff simple_cmdline.f90 $DIFF_PATH
echo "::::::::::::::"
echo ""
echo "::::::::::::::"
diff simple_comlin_corr.f90 $DIFF_PATH
echo "::::::::::::::"
echo ""
echo "::::::::::::::"
diff simple_comlin.f90 $DIFF_PATH
echo "::::::::::::::"
echo ""
echo "::::::::::::::"
diff simple_comlin_srch.f90 $DIFF_PATH
echo "::::::::::::::"
echo ""
echo "::::::::::::::"
diff simple_comlin_sym.f90 $DIFF_PATH
echo "::::::::::::::"
echo ""
echo "::::::::::::::"
diff simple_comlin_symsrch.f90 $DIFF_PATH
echo "::::::::::::::"
echo ""
echo "::::::::::::::"
diff simple_corrmat.f90 $DIFF_PATH
echo "::::::::::::::"
echo ""
echo "::::::::::::::"
diff simple_ctf.f90 $DIFF_PATH
echo "::::::::::::::"
echo ""
echo "::::::::::::::"
diff simple_de_opt.f90 $DIFF_PATH
echo "::::::::::::::"
echo ""
echo "::::::::::::::"
diff simple_eo_reconstructor.f90 $DIFF_PATH
echo "::::::::::::::"
echo ""
echo "::::::::::::::"
diff simple_filterer.f90 $DIFF_PATH
echo "::::::::::::::"
echo ""
echo "::::::::::::::"
diff simple_ft_expanded.f90 $DIFF_PATH
echo "::::::::::::::"
echo ""
echo "::::::::::::::"
diff simple_ftiter.f90 $DIFF_PATH
echo "::::::::::::::"
echo ""
echo "::::::::::::::"
diff simple_ft_shsrch_devel.f90 $DIFF_PATH
echo "::::::::::::::"
echo ""
echo "::::::::::::::"
diff simple_ft_shsrch.f90 $DIFF_PATH
echo "::::::::::::::"
echo ""
echo "::::::::::::::"
diff simple_gridding.f90 $DIFF_PATH
echo "::::::::::::::"
echo ""
echo "::::::::::::::"
diff simple_hac.f90 $DIFF_PATH
echo "::::::::::::::"
echo "simple_hadamard2D_matcher.f90"
echo "::::::::::::::"
diff simple_hadamard2D_matcher.f90 $DIFF_PATH
echo "::::::::::::::"
echo "simple_hadamard3D_matcher.f90"
echo "::::::::::::::"
diff simple_hadamard3D_matcher.f90 $DIFF_PATH
echo "::::::::::::::"
echo ""
echo "::::::::::::::"
diff simple_hash.f90 $DIFF_PATH
echo "::::::::::::::"
echo ""
echo "::::::::::::::"
diff simple_image.f90 $DIFF_PATH
echo "::::::::::::::"
echo ""
echo "::::::::::::::"
diff simple_imgfile.f90 $DIFF_PATH
echo "::::::::::::::"
echo ""
echo "::::::::::::::"
diff simple_imghead.f90 $DIFF_PATH
echo "::::::::::::::"
echo ""
echo "::::::::::::::"
diff simple_imgheadrec.f90 $DIFF_PATH
echo "::::::::::::::"
echo ""
echo "::::::::::::::"
diff simple_kmeans.f90 $DIFF_PATH
echo "::::::::::::::"
echo ""
echo "::::::::::::::"
diff simple_map_reduce.f90 $DIFF_PATH
echo "::::::::::::::"
echo ""
echo "::::::::::::::"
diff simple_masker.f90 $DIFF_PATH
echo "::::::::::::::"
echo ""
echo "::::::::::::::"
diff simple_matchpursuit.f90 $DIFF_PATH
echo "::::::::::::::"
echo ""
echo "::::::::::::::"
diff simple_nrtxtfile.f90 $DIFF_PATH
echo "::::::::::::::"
echo ""
echo "::::::::::::::"
diff simple_oasis_opt.f90 $DIFF_PATH
echo "::::::::::::::"
echo ""
echo "::::::::::::::"
diff simple_opt_factory.f90 $DIFF_PATH
echo "::::::::::::::"
echo ""
echo "::::::::::::::"
diff simple_optimizer.f90 $DIFF_PATH
echo "::::::::::::::"
echo ""
echo "::::::::::::::"
diff simple_opt_spec.f90 $DIFF_PATH
echo "::::::::::::::"
echo ""
echo "::::::::::::::"
diff simple_opt_subs.f90 $DIFF_PATH
echo "::::::::::::::"
echo ""
echo "::::::::::::::"
diff simple_ori.f90 $DIFF_PATH
echo "::::::::::::::"
echo ""
echo "::::::::::::::"
diff simple_oris.f90 $DIFF_PATH
echo "::::::::::::::"
echo ""
echo "::::::::::::::"
diff simple_pair_dtab.f90 $DIFF_PATH
echo "::::::::::::::"
echo "simple_params.f90"
echo "::::::::::::::"
diff simple_params.f90 $DIFF_PATH
echo "::::::::::::::"
echo ""
echo "::::::::::::::"
diff simple_particle_swarm_opt.f90 $DIFF_PATH
echo "::::::::::::::"
echo ""
echo "::::::::::::::"
diff simple_pftcc_shsrch.f90 $DIFF_PATH
echo "::::::::::::::"
echo ""
echo "::::::::::::::"
diff simple_polarft_corrcalc.f90 $DIFF_PATH
echo "::::::::::::::"
echo ""
echo "::::::::::::::"
diff simple_polarft.f90 $DIFF_PATH
echo "::::::::::::::"
echo ""
echo "::::::::::::::"
diff simple_polarft_shsrch.f90 $DIFF_PATH
echo "::::::::::::::"
echo ""
echo "::::::::::::::"
diff simple_powell_opt.f90 $DIFF_PATH
echo "::::::::::::::"
echo ""
echo "::::::::::::::"
diff simple_ppca.f90 $DIFF_PATH
echo "::::::::::::::"
echo ""
echo "::::::::::::::"
diff simple_prime2D_srch.f90 $DIFF_PATH
echo "::::::::::::::"
echo "simple_prime3D_srch.f90"
echo "::::::::::::::"
diff simple_prime3D_srch.f90 $DIFF_PATH
echo "::::::::::::::"
echo ""
echo "::::::::::::::"
diff simple_procimgfile.f90 $DIFF_PATH
echo "::::::::::::::"
echo ""
echo "::::::::::::::"
diff simple_projector.f90 $DIFF_PATH
echo "::::::::::::::"
echo ""
echo "::::::::::::::"
diff simple_rec_master.f90 $DIFF_PATH
echo "::::::::::::::"
echo ""
echo "::::::::::::::"
diff simple_reconstructor.f90 $DIFF_PATH
echo "::::::::::::::"
echo ""
echo "::::::::::::::"
diff simple_refinecluster.f90 $DIFF_PATH
echo "::::::::::::::"
echo ""
echo "::::::::::::::"
diff simple_sa_opt.f90 $DIFF_PATH
echo "::::::::::::::"
echo ""
echo "::::::::::::::"
diff simple_shalgn.f90 $DIFF_PATH
echo "::::::::::::::"
echo ""
echo "::::::::::::::"
diff simple_shc_cluster.f90 $DIFF_PATH
echo "::::::::::::::"
echo ""
echo "::::::::::::::"
diff simple_simplex_opt.f90 $DIFF_PATH
echo "::::::::::::::"
echo ""
echo "::::::::::::::"
diff simple_sll.f90 $DIFF_PATH
echo "::::::::::::::"
echo ""
echo "::::::::::::::"
diff simple_softmax_weights.f90 $DIFF_PATH
echo "::::::::::::::"
echo ""
echo "::::::::::::::"
diff simple_spatial_median.f90 $DIFF_PATH
echo "::::::::::::::"
echo ""
echo "::::::::::::::"
diff simple_sym.f90 $DIFF_PATH
echo "::::::::::::::"
echo ""
echo "::::::::::::::"
diff simple_unblur_devel.f90 $DIFF_PATH
echo "::::::::::::::"
echo ""
echo "::::::::::::::"
diff simple_unblur.f90 $DIFF_PATH
echo "::::::::::::::"
echo ""
echo "::::::::::::::"
diff simple_winfuns.f90 $DIFF_PATH
