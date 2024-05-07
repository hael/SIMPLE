# CTest - Simple Test Set
# emulate GNU Autotools `make check`
#add_custom_target(check COMMAND ${CMAKE_CTEST_COMMAND})

#add_test(
#NAME Test_2D_core_finder
#COMMAND simple_2D_core_finder
#)
#add_test(
#NAME Test_2D_opt_filt
#COMMAND simple_2D_opt_filt
#)
#add_test(
#NAME Test_3D_core_finder
#COMMAND simple_3D_core_finder
#)
#add_test(
#NAME Test_3D_lpbasic
#COMMAND simple_3D_lpbasic
#)
#add_test(
#NAME Test_3D_opt_filt
#COMMAND simple_3D_opt_filt
#)
add_test(
NAME Test_angres
COMMAND simple_test_angres
)
add_test(
NAME Test_ansi_colors
COMMAND simple_test_ansi_colors
)
add_test(
NAME Test_binimage
COMMAND simple_test_binimage
)
add_test(
NAME Test_binoris
COMMAND simple_test_binoris
)
add_test(
NAME Test_bounds_from_mask3D
COMMAND simple_test_bounds_from_mask3D
)
add_test(
NAME Test_cartcorr_shifted
COMMAND simple_test_cartcorr_shifted
)
add_test(
NAME Test_cavg_kpca
COMMAND simple_test_cavg_kpca
)
add_test(
NAME Test_cc_gradient
COMMAND simple_test_cc_gradient
)
add_test(
NAME Test_cont_inplane
COMMAND simple_test_cont_inplane
)
add_test(
NAME Test_corrs2weights
COMMAND simple_test_corrs2weights
)
#add_test(
#NAME Test_corrweights
#COMMAND simple_test_corrweights
#)
add_test(
NAME Test_CTF
COMMAND simple_test_ctf
)
if(USE_CUDA)
    add_test(
    NAME Test_CUDA
    COMMAND simple_test_cuda
    )
endif()
#add_test(
#NAME Test_eo_diff
#COMMAND simple_eo_diff
#)
add_test(
NAME Test_Euler
COMMAND simple_test_euler
)
add_test(
NAME Test_eval_polarftcc
COMMAND simple_test_eval_polarftcc
)
add_test(
NAME Test_extr_frac
COMMAND simple_test_extr_frac
)
add_test(
NAME Test_fast_corrcalc
COMMAND simple_test_fast_corrcalc
)
add_test(
NAME Test_fftw_plan_many_1D
COMMAND simple_test_fftw_plan_many_1D
)
add_test(
NAME Test_fileio
COMMAND simple_test_fileio
)
add_test(
NAME Test_find_boundaries
COMMAND simple_test_find_boundaries
)
add_test(
NAME Test_fit_lattice
COMMAND simple_test_fit_lattice
)
add_test(
NAME Test_fplane_proj
COMMAND simple_test_fplane_proj
)
#add_test(
#NAME Test_ft_expanded
#COMMAND simple_test_ft_expanded
#)
#add_test(
#NAME Test_gauss2Dfit
#COMMAND simple_test_gauss2Dfit
#)
add_test(
NAME Test_gaussian1D
COMMAND simple_test_gaussian1D
)
add_test(
NAME Test_gencorrs_fft
COMMAND simple_test_gencorrs_fft
)
add_test(
NAME Test_graphene_mask
COMMAND simple_test_graphene_mask
)
add_test(
NAME Test_imgfile
COMMAND simple_test_imgfile
)
add_test(
NAME Test_inside_write
COMMAND simple_test_inside_write
)
add_test(
NAME Test_install
COMMAND simple_test_install
)
add_test(
NAME Test_io
COMMAND simple_test_io
)
add_test(
NAME Test_io_parallel
COMMAND simple_test_io_parallel
)
add_test(
NAME Test_lplims
COMMAND simple_test_lplims
)
add_test(
NAME Test_mask
COMMAND simple_test_mask
)
add_test(
NAME Test_masscen_nano
COMMAND simple_test_masscen_nano
)
add_test(
NAME Test_maxnloc
COMMAND simple_test_maxnloc
)
add_test(
NAME Test_memuse
COMMAND simple_test_memuse
)
if(USE_MPI)
   add_test(
   NAME Test_MPI
   COMMAND simple_test_mpi
   )
endif()
add_test(
NAME Test_multinomal
COMMAND simple_test_multinomal
)
#add_test(
#NAME Test_nano_detect_atoms
#COMMAND simple_nano_detect_atoms
#)
add_test(
NAME Test_neigh
COMMAND simple_test_neigh
)
if(USE_OPENACC)
   add_test(
   NAME Test_OpenACC
   COMMAND simple_test_openacc
   )
endif()
add_test(
NAME Test_OpenMP
COMMAND simple_test_openmp
)
add_test(
NAME Test_opt_genetic
COMMAND simple_test_opt_genetic
)
add_test(
NAME Test_opt_lp
COMMAND simple_test_opt_lp
)
#add_test(
#NAME Test_order2D_byresol
#COMMAND simple_test_order2D_byresol
#)
add_test(
NAME Test_order_corr
COMMAND simple_test_order_corr
)
add_test(
NAME Test_ori
COMMAND simple_test_ori
)
add_test(
NAME Test_Otsu
COMMAND simple_test_otsu
)
add_test(
NAME Test_phasecorr
COMMAND simple_test_phasecorr
)
add_test(
NAME Test_phaseplate_correct_FSC
COMMAND simple_test_phaseplate_correct_fsc
)
#add_test(
#NAME Test_picker_comp
#COMMAND simple_test_picker_comp
#)
#add_test(
#NAME Test_polar2cartesian
#COMMAND simple_test_polar2cartesian
#)
add_test(
NAME Test_ppca_kpca
COMMAND simple_test_ppca_kpca
)
add_test(
NAME Test_radial_cc
COMMAND simple_test_radial_cc_
)
add_test(
NAME Test_rad_med
COMMAND simple_test_rad_med
)
add_test(
NAME Test_rank_weights
COMMAND simple_test_rank_weights
)
add_test(
NAME Test_ref_assign
COMMAND simple_test_ref_assign
)
add_test(
NAME Test_replace_substr
COMMAND simple_test_replace_substr
)
add_test(
NAME Test_rotate_ref
COMMAND simple_test_rotate_ref
)
#add_test(
#NAME Test_segpicker
#COMMAND simple_test_segpicker
#)
add_test(
NAME Test_serialize
COMMAND simple_test_serialize
)
#add_test(
#NAME Test_shiftsrch
#COMMAND simple_test_shiftsrch
#)
add_test(
NAME Test_simd
COMMAND simple_test_simd
)
#add_test(
#NAME Test_speed
#COMMAND simple_test_speed
#)
add_test(
NAME Test_sp_project
COMMAND simple_test_sp_project
)
add_test(
NAME Test_srch
COMMAND simple_test_srch
)
add_test(
NAME Test_stack
COMMAND simple_test_stack
)
#add_test(
#NAME Test_stack_io
#COMMAND simple_test_stack_io
#)
add_test(
NAME Test_starfile
COMMAND simple_test_starfile
)
add_test(
NAME Test_stringmatch
COMMAND simple_test_stringmatch
)
add_test(
NAME Test_sym
COMMAND simple_test_sym
)
add_test(
NAME Test_tseries_neigh
COMMAND simple_test_tseries_neigh
)
add_test(
NAME Test_uniform_cluster
COMMAND simple_test_uniform_cluster
)
add_test(
NAME Test_uniform_euler
COMMAND simple_test_uniform_euler
)
#add_test(
#NAME Test_uniform_filter_2D
#COMMAND simple_test_uniform_filter_2D
#)
#add_test(
#NAME Test_uniform_filter_3D
#COMMAND simple_test_uniform_filter_3D
#)
add_test(
NAME Test_uniform_rot
COMMAND simple_test_uniform_rot
)
add_test(
NAME Test_units 
COMMAND simple_test_units
)
