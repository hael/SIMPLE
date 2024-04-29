# CTest - Simple Test Set
# emulate GNU Autotools `make check`
#add_custom_target(check COMMAND ${CMAKE_CTEST_COMMAND})

add_test(
    NAME Test_install 
    COMMAND simple_test_install 
)

add_test(
NAME Test_angres
COMMAND simple_test_angres
)

add_test(
NAME Test_ansi_colors
COMMAND simple_test_ansi_colors
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
NAME Test_cont_inplane
COMMAND simple_test_cont_inplane
)

add_test(
NAME Test_corrs2weights
COMMAND simple_test_corrs2weights
)

add_test(
NAME Test_CTF
COMMAND simple_test_ctf
)

add_test(
NAME Test_Euler
COMMAND simple_test_euler
)

add_test(
NAME Test_extr_frac
COMMAND simple_test_extr_frac
)

add_test(
NAME Test_fileio
COMMAND simple_test_fileio
)

add_test(
NAME Test_io_parallel
COMMAND simple_test_io_parallel
)

add_test(
NAME Test_fplane_proj
COMMAND simple_test_fplane_proj
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
NAME Test_lplims
COMMAND simple_test_lplims
)

add_test(
NAME Test_mask
COMMAND simple_test_mask
)

add_test(
NAME Test_maxnloc
COMMAND simple_test_maxnloc
)

add_test(
NAME Test_memuse
COMMAND simple_test_memuse
)

add_test(
NAME Test_multinomal
COMMAND simple_test_multinomal
)

add_test(
NAME Test_neigh
COMMAND simple_test_neigh
)

add_test(
NAME Test_OpenMP
COMMAND simple_test_openmp
)

add_test(
NAME Test_opt_genetic
COMMAND simple_test_opt_genetic
)

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

add_test(
NAME Test_rad_med
COMMAND simple_test_rad_med
)

add_test(
NAME Test_rank_weights
COMMAND simple_test_rank_weights
)

add_test(
NAME Test_replace_substr
COMMAND simple_test_replace_substr
)

add_test(
NAME Test_serialize
COMMAND simple_test_serialize
)

add_test(
NAME Test_simd
COMMAND simple_test_simd
)

add_test(
NAME Test_srch
COMMAND simple_test_srch
)

add_test(
NAME Test_stack
COMMAND simple_test_stack
)

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
NAME Test_uniform_euler
COMMAND simple_test_uniform_euler
)

add_test(
NAME Test_uniform_rot
COMMAND simple_test_uniform_rot
)

add_test(
NAME Test_units 
COMMAND simple_test_units
)

