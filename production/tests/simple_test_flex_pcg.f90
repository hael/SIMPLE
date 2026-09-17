!@descr: validates the flex_pca PCG M-step operator: the pair Gram kernel on the 2x lattice against the
!  direct KB gather/scatter Gram of a synthetic sample set (scale and residual)
program simple_test_flex_pcg
use simple_core_module_api
use simple_flex_pca_pcg, only: test_flex_pcg_operator
implicit none
#include "simple_local_flags.inc"
logical :: l_pass1, l_pass2, l_dbg
! debug: single samples first (profiles), then the random sets
call test_flex_pcg_operator(32, 8, l_dbg, loc_fixed=[0.0, 0.0, 0.0])
call test_flex_pcg_operator(32, 8, l_dbg, loc_fixed=[6.0, 0.0, 0.0])
call test_flex_pcg_operator(32, 8, l_dbg, loc_fixed=[3.3, 2.1, -1.7])
call test_flex_pcg_operator(32, 200, l_pass1)
call test_flex_pcg_operator(64, 400, l_pass2)
if( .not. (l_pass1 .and. l_pass2) ) THROW_HARD('flex PCG operator test failed')
call simple_end('**** SIMPLE_TEST_FLEX_PCG NORMAL STOP ****')
end program simple_test_flex_pcg
