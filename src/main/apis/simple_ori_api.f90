module simple_ori_api
!$ use omp_lib
!$ use omp_lib_kinds
use json_kinds
use json_module
use simple_defs
use simple_defs_ori
use simple_type_defs
use simple_ori_utils        ! appropriate to use everything here
use simple_chash,           only: chash
use simple_error,           only: simple_exception
use simple_fileio,          only: fopen, fname2ext, fclose, fileiochk, file_exists
use simple_hash,            only: hash
use simple_is_check_assert, only: is_a_number, is_even, is_even, is_odd, is_equal, is_zero
use simple_linalg,          only: myacos, deg2rad, rad2deg, trace
use simple_math_ft,         only: calc_fourier_index
use simple_math,            only: otsu
use simple_nrtxtfile,       only: nrtxtfile
use simple_ran_tabu,        only: ran_tabu
use simple_rnd,             only: gasdev, ran3, irnd_uni
use simple_srch_sort_loc,   only: hpsort, min3, reverse
use simple_stat,            only: z_scores, median_nocopy, pearsn, moment, robust_scaling
use simple_string_utils,    only: real2str, int2str, uppercase
use simple_string,          only: string
end module simple_ori_api
