module simple_core_module_api

! DEFINE API
use simple_binoris,      only: binoris
use simple_chash,        only: chash
use simple_cmdline,      only: cmdline, cmdline_err
use simple_defs,         only: STDLEN, XLONGSTRLEN, logfhandle, C_FALSE, C_TRUE, SP, DP, MAXS, JPEG_DIM, PI, PIO2, TWOPI, TINY,&
                         &SMALL, CMPLX_ZERO, DCMPLX_ZERO, GRAPHENE_BAND1, GRAPHENE_BAND2
use simple_error,        only: simple_exception
use simple_fileio,       only: fopen, fclose, basename, file_exists, fname_new_ext, del_file, del_files, simple_list_files,&
                              &move_files2dir, append2basename
use simple_hash,         only: hash
use simple_imghead,      only: find_ldim_nptcls, update_stack_nimgs
use simple_linalg,       only: euclid
use simple_math_ft,      only: csq_fast, csq, calc_fourier_index, calc_graphene_mask, calc_lowpass_lim, cyci_1d, cyci_1d_static,&
                              &fdim, get_find_at_res, get_find_at_crit, get_resarr, mycabs, phase_angle, resang
use simple_ori,          only: ori
use simple_oris,         only: oris
use simple_rnd,          only: ran3, seed_rnd, irnd_uni
use simple_stat,         only: avg_sdev, moment, skewness, kurtosis, pearsn, normalize, normalize_minmax, merge_dmats,&
                              &avg_frac_smallest, pearsn_serial, std_mean_diff, calc_stats, corrs2weights, kstwo, analyze_smat,&
                              &dmat2smat, smat2dmat, scores2scores_percen, dists2scores_percen, merge_smats, medoid_from_smat,&
                              &medoid_from_dmat, median, median_nocopy
use simple_string,       only: string
use simple_string_utils, only: str_has_substr, str2real, real2str, int2str
use simple_sym,          only: sym
use simple_syslib,       only: is_open, syslib_symlink, exec_cmdline, simple_mkdir
use simple_timer,        only: timer_int_kind, tic, toc
use simple_type_defs,    only: ENUM_CTFFLAG, CTFFLAG_NO, CTFFLAG_YES, CTFFLAG_FLIP, ENUM_OBJFUN, OBJFUN_CC, OBJFUN_EUCLID,&
                              &class_sample, CORRW_CRIT, ctfparams

use, intrinsic :: iso_fortran_env, only: output_unit

! PUBLICIZE API
! binoris
public :: binoris
! chash
public :: chash
! cmdline
public :: cmdline, cmdline_err
! defs
public :: STDLEN, XLONGSTRLEN, logfhandle, C_FALSE, C_TRUE, SP, DP, MAXS, JPEG_DIM, PI, PIO2, TWOPI, TINY,&
         &SMALL, CMPLX_ZERO, DCMPLX_ZERO
! error
public :: simple_exception
! fileio
public :: fopen, fclose, basename, file_exists, fname_new_ext, del_file, del_files, simple_list_files,&
         &move_files2dir, append2basename
! hash
public :: hash
! imghead
public :: find_ldim_nptcls, update_stack_nimgs
! linalg
public :: euclid
! math_ft
public :: csq_fast, csq, calc_fourier_index, calc_graphene_mask, calc_lowpass_lim, cyci_1d, cyci_1d_static,&
         &fdim, get_find_at_res, get_find_at_crit, get_resarr, mycabs, phase_angle, resang
! ori
public :: ori
! oris
public :: oris
! rnd
public :: ran3, seed_rnd, irnd_uni
! stat
public :: avg_sdev, moment, skewness, kurtosis, pearsn, normalize, normalize_minmax, merge_dmats,&
         &avg_frac_smallest, pearsn_serial, std_mean_diff, calc_stats, corrs2weights, kstwo, analyze_smat,&
         &dmat2smat, smat2dmat, scores2scores_percen, dists2scores_percen, merge_smats, medoid_from_smat,&
         &medoid_from_dmat, median, median_nocopy
! string
public :: string
! string_utils
public :: str_has_substr, str2real, real2str, int2str
! sym
public :: sym
! syslib
public :: is_open, syslib_symlink, exec_cmdline, simple_mkdir
! timer
public :: timer_int_kind, tic, toc
! type_defs
public :: ENUM_CTFFLAG, CTFFLAG_NO, CTFFLAG_YES, CTFFLAG_FLIP, ENUM_OBJFUN, OBJFUN_CC, OBJFUN_EUCLID,&
         &class_sample, CORRW_CRIT, ctfparams

! intrisic :: iso_fortran_env
public :: output_unit

end module simple_core_module_api