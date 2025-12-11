module simple_core_module_api

! DEFINE API
use simple_chash,          only: chash
use simple_cmdline,        only: cmdline, cmdline_err
use simple_defs,           only: STDLEN, XLONGSTRLEN, logfhandle, class_sample, C_FALSE, C_TRUE, SP, DP
use simple_error,          only: simple_exception
use simple_fileio,         only: fopen, fclose, basename, file_exists, fname_new_ext, del_file, del_files
use simple_hash,           only: hash
use simple_linalg,         only: euclid
use simple_ori,            only: ori
use simple_oris,           only: oris
use simple_rnd,            only: ran3, seed_rnd
use simple_string,         only: string
use simple_string_utils,   only: str_has_substr, str2real, real2str, int2str
use simple_syslib,         only: is_open, syslib_symlink, exec_cmdline, simple_mkdir
use simple_timer,          only: timer_int_kind, tic, toc

use, intrinsic :: iso_fortran_env, only: output_unit

! PUBLICIZE API
! chash
public :: chash
! cmdline
public :: cmdline, cmdline_err
! defs
public :: STDLEN, XLONGSTRLEN, logfhandle, class_sample, C_FALSE, C_TRUE, SP, DP
! error
public :: simple_exception
! fileio
public :: fopen, fclose, basename, file_exists, fname_new_ext, del_file, del_files
! hash
public :: hash
! linalg
public :: euclid
! ori
public :: ori
! oris
public :: oris
! rnd
public :: ran3, seed_rnd
! string
public :: string
! string_utils
public :: str_has_substr, str2real, real2str, int2str
! syslib
public :: is_open, syslib_symlink, exec_cmdline, simple_mkdir
! timer
public :: timer_int_kind, tic, toc

! intrisic :: iso_fortran_env
public :: output_unit

end module simple_core_module_api