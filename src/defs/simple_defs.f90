! module with global variables
! DEBUG and VERBOSE modes can be enabled at compile-time or run-time
! BUILD_DESC uses the cmake build to describe <GNU/Intel/PGI>_<RELEASE/DEBUG>_<FFTW/MKL>
! SimpleGitVersion include file contains the string definitions for SIMPLE_PATH and version info
module simple_defs
use, intrinsic :: iso_c_binding
implicit none
character(len=1), parameter :: default_file_format = 'M' !<  'I', 'M' or 'S' for imagic, mrc, spider
integer, parameter  :: IMPORTANT=10 !< number of solutions considered important
integer, parameter  :: MAXS=99      !< maximum number of states
integer, parameter  :: STDLEN     = 256  !< standard string length
integer, parameter  :: LONGSTRLEN = 2048 !< longer string length
integer, parameter  :: short = selected_int_kind(4)
integer, parameter  :: long  = selected_int_kind(9)
integer, parameter  :: longer  = selected_int_kind(16)
integer, parameter  :: I4B = SELECTED_INT_KIND(9)
integer, parameter  :: I2B = SELECTED_INT_KIND(4)
integer, parameter  :: I1B = SELECTED_INT_KIND(2)
integer, parameter  :: SP = KIND(1.0)        !< single float precision
integer, parameter  :: DP = KIND(1.0D0)      !< double float precision
integer, parameter  :: DOUBLE = KIND(1.0D0)
integer, parameter  :: SPC = KIND((1.0,1.0))
integer, parameter  :: DPC = KIND((1.0D0,1.0D0))
integer, parameter  :: LGT = KIND(.true.)
integer, parameter  :: line_max_len = 8192 !< Max number of characters on line
real(sp), parameter :: PI=acos(-1.)
real(dp), parameter :: DPI=acos(-1.d0)
real(sp), parameter :: PIO2=acos(-1.)/2.
real(sp), parameter :: TWOPI=2.*acos(-1.)
real(sp), parameter :: DTWOPI=2.d0*acos(-1.d0)
real(sp), parameter :: FOURPI=4.*acos(-1.)
real(sp), parameter :: SQRT2=sqrt(2.)
real(sp), parameter :: EUL=0.5772156649015328606065120900824024310422_sp
real(sp), parameter :: um2a = 10000.
real(sp), parameter :: TINY=1e-10
real(dp), parameter :: DTINY=1e-10
real(sp), parameter :: SMALL=1e-6
real(sp), parameter :: MINEULSDEV=3.
real(sp), parameter :: MINTRSSDEV=0.5
real(sp), parameter :: FTOL=1e-4
real(dp), parameter :: DSMALL=1e-6
real(dp), parameter :: pisqr = PI*PI   ! PI^2.
real(sp), parameter :: ATHRES_LIM = 5.

! constants for masks
real, parameter :: COSMSKHALFWIDTH = 3.0

! plan for the CTF
type :: ctfplan
    character(len=STDLEN) :: mode='' !< astig/noastig
    character(len=STDLEN) :: flag='' !< flag: <mul|flip|no>
end type ctfplan

! constants for picker
real,    parameter :: PICKER_SHRINK        = 4.
real,    parameter :: PICKER_SHRINK_REFINE = 2.
integer, parameter :: PICKER_OFFSET        = 3

! constants for interpolation
real, parameter :: KBWINSZ = 1.5    !< interpolation window size
real, parameter :: KBALPHA = 2.0    !< interpolation alpha (smoothing constant)

integer, parameter :: SPECWMINPOP=2000 !< minimum population for spectral weighting

! SNHC-related global vars
character(len=32), parameter :: SNHCDOC = 'snhc_oris.txt'
character(len=32), parameter :: SNHCVOL = 'snhc_recvol_state'
integer,           parameter :: SZSN_INIT = 5
integer,           parameter :: SZSN_STEP = 3
integer,           parameter :: SZSN_MAX  = 20

! constants that control search and convergence
real,    parameter :: FRAC_SH_LIM      = 80.0 !< at what frac to turn on the shift search
real,    parameter :: FRAC_INTERPOL    = 60.0 ! at what frac to turn on the gridding interpolation (2D)
real,    parameter :: EXTRINITHRESH    = 0.5
real,    parameter :: EXTRTHRESH_CONST = 0.2
real,    parameter :: LP2SMPDFAC       = 0.4125
real,    parameter :: NPEAKSATHRES     = 12.0
integer, parameter :: LPLIM1ITERBOUND  = 5
integer, parameter :: LPLIM3ITERBOUND  = 7
integer, parameter :: GRIDNPEAKS       = 3
integer, parameter :: MAXNPEAKS        = 40

character(len=:), allocatable :: endconv           !< endianness conversion
integer(kind=c_int)           :: nthr_glob         !< number of threads global variable
logical                       :: l_distr_exec_glob !< global distributed execution flag
character(len=STDLEN)         :: exec_abspath_glob !< global executable absolute path
character(len=LONGSTRLEN)     :: cmdline_glob      !< global command line string

#ifndef IMAGE_SINGLE_PRECISION
  integer, parameter :: img_kind = DP
#else
  integer, parameter :: img_kind = SP
#endif
integer, parameter :: fp_kind = DP

! Debugging and print verbosity flags
#ifdef _DEBUG
  logical :: global_debug=.true.          !< global debugging flag
  logical :: global_verbose=.true.        !< global flag for verbosity set to TRUE in debug mode
#else
  logical :: global_debug=.false.         !< global flag for debugging disabled
#ifdef VERBOSE
  logical :: global_verbose=.true.        !< global flag for verbosity TRUE with VERBOSE compilation flag
#else
  logical :: global_verbose=.false.       !< global flag for verbosity FALSE by default
#endif
#endif
  logical :: global_warn=.false.          !< warning flag

! append SIMPLE_VERSION and SIMPLE_GIT_VERSION strings to simple_defs
#include "SimpleGitVersion.h"

end module simple_defs
