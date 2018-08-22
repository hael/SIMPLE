module simple_defs
use, intrinsic :: iso_c_binding
use, intrinsic :: iso_fortran_env
use simple_defs_fname
implicit none
integer,  parameter :: ascii        = selected_char_kind ("ascii")
integer,  parameter :: ucs4         = selected_char_kind ('ISO_10646')
integer,  parameter :: MAXS         = 99   !< maximum number of states
integer,  parameter :: short        = selected_int_kind(4)
integer,  parameter :: long         = selected_int_kind(9)
integer,  parameter :: longer       = selected_int_kind(16)
integer,  parameter :: I4B          = selected_int_kind(9)
integer,  parameter :: I2B          = selected_int_kind(4)
integer,  parameter :: I1B          = selected_int_kind(2)
integer,  parameter :: SP           = kind(1.0)
integer,  parameter :: DP           = selected_real_kind(8) !kind(1.0d0)
integer,  parameter :: DOUBLE       = kind(1.0d0)
integer,  parameter :: SPC          = kind((1.0,1.0))
integer,  parameter :: DPC          = kind((1.0d0,1.0d0))
integer,  parameter :: LGT          = kind(.true.)
integer,  parameter :: LINE_MAX_LEN = 8192
real(sp), parameter :: PI           = acos(-1.)
real(dp), parameter :: DPI          = acos(-1.d0)
real(sp), parameter :: PIO2         = acos(-1.)/2.
real(sp), parameter :: TWOPI        = 2.*acos(-1.)
real(dp), parameter :: DTWOPI       = 2.d0*acos(-1.d0)
real(sp), parameter :: FOURPI       = 4.*acos(-1.)
real(sp), parameter :: SQRT2        = sqrt(2.)
real(sp), parameter :: EUL          = 0.5772156649015328606065120900824024310422_sp
real(sp), parameter :: TINY         = 1e-10
real(dp), parameter :: DTINY        = 1e-10
real(sp), parameter :: SMALL        = 1e-6
real(sp), parameter :: FTOL         = 1e-4
real(dp), parameter :: DSMALL       = 1e-6
real(dp), parameter :: PISQR        = DPI*DPI

! directory-based execution model
character(len=:), allocatable :: CWD_GLOB_ORIGINAL, CWD_GLOB

! other global variables
integer(kind=c_int)           :: nthr_glob         !< number of threads global variable
logical                       :: l_distr_exec_glob !< global distributed execution flag
integer                       :: part_glob         !< global part index
character(len=:), allocatable :: cmdline_glob      !< global command line string

! plan for the CTF
type :: ctfplan
    character(len=STDLEN) :: mode=''                !< astig/noastig
    character(len=STDLEN) :: flag=''                !< flag: <mul|flip|no>
    logical               :: l_phaseplate = .false. !< image obtained with Volta phaseplate
end type ctfplan

! type for arrays of allocatable strings
type str4arr
    character(len=:), allocatable :: str
end type str4arr

! CTF flag type
enum, bind(c)
    enumerator :: CTFFLAG_NO = 0, CTFFLAG_YES = 1,  CTFFLAG_FLIP = 2
end enum

! Objective function
enum, bind(c)
    enumerator :: OBJFUN_CC = 0, OBJFUN_RES = 1, OBJFUN_EUCLID = 2
end enum

! type for CTF parameters
type :: ctfparams
    integer :: ctfflag = 0
    real    :: smpd    = 0.
    real    :: kv      = 0.
    real    :: cs      = 0.
    real    :: fraca   = 0.
    real    :: dfx     = 0.
    real    :: dfy     = 0.
    real    :: angast  = 0.
    real    :: phshift = 0.
    logical :: l_phaseplate = .false. !< image obtained with Volta phaseplate
end type ctfparams

! character constants
character(len=*), parameter :: NEWLINE = new_line('a')
character(len=*), parameter :: SUPPRESS_MSG='2>/dev/null'

! constants for picker
real,    parameter :: PICKER_SHRINK        = 4.        !< picker shrink factor
real,    parameter :: PICKER_SHRINK_REFINE = 2.        !< picker shrink factor, peak refine step
real,    parameter :: GAUPICK_SIGMA_SHRINK = 13.0      !< picker shrink factor, Gaussian convolution
integer, parameter :: PICKER_OFFSET        = 3         !< picker offset for grid search

! constants for masking/interpolation
real, parameter :: COSMSKHALFWIDTH         = 6.0       !< spherical soft masking
real, parameter :: KBWINSZ                 = 1.5       !< interpolation window size for 2D
real, parameter :: KBALPHA                 = sqrt(2.0) !< interpolation alpha (oversampling constant)
real, parameter :: RECWINSZ                = 1.5       !< half-window size for 3D reconstruction

! real constants that control search and convergence
real, parameter    :: FRAC_SH_LIM          = 80.0      !< at what frac to turn on the shift search
real, parameter    :: EXTRINITHRESH        = 0.5       !< initial randomization threshold for extremal search
real, parameter    :: EXTRTHRESH_CONST     = 0.2       !< threshold for factorial decay in extremal search
real, parameter    :: SNHC2D_INITFRAC      = 0.5       !< initial neighbourhood fraction for 2D SNHC
real, parameter    :: SNHC2D_DECAY         = 0.2       !< factorial decay in 2D SNHC
real, parameter    :: GREEDY_FREQ          = 0.2       !< frequency of greedy search in refine3D with refine=single
real, parameter    :: LP2SMPDFAC           = 0.4125    !< low-pass limit scaling constant
real, parameter    :: LP2SMPDFAC2D         = 0.4       !< low-pass limit scaling constant
real, parameter    :: NPEAKSATHRES         = 12.0      !< angular threshold for determining npeaks (PRIME3D)
real, parameter    :: SHC_INPL_TRSHWDTH    = 2.0       !< shift search halfwidht (pixels)
real, parameter    :: TAU_DEFAULT          = 0.01      !< controls the sharpeness of the orientation weight distribution
                                                       !! smaller number means sharper distribution
integer, parameter :: MAX_EXTRLIM2D        = 15        !< maximum # of iterations for which 2D extremal opt is performed
real,    parameter :: SOFTMAXW_THRESH      = 1.0       !< threshold for orientations softmax weights, # sigmas to the right of mean
real,    parameter :: BSC3D                = 20.       !< for shell-weighted 3D reconstruction (shellw), used in B-factor calculation
real,    parameter :: BSC2D                = 15.       !< for shell-weighted 2D reconstruction (shellw), used in B-factor calculation
real,    parameter :: HP_CORR_VALID        = 20.       !< high-pass limit for validation corr calculation
                                                       !! 20 A is where domain structure starts
real,    parameter :: LP_CORR_VALID        = 8.        !< low-pass limit for validation corr calculation
                                                       !! signal is strong out to A with DDDs and 2ndary structure appears here

! integer #/threshold constants
integer, parameter :: LPLIM1ITERBOUND      = 5         !< # iteration bound lplim stage 1 (PRIME2D)
integer, parameter :: LPLIM3ITERBOUND      = 7         !< # iteration bound lplim stage 2 (PRIME2D)
integer, parameter :: MINCLSPOPLIM         = 5         !< limit for adaptive cluster splitting/spreading (PRIME2D)
integer, parameter :: GRIDNPEAKS           = 3         !< # peaks to consider in angular grid search (PRIME3D)
integer, parameter :: CONTNPEAKS           = 5         !< # peaks to refine continuously
integer, parameter :: NPEAKS2REFINE        = 200       !< # peaks to be further optimised
integer, parameter :: NINPLPEAKS2SORT      = 5         !< maximum # in-plane peaks to be considered for sorting
integer, parameter :: MINNPEAKS            = 6         !< minimum # of peaks in soft assignment
integer, parameter :: NSPACE_REDUCED       = 600       !< # projection directions for the balancing constraint (PRIME3D)
integer, parameter :: GRIDCORR_MAXITS      = 5         !< # iterations for reconstruction gridding correction
integer, parameter :: MAXIMGBATCHSZ        = 500       !< max # images in batch

! criterion for even/odd averaging in gold-FSC
real,    parameter :: FSC4EOAVG3D = 0.9                !< corr criterium for eo-averaging in 3D
real,    parameter :: FSC4EOAVG2D = 0.7                !< corr criterium for eo-averaging in 2D
integer, parameter :: K4EOAVGLB   = 4                  !< Fourier index lower-bound

! SNHC-related global constants, PRIME3D, refine=snhc
integer, parameter :: SZSN_INIT  = 5
integer, parameter :: SZSN_STEP  = 3
integer, parameter :: SZSN_MAX   = 20

! computer related
integer, parameter :: JOB_MEMORY_PER_TASK_DEFAULT = 16000
integer, parameter :: TIME_PER_IMAGE_DEFAULT      = 100

! precision constants
#ifndef IMAGE_SINGLE_PRECISION
integer, parameter :: img_kind = DP
#else
integer, parameter :: img_kind = SP
#endif
integer, parameter :: fp_kind = DP

! debugging and print verbosity flags
#ifdef _DEBUG
logical :: global_debug   = .true.  !< global debugging flag
logical :: global_verbose = .true.  !< global flag for verbosity set to TRUE in debug mode
#else
logical :: global_debug   = .false. !< global flag for debugging disabled
#ifdef VERBOSE
logical :: global_verbose = .true.  !< global flag for verbosity TRUE with VERBOSE compilation flag
#else
logical :: global_verbose = .false. !< global flag for verbosity FALSE by default
#endif
#endif
logical :: global_warn    = .false. !< warning flag
integer :: alloc_stat
! append SIMPLE_VERSION and SIMPLE_GIT_VERSION strings to simple_defs
#include "SimpleGitVersion.h"

end module simple_defs
