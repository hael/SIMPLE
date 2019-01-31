module simple_defs
use, intrinsic :: iso_fortran_env
use, intrinsic :: iso_c_binding
use simple_defs_fname
implicit none
integer,  parameter :: MAXS         = 99 !< maximum number of states
integer,  parameter :: short        = selected_int_kind(4)
integer,  parameter :: long         = selected_int_kind(9)
integer,  parameter :: longer       = selected_int_kind(16)
integer,  parameter :: I4B          = selected_int_kind(9)
integer,  parameter :: I2B          = selected_int_kind(4)
integer,  parameter :: I1B          = selected_int_kind(2)
integer,  parameter :: SP           = kind(1.0)
integer,  parameter :: DP           = selected_real_kind(8) ! kind(1.0d0)
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
character(len=:), allocatable :: cwd_glob_orig, cwd_glob
character(len=*), parameter   :: LOGFNAME   = 'simple.log' !< log file name
integer                       :: logfhandle = OUTPUT_UNIT  !< log file handle, default to STDOUT
logical, parameter            :: STDOUT2LOG = .false.

! other global variables
integer                       :: nthr_glob         !< number of threads global variable
logical                       :: l_distr_exec_glob !< global distributed execution flag
integer                       :: part_glob         !< global part index
character(len=:), allocatable :: cmdline_glob      !< global command line string

! type for arrays of allocatable strings
type str4arr
    character(len=:), allocatable :: str
end type str4arr

! CTF flag type
enum, bind(c)
    enumerator :: ENUM_CTFFLAG = -1
    enumerator :: CTFFLAG_NO = 0, CTFFLAG_YES = 1,  CTFFLAG_FLIP = 2
end enum

! Objective function type
enum, bind(c)
    enumerator :: ENUM_OBJFUN= -1
    enumerator :: OBJFUN_CC = 0, OBJFUN_RES = 1, OBJFUN_EUCLID = 2
end enum

! type for CTF parameters
type ctfparams
    integer(kind(ENUM_CTFFLAG)) :: ctfflag = CTFFLAG_NO
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

type stats_struct
    real :: avg  = 0.
    real :: med  = 0.
    real :: sdev = 0.
    real :: maxv = 0
    real :: minv = 0.
end type stats_struct

! oritype enumeration
enum, bind(c)
    enumerator :: ENUM_ORISEG = 0
    enumerator :: MIC_SEG = 1, STK_SEG = 2, PTCL2D_SEG = 3, CLS2D_SEG = 4
    enumerator :: CLS3D_SEG = 5, PTCL3D_SEG = 6, OUT_SEG = 7
    enumerator :: PROJINFO_SEG = 11, JOBPROC_SEG = 12, COMPENV_SEG = 13
end enum
integer(kind=kind(ENUM_ORISEG)), parameter :: GENERIC_SEG = PTCL3D_SEG

! export (to STAR) type enumeration
enum, bind(c)
    enumerator :: ENUM_STARTYPE = 0
    enumerator :: MOV_STAR = 1, MIC_STAR = 2, STK_STAR = 3, PTCL_STAR = 4, CLSAVG_STAR = 5, OUT_STAR = 6
    enumerator :: PROJINFO_STAR = 7
end enum
integer(kind=kind(ENUM_STARTYPE)), parameter :: GENERIC_STAR = PTCL_STAR

! power spectrum related stuff
real,    parameter :: SMPD4VIZ              = 1.25    !< default sampling distance for powerspectrum visualisation
real,    parameter :: LP_PSPEC_BACKGR_SUBTR = 20.     !< default low-pass limit for power spectrum background subtraction

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
real, parameter    :: TAU_DEFAULT          = 0.01      !< controls the sharpeness of the orientation weight distribution when objfun .ne. euclid
                                                       !! smaller number means sharper distribution
real, parameter    :: SIGMA2_FUDGE_DEFAULT = 100.      !< controls the sharpeness of the orientation weight distribution when objfun .eq. euclid
                                                       !! smaller number means sharper distribution, but with sigma2_fudge == 1 it turns into a delta
integer, parameter :: MAX_EXTRLIM2D        = 15        !< maximum # of iterations for which 2D extremal opt is performed
real,    parameter :: SOFTMAXW_THRESH      = 1.0       !< threshold for orientations softmax weights, # sigmas to the right of mean

! integer #/threshold constants
integer, parameter :: LPLIM1ITERBOUND      = 5         !< # iteration bound lplim stage 1 (PRIME2D)
integer, parameter :: LPLIM3ITERBOUND      = 7         !< # iteration bound lplim stage 2 (PRIME2D)
integer, parameter :: MINCLSPOPLIM         = 5         !< limit for adaptive cluster splitting/spreading (PRIME2D)
integer, parameter :: GRIDNPEAKS           = 3         !< # peaks to consider in angular grid search (PRIME3D)
integer, parameter :: CONTNPEAKS           = 5         !< # peaks to refine continuously
integer, parameter :: NPROJPEAKS           = 40        !< minimum # projection directions evaluated in more random stochastic search
integer, parameter :: NPEAKS2REFINE        = 200       !< # peaks to be further optimised
integer, parameter :: NINPLPEAKS2SORT      = 5         !< maximum # in-plane peaks to be considered for sorting
integer, parameter :: MINNPEAKS            = 6         !< minimum # of peaks in soft assignment
integer, parameter :: NSPACE_REDUCED       = 600       !< # projection directions for the balancing constraint (PRIME3D)
integer, parameter :: GRIDCORR_MAXITS      = 2         !< # iterations for reconstruction gridding correction
integer, parameter :: MAXIMGBATCHSZ        = 500       !< max # images in batch
integer, parameter :: RANDOMNESS_FAC       = 3         !< controls randomness of stochastic search, 1 is most random, 6 is least

! weighting scheme
logical, parameter :: SOFT_PTCL_WEIGHTS    = .true.    !< soft (as opposed to hard) particles weights based on B-factor
logical, parameter :: WEIGHT_SCHEME_GLOBAL = .false.    !< per-particle or global orientation weighting scheme
real,    parameter :: GLOBAL_WEIGHT_FRAC   = 0.16      !< corresponds to threshold of mean + one sigma
real,    parameter :: BFAC_SDEV            = 50.       !< B-factor standard deviation for shell-weighting in 2D/3D
! real,    parameter :: BFAC_LBOUND        = -50.      !< lower B-factor bound for shell-weighting in 2D/3D based on specscore


! criterion for even/odd averaging in gold-FSC
real,    parameter :: FSC4EOAVG3D = 0.95               !< corr criterium for eo-averaging in 3D
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

!! Control exit status
!! STOP is standard, CALL EXIT is an extension. STOP lets you display a text message where EXIT
!! doesn't. Both let you set an exit status value. Otherwise they are pretty much the same.
!! use call exit(EXIT_FAILURE)
integer, parameter :: EXIT_SUCCESS = 0
integer, parameter :: EXIT_FAILURE = 1
integer, parameter :: EXIT_FAILURE2 = 2
integer, parameter :: EXIT_FAILURE3 = 3
integer, parameter :: EXIT_FAILURE4 = 4
integer, parameter :: EXIT_FAILURE5 = 5
integer, parameter :: EXIT_FAILURE6 = 6
integer, parameter :: EXIT_FAILURE7 = 7

logical :: mir_projns = .false.    ! flag indicating if mirrored projection should be computed simultaneously
                                   ! (remove this flag once integration is done)

! append SIMPLE_VERSION and SIMPLE_GIT_VERSION strings to simple_defs
#include "SimpleGitVersion.h"
end module simple_defs
