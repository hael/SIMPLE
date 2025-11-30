module simple_defs
use, intrinsic :: iso_fortran_env
use, intrinsic :: iso_c_binding
use simple_defs_fname
use simple_defs_environment
implicit none
integer,     parameter :: MAXS         = 99 !< maximum number of states
integer,     parameter :: short        = selected_int_kind(4)
integer,     parameter :: long         = selected_int_kind(9)
integer,     parameter :: longer       = selected_int_kind(16)
integer,     parameter :: I4B          = selected_int_kind(9)
integer,     parameter :: I2B          = selected_int_kind(4)
integer,     parameter :: I1B          = selected_int_kind(2)
integer,     parameter :: SP           = kind(1.0)
integer,     parameter :: DP           = selected_real_kind(8) ! kind(1.0d0)
integer,     parameter :: QP           = selected_real_kind(4*precision(1.0_sp))
integer,     parameter :: DOUBLE       = kind(1.0d0)
integer,     parameter :: SPC          = kind((1.0,1.0))
integer,     parameter :: DPC          = kind((1.0d0,1.0d0))
integer,     parameter :: LGT          = kind(.true.)
integer,     parameter :: LINE_MAX_LEN = 8192
integer,     parameter :: JPEG_DIM     = 100                 ! dimension of jpeg thumbnails
real(sp),    parameter :: PI           = acos(-1.)
real(dp),    parameter :: DPI          = acos(-1.d0)
real(sp),    parameter :: PIO2         = acos(-1.)/2.
real(dp),    parameter :: DPIO2        = acos(-1.d0)/2.d0
real(sp),    parameter :: TWOPI        = 2.*acos(-1.)
real(dp),    parameter :: DTWOPI       = 2.d0*acos(-1.d0)
real(sp),    parameter :: FOURPI       = 4.*acos(-1.)
real(sp),    parameter :: SQRT2        = sqrt(2.)
real(sp),    parameter :: EUL          = 0.5772156649015328606065120900824024310422_sp
real(sp),    parameter :: TINY         = 1e-10
real(dp),    parameter :: DTINY        = 1e-10
real(sp),    parameter :: SMALL        = 1e-6
real(sp),    parameter :: FTOL         = 1e-4
real(dp),    parameter :: DSMALL       = 1.d-6
real(dp),    parameter :: PISQR        = DPI*DPI
complex(sp), parameter :: CMPLX_ZERO   = cmplx(0.,0.)
complex(dp), parameter :: DCMPLX_ZERO  = cmplx(0.d0,0.d0, kind=dp)

! directory-based execution model
character(len=:), allocatable :: CWD_GLOB_ORIG, CWD_GLOB
character(len=*), parameter   :: LOGFNAME   = 'simple.log'     !< log file name
integer                       :: logfhandle = OUTPUT_UNIT      !< log file handle, default to STDOUT
logical, parameter            :: STDOUT2LOG = .false.

! other global variables
integer                       :: nthr_glob = 1                 !< number of threads global variable
logical                       :: l_distr_exec_glob             !< global distributed execution flag
integer                       :: part_glob                     !< global part index
character(len=:), allocatable :: cmdline_glob                  !< global command line string
integer,          parameter   :: NTHR_SHMEM_MAX     = 20       !< maximum number of shared-memory threads used by master process
logical,          parameter   :: L_BENCH_GLOB       = .false.  !< global benchmarking flag
logical,          parameter   :: L_DO_GRIDCORR_GLOB = .false.  !< global gridding correction flag
logical,          parameter   :: L_USE_SLURM_ARR    = .false.  !< use SLURM arrays for jobs where we know nparts
logical,          parameter   :: L_USE_AUTO_MEM     = .false.  !< auto estmate memory usage for parts
logical,          parameter   :: L_DEV_GLOB         = .false.  !< global development flag
logical,          parameter   :: L_VERBOSE_GLOB     = .false.  !< verbose output or not
real,             parameter   :: HPLIM_GUINIER      = 20.      !< high-pass limit for Guinier plot

! CTF flag type
enum, bind(c)
    enumerator :: ENUM_CTFFLAG = -1
    enumerator :: CTFFLAG_NO = 0, CTFFLAG_YES = 1, CTFFLAG_FLIP = 2
end enum

! CTF treatment flag type
enum, bind(c)
    enumerator :: ENUM_CTFLIMFLAG = -1
    enumerator :: CTFLIMFLAG_FULL = 0, CTFLIMFLAG_PI = 1, CTFLIMFLAG_PIO2 = 2
end enum

! Objective function type
enum, bind(c)
    enumerator :: ENUM_OBJFUN= -1
    enumerator :: OBJFUN_CC = 0, OBJFUN_EUCLID = 1
end enum

! type for CTF parameters
type ctfparams
    integer(kind(ENUM_CTFFLAG)) :: ctfflag = CTFFLAG_YES
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
    real :: maxv = 0.
    real :: minv = 0.
end type stats_struct

type inpl_struct
    real    :: e3 = 0., x = 0., y = 0., corr = 0.
    logical :: l_mirr = .false.
    integer :: find_fsc05 = 0
end type inpl_struct

type clust_inpl
    type(inpl_struct), allocatable :: params(:)
end type clust_inpl

 type clust_info
    integer          :: nptcls      = 0
    integer          :: pop         = 0
    integer          :: good_bad    = 0
    integer          :: scoreclust  = 0
    type(clust_inpl) :: algninfo
    real             :: res         = 0.
    real             :: euclid      = 0.
    real             :: homogeneity = 0.
    real             :: resscore    = 0.
    real             :: clustscore  = 0.
    real             :: jointscore  = 0.
    real             :: icescore    = 0.
end type clust_info

! type for particle reference relation in eul_prob_tab
type ptcl_ref
    integer :: pind = 0, iproj = 0, inpl = 0, istate=0
    real    :: dist = 0., x = 0., y = 0.
    logical :: has_sh = .false.
end type ptcl_ref

type lp_crop_inf
    real    :: lp=0., smpd_crop=0., scale=1., trslim=0., frc_crit=0.
    integer :: box_crop=0
    logical :: l_autoscale=.false., l_lpset=.false.
end type lp_crop_inf

type class_sample
    integer :: clsind = 0, pop = 0, nsample = 0
    integer, allocatable :: pinds(:)
    real,    allocatable :: ccs(:)
end type class_sample

! oritype enumeration
enum, bind(c)
    enumerator :: ENUM_ORISEG  = 0
    enumerator :: MIC_SEG      = 1,  STK_SEG     = 2,  PTCL2D_SEG  = 3, CLS2D_SEG  = 4
    enumerator :: CLS3D_SEG    = 5,  PTCL3D_SEG  = 6,  OUT_SEG     = 7, OPTICS_SEG = 8
    enumerator :: PROJINFO_SEG = 11, JOBPROC_SEG = 12, COMPENV_SEG = 13
end enum
integer(kind=kind(ENUM_ORISEG)), parameter :: GENERIC_SEG = PTCL3D_SEG

! weighting methods enumeration
enum, bind(c)
    enumerator :: ENUM_WCRIT        = 0
    enumerator :: RANK_SUM_CRIT     = 1
    enumerator :: RANK_CEN_CRIT     = 2
    enumerator :: RANK_EXP_CRIT     = 3
    enumerator :: RANK_INV_CRIT     = 4
    enumerator :: CORRW_CRIT        = 5
    enumerator :: CORRW_ZSCORE_CRIT = 6
end enum

! export (to STAR) type enumeration
enum, bind(c)
    enumerator :: ENUM_STARTYPE = 0
    enumerator :: MOV_STAR = 1, MIC_STAR = 2, STK_STAR = 3, PTCL_STAR = 4, CLSAVG_STAR = 5, OUT_STAR = 6
    enumerator :: PROJINFO_STAR = 7
end enum
integer(kind=kind(ENUM_STARTYPE)), parameter :: GENERIC_STAR = PTCL_STAR

! general parameters
real,    parameter :: PRUNE_FRAC                = 0.3            !< fraction of particles after which a project is automatically pruned
integer, parameter :: BUFSZ_DEFAULT             = 1024           !< Default stack_io buffer size

! power spectrum related stuff
integer, parameter :: GUI_PSPECSZ               = 512            !< hard-coded image size for gui
real,    parameter :: SMPD4VIZ                  = 1.25           !< default sampling distance for powerspectrum visualisation
real,    parameter :: SMPD4DOWNSCALE            = 1.3            !< default sampling distance for downscaling in motion correction 
real,    parameter :: LP_PSPEC_BACKGR_SUBTR     = 20.            !< default low-pass limit for power spectrum background subtraction

! constants for picker & extraction
real,    parameter :: PICKER_SHRINK             = 4.             !< picker shrink factor
real,    parameter :: PICKER_SHRINK_REFINE      = 2.             !< picker shrink factor, peak refine step
real,    parameter :: GAUPICK_SIGMA_SHRINK      = 13.0           !< picker shrink factor, Gaussian convolution
real,    parameter :: RADFRAC_NORM_EXTRACT      = 0.75           !< radius fraction for extraction background normalization
integer, parameter :: PICKER_OFFSET             = 3              !< picker offset for grid search

! constants for masking/interpolation
real, parameter    :: COSMSKHALFWIDTH           = 6.0            !< spherical soft masking
real, parameter    :: KBWINSZ                   = 1.5            !< interpolation window size for 2D
real, parameter    :: KBALPHA                   = 2.             !< interpolation alpha (oversampling constant), previously sqrt(2)
real, parameter    :: RECWINSZ                  = 1.5            !< half-window size for 3D reconstruction

! real constants that control search and convergence
real, parameter    :: FRAC_SH_LIM               = 75.0           !< at what frac to turn on the shift search
real, parameter    :: NEIGH_MINFRAC             = 0.3            !< minimum fraction of search space scanned in refine=neigh
real, parameter    :: FRAC_GREEDY_LIM           = 99.0           !< at what frac to turn to greedy search
real, parameter    :: EXTRINITHRES              = 0.5            !< initial randomization threshold for extremal search
real, parameter    :: EXTRTHRESH_CONST          = 0.2            !< threshold for factorial decay in extremal search
real, parameter    :: SNHC2D_INITFRAC           = 0.5            !< initial neighbourhood fraction for 2D SNHC
real, parameter    :: SNHC2D_DECAY              = 0.2            !< factorial decay in 2D SNHC
real, parameter    :: GREEDY_FREQ               = 0.1            !< frequency of greedy search in refine3D
real, parameter    :: GLOB_FREQ                 = 0.2            !< frequency of global stoachastic search in  with refine=neigh
real, parameter    :: LP2SMPDFAC                = 0.4125         !< low-pass limit scaling constant
real, parameter    :: LP2SMPDFAC2D              = 0.4            !< low-pass limit scaling constant
real, parameter    :: SHC_INPL_TRSHWDTH         = 2.0            !< shift search halfwidht (pixels)ch
real, parameter    :: MC_PATCHSZ                = 200.           !< recommended patch size (in Angstroms) for motion correction
real, parameter    :: ENVMSK_FSC_THRESH         = 0.8            !< FSC value after which phase-randomization and FSC correction is applied in enveloppe masking
real, parameter    :: MAX_SMPD                  = 2.67           !< maximum sampling distance in scaling
real, parameter    :: TAU_DEFAULT               = 1.0            !< TAU fudge factor to control strength or regularization [0.5,5] more -> less low-pass effect
real, parameter    :: CENTHRESH                 = 0.5            ! threshold for performing volume/cavg centering in pixels
real, parameter    :: MAXCENTHRESH2D            = 3.0            ! max threshold for performing cavg centering in pixels
real, parameter    :: EXTR_POWER                = 2.0            ! Exponent of the sampling function during extremal stochastic phase of 2D analysis
real, parameter    :: POST_EXTR_POWER           = 4.0            ! Exponent of the sampling function after the extremal stochastic phase
integer, parameter :: MAXPOP_CLS                = 5000
integer, parameter :: MAXPOP_PTCLS              = 1200000

! preprocessing constants
real, parameter    :: FRACTION_DOSE_TARGET_DEFAULT=1.0           !< EER target fraction dose in e/A2
real, parameter    :: DFMAX_DEFAULT             = 5.0            !< Default maximum bound for defocus search (microns)
real, parameter    :: DFMIN_DEFAULT             = 0.2            !< Default minimum bound for defocus search (microns)
real, parameter    :: LP_CTF_ESTIMATE           = 5.0            !< Default low-pass limit for defocus search (Angstroms)
real, parameter    :: HP_CTF_ESTIMATE           = 30.0           !< Default high-pass limit for defocus search (Angstroms)
real, parameter    :: HP_BACKGR_SUBTR           = 400.0          !< High-pass frequency for micrograph background subtraction (Angstroms)
real, parameter    :: CTFRES_THRESHOLD          = 50.0           !< Ctfres rejection threshold (Angstroms)
real, parameter    :: ICEFRAC_THRESHOLD         = 1.0            !< Icefrac rejection threshold
real, parameter    :: ASTIG_THRESHOLD           = 10.0           !< Astigmatism rejection threshold
real, parameter    :: PICK_LP_DEFAULT           = 20.            !< Picking resolution limit
real, parameter    :: BOX_EXP_FAC               = 1.0            !< Multiplication factor box = (moldiam + moldiam * BOXFAC)/smpd
real, parameter    :: MSK_EXP_FAC               = 1.2            !< Multiplication factor mskdiam  = smpd * box_for_pick * MSK_EXP_FAC

! integer #/threshold constants
integer, parameter :: LPLIM1ITERBOUND           = 5              !< # iteration bound lplim stage 1 (PRIME2D)
integer, parameter :: LPLIM3ITERBOUND           = 7              !< # iteration bound lplim stage 2 (PRIME2D)
integer, parameter :: MINCLSPOPLIM              = 5              !< limit for adaptive cluster splitting/spreading (PRIME2D)
integer, parameter :: GRIDCORR_MAXITS           = 2              !< # iterations for reconstruction gridding correction
integer, parameter :: MAXIMGBATCHSZ             = 500            !< max # images in batch
integer, parameter :: MAX_EXTRLIM2D             = 15             !< maximum # of iterations for which 2D extremal opt is performed
integer, parameter :: NPEAKS_DEFAULT            = 3              !< # of greedy subspace peaks to construct multi-neighborhood search spaces from
integer, parameter :: NPEAKS_INPL_DEFAULT       = 10             !< # neighborhood search peaks to refine with L-BFGS
integer, parameter :: NSAMPLE_MINMAX_DEFAULT(2) = [10000,25000]  !< default minimum and maximum particle sampling size
integer, parameter :: MC_NPATCH                 = 5              !< number of patches in x/y-direction for motion correction
integer, parameter :: MC_MINPATCHSZ             = 200            !< Minimum patch size in pixels for motion correction
integer, parameter :: MIN_ITERS_SHC             = 5              !< minimum number of iterations of stochastic search
integer, parameter :: BATCHTHRSZ                = 50             !< # of images per thread
integer, parameter :: AMSK_FREQ                 = 3              !< automasking every third iteration

! stream-related constants & thresholds
real,    parameter :: CTFRES_THRESHOLD_STREAM   = 10.0           !< preprocessing: Stream ctfres rejection threshold (Angstroms)
real,    parameter :: ICEFRAC_THRESHOLD_STREAM  = 1.0            !< preprocessing: Stream icefrac rejection threshold
real,    parameter :: ASTIG_THRESHOLD_STREAM    = 10.0           !< preprocessing: Stream astigmatism rejection threshold
real,    parameter :: RES_THRESHOLD_STREAM      = 35.0           !< class rejection: Default streaming best resolution rejection threshold
real,    parameter :: LOWRES_REJECT_THRESHOLD   = 199.           !< class rejection: Deactivates resolution-based rejection when lpthres > LOWRES_REJECT_THRESHOLD
real,    parameter :: CLS_REJECT_STD            = 2.5            !< class rejection: # deviations for 2D class selection/rejection
real,    parameter :: MEAN_THRESHOLD            = -8.0           !< class rejection: image mean     threshold (Mahalabonis distance)
real,    parameter :: REL_VAR_THRESHOLD         = 6.0            !< class rejection: image variance threshold (Mahalabonis distance)
real,    parameter :: ABS_VAR_THRESHOLD         = 1.5            !< class rejection: image variance threshold (absolute value)
real,    parameter :: TVD_THRESHOLD             = 0.55           !< class rejection: Total Variation Distance of image distributions
real,    parameter :: MINMAX_THRESHOLD          = 2.0            !< class rejection: image min & max threshold (absolute value)
real,    parameter :: FRAC_SKIP_REJECTION       = 0.7            !< 2D analysis: When the number of classes to reject is too high rejection is skipped
integer, parameter :: STREAM_SRCHLIM            = 5              !< 2D analysis: maximum # of systematic iterations for streaming 2D pool
integer, parameter :: MAX_STREAM_NPTCLS         = 500000         !< 2D analysis: cap for adjusting update_frac in 2D streaming
integer, parameter :: STREAM_NMOVS_SET          = 5              !< number of movies processed at once (>1)
integer, parameter :: STREAM_NMOVS_SET_TIFF     = 3              !< number of TIFF movies processed at once (>1)

! nanoparticles
real,    parameter :: AMSKLP_NANO               = 5.

! Graphene
real, parameter    :: GRAPHENE_BAND1            = 2.14           !< graphene band 1 for omission in score function
real, parameter    :: GRAPHENE_BAND2            = 1.23           !< graphene band 2 for omission in score function

! Ice
real, parameter    :: ICE_BAND1                 = 3.7
real, parameter    :: ICE_BAND2                 = 1.23

! criterion for even/odd averaging in gold-FSC
real,    parameter :: FREQ4EOAVG3D              = 20.            !< Frequencry criterion for eo-averaging in 3D
real,    parameter :: FSC4EOAVG3D               = 0.95           !< corr criterion for eo-averaging in 3D
real,    parameter :: FSC4EOAVG2D               = 0.7            !< corr criterion for eo-averaging in 2D
integer, parameter :: K4EOAVGLB                 = 4              !< Fourier index lower-bound

! qsys related
integer, parameter :: QSYS_SUBMISSION_RETRY_LIMIT = 5
integer, parameter :: QSYS_SUBMISSION_RETRY_SLEEP = 5
integer, parameter :: QSYS_SUBMISSION_RETRY_MULTI = 3

! computer related
integer, parameter :: JOB_MEMORY_PER_TASK_DEFAULT = 16000
integer, parameter :: TIME_PER_IMAGE_DEFAULT      = 100          !< seconds
integer, parameter :: WALLTIME_DEFAULT            = 172740       !< seconds, 47h59mins

! C-compatible boolean constants
logical(c_bool), parameter :: C_FALSE = logical(.false.,kind=c_bool)
logical(c_bool), parameter :: C_TRUE  = logical(.true. ,kind=c_bool)

! append SIMPLE_VERSION and SIMPLE_GIT_VERSION strings to simple_defs
#include "SimpleGitVersion.h"
end module simple_defs
