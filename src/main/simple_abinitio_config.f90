!@descr: singleton for common state variables used in abinitio_utils and commanders_abinitio
module simple_abinitio_config
use simple_type_defs, only: lp_crop_inf
use simple_sym,       only: sym
use simple_cmdline,   only: cmdline
implicit none

! singleton constants
character(len=*), parameter :: REC_FBODY             = 'rec_final_state'
character(len=*), parameter :: STR_STATE_GLOB        = '01'
real,             parameter :: LPSTOP_BOUNDS(2)      = [4.5,6.0]
real,             parameter :: LPSTART_BOUNDS(2)     = [10.,20.]
real,             parameter :: CENLP_DEFAULT         = 30.
real,             parameter :: LPSYMSRCH_LB          = 12.
real,             parameter :: UPDATE_FRAC_MAX       = 0.9                  ! to ensure fractional update is always on
real,             parameter :: UPDATE_FRAC_MIN       = 0.1                  ! 10% of the particles updated each iteration
real,             parameter :: LPSTART_INI3D         = 20.                  ! Default lpstart for abinitio3D_cavgs/cavgs_ini
real,             parameter :: LPSTOP_INI3D          = 8.                   ! Default lpstop for abinitio3D_cavgs/cavgs_ini
integer,          parameter :: NSTAGES               = 8
integer,          parameter :: NSTAGES_INI3D         = 4 ! # of ini3D stages used for initialization
integer,          parameter :: NSTAGES_INI3D_MAX     = 7
integer,          parameter :: PHASES(3)             = [2,6,NSTAGES]
integer,          parameter :: MAXITS(8)             = [20,20,17,17,17,17,15,30]
integer,          parameter :: MAXITS_GLOB           = SUM(MAXITS(1:7))  ! the last 30 iterations are not included in this estimate since the sampling method changes
integer,          parameter :: NSPACE(3)             = [500,1000,2500]
integer,          parameter :: SYMSRCH_STAGE         = 3
integer,          parameter :: PROBREFINE_STAGE      = 5
integer,          parameter :: ICM_STAGE             = PROBREFINE_STAGE     ! we switch from ML regularization when prob is switched on
integer,          parameter :: STOCH_SAMPL_STAGE     = PROBREFINE_STAGE     ! we switch from greedy to stochastic blanced class sampling when prob is switched on
integer,          parameter :: TRAILREC_STAGE_SINGLE = STOCH_SAMPL_STAGE    ! we start trailing when we start sampling particles randomly
integer,          parameter :: TRAILREC_STAGE_MULTI  = NSTAGES              ! we start trailing in the last stage
integer,          parameter :: LPAUTO_STAGE          = NSTAGES - 1          ! cannot be switched on too early
integer,          parameter :: AUTOMSK_STAGE         = LPAUTO_STAGE         ! swith on automasking when lpauto is switched on
integer,          parameter :: HET_DOCKED_STAGE      = NSTAGES              ! stage at which state splitting is done when multivol_mode==docked
integer,          parameter :: STREAM_ANALYSIS_STAGE = 5                    ! when streaming on some analysis will be performed
integer,          parameter :: CAVGWEIGHTS_STAGE     = 3                    ! when to activate optional cavg weighing in abinitio3D_cavgs/cavgs_fast
integer,          parameter :: GAUREF_LAST_STAGE     = PHASES(1)            ! When to stop using gaussian filtering of the references with polar=yes
integer,          parameter :: MAXITS_BETWEEN        = 10                   ! Development

! singleton variables
type(lp_crop_inf), allocatable :: lpinfo(:)
logical          :: l_srch4symaxis=.false., l_symran=.false., l_sym=.false., l_update_frac_dyn=.false., l_polar=.false.
logical          :: l_ini3D=.false., l_lpauto=.false., l_nsample_given=.false., l_nsample_stop_given=.false., l_automsk=.false.
type(sym)        :: se1, se2
type(cmdline)    :: cline_refine3D, cline_symmap, cline_reconstruct3D, cline_postprocess, cline_reproject
real             :: update_frac  = 1.0, update_frac_dyn  = 1.0
integer          :: nstates_glob = 1, nptcls_eff = 0, nsample_minmax(2), maxits_dyn=0

end module simple_abinitio_config