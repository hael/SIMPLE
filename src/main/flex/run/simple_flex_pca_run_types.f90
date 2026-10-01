!@descr: flex_pca run settings: every process-wide switch of a run, resolved once
!!
!! The one place a flex_pca run reads its environment. `flex_run_settings` is an immutable
!! projection of facts the typed parameters object does not retain (whether a key was given on
!! the command line) together with the SIMPLE_COV_* overrides, each read exactly once by `new`.
!! Every consumer receives the resolved value; no numerical module reads the environment or the
!! command line itself. Defaults are the values the readers used before the switches were
!! centralised, so an environment without any SIMPLE_COV_* key resolves to the same run.
module simple_flex_pca_run_types
use simple_core_module_api
use simple_parameters,        only: parameters
use simple_cmdline,           only: cmdline
use simple_reconstructor,     only: reconstructor
use simple_flex_pca_records,  only: flex_selection, flex_fit_model, flex_latent, flex_state_set, flex_latent_readout
implicit none
private
#include "simple_local_flags.inc"

public :: flex_run_settings, flex_run_session

!> stage particle caps: 0 = off (every particle); a positive SIMPLE_COV_PROBE_MAX / SIMPLE_COV_CALIB_MAX
!! overrides. The probe cap is off because capping traded a recovered conformation for speed
!! (doc/policies/flex_pca_policy.md); the calibration pass estimates two scalars whose precision
!! improves as 1/sqrt(N), so 20k particles suffice.
integer,  parameter :: COV_PROBE_MAX_PTCLS = 0
integer,  parameter :: COV_CALIB_MAX_PTCLS = 20000

type :: flex_run_settings
    ! ---- command-line facts the parameters object does not keep ----
    logical :: l_npreimages_explicit = .false. !< npreimages was given (the auto ceiling may not raise it)
    logical :: l_resume              = .false. !< infile was given: states-only resume from an embedding cache
    logical :: l_vol1_explicit       = .false. !< vol1 was given (or picked up from the project by the commander)
    logical :: l_pindfile            = .false. !< pindfile was given (a distributed worker's partition)
    ! ---- state-map delivery (PCG state backend) ----
    logical :: l_state_eofilt = .false.        !< SIMPLE_COV_STATE_EOFILT=1: per-state eo-FSC optimal filter
    logical :: l_state_filt   = .true.         !< SIMPLE_COV_STATE_FILT=0: no filter on the delivered maps
    ! ---- two-gate state merge ----
    logical  :: l_merge              = .false. !< SIMPLE_COV_MERGE (0 always wins), else on under preimage_auto
    logical  :: l_merge_r_set        = .false. !< SIMPLE_COV_MERGE_R present: fixed gate, adaptive cut off
    logical  :: l_merge_r_value      = .false. !< ... and it parsed: merge_r carries the gate
    real(dp) :: merge_r              = 0.d0
    real(dp) :: merge_margin         = 0.d0    !< SIMPLE_COV_MERGE_MARGIN
    logical  :: l_merge_eo           = .false. !< SIMPLE_COV_MERGE_EO=1
    logical  :: l_merge_madk_set     = .false. !< SIMPLE_COV_MERGE_MADK given
    real(dp) :: merge_madk           = 0.d0
    real(dp) :: merge_lp             = 0.d0    !< SIMPLE_COV_MERGE_LP (A); 0 = no band cap
    logical  :: l_merge_link_single  = .false. !< SIMPLE_COV_MERGE_LINK=single
    real(dp) :: auto_k               = 0.d0    !< SIMPLE_COV_AUTO_K (>0: gate 1 reports, never folds)
    ! ---- basis composition from finished runs ----
    logical :: l_compose     = .false.         !< SIMPLE_COV_COMPOSE set
    character(len=:), allocatable :: compose_dirs !< its value: comma-separated run directories
    logical  :: l_compose_cut = .true.         !< SIMPLE_COV_COMPOSE_CUT=0 opts out of the signal-subspace cut
    real(dp) :: cut_snr       = 1.d0           !< SIMPLE_COV_CUT_SNR
    ! ---- fit policy ----
    integer :: mod4_pairing = 1                !< SIMPLE_COV_MOD4_PAIRING (1 or 3; validated by the paired driver)
    integer :: probe_max    = COV_PROBE_MAX_PTCLS !< probe-stage particle cap, total across processes
    integer :: calib_max    = COV_CALIB_MAX_PTCLS !< noise/prior calibration particle cap
    logical :: l_deflate_bg       = .true.     !< SIMPLE_COV_DEFLATE_BG=0 opts out
    logical :: l_deflate_dilation = .true.     !< SIMPLE_COV_DEFLATE_DILATION=0 opts out
    ! ---- coupled PCG M-step ----
    logical :: l_pcg_lambda_set = .false.      !< SIMPLE_COV_PCG_LAMBDA given
    real    :: pcg_lambda_rel   = 0.0
    integer :: pcg_verbose      = 0            !< SIMPLE_COV_PCG_VERBOSE
  contains
    procedure :: new  => settings_new
    procedure :: kill => settings_kill
end type flex_run_settings

!> Everything a run owns between its phases. Allocated by the phases as they always were; freed
!! once, in one place, by `kill` (safe on a partially built session: every free is guarded).
type :: flex_run_session
    type(flex_run_settings)   :: cfg
    type(flex_selection)      :: sel      !< the particle selection
    type(flex_fit_model)      :: model    !< mean, basis, prior variances, rank, noise level
    type(flex_latent)         :: latent   !< the embedding and its statistics
    type(flex_state_set)      :: states   !< the delivered state set
    type(flex_latent_readout) :: readout  !< UMAP readout coordinates for the state figure
    logical :: sigma_loaded = .false., l_resume = .false., l_compose = .false., l_paired_states = .false.
    integer, allocatable :: deconv_labels(:)
    logical :: l_deconv_applied = .false., l_deconv_adopted = .false.
    integer :: min_neff = 0, state_axis = 0, col_sep = 1, neigs_req = 0, nkern = 0
    real(dp), allocatable :: pviews(:,:)
    logical :: l_pop_floor = .false., l_merged = .false., l_state_rec = .true.
  contains
    procedure :: clamp_state_axis => session_clamp_state_axis
    procedure :: kill => session_kill
end type flex_run_session

contains

    !> Free everything the session may hold: model handles first (a resume never built them), then
    !! the embedding, the state table, the readout and the settings.
    !> the state axis can never exceed the rank or the kernel dimension; applied after every rank change
    subroutine session_clamp_state_axis( self )
        class(flex_run_session), intent(inout) :: self
        if( self%state_axis > 0 ) self%state_axis = min(self%state_axis, min(self%model%ncomp, self%nkern))
    end subroutine session_clamp_state_axis

    subroutine session_kill( self )
        class(flex_run_session), intent(inout) :: self
        call self%sel%kill
        call self%model%kill
        call self%latent%kill
        call self%states%kill
        call self%readout%kill
        if( allocated(self%deconv_labels) ) deallocate(self%deconv_labels)
        if( allocated(self%pviews) )        deallocate(self%pviews)
        call self%cfg%kill
    end subroutine session_kill

    !> Resolve every switch of the run. Reads the command line for the facts `parameters` drops
    !! and the environment for the SIMPLE_COV_* overrides, once.
    subroutine settings_new( self, params, cline )
        class(flex_run_settings), intent(inout) :: self
        class(parameters),        intent(in)    :: params
        class(cmdline),           intent(in)    :: cline
        character(len=XLONGSTRLEN) :: envc
        integer  :: ival, envlen
        real(dp) :: dval
        real     :: rval
        logical  :: l_set
        call self%kill
        self%l_npreimages_explicit = cline%defined('npreimages')
        self%l_resume              = cline%defined('infile')
        self%l_vol1_explicit       = cline%defined('vol1')
        self%l_pindfile            = cline%defined('pindfile')
        ! state-map delivery
        self%l_state_eofilt = env_is('SIMPLE_COV_STATE_EOFILT', '1')
        self%l_state_filt   = .not. env_is('SIMPLE_COV_STATE_FILT', '0')
        ! merge: an explicit SIMPLE_COV_MERGE wins (0 = off); unset, the auto state ceiling turns it on,
        ! since a ceiling without the collapse is just a large state count
        call env_get('SIMPLE_COV_MERGE', envc, envlen, l_set)
        if( l_set )then
            self%l_merge = trim(adjustl(envc(:envlen))) /= '0'
        else
            self%l_merge = params%l_preimage_auto
        endif
        self%l_merge_r_set = env_present('SIMPLE_COV_MERGE_R')
        call env_dp_silent('SIMPLE_COV_MERGE_R',      self%merge_r,      self%l_merge_r_value)
        call env_dp_silent('SIMPLE_COV_MERGE_MARGIN', self%merge_margin, l_set)
        self%l_merge_eo = env_is('SIMPLE_COV_MERGE_EO', '1')
        call env_dp_silent('SIMPLE_COV_MERGE_MADK',   self%merge_madk,   self%l_merge_madk_set)
        call env_dp_silent('SIMPLE_COV_MERGE_LP',     self%merge_lp,     l_set)
        self%l_merge_link_single = env_is('SIMPLE_COV_MERGE_LINK', 'single')
        call env_dp_silent('SIMPLE_COV_AUTO_K',       self%auto_k,       l_set)
        ! composition
        call env_get('SIMPLE_COV_COMPOSE', envc, envlen, self%l_compose)
        if( self%l_compose ) self%compose_dirs = envc(:envlen)
        self%l_compose_cut = .not. env_is('SIMPLE_COV_COMPOSE_CUT', '0')
        dval = self%cut_snr
        call env_dp('SIMPLE_COV_CUT_SNR', dval)
        self%cut_snr = dval
        ! fit policy (integer overrides take values > 0 only, as cov_env_int always did)
        ival = self%mod4_pairing; call env_int('SIMPLE_COV_MOD4_PAIRING', ival); self%mod4_pairing = ival
        ival = self%probe_max;    call env_int('SIMPLE_COV_PROBE_MAX',    ival); self%probe_max    = ival
        ival = self%calib_max;    call env_int('SIMPLE_COV_CALIB_MAX',    ival); self%calib_max    = ival
        self%l_deflate_bg       = .not. env_is('SIMPLE_COV_DEFLATE_BG',       '0')
        self%l_deflate_dilation = .not. env_is('SIMPLE_COV_DEFLATE_DILATION', '0')
        ! PCG
        rval = 0.0
        call env_real_nonneg('SIMPLE_COV_PCG_LAMBDA', rval, self%l_pcg_lambda_set)
        if( self%l_pcg_lambda_set ) self%pcg_lambda_rel = rval
        ival = self%pcg_verbose;  call env_int('SIMPLE_COV_PCG_VERBOSE',  ival); self%pcg_verbose  = ival
    end subroutine settings_new

    subroutine settings_kill( self )
        class(flex_run_settings), intent(inout) :: self
        if( allocated(self%compose_dirs) ) deallocate(self%compose_dirs)
        self%l_compose = .false.
    end subroutine settings_kill

    ! ---- environment readers (the only ones in the flex tree) ----

    !> raw value of a variable; l_set is .true. when it is present and non-empty
    subroutine env_get( name, val, ln, l_set )
        character(len=*), intent(in)  :: name
        character(len=*), intent(out) :: val
        integer,          intent(out) :: ln
        logical,          intent(out) :: l_set
        integer :: stat
        val = ''
        call get_environment_variable(name, val, ln, stat)
        l_set = stat == 0 .and. ln >= 1
        if( .not. l_set ) ln = 0
    end subroutine env_get

    !> .true. when the variable is present and non-empty, whatever its value
    logical function env_present( name )
        character(len=*), intent(in) :: name
        character(len=32) :: envval
        integer :: ln
        call env_get(name, envval, ln, env_present)
    end function env_present

    !> .true. when the variable is set to exactly `want` (trimmed, case-sensitive)
    logical function env_is( name, want ) result( yes )
        character(len=*), intent(in) :: name, want
        character(len=32) :: envval
        integer :: ln
        logical :: l_set
        yes = .false.
        call env_get(name, envval, ln, l_set)
        if( .not. l_set ) return
        yes = trim(adjustl(envval)) == want
    end function env_is

    !> integer override, values > 0 only; logs the override
    subroutine env_int( name, val )
        character(len=*), intent(in)    :: name
        integer,          intent(inout) :: val
        character(len=32) :: envval
        integer :: ln, stat, ival
        logical :: l_set
        call env_get(name, envval, ln, l_set)
        if( .not. l_set ) return
        read(envval(:ln), *, iostat=stat) ival
        if( stat == 0 .and. ival > 0 )then
            val = ival
            write(logfhandle,'(A,A,A,I0)') '>>> FLEX_PCA ',trim(name),' override: ',ival
            call flush(logfhandle)
        endif
    end subroutine env_int

    !> double override; logs the override
    subroutine env_dp( name, val )
        character(len=*), intent(in)    :: name
        real(dp),         intent(inout) :: val
        character(len=32) :: envval
        integer  :: ln, stat
        real(dp) :: rval
        logical  :: l_set
        call env_get(name, envval, ln, l_set)
        if( .not. l_set ) return
        read(envval(:ln), *, iostat=stat) rval
        if( stat == 0 )then
            val = rval
            write(logfhandle,'(A,A,A,ES12.4)') '>>> FLEX_PCA ',trim(name),' override: ',rval
            call flush(logfhandle)
        endif
    end subroutine env_dp

    !> double override without a log line (the merge gates report their own values); l_parsed
    !! is .true. only when the value was present and readable
    subroutine env_dp_silent( name, val, l_parsed )
        character(len=*), intent(in)    :: name
        real(dp),         intent(inout) :: val
        logical,          intent(out)   :: l_parsed
        character(len=32) :: envval
        integer  :: ln, stat
        real(dp) :: rval
        logical  :: l_set
        l_parsed = .false.
        call env_get(name, envval, ln, l_set)
        if( .not. l_set ) return
        read(envval, *, iostat=stat) rval
        if( stat == 0 )then
            val      = rval
            l_parsed = .true.
        endif
    end subroutine env_dp_silent

    !> finite non-negative real override (the PCG lambda reader)
    subroutine env_real_nonneg( name, val, l_set )
        use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
        character(len=*), intent(in)    :: name
        real,             intent(inout) :: val
        logical,          intent(out)   :: l_set
        character(len=64) :: envval
        integer :: ln, stat
        real    :: v
        call env_get(name, envval, ln, l_set)
        if( .not. l_set ) return
        read(envval(1:ln), *, iostat=stat) v
        l_set = .false.
        if( stat /= 0 ) return
        if( ieee_is_finite(v) .and. v >= 0.0 )then
            val   = v
            l_set = .true.
        endif
    end subroutine env_real_nonneg

end module simple_flex_pca_run_types
