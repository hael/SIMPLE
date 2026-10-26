!@descr: reconstruct3D_wave: half maps under a single-scatter wave model, solved by nonlinear conjugate gradients from the flat PCG maps
! (the program: doc/policies/3D/wave_model_policy.md, section 6)
module simple_rec3D_wave_strategy
use simple_defs,                only: dp, MAXIMGBATCHSZ, OSMPL_PAD_FAC, logfhandle
use simple_type_defs,           only: ctfparams, OBJFUN_EUCLID
use simple_error,               only: simple_exception
use simple_string,              only: string
use simple_string_utils,        only: int2str_pad
use simple_fileio,              only: arr2file, fopen, fclose
use simple_timer,               only: timer_int_kind
use simple_parameters,          only: parameters
use simple_builder,             only: builder
use simple_image,               only: image
use simple_imgarr_utils,        only: alloc_imgarr, dealloc_imgarr
use simple_oris,                only: oris
use simple_math_ctf,            only: electron_scale_link
use simple_eo_partition,        only: ensure_field_halves
use simple_sigma2_files,        only: load_sigma2
use simple_nlcg_solver,         only: nlcg_objective, nlcg_solver_options, nlcg_solver_outcome, nlcg_solve
use simple_wave_imaging,        only: WAVE_PROJECTION, WAVE_FIRST_BORN, WAVE_RYTOV, WAVE_MSLICE
use simple_wave_projector,      only: wave_projector
use simple_reconstructor_pcg,   only: reconstructor_pcg, pcg_solver_outcome, PCG_LAMBDA, PCG_OP_KERNEL, &
    &handle_cold_restart_outcome, validate_solved_map, report_pcg_solve
use simple_matcher_ptcl_io,     only: prepimgbatch, killimgbatch, ptcl_read_session, rec_state_half, rec_obs_oris, &
    &rec_obs_inv_sig2, read_rec_obs_batch, ptcl_buffer_budget, fill_pass_buffer, release_pass_buffer, &
    &pass_buffer_status_name
use simple_rec3D_pcg_strategy,  only: set_pcg_solve_support, calculate_pcg_state_diagnostics, accumulate_pcg_state_half
use simple_halfmap_diagnostics, only: halfmap_diagnostics_result, write_halfmap_diagnostics
use simple_refine3D_fnames,     only: refine3D_state_halfvol_fname, refine3D_state_vol_fname, refine3D_fsc_fname, &
    &refine3D_resolution_txt_fbody, refine3D_state_vol_suffix_fname
!$ use omp_lib,                 only: omp_get_max_threads
implicit none

public :: exec_rec3D_wave, rec3D_wave_objective
private
#include "simple_local_flags.inc"

! the step subset of the line search (curvature and parabola), as specter's 5%, at least STEP_SUBSET_MIN
! particles and at most STEP_SUBSET_MAX (its planes stay in memory for the solve)
real,    parameter :: STEP_SUBSET_FRAC = 0.05
integer, parameter :: STEP_SUBSET_MIN  = 256
integer, parameter :: STEP_SUBSET_MAX  = 1000

!> the raw density of a half's flat operator, kept from end_accum for its regularization
type :: raw_density
    real, allocatable :: d(:,:,:)
end type raw_density

!> One half's objective in volts, L(V) = 1/2 sum w |y - p(V)|^2 + 1/2 kappa^2 V^T P V / boxpd^3 (P the flat PCG
!! operator's ML prior, kappa the scale link); preconditioner and support the PCG operator's, curvature and line
!! losses on the step subset (wave_model_policy.md, sections 5 and 6)
type, extends(nlcg_objective) :: rec3D_wave_objective
    private
    class(parameters),        pointer :: params => null()
    class(builder),           pointer :: build  => null()
    class(reconstructor_pcg), pointer :: pcgop  => null()
    type(wave_projector),    allocatable :: wp
    type(ptcl_read_session), allocatable :: rsess
    type(image), allocatable :: obs(:)          !< per-thread observation work images
    integer, allocatable :: pinds(:)            !< the half's particles, the step subset first
    complex, allocatable :: ybatch(:,:,:)       !< one batch of observation planes
    complex, allocatable :: ysub(:,:,:)         !< the step subset's observation planes
    real,    allocatable :: supp(:,:,:)         !< the solve support, 0 or 1
    integer  :: n = 0, nsub = 0, box = 0, lims2(2,2) = 0, sqlp = 0
    real(dp) :: kappa2 = 0.0_dp, subscale = 1.0_dp
    real     :: gain = 1.0                      !< g of the scale link
    integer  :: npasses = 0, nsubpasses = 0
    real(timer_int_kind) :: rt_prep = 0.
  contains
    procedure :: new         => objective_new
    procedure :: size        => objective_size
    procedure :: eval        => objective_eval
    procedure :: curvature   => objective_curvature
    procedure :: line_losses => objective_line_losses
    procedure :: precond     => objective_precond
    procedure :: dot         => objective_dot
    procedure :: project     => objective_project
    procedure :: fitted_gain => objective_fitted_gain
    procedure :: kill        => objective_kill
end type rec3D_wave_objective

contains

    !> The program body (shared memory): per state, the flat PCG start of each half, the ML prior from the
    !! FSC of the two start maps, the wave solve of each half, the maps at SIMPLE's scale, the FSC of the
    !! unfiltered pair; the maps and FSC files registered in the out segment, which is written once
    subroutine exec_rec3D_wave( params, build )
        type(parameters), target, intent(inout) :: params
        type(builder),    target, intent(inout) :: build
        type(reconstructor_pcg),    allocatable, target :: pcgops(:)
        type(raw_density)                       :: d_raws(2)
        type(rec3D_wave_objective), allocatable :: obj
        type(ptcl_read_session),    allocatable :: rsess
        type(halfmap_diagnostics_result) :: diag
        type(nlcg_solver_options) :: opts
        type(nlcg_solver_outcome) :: outcome
        type(image)               :: vols(2), unfil(2), avg
        type(string)              :: fname, support_kind
        integer, allocatable :: all_pinds(:), halves(:,:), even_pinds(:), odd_pinds(:), state_pinds(:)
        real,    allocatable :: x(:,:,:), xreg(:,:,:), supp(:,:,:)
        real,    allocatable, target :: v(:,:,:)
        real,    pointer :: vp(:)
        integer  :: model, kappa, state, ihalf, nhalf(2), i, buffer_status, nstates_done
        real     :: gain, kv, kv2, scale_volts
        real(dp) :: gain_fit, gain_num, gain_den
        logical  :: l_sigma_loaded
        character(len=4), parameter :: HALF_NAMES(2) = ['even', 'odd ']
        select case(trim(params%wave_model))
            case('projection')
                model = WAVE_PROJECTION
            case('first_born')
                model = WAVE_FIRST_BORN
            case('rytov')
                model = WAVE_RYTOV
            case('multislice')
                model = WAVE_MSLICE
            case DEFAULT
                THROW_HARD('unknown wave_model; reconstruct3D_wave')
        end select
        select case(trim(params%ewald))
            case('neg')
                kappa = -1
            case('pos')
                kappa = 1
            case DEFAULT
                kappa = 0
        end select
        if( model == WAVE_PROJECTION .and. kappa /= 0 ) THROW_HARD('wave_model=projection samples the central section: ewald=no')
        if( model == WAVE_MSLICE .and. kappa == 0 ) THROW_HARD('wave_model=multislice propagates through the box: ewald=neg or pos')
        if( params%l_pcg_solvent ) THROW_HARD('reconstruct3D_wave has no solvent prior: pcg_solvent=no')
        if( params%total_dose <= 0. ) THROW_HARD('total_dose must be positive; reconstruct3D_wave')
        if( trim(params%automsk) /= 'no' ) THROW_HARD('reconstruct3D_wave solves on the mskdiam sphere or pcg_mskfile: automsk=no')
        ! the gain from relative intensity to the particles, normalised to unit noise at the native box
        gain = sqrt(params%total_dose) * params%smpd
        call ensure_field_halves(build%spproj, params%oritype, params%projfile)
        if( params%cc_objfun == OBJFUN_EUCLID )then
            call load_sigma2(params, build%pftc, build%esig, build%spproj, build%spproj_field, l_sigma_loaded)
            if( .not. l_sigma_loaded ) THROW_HARD('reconstruct3D_wave with objfun=euclid needs the project''s sigma2')
        endif
        allocate(all_pinds(build%spproj_field%get_noris()))
        all_pinds = [(i, i = 1, size(all_pinds))]
        call prepimgbatch(params, build, MAXIMGBATCHSZ)
        allocate(rsess)
        call rsess%new(params)
        opts%maxits      = params%maxits
        opts%refine_step = .true.
        opts%verbose     = 1
        nstates_done = 0
        do state = 1, params%nstates
            ! the state's halves
            call rec_state_half(build, state, 0, all_pinds, even_pinds)
            nhalf(1) = size(even_pinds)
            call rec_state_half(build, state, 1, all_pinds, odd_pinds)
            nhalf(2) = size(odd_pinds)
            if( all(nhalf == 0) )then
                write(logfhandle,'(A,I0,A)') '>>> RECONSTRUCT3D_WAVE: STATE ', state, ' has no particles; skipped'
                cycle
            endif
            if( any(nhalf == 0) ) THROW_HARD('a state has particles in one half only; reconstruct3D_wave')
            if( allocated(halves) ) deallocate(halves)
            allocate(halves(maxval(nhalf),2), source=0)
            halves(1:nhalf(1),1) = even_pinds
            halves(1:nhalf(2),2) = odd_pinds
            state_pinds = [even_pinds, odd_pinds]
            ! the state's particles read once into the pass buffer: every later pass copies from it
            buffer_status = fill_pass_buffer(params, build, rsess, size(state_pinds), state_pinds, &
                &ptcl_buffer_budget(params, build))
            write(logfhandle,'(A,I0,A)') '>>> RECONSTRUCT3D_WAVE: STATE ', state, ' particle buffer: '// &
                &pass_buffer_status_name(buffer_status)
            ! the flat PCG start of each half, and the scale of its map in volts
            allocate(pcgops(2))
            call state_scale(halves(1:nhalf(1),1), kv)
            call state_scale(halves(1:nhalf(2),2), kv2)
            if( abs(kv - kv2) > 1.e-3 ) THROW_HARD('the halves of a state differ in voltage; reconstruct3D_wave')
            scale_volts = electron_scale_link(params%total_dose, kv, params%smpd, params%box)
            do ihalf = 1, 2
                call base_solve(ihalf, halves(1:nhalf(ihalf),ihalf), vols(ihalf))
            end do
            support_kind = 'sphere'
            if( params%pcg_mskfile%is_allocated() )then
                if( len_trim(params%pcg_mskfile%to_char()) > 0 ) support_kind = 'explicit'
            endif
            call avg%copy(vols(1))
            call avg%add(vols(2))
            call avg%mul(0.5)
            call calculate_pcg_state_diagnostics(params, state, 'WAVE START', vols(1), vols(2), avg, diag, &
                &.false., support_kind%to_char())
            ! each half: the ML prior from the pair's FSC, the regularized flat map, then the wave solve
            do ihalf = 1, 2
                call pcgops(ihalf)%set_ml_prior(diag%fsc, params%tau, params%hp)
                call regularize_start(ihalf, vols(ihalf), xreg)
                ! the flat start beside the wave maps: the flat model under the same prior and support
                call vols(ihalf)%set_rmat(xreg, .false.)
                fname = refine3D_state_vol_suffix_fname(state, '_'//trim(HALF_NAMES(ihalf))//'_flat')
                call vols(ihalf)%write(fname)
                allocate(obj)
                call obj%new(params, build, pcgops(ihalf), halves(1:nhalf(ihalf),ihalf), model, kappa, gain, &
                    &scale_volts)
                v = xreg / scale_volts
                opts%tag = 'WAVE STATE '//int2str_pad(state,2)//' '//trim(HALF_NAMES(ihalf))
                vp(1:size(v)) => v
                call nlcg_solve(obj, vp, opts, outcome)
                nullify(vp)
                gain_fit = obj%fitted_gain(v, gain_num, gain_den)
                call write_wave_diagnostics(state, HALF_NAMES(ihalf), nhalf(ihalf), outcome, gain, gain_fit, gain_num, &
                    &gain_den, obj)
                call obj%kill
                deallocate(obj)
                ! the regularized map and the unfiltered one (the Wiener factor undone), at SIMPLE's scale
                call vols(ihalf)%set_rmat(v * scale_volts, .false.)
                x = v
                call pcgops(ihalf)%unshrink_by_ml_prior(x)
                call pcgops(ihalf)%get_solve_support(supp)
                x = x * supp
                call unfil(ihalf)%new([params%box_crop,params%box_crop,params%box_crop], params%smpd_crop)
                call unfil(ihalf)%set_rmat(x * scale_volts, .false.)
                call pcgops(ihalf)%kill
            end do
            deallocate(pcgops)
            ! the maps, and the FSC of the unfiltered pair
            do ihalf = 1, 2
                fname = refine3D_state_halfvol_fname(state, trim(HALF_NAMES(ihalf)))
                call vols(ihalf)%write(fname)
                fname = refine3D_state_halfvol_fname(state, trim(HALF_NAMES(ihalf)), unfil=.true.)
                call unfil(ihalf)%write(fname)
            end do
            call avg%copy(vols(1))
            call avg%add(vols(2))
            call avg%mul(0.5)
            fname = refine3D_state_vol_fname(state)
            call avg%write(fname)
            call avg%copy(unfil(1))
            call avg%add(unfil(2))
            call avg%mul(0.5)
            call calculate_pcg_state_diagnostics(params, state, 'WAVE', unfil(1), unfil(2), avg, diag, &
                &.false., support_kind%to_char())
            fname = refine3D_fsc_fname(state)
            call arr2file(diag%fsc, fname)
            call write_halfmap_diagnostics(diag, params%box_crop, params%smpd_crop, refine3D_resolution_txt_fbody(state))
            call build%spproj%add_fsc2os_out(fname, state, params%box_crop)
            call build%spproj%add_vol2os_out(refine3D_state_vol_fname(state), params%smpd_crop, state, 'vol')
            call release_pass_buffer
            do ihalf = 1, 2
                call vols(ihalf)%kill
                call unfil(ihalf)%kill
            end do
            call avg%kill
            nstates_done = nstates_done + 1
        end do
        if( nstates_done == 0 ) THROW_HARD('no state has particles; reconstruct3D_wave')
        call build%spproj%write_segment_inside('out', params%projfile)
        call rsess%kill
        call killimgbatch(build)

    contains

        !> the voltage of the particles pinds_here, which must share it (the PCG map's scale in volts is one number)
        subroutine state_scale( pinds_here, kv_here )
            integer, intent(in)  :: pinds_here(:)
            real,    intent(out) :: kv_here
            type(ctfparams) :: ctfparms
            integer :: ii
            kv_here = 0.
            do ii = 1, size(pinds_here)
                ctfparms = build%spproj%get_ctfparams(params%oritype, pinds_here(ii))
                if( ii == 1 ) kv_here = ctfparms%kv
                if( abs(ctfparms%kv - kv_here) > 1.e-3 ) THROW_HARD('the particles of a state differ in voltage; reconstruct3D_wave')
            end do
        end subroutine state_scale

        !> the flat PCG solve of one half from its particle pass (the PCG worker's accumulation and the master's
        !! base solve in one process); pcgops(ihalf_here) keeps its finalized operator for the regularization
        subroutine base_solve( ihalf_here, pinds_here, vol )
            integer,     intent(in)    :: ihalf_here, pinds_here(:)
            type(image), intent(inout) :: vol
            type(pcg_solver_outcome) :: pcg_outcome
            real, allocatable :: rel_res_hist(:), xbase(:,:,:)
            real(timer_int_kind) :: rt_prep
            integer(timer_int_kind) :: t0, t1, rate
            integer :: niters
            call system_clock(t0, rate)
            associate( op => pcgops(ihalf_here) )
            call op%new(params%box_crop, params%smpd_crop, PCG_LAMBDA, fft_nthreads=params%nthr)
            call set_pcg_solve_support(op, params)
            rt_prep = 0.
            call accumulate_pcg_state_half(params, build, rsess, ihalf_here-1, pinds_here, op, rt_prep)
            call op%end_accum(.true., d_raws(ihalf_here)%d)
            call op%set_op_mode(PCG_OP_KERNEL)
            allocate(xbase(params%box_crop,params%box_crop,params%box_crop), source=0.0)
            call op%solve_with_cold_restart(xbase, .false., params%maxits_pcg, params%rtol, rel_res_hist, niters, pcg_outcome)
            call handle_cold_restart_outcome(pcg_outcome, 'wave start', trim(HALF_NAMES(ihalf_here)), 'base')
            call validate_solved_map(xbase, 'wave start', state, trim(HALF_NAMES(ihalf_here)), 'base')
            call system_clock(t1)
            call report_pcg_solve('WAVE START', state, trim(HALF_NAMES(ihalf_here)), 'base', size(pinds_here), niters, &
                &real(t1-t0,dp)/real(rate,dp), pcg_outcome)
            end associate
            call vol%new([params%box_crop,params%box_crop,params%box_crop], params%smpd_crop)
            call vol%set_rmat(xbase, .false.)
            deallocate(xbase)
        end subroutine base_solve

        !> the regularized flat map of a half: its operator under the ML prior (set_ml_prior first), the
        !! closed-form shrink of the base map and maxits_ml coupled iterations, as the PCG master's
        subroutine regularize_start( ihalf_here, vol, xout )
            integer,           intent(in)    :: ihalf_here
            type(image),       intent(in)    :: vol
            real, allocatable, intent(inout) :: xout(:,:,:)
            type(pcg_solver_outcome) :: pcg_outcome
            real, allocatable :: rel_res_hist(:), x_cf(:,:,:)
            integer(timer_int_kind) :: t0, t1, rate
            integer :: niters
            call system_clock(t0, rate)
            call pcgops(ihalf_here)%regularize_finalized(d_raws(ihalf_here)%d)
            call pcgops(ihalf_here)%set_op_mode(PCG_OP_KERNEL)
            xout = vol%get_rmat()
            call pcgops(ihalf_here)%solve_regularized(xout, params%maxits_ml, rel_res_hist, niters, pcg_outcome, x_cf)
            call handle_cold_restart_outcome(pcg_outcome, 'wave start', trim(HALF_NAMES(ihalf_here)), 'ml')
            call validate_solved_map(xout, 'wave start', state, trim(HALF_NAMES(ihalf_here)), 'ml')
            call system_clock(t1)
            call report_pcg_solve('WAVE START', state, trim(HALF_NAMES(ihalf_here)), 'ml', nhalf(ihalf_here), niters, &
                &real(t1-t0,dp)/real(rate,dp), pcg_outcome)
            if( allocated(x_cf) ) deallocate(x_cf)
        end subroutine regularize_start

        subroutine write_wave_diagnostics( state_here, half, nptcls, res, g, g_fit, g_num, g_den, objective )
            integer,                    intent(in) :: state_here, nptcls
            character(len=*),           intent(in) :: half
            type(nlcg_solver_outcome),  intent(in) :: res
            real,                       intent(in) :: g
            real(dp),                   intent(in) :: g_fit, g_num, g_den
            type(rec3D_wave_objective), intent(in) :: objective
            type(string) :: dname
            integer :: funit, it
            dname = 'reconstruct3D_wave_state'//int2str_pad(state_here,2)//'_'//trim(half)//'_diagnostics.txt'
            call fopen(funit, file=dname, status='replace', action='write')
            write(funit,'(A,A)')      'wave_model=',        trim(params%wave_model)
            write(funit,'(A,A)')      'ewald=',             trim(params%ewald)
            write(funit,'(A,I0)')     'nptcls=',            nptcls
            write(funit,'(A,I0)')     'step_subset=',       objective%nsub
            write(funit,'(A,ES14.6)') 'gain_dose=',         g
            write(funit,'(A,ES14.6)') 'gain_fitted=',       g_fit
            write(funit,'(A,ES14.6)') 'gain_fit_num=',      g_num
            write(funit,'(A,ES14.6)') 'gain_fit_den=',      g_den
            write(funit,'(A,ES14.6)') 'map_scale_per_volt=', scale_volts
            write(funit,'(A,A)')      'stop_reason=',       trim(res%stop_reason)
            write(funit,'(A,I0)')     'iteration_count=',   res%iteration_count
            write(funit,'(A,I0)')     'restarts=',          res%restarts
            write(funit,'(A,I0)')     'backtracks=',        res%backtracks
            write(funit,'(A,I0)')     'full_passes=',       objective%npasses
            write(funit,'(A,I0)')     'subset_passes=',     objective%nsubpasses
            write(funit,'(A,ES20.12)') 'initial_loss=',     res%initial_loss
            do it = 1, res%iteration_count
                write(funit,'(A,I0,A,ES20.12,A,ES12.4,A,ES12.4,A,ES12.4)') 'iter=', it, ' loss=', res%loss_history(it), &
                    &' step=', res%step_history(it), ' gn_step=', res%gn_step_history(it), &
                    &' grad_rel=', res%grad_history(it)
            end do
            call fclose(funit)
            write(logfhandle,'(A,I0,A,A,A,F9.4,A,F9.4)') '>>> RECONSTRUCT3D_WAVE: STATE ', state_here, ' HALF ', trim(half), &
                &' GAIN FROM DOSE ', g, '  GAIN FITTED ', real(g_fit)
        end subroutine write_wave_diagnostics

    end subroutine exec_rec3D_wave

    ! THE OBJECTIVE

    !> the objective of the particles pinds_here (one half of a state) on the operator pcgop_here, which
    !! holds the half's ML prior and preconditioner; model, kappa_here, gain_here as wave_projector%new;
    !! scale_volts the PCG map's scale per volt
    subroutine objective_new( self, params, build, pcgop_here, pinds_here, model, kappa_here, gain_here, scale_volts )
        class(rec3D_wave_objective),      intent(inout) :: self
        class(parameters),        target, intent(inout) :: params
        class(builder),           target, intent(inout) :: build
        class(reconstructor_pcg), target, intent(inout) :: pcgop_here
        integer,                          intent(in)    :: pinds_here(:), model, kappa_here
        real,                             intent(in)    :: gain_here, scale_volts
        type(oris) :: selection
        real, allocatable :: inv_sig2(:)
        integer :: i, j, stride, batchlims(2), eo, nthr
        logical, allocatable :: l_sub(:)
        call self%kill
        self%params => params
        self%build  => build
        self%pcgop  => pcgop_here
        self%n      =  size(pinds_here)
        self%box    =  params%box_crop
        ! PCG's normal operator and prior carry boxpd^3 beside kappa^2: its adjoint folds without the 1/boxpd^3 of
        ! the transform's transpose (the wave projector's tester checks the link)
        self%kappa2 =  real(scale_volts,dp)**2 / real(OSMPL_PAD_FAC * params%box_crop,dp)**3
        self%gain   =  gain_here
        ! the step subset: every stride-th particle, first in the objective's order
        self%nsub = min(self%n, max(STEP_SUBSET_MIN, nint(STEP_SUBSET_FRAC * real(self%n))))
        self%nsub = min(self%nsub, STEP_SUBSET_MAX)
        self%subscale = real(self%n,dp) / real(self%nsub,dp)
        stride = max(1, self%n / self%nsub)
        allocate(self%pinds(self%n), source=0)
        allocate(l_sub(self%n), source=.false.)
        j = 0
        do i = 1, self%n, stride
            if( j == self%nsub ) exit
            j = j + 1
            l_sub(i) = .true.
            self%pinds(j) = pinds_here(i)
        end do
        do i = 1, self%n
            if( l_sub(i) ) cycle
            j = j + 1
            self%pinds(j) = pinds_here(i)
        end do
        ! the half's particles at the cropped box, the same selection the start solve used
        allocate(self%wp)
        ! winsz: the KB window half-width of the operator's lattice (KBWINSZ by default, PCG's three taps)
        call self%wp%new(params%box_crop, params%smpd_crop, model, kappa_here, gain_here, fft_nthreads=params%nthr, &
            &winsz=params%winsz)
        call self%wp%set_sym(build%pgrpsyms)
        call rec_obs_oris(params, build, self%pinds, selection)
        if( params%cc_objfun == OBJFUN_EUCLID )then
            eo = build%spproj_field%get_eo(self%pinds(1))
            call rec_obs_inv_sig2(params, build, eo, params%box_crop/2, inv_sig2)
            call self%wp%prep_particles(selection, inv_sig2=inv_sig2)
        else
            call self%wp%prep_particles(selection)
        endif
        call selection%kill
        self%lims2 = pcgop_here%get_lims2()
        self%sqlp  = self%lims2(1,2)**2
        allocate(self%ybatch(self%lims2(1,1):self%lims2(1,2), self%lims2(2,1):self%lims2(2,2), MAXIMGBATCHSZ))
        allocate(self%ysub(self%lims2(1,1):self%lims2(1,2), self%lims2(2,1):self%lims2(2,2), self%nsub))
        nthr = 1
        !$ nthr = omp_get_max_threads()
        call alloc_imgarr(nthr, [params%box_crop,params%box_crop,1], params%smpd_crop, self%obs)
        allocate(self%rsess)
        call self%rsess%new(params)
        ! the step subset's planes, read once
        do i = 1, self%nsub, MAXIMGBATCHSZ
            batchlims = [i, min(self%nsub, i+MAXIMGBATCHSZ-1)]
            call read_rec_obs_batch(params, build, self%rsess, self%pinds, batchlims, self%obs, self%lims2, self%sqlp, &
                &self%ysub(:,:,batchlims(1):batchlims(2)), self%rt_prep)
        end do
        call pcgop_here%get_solve_support(self%supp)
    end subroutine objective_new

    integer function objective_size( self ) result( n )
        class(rec3D_wave_objective), intent(in) :: self
        n = self%box**3
    end function objective_size

    !> L(x) and its gradient: one pass over the half's particles
    subroutine objective_eval( self, x, loss, g )
        class(rec3D_wave_objective), intent(inout) :: self
        real, contiguous, target,    intent(in)    :: x(:)
        real(dp),                    intent(out)   :: loss
        real, contiguous, target,    intent(out)   :: g(:)
        real, pointer     :: x3(:,:,:), g3(:,:,:)
        real, allocatable :: q(:,:,:)
        integer :: i, nb, batchlims(2)
        x3(1:self%box,1:self%box,1:self%box) => x
        g3(1:self%box,1:self%box,1:self%box) => g
        call self%wp%set_volume(x3)
        call self%wp%begin_grad
        loss = 0.0_dp
        do i = 1, self%n, MAXIMGBATCHSZ
            batchlims = [i, min(self%n, i+MAXIMGBATCHSZ-1)]
            nb = batchlims(2) - batchlims(1) + 1
            call read_rec_obs_batch(self%params, self%build, self%rsess, self%pinds, batchlims, self%obs, self%lims2, &
                &self%sqlp, self%ybatch, self%rt_prep)
            loss = loss + self%wp%loss_grad(self%ybatch(:,:,1:nb), nb, i)
        end do
        call self%wp%get_grad(g3)
        q    = self%pcgop%apply_ml_prior(x3)
        loss = loss + 0.5_dp * self%kappa2 * sum(real(x3,dp) * real(q,dp))
        g3   = g3 + real(self%kappa2) * q
        self%npasses = self%npasses + 1
    end subroutine objective_eval

    !> the Gauss-Newton curvature along d at x on the step subset, scaled to the half, plus the prior's
    real(dp) function objective_curvature( self, x, d ) result( curv )
        class(rec3D_wave_objective), intent(inout) :: self
        real, contiguous, target,    intent(in)    :: x(:), d(:)
        real, pointer     :: x3(:,:,:), d3(:,:,:)
        real, allocatable :: q(:,:,:)
        x3(1:self%box,1:self%box,1:self%box) => x
        d3(1:self%box,1:self%box,1:self%box) => d
        call self%wp%set_volume(x3)
        call self%wp%set_direction(d3)
        curv = self%subscale * self%wp%curvature(self%nsub, 1)
        q    = self%pcgop%apply_ml_prior(d3)
        curv = curv + self%kappa2 * sum(real(d3,dp) * real(q,dp))
        self%nsubpasses = self%nsubpasses + 1
    end function objective_curvature

    !> L(x + t(j) d) on the step subset, its data term scaled to the half
    subroutine objective_line_losses( self, x, d, t, losses )
        class(rec3D_wave_objective), intent(inout) :: self
        real, contiguous, target,    intent(in)    :: x(:), d(:)
        real(dp),                    intent(in)    :: t(:)
        real(dp),                    intent(out)   :: losses(:)
        real, pointer     :: x3(:,:,:), d3(:,:,:)
        real, allocatable :: qx(:,:,:), qd(:,:,:)
        real(dp) :: pxx, pxd, pdd
        x3(1:self%box,1:self%box,1:self%box) => x
        d3(1:self%box,1:self%box,1:self%box) => d
        call self%wp%set_volume(x3)
        call self%wp%set_direction(d3)
        call self%wp%line_losses(self%ysub, self%nsub, 1, t, losses)
        qx  = self%pcgop%apply_ml_prior(x3)
        qd  = self%pcgop%apply_ml_prior(d3)
        pxx = sum(real(x3,dp) * real(qx,dp))
        pxd = sum(real(x3,dp) * real(qd,dp))
        pdd = sum(real(d3,dp) * real(qd,dp))
        losses = self%subscale * losses + 0.5_dp * self%kappa2 * (pxx + 2.0_dp * t * pxd + t**2 * pdd)
        self%nsubpasses = self%nsubpasses + 1
    end subroutine objective_line_losses

    subroutine objective_precond( self, r, z )
        class(rec3D_wave_objective), intent(inout) :: self
        real, contiguous, target,    intent(in)    :: r(:)
        real, contiguous, target,    intent(out)   :: z(:)
        real, pointer :: r3(:,:,:), z3(:,:,:)
        r3(1:self%box,1:self%box,1:self%box) => r
        z3(1:self%box,1:self%box,1:self%box) => z
        z3 = self%pcgop%apply_preconditioner(r3)
    end subroutine objective_precond

    real(dp) function objective_dot( self, a, b ) result( value )
        class(rec3D_wave_objective), intent(in) :: self
        real, contiguous, target,    intent(in) :: a(:), b(:)
        value = sum(real(a,dp) * real(b,dp))
    end function objective_dot

    subroutine objective_project( self, v )
        class(rec3D_wave_objective), intent(inout) :: self
        real, contiguous, target,    intent(inout) :: v(:)
        real, pointer :: v3(:,:,:)
        v3(1:self%box,1:self%box,1:self%box) => v
        v3 = v3 * self%supp
    end subroutine objective_project

    !> The gain the data prefer (a diagnostic, wave_model_policy.md section 5): the one-scalar fit of the data's nonlinear part
    !! at v over the half (wave_projector%gain_terms), with its two sums; zero when the fit of the nonlinear
    !! part is not positive (the data hold no measurable nonlinear part)
    real(dp) function objective_fitted_gain( self, v, snum, sden ) result( g_fit )
        class(rec3D_wave_objective), intent(inout) :: self
        real,                        intent(in)    :: v(self%box,self%box,self%box)
        real(dp),                    intent(out)   :: snum, sden
        real(dp), allocatable :: num(:), den(:)
        integer  :: i, nb, ib, batchlims(2)
        call self%wp%set_volume(v)
        allocate(num(MAXIMGBATCHSZ), den(MAXIMGBATCHSZ))
        snum = 0.0_dp
        sden = 0.0_dp
        do i = 1, self%n, MAXIMGBATCHSZ
            batchlims = [i, min(self%n, i+MAXIMGBATCHSZ-1)]
            nb = batchlims(2) - batchlims(1) + 1
            call read_rec_obs_batch(self%params, self%build, self%rsess, self%pinds, batchlims, self%obs, self%lims2, &
                &self%sqlp, self%ybatch, self%rt_prep)
            call self%wp%gain_terms(self%ybatch(:,:,1:nb), nb, i, num(1:nb), den(1:nb))
            do ib = 1, nb
                snum = snum + num(ib)
                sden = sden + den(ib)
            end do
        end do
        g_fit = 0.0_dp
        if( snum > 0.0_dp ) g_fit = real(self%gain,dp) * sden / snum
    end function objective_fitted_gain

    !> idempotent
    subroutine objective_kill( self )
        class(rec3D_wave_objective), intent(inout) :: self
        if( allocated(self%wp) )then
            call self%wp%kill
            deallocate(self%wp)
        endif
        if( allocated(self%rsess) )then
            call self%rsess%kill
            deallocate(self%rsess)
        endif
        if( allocated(self%obs)    ) call dealloc_imgarr(self%obs)
        if( allocated(self%pinds)  ) deallocate(self%pinds)
        if( allocated(self%ybatch) ) deallocate(self%ybatch)
        if( allocated(self%ysub)   ) deallocate(self%ysub)
        if( allocated(self%supp)   ) deallocate(self%supp)
        nullify(self%params, self%build, self%pcgop)
        self%n = 0; self%nsub = 0; self%box = 0
        self%kappa2 = 0.0_dp; self%subscale = 1.0_dp; self%gain = 1.0
        self%npasses = 0; self%nsubpasses = 0
        self%rt_prep = 0.
    end subroutine objective_kill

end module simple_rec3D_wave_strategy
