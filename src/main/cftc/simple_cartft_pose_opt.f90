!@descr: bounded Levenberg-Marquardt optimizer of one five-parameter Cartesian particle pose (refine=cont)
!
! The role pftc_shsrch_grad has for the polar branch. One optimizer belongs to one particle
! search: built in new, released in kill, used by the one thread running that search. It holds
! fixed-size state only (policy, the last solve's record), allocates nothing during a solve,
! and evaluates through a cartft_calc passed intent(in) to every solve. It imports no
! matcher, commander, UI, project or parallelization module.
!
! Coordinates: a pose is the particle rotation matrix and a 2D shift. Project orientations
! (ori) carry native-box pixels; the calculator works on the cropped box, so refine converts
! by box_crop/box on the way out. The rotation is updated on the right, R exp([omega]x)
! (right_increment_rotation); omega in radians. As in the polar branch (prepimg4align and
! assign_ori, O4), refine solves for the shift increment from the stored shift: the particle
! slot holds the observation centred on the stored shift, the solve starts from a zero
! increment and the committed shift is the stored shift plus the increment.
!
! Transaction (refine, refine_pose): evaluate the seed, run the joint five-parameter stage
! (preceded by the shift-only stage on the shift-then-joint route), bounded per step and
! cumulatively from the seed, and commit the staged pose only if its objective is finite and below the
! seed's; every other outcome returns the seed bit for bit (rollback). The objective is the
! one the calculator's particle slot was prepared for.
!
! Bounds (O6): the total shift from the seed is bounded by trs and the total rotation by
! athres_cont (the capture range for wrongly assigned orientations, independent of the polar
! searches' athres and prob_athres); trs = 0 freezes the shifts as l_doshift=.false. does in
! the polar branch. Route (C16, ruling of 2026-10-02): the joint stage alone in production,
! about a third faster per transaction; shift then joint (the route of the E25 floors) when
! new is asked for it. Step and iteration caps are private defaults the tester may override.
module simple_cartft_pose_opt
use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
use simple_core_module_api, only: dp, PI, ori, dm2euler, simple_exception
use simple_cartft_calc,     only: cartft_calc
implicit none
private
public :: cartft_pose_opt, right_increment_rotation
! typed per-particle outcomes of a solve or a stage
public :: CARTFT_NOT_ATTEMPTED, CARTFT_INVALID_PREPARATION, CARTFT_ACCEPTED, CARTFT_NO_IMPROVEMENT
public :: CARTFT_NO_RELIABLE_UPDATE, CARTFT_BOUND_REJECTED, CARTFT_INVALID_NUMERICS, CARTFT_ITERATION_LIMIT
! stage selectors of get_stage
public :: CARTFT_STAGE_SHIFT, CARTFT_STAGE_JOINT

#include "simple_local_flags.inc"

integer, parameter :: CARTFT_NOT_ATTEMPTED       =  0 !< stage or solve not run
integer, parameter :: CARTFT_INVALID_PREPARATION = -1 !< the particle slot holds no valid preparation
integer, parameter :: CARTFT_ACCEPTED            =  1 !< an objective-reducing pose was committed
integer, parameter :: CARTFT_NO_IMPROVEMENT      =  2 !< finite, no reduction; seed kept
integer, parameter :: CARTFT_NO_RELIABLE_UPDATE  =  3 !< unobservable or ill-conditioned system; seed kept
integer, parameter :: CARTFT_BOUND_REJECTED      =  4 !< every proposal left the step or cumulative bound; seed kept
integer, parameter :: CARTFT_INVALID_NUMERICS    =  5 !< non-finite objective or step; seed kept
integer, parameter :: CARTFT_ITERATION_LIMIT     =  6 !< iteration cap reached without a decision
integer, parameter :: CARTFT_STAGE_SHIFT = 1, CARTFT_STAGE_JOINT = 2

! LM numerics
real(dp), parameter :: NUMERIC_FLOOR                = epsilon(1._dp)**2
real(dp), parameter :: LM_INITIAL_DAMPING           = 1.e-3_dp
real(dp), parameter :: LM_INITIAL_REJECTION_MULT    = 4._dp
real(dp), parameter :: LM_MAX_DAMPING               = 1._dp/epsilon(1._dp)
real(dp), parameter :: LM_ACCEPTED_RELATIVE_TOL     = 1.e-6_dp
integer,  parameter :: LM_MAX_CONSECUTIVE_REJECTIONS = 8
! default policy
integer,  parameter :: DEFAULT_MAXITS             = 40
real(dp), parameter :: DEFAULT_ROTATION_SCALE     = 0.1_dp
real(dp), parameter :: DEFAULT_SHIFT_STEP_NATIVE  = 1._dp  !< native pixels per proposal

!> Record of one LM stage.
type :: cartft_stage
    integer  :: status         = CARTFT_NOT_ATTEMPTED
    integer  :: niter          = 0
    integer  :: nattempted     = 0     !< trial poses evaluated
    integer  :: naccepted      = 0     !< trial poses accepted
    integer  :: nbound_hits    = 0     !< proposals capped or rejected by a bound
    real(dp) :: max_rotation_step = 0._dp
    real(dp) :: max_shift_step    = 0._dp
    real(dp) :: objective_before  = -1._dp !< set by a transaction
    real(dp) :: objective_after   = -1._dp !< set by a transaction
end type cartft_stage

type :: cartft_pose_opt
    private
    ! geometry
    integer  :: box      = 0               !< native box of the project shifts
    integer  :: box_crop = 0               !< box of the calculator
    ! policy
    logical  :: refine_shifts = .true.     !< false when trs = 0: shifts frozen
    logical  :: shift_first   = .false.    !< the shift-only stage before the joint stage
    integer  :: maxits   = DEFAULT_MAXITS  !< iteration cap of each stage
    real(dp) :: rotation_scale     = 0._dp !< radians per joint proposal (scaling and cap)
    real(dp) :: max_total_rotation = 0._dp !< radians from the seed (athres_cont)
    real(dp) :: shift_step_bound   = 0._dp !< cropped pixels per proposal
    real(dp) :: max_total_shift    = 0._dp !< cropped pixels from the seed (trs)
    ! record of the last solve
    integer  :: status = CARTFT_NOT_ATTEMPTED
    real(dp) :: objective_before = -1._dp, objective_after = -1._dp
    real(dp) :: cumulative_rotation = 0._dp, cumulative_shift = 0._dp
    type(cartft_stage) :: stages(2)
    logical  :: exists = .false.
contains
    ! lifecycle
    procedure          :: new
    procedure          :: kill
    ! solves
    procedure          :: refine
    procedure          :: refine_pose
    procedure          :: refine_shift
    procedure          :: refine_joint
    ! record of the last solve
    procedure          :: get_status
    procedure          :: get_objectives
    procedure          :: get_motion
    procedure          :: get_stage
    ! private helpers
    procedure, private :: reset_record
end type cartft_pose_opt

contains

    ! LIFECYCLE

    !> Geometry (box of the project shifts, box_crop of the calculator), bounds (trs the total
    !! shift in native pixels, 0 freezes shifts; rotation_bound the total rotation in degrees) and
    !! route (shift_first: the shift-only stage first; default the joint stage alone).
    subroutine new( self, box, box_crop, trs, rotation_bound, shift_first, maxits, rotation_scale, shift_step_bound )
        class(cartft_pose_opt), intent(inout) :: self
        integer,                intent(in)    :: box, box_crop
        real,                   intent(in)    :: trs, rotation_bound
        logical,  optional,     intent(in)    :: shift_first
        integer,  optional,     intent(in)    :: maxits
        real(dp), optional,     intent(in)    :: rotation_scale, shift_step_bound
        real(dp) :: crop_scale
        call self%kill
        if( box < 2 .or. box_crop < 2 .or. mod(box,2) /= 0 .or. mod(box_crop,2) /= 0 ) &
            &THROW_HARD('cartft_pose_opt requires positive even boxes')
        if( .not. ieee_is_finite(trs) .or. trs < 0. ) THROW_HARD('cartft_pose_opt shift bound trs must be finite and >= 0')
        if( .not. positive_finite(real(rotation_bound, dp)) ) THROW_HARD('cartft_pose_opt rotation bound must be positive and finite')
        self%box      = box
        self%box_crop = box_crop
        crop_scale    = real(box_crop, dp)/real(box, dp)
        self%max_total_rotation = real(rotation_bound, dp)*real(PI, dp)/180._dp
        self%max_total_shift    = real(trs, dp)*crop_scale
        self%refine_shifts      = trs > 0.
        self%shift_first        = .false.
        if( present(shift_first) ) self%shift_first = shift_first
        self%maxits             = DEFAULT_MAXITS
        self%rotation_scale     = DEFAULT_ROTATION_SCALE
        self%shift_step_bound   = DEFAULT_SHIFT_STEP_NATIVE*crop_scale
        if( present(maxits)           ) self%maxits           = maxits
        if( present(rotation_scale)   ) self%rotation_scale   = rotation_scale
        if( present(shift_step_bound) ) self%shift_step_bound = shift_step_bound
        if( self%maxits < 1 ) THROW_HARD('cartft_pose_opt requires at least one LM iteration')
        if( .not. positive_finite(self%rotation_scale) )   THROW_HARD('cartft_pose_opt rotation scale must be positive and finite')
        if( .not. positive_finite(self%shift_step_bound) ) THROW_HARD('cartft_pose_opt shift step bound must be positive and finite')
        self%exists = .true.
    end subroutine new

    !> Reset policy and record.
    subroutine kill( self )
        class(cartft_pose_opt), intent(inout) :: self
        self%box                = 0
        self%box_crop           = 0
        self%refine_shifts      = .true.
        self%shift_first        = .false.
        self%maxits             = DEFAULT_MAXITS
        self%rotation_scale     = 0._dp
        self%max_total_rotation = 0._dp
        self%shift_step_bound   = 0._dp
        self%max_total_shift    = 0._dp
        call self%reset_record
        self%exists = .false.
    end subroutine kill

    ! SOLVES

    !> Refine the pose stored in o against particle slot iptcl and the reference of o's state
    !! and half (eo 0 even, 1 odd). The slot holds the observation centred on o's stored shift
    !! (O4), so the solve starts from o's rotation and a zero shift increment. On
    !! CARTFT_ACCEPTED only the Euler angles and the native shift (the stored shift plus the
    !! increment) of o change; state, half, proj and every other field stay as they are. On any
    !! other status o is unchanged. THROW_HARD for a state below 1 or an eo outside 0/1.
    subroutine refine( self, calc, iptcl, o )
        class(cartft_pose_opt), intent(inout) :: self
        class(cartft_calc),     intent(in)    :: calc
        integer,                intent(in)    :: iptcl
        class(ori),             intent(inout) :: o
        real(dp) :: rotmat(3,3), shift_increment(2)
        integer  :: state
        logical  :: iseven
        if( .not. self%exists ) THROW_HARD('cartft_pose_opt used before new')
        state = o%get_state()
        if( state < 1 ) THROW_HARD('cartft_pose_opt requires a positive state')
        select case(o%get_eo())
            case(0)
                iseven = .true.
            case(1)
                iseven = .false.
            case DEFAULT
                THROW_HARD('cartft_pose_opt requires an even/odd assignment')
        end select
        ! the rotation matrix is copied directly, avoiding an Euler round trip
        rotmat          = real(o%get_mat(), dp)
        shift_increment = 0._dp
        call self%refine_pose(calc, state, iseven, iptcl, rotmat, shift_increment)
        if( self%status /= CARTFT_ACCEPTED ) return
        call o%set_euler(real(dm2euler(rotmat)))
        call o%set_shift(o%get_2Dshift() + real(shift_increment)*real(self%box)/real(self%box_crop))
    end subroutine refine

    !> The transaction of refine in calculator coordinates: rotmat and shift (cropped pixels,
    !! the model shift relative to the slot's observation) are the seed on entry and the
    !! committed pose on return (the seed unless accepted). With trs = 0 only the rotation is
    !! refined. The Phase 2 forwarders and the testers use it; Phase 7 decides whether it stays
    !! public.
    subroutine refine_pose( self, calc, state, iseven, iptcl, rotmat, shift )
        class(cartft_pose_opt), intent(inout) :: self
        class(cartft_calc),     intent(in)    :: calc
        integer,                intent(in)    :: state, iptcl
        logical,                intent(in)    :: iseven
        real(dp),               intent(inout) :: rotmat(3,3), shift(2)
        real(dp) :: staged_rotmat(3,3), staged_shift(2), gradient(5)
        real(dp) :: shift_objective_after, joint_objective_after, sine_half
        if( .not. self%exists ) THROW_HARD('cartft_pose_opt used before new')
        ! every early exit leaves rotmat and shift at the seed (rollback)
        call self%reset_record
        if( .not. calc%ptcl_is_valid(iptcl) )then
            self%status = CARTFT_INVALID_PREPARATION
            return
        endif
        ! the seed objective is the single acceptance baseline of the transaction
        call calc%objective_gradient(state, iseven, iptcl, rotmat, shift, self%objective_before, gradient)
        if( .not. ieee_is_finite(self%objective_before) )then
            self%status = CARTFT_INVALID_NUMERICS
            return
        endif
        staged_rotmat = rotmat
        staged_shift  = shift
        shift_objective_after = self%objective_before
        if( self%refine_shifts .and. self%shift_first )then
            ! stage 1 of the shift-then-joint route: translation only, the seed rotation fixed
            call self%refine_shift(calc, state, iseven, iptcl, staged_rotmat, staged_shift)
            call calc%objective_gradient(state, iseven, iptcl, staged_rotmat, staged_shift, &
                &shift_objective_after, gradient)
            associate( st => self%stages(CARTFT_STAGE_SHIFT) )
                st%objective_before = self%objective_before
                st%objective_after  = shift_objective_after
                ! the cumulative shift bound from the seed, explicitly
                if( sqrt(sum((staged_shift - shift)**2)) > self%max_total_shift + 10._dp*epsilon(1._dp) )then
                    st%status      = CARTFT_BOUND_REJECTED
                    st%nbound_hits = st%nbound_hits + 1
                endif
                if( aborts_route(st%status) )then
                    ! a failed shift stage cannot supply a trustworthy joint-stage seed
                    self%status          = st%status
                    self%objective_after = self%objective_before
                    return
                endif
            end associate
        endif
        ! joint stage from the route's current endpoint (the shift result or the seed);
        ! both cumulative guards stay anchored at the transaction seed
        if( self%refine_shifts )then
            call self%refine_joint(calc, state, iseven, iptcl, staged_rotmat, staged_shift, &
                &anchor_rotmat=rotmat, anchor_shift=shift)
        else
            call self%refine_joint(calc, state, iseven, iptcl, staged_rotmat, staged_shift, &
                &active=[.true., .true., .true., .false., .false.], anchor_rotmat=rotmat, anchor_shift=shift)
        endif
        call calc%objective_gradient(state, iseven, iptcl, staged_rotmat, staged_shift, &
            &joint_objective_after, gradient)
        associate( st => self%stages(CARTFT_STAGE_JOINT) )
            st%objective_before = shift_objective_after
            st%objective_after  = joint_objective_after
            if( aborts_route(st%status) )then
                self%status          = st%status
                self%objective_after = self%objective_before
                return
            endif
        end associate
        ! commit the staged pose only if its final objective beats the seed
        self%objective_after = joint_objective_after
        if( .not. ieee_is_finite(self%objective_before) .or. .not. ieee_is_finite(self%objective_after) )then
            self%status          = CARTFT_INVALID_NUMERICS
            self%objective_after = self%objective_before
            return
        endif
        if( .not.(self%objective_after < self%objective_before) )then
            self%status = CARTFT_NO_IMPROVEMENT
            return
        endif
        self%status = CARTFT_ACCEPTED
        ! cumulative motion of the committed transaction
        sine_half = sqrt(sum((staged_rotmat - rotmat)**2))/(2._dp*sqrt(2._dp))
        self%cumulative_rotation = 2._dp*asin(max(0._dp, min(1._dp, sine_half)))
        self%cumulative_shift    = sqrt(sum((staged_shift - shift)**2))
        rotmat = staged_rotmat
        shift  = staged_shift
    end subroutine refine_pose

    !> The shift-only LM stage alone, with the rotation fixed: damped 2x2 Gauss-Newton steps
    !! capped at shift_step_bound. shift is updated through accepted steps only. Overwrites
    !! the CARTFT_STAGE_SHIFT record; an invalid particle slot gives CARTFT_INVALID_NUMERICS.
    subroutine refine_shift( self, calc, state, iseven, iptcl, rotmat, shift )
        class(cartft_pose_opt), intent(inout) :: self
        class(cartft_calc),     intent(in)    :: calc
        integer,                intent(in)    :: state, iptcl
        logical,                intent(in)    :: iseven
        real(dp),               intent(in)    :: rotmat(3,3)
        real(dp),               intent(inout) :: shift(2)
        real(dp) :: gradient(2), hessian(2,2), trial_gradient(2), trial_hessian(2,2)
        real(dp) :: solve_matrix(2,2), diagonal(2), direction(2), trial_shift(2)
        real(dp) :: objective, trial_objective, mu, rejection_multiplier
        real(dp) :: det, predicted, actual, ratio, maxdiag
        real(dp) :: discriminant, lambda_max, lambda_min, step_norm, relative_reduction
        integer  :: axis, iteration, naccepted, consecutive_rejections
        logical  :: bounded_trial
        if( .not. self%exists ) THROW_HARD('cartft_pose_opt used before new')
        self%stages(CARTFT_STAGE_SHIFT) = cartft_stage(status=CARTFT_ITERATION_LIMIT)
        associate( st => self%stages(CARTFT_STAGE_SHIFT) )
            if( .not. calc%ptcl_is_valid(iptcl) )then
                st%status = CARTFT_INVALID_NUMERICS
                return
            endif
            call calc%shift_normal_terms(state, iseven, iptcl, rotmat, shift, objective, gradient, hessian)
            naccepted = 0
            consecutive_rejections = 0
            bounded_trial = .false.
            if( .not. ieee_is_finite(objective) .or. any(.not. ieee_is_finite(gradient)) .or. &
                &any(.not. ieee_is_finite(hessian)) )then
                st%status = CARTFT_INVALID_NUMERICS
                return
            endif
            mu = LM_INITIAL_DAMPING
            rejection_multiplier = LM_INITIAL_REJECTION_MULT
            do iteration = 1, self%maxits
                st%niter = iteration
                maxdiag  = max(maxval([(hessian(axis,axis), axis=1,2)]), 1._dp)
                ! eigenvalues diagnose whether both shift directions are observable
                discriminant = sqrt(max(0._dp, (hessian(1,1) - hessian(2,2))**2 + 4._dp*hessian(1,2)*hessian(2,1)))
                lambda_max   = 0.5_dp*(hessian(1,1) + hessian(2,2) + discriminant)
                lambda_min   = 0.5_dp*(hessian(1,1) + hessian(2,2) - discriminant)
                if( lambda_max <= sqrt(epsilon(1._dp))*maxdiag .or. &
                    &lambda_min <= sqrt(epsilon(1._dp))*max(lambda_max, 1._dp) )then
                    st%status = CARTFT_NO_RELIABLE_UPDATE
                    exit
                endif
                if( sqrt(dot_product(gradient, gradient)) < 1.e-8_dp )then
                    st%status = merge(CARTFT_ACCEPTED, CARTFT_NO_IMPROVEMENT, naccepted > 0)
                    exit
                endif
                do axis = 1, 2
                    diagonal(axis) = max(hessian(axis,axis), sqrt(epsilon(1._dp))*maxdiag, epsilon(1._dp))
                end do
                solve_matrix      = hessian
                solve_matrix(1,1) = solve_matrix(1,1) + mu*diagonal(1)
                solve_matrix(2,2) = solve_matrix(2,2) + mu*diagonal(2)
                det = solve_matrix(1,1)*solve_matrix(2,2) - solve_matrix(1,2)*solve_matrix(2,1)
                if( abs(det) <= epsilon(1._dp)*maxdiag*maxdiag )then
                    st%status = CARTFT_NO_RELIABLE_UPDATE
                    exit
                endif
                direction(1) = (-solve_matrix(2,2)*gradient(1) + solve_matrix(1,2)*gradient(2))/det
                direction(2) = ( solve_matrix(2,1)*gradient(1) - solve_matrix(1,1)*gradient(2))/det
                if( any(.not. ieee_is_finite(direction)) )then
                    st%status = CARTFT_INVALID_NUMERICS
                    exit
                endif
                ! shift coordinates are pixels; cap every trial displacement at the step bound
                step_norm = sqrt(dot_product(direction, direction))
                if( step_norm > self%shift_step_bound )then
                    direction     = direction*(self%shift_step_bound/step_norm)
                    bounded_trial = .true.
                    st%nbound_hits = st%nbound_hits + 1
                endif
                step_norm = min(step_norm, self%shift_step_bound)
                st%max_shift_step = max(st%max_shift_step, step_norm)
                predicted = -dot_product(gradient, direction) - 0.5_dp*dot_product(direction, matmul(hessian, direction))
                if( .not. ieee_is_finite(predicted) )then
                    st%status = CARTFT_INVALID_NUMERICS
                    exit
                else if( predicted <= 0._dp )then
                    consecutive_rejections = consecutive_rejections + 1
                    call increase_lm_damping(mu, rejection_multiplier)
                    if( consecutive_rejections >= LM_MAX_CONSECUTIVE_REJECTIONS )then
                        st%status = rejection_terminal_status(naccepted, bounded_trial)
                        exit
                    endif
                    cycle
                endif
                trial_shift   = shift + direction
                st%nattempted = st%nattempted + 1
                call calc%shift_normal_terms(state, iseven, iptcl, rotmat, trial_shift, &
                    &trial_objective, trial_gradient, trial_hessian)
                if( .not. ieee_is_finite(trial_objective) .or. any(.not. ieee_is_finite(trial_gradient)) .or. &
                    &any(.not. ieee_is_finite(trial_hessian)) )then
                    consecutive_rejections = consecutive_rejections + 1
                    call increase_lm_damping(mu, rejection_multiplier)
                    st%status = CARTFT_INVALID_NUMERICS
                    if( consecutive_rejections >= LM_MAX_CONSECUTIVE_REJECTIONS ) exit
                    cycle
                endif
                actual = objective - trial_objective
                ratio  = actual/predicted
                if( actual > 0._dp .and. ratio >= 0.25_dp )then
                    relative_reduction = actual/max(abs(objective), NUMERIC_FLOOR)
                    shift     = trial_shift
                    objective = trial_objective
                    gradient  = trial_gradient
                    hessian   = trial_hessian
                    naccepted = naccepted + 1
                    consecutive_rejections = 0
                    st%naccepted = naccepted
                    if( ratio > 0.75_dp ) mu = max(mu/2._dp, epsilon(1._dp))
                    rejection_multiplier = LM_INITIAL_REJECTION_MULT
                    st%status = CARTFT_ACCEPTED
                    if( step_norm < 1.e-8_dp .or. relative_reduction < LM_ACCEPTED_RELATIVE_TOL ) exit
                else
                    consecutive_rejections = consecutive_rejections + 1
                    call increase_lm_damping(mu, rejection_multiplier)
                    if( consecutive_rejections >= LM_MAX_CONSECUTIVE_REJECTIONS )then
                        st%status = rejection_terminal_status(naccepted, bounded_trial)
                        exit
                    endif
                endif
            end do
            if( st%status == CARTFT_ITERATION_LIMIT .and. naccepted == 0 .and. bounded_trial ) &
                &st%status = CARTFT_BOUND_REJECTED
        end associate
    end subroutine refine_shift

    !> The joint five-parameter LM stage alone: scaled, damped steps with the rotation capped
    !! at rotation_scale and the shift at shift_step_bound. active masks parameters (1:3
    !! rotation, 4:5 shift; default all). With anchor_rotmat and anchor_shift (both or
    !! neither) every trial pose must stay within max_total_rotation and max_total_shift of
    !! the anchor. rotmat and shift are updated through accepted steps only. Overwrites the
    !! CARTFT_STAGE_JOINT record; an invalid particle slot gives CARTFT_INVALID_NUMERICS.
    !! THROW_HARD for no active parameter, one anchor without the other, or a non-finite anchor.
    subroutine refine_joint( self, calc, state, iseven, iptcl, rotmat, shift, active, anchor_rotmat, anchor_shift )
        class(cartft_pose_opt), intent(inout) :: self
        class(cartft_calc),     intent(in)    :: calc
        integer,                intent(in)    :: state, iptcl
        logical,                intent(in)    :: iseven
        real(dp),               intent(inout) :: rotmat(3,3), shift(2)
        logical,  optional,     intent(in)    :: active(5)
        real(dp), optional,     intent(in)    :: anchor_rotmat(3,3), anchor_shift(2)
        real(dp) :: gradient(5), hessian(5,5), trial_gradient(5), trial_hessian(5,5)
        real(dp) :: scaled_gradient(5), scaled_hessian(5,5), solve_matrix(5,5)
        real(dp) :: diagonal(5), scaled_direction(5), direction(5)
        real(dp) :: trial_rotmat(3,3), trial_shift(2)
        real(dp) :: objective, trial_objective, mu, rejection_multiplier
        real(dp) :: predicted, actual, ratio, rotation_norm, shift_norm
        real(dp) :: relative_reduction, cumulative_rotation, cumulative_shift, sine_half
        integer  :: iteration, naccepted, consecutive_rejections
        logical  :: l_active(5), bounded_trial, bounded_step, cumulative_guard
        logical  :: accept_trial, identifiable, reliable, stationary
        if( .not. self%exists ) THROW_HARD('cartft_pose_opt used before new')
        l_active = .true.
        if( present(active) ) l_active = active
        if( .not. any(l_active) ) THROW_HARD('cartft_pose_opt joint stage requires one active parameter')
        if( present(anchor_rotmat) .neqv. present(anchor_shift) ) &
            &THROW_HARD('cartft_pose_opt joint stage requires both anchors or none')
        cumulative_guard = present(anchor_rotmat)
        if( cumulative_guard )then
            if( any(.not. ieee_is_finite(anchor_rotmat)) .or. any(.not. ieee_is_finite(anchor_shift)) ) &
                &THROW_HARD('cartft_pose_opt cumulative anchor must be finite')
        endif
        self%stages(CARTFT_STAGE_JOINT) = cartft_stage(status=CARTFT_ITERATION_LIMIT)
        associate( st => self%stages(CARTFT_STAGE_JOINT) )
            if( .not. calc%ptcl_is_valid(iptcl) )then
                st%status = CARTFT_INVALID_NUMERICS
                return
            endif
            naccepted = 0
            consecutive_rejections = 0
            ! evaluate the seed once: objective, five-vector gradient and 5x5 Gauss-Newton matrix
            call calc%pose_normal_terms(state, iseven, iptcl, rotmat, shift, objective, gradient, hessian)
            bounded_trial = .false.
            if( .not. ieee_is_finite(objective) .or. any(.not. ieee_is_finite(gradient)) .or. &
                &any(.not. ieee_is_finite(hessian)) )then
                st%status = CARTFT_INVALID_NUMERICS
                return
            endif
            mu = LM_INITIAL_DAMPING
            rejection_multiplier = LM_INITIAL_REJECTION_MULT
            do iteration = 1, self%maxits
                st%niter = iteration
                ! form and solve the small damped LM system, inexpensive next to a Fourier-plane evaluation
                call build_pose_lm_system(gradient, hessian, self%rotation_scale, mu, l_active, scaled_gradient, &
                    &scaled_hessian, diagonal, solve_matrix, scaled_direction, direction, identifiable, &
                    &stationary, reliable, bounded_step, self%shift_step_bound)
                ! stop when the local system is stationary, unidentifiable or numerically
                ! unreliable; retain any improvement already accepted
                if( .not. identifiable )then
                    st%status = merge(CARTFT_ACCEPTED, CARTFT_NO_RELIABLE_UPDATE, naccepted > 0)
                    exit
                endif
                if( stationary )then
                    st%status = merge(CARTFT_ACCEPTED, CARTFT_NO_IMPROVEMENT, naccepted > 0)
                    exit
                endif
                if( .not. reliable )then
                    st%status = CARTFT_NO_RELIABLE_UPDATE
                    exit
                endif
                if( any(.not. ieee_is_finite(direction)) )then
                    st%status = CARTFT_INVALID_NUMERICS
                    exit
                endif
                bounded_trial = bounded_trial .or. bounded_step
                if( bounded_step ) st%nbound_hits = st%nbound_hits + 1
                rotation_norm = sqrt(dot_product(direction(1:3), direction(1:3)))
                shift_norm    = sqrt(dot_product(direction(4:5), direction(4:5)))
                st%max_rotation_step = max(st%max_rotation_step, rotation_norm)
                st%max_shift_step    = max(st%max_shift_step, shift_norm)
                ! reduction predicted by the local quadratic model: -g^T d - 1/2 d^T H d
                predicted = -dot_product(gradient, direction) - 0.5_dp*dot_product(direction, matmul(hessian, direction))
                if( .not. ieee_is_finite(predicted) )then
                    st%status = CARTFT_INVALID_NUMERICS
                    exit
                else if( predicted <= 0._dp )then
                    consecutive_rejections = consecutive_rejections + 1
                    call increase_lm_damping(mu, rejection_multiplier)
                    if( consecutive_rejections >= LM_MAX_CONSECUTIVE_REJECTIONS )then
                        st%status = rejection_terminal_status(naccepted, bounded_trial)
                        exit
                    endif
                    cycle
                endif
                ! the proposed SO(3) and shift increments, without changing the accepted pose yet
                trial_rotmat  = right_increment_rotation(rotmat, direction(1:3))
                trial_shift   = shift + direction(4:5)
                st%nattempted = st%nattempted + 1
                if( cumulative_guard )then
                    ! reject proposals that leave the capture basin of the anchor; increase
                    ! damping and try a smaller step
                    sine_half = sqrt(sum((trial_rotmat - anchor_rotmat)**2))/(2._dp*sqrt(2._dp))
                    cumulative_rotation = 2._dp*asin(max(0._dp, min(1._dp, sine_half)))
                    cumulative_shift    = sqrt(sum((trial_shift - anchor_shift)**2))
                    if( cumulative_rotation > self%max_total_rotation + 10._dp*epsilon(1._dp) .or. &
                        &cumulative_shift > self%max_total_shift + 10._dp*epsilon(1._dp) )then
                        consecutive_rejections = consecutive_rejections + 1
                        call increase_lm_damping(mu, rejection_multiplier)
                        if( .not. bounded_step ) st%nbound_hits = st%nbound_hits + 1
                        bounded_trial = .true.
                        if( consecutive_rejections >= LM_MAX_CONSECUTIVE_REJECTIONS )then
                            st%status = rejection_terminal_status(naccepted, bounded_trial)
                            exit
                        endif
                        cycle
                    endif
                endif
                ! the trial objective, gradient and normal matrix over the full active disk:
                ! the dominant cost of an LM iteration
                call calc%pose_normal_terms(state, iseven, iptcl, trial_rotmat, trial_shift, &
                    &trial_objective, trial_gradient, trial_hessian)
                if( .not. ieee_is_finite(trial_objective) .or. any(.not. ieee_is_finite(trial_gradient)) .or. &
                    &any(.not. ieee_is_finite(trial_hessian)) )then
                    consecutive_rejections = consecutive_rejections + 1
                    call increase_lm_damping(mu, rejection_multiplier)
                    st%status = CARTFT_INVALID_NUMERICS
                    if( consecutive_rejections >= LM_MAX_CONSECUTIVE_REJECTIONS ) exit
                    cycle
                endif
                ! accept a trustworthy improvement; otherwise damp and retry from the last accepted pose
                actual = objective - trial_objective
                ratio  = actual/predicted
                accept_trial = actual > 0._dp .and. ratio >= 0.25_dp
                if( accept_trial )then
                    relative_reduction = actual/max(abs(objective), NUMERIC_FLOOR)
                    rotmat    = trial_rotmat
                    shift     = trial_shift
                    objective = trial_objective
                    gradient  = trial_gradient
                    hessian   = trial_hessian
                    naccepted = naccepted + 1
                    consecutive_rejections = 0
                    st%naccepted = naccepted
                    if( ratio > 0.75_dp ) mu = max(mu/2._dp, epsilon(1._dp))
                    rejection_multiplier = LM_INITIAL_REJECTION_MULT
                    st%status = CARTFT_ACCEPTED
                    ! stop after an accepted but negligible step or reduction
                    if( max(rotation_norm, shift_norm) < 1.e-8_dp .or. &
                        &relative_reduction < LM_ACCEPTED_RELATIVE_TOL ) exit
                else
                    consecutive_rejections = consecutive_rejections + 1
                    call increase_lm_damping(mu, rejection_multiplier)
                    if( consecutive_rejections >= LM_MAX_CONSECUTIVE_REJECTIONS )then
                        st%status = rejection_terminal_status(naccepted, bounded_trial)
                        exit
                    endif
                endif
            end do
            if( st%status == CARTFT_ITERATION_LIMIT .and. naccepted == 0 .and. bounded_trial ) &
                &st%status = CARTFT_BOUND_REJECTED
        end associate
    end subroutine refine_joint

    ! RECORD OF THE LAST SOLVE

    !> Terminal status of the last transaction (CARTFT_* outcome).
    pure integer function get_status( self )
        class(cartft_pose_opt), intent(in) :: self
        get_status = self%status
    end function get_status

    !> Objective at the seed and at the returned pose of the last transaction (the seed's
    !! value when nothing was committed; -1 when the seed was never evaluated).
    pure subroutine get_objectives( self, objective_before, objective_after )
        class(cartft_pose_opt), intent(in)  :: self
        real(dp),               intent(out) :: objective_before, objective_after
        objective_before = self%objective_before
        objective_after  = self%objective_after
    end subroutine get_objectives

    !> Seed-to-result motion of the last transaction: geodesic rotation (radians) and shift
    !! (cropped pixels); zero unless a pose was committed.
    pure subroutine get_motion( self, rotation, shift )
        class(cartft_pose_opt), intent(in)  :: self
        real(dp),               intent(out) :: rotation, shift
        rotation = self%cumulative_rotation
        shift    = self%cumulative_shift
    end subroutine get_motion

    !> Record of one stage (CARTFT_STAGE_SHIFT or CARTFT_STAGE_JOINT) of the last solve.
    !! objective_before/after are set by a transaction (-1 after a stage-only call).
    pure subroutine get_stage( self, istage, status, niter, nattempted, naccepted, nbound_hits, &
        &max_rotation_step, max_shift_step, objective_before, objective_after )
        class(cartft_pose_opt), intent(in)  :: self
        integer,                intent(in)  :: istage
        integer,                intent(out) :: status, niter, nattempted, naccepted, nbound_hits
        real(dp),               intent(out) :: max_rotation_step, max_shift_step
        real(dp),               intent(out) :: objective_before, objective_after
        type(cartft_stage) :: st
        if( istage == CARTFT_STAGE_SHIFT .or. istage == CARTFT_STAGE_JOINT ) st = self%stages(istage)
        status            = st%status
        niter             = st%niter
        nattempted        = st%nattempted
        naccepted         = st%naccepted
        nbound_hits       = st%nbound_hits
        max_rotation_step = st%max_rotation_step
        max_shift_step    = st%max_shift_step
        objective_before  = st%objective_before
        objective_after   = st%objective_after
    end subroutine get_stage

    subroutine reset_record( self )
        class(cartft_pose_opt), intent(inout) :: self
        self%status              = CARTFT_NOT_ATTEMPTED
        self%objective_before    = -1._dp
        self%objective_after     = -1._dp
        self%cumulative_rotation = 0._dp
        self%cumulative_shift    = 0._dp
        self%stages              = cartft_stage()
    end subroutine reset_record

    ! FREE PROCEDURE

    !> R exp([omega]x): apply a tangent-space rotation increment omega (radians) on the right,
    !! the convention of the rotation derivatives in the joint stage and of the fixtures that
    !! perturb poses for it.
    pure function right_increment_rotation( rotmat, omega ) result( updated_rotmat )
        real(dp), intent(in) :: rotmat(3,3), omega(3)
        real(dp) :: updated_rotmat(3,3), skew(3,3), exp_skew(3,3)
        real(dp) :: identity(3,3), theta2, theta4, sinc_theta, cosc_theta
        identity      = 0._dp
        identity(1,1) = 1._dp
        identity(2,2) = 1._dp
        identity(3,3) = 1._dp
        ! [omega]x u = omega x u
        skew   = reshape([0._dp, omega(3), -omega(2), -omega(3), 0._dp, omega(1), omega(2), -omega(1), 0._dp], [3,3])
        theta2 = dot_product(omega, omega)
        if( theta2 < 1.e-8_dp )then
            ! Taylor forms avoid cancellation as the rotation angle approaches zero
            theta4     = theta2*theta2
            sinc_theta = 1._dp - theta2/6._dp + theta4/120._dp
            cosc_theta = 0.5_dp - theta2/24._dp + theta4/720._dp
        else
            sinc_theta = sin(sqrt(theta2))/sqrt(theta2)
            cosc_theta = (1._dp - cos(sqrt(theta2)))/theta2
        endif
        ! Rodrigues' formula of the SO(3) exponential map; right multiplication keeps omega
        ! in the current particle-pose frame
        exp_skew       = identity + sinc_theta*skew + cosc_theta*matmul(skew, skew)
        updated_rotmat = matmul(rotmat, exp_skew)
    end function right_increment_rotation

    ! PRIVATE

    pure logical function positive_finite( val )
        real(dp), intent(in) :: val
        positive_finite = ieee_is_finite(val) .and. val > 0._dp
    end function positive_finite

    !> A stage outcome other than an accepted or finite non-improving result ends the route.
    pure logical function aborts_route( status )
        integer, intent(in) :: status
        select case(status)
            case(CARTFT_ACCEPTED, CARTFT_NO_IMPROVEMENT)
                aborts_route = .false.
            case DEFAULT
                aborts_route = .true.
        end select
    end function aborts_route

    !> Classify a finite rejection limit without discarding an earlier accepted endpoint.
    pure integer function rejection_terminal_status( naccepted, bounded_trial ) result( status )
        integer, intent(in) :: naccepted
        logical, intent(in) :: bounded_trial
        if( naccepted > 0 )then
            status = CARTFT_ACCEPTED
        else if( bounded_trial )then
            status = CARTFT_BOUND_REJECTED
        else
            status = CARTFT_NO_IMPROVEMENT
        endif
    end function rejection_terminal_status

    !> Increase damping aggressively across consecutive rejected proposals; the solver resets
    !! the multiplier after an accepted proposal.
    pure subroutine increase_lm_damping( mu, rejection_multiplier )
        real(dp), intent(inout) :: mu, rejection_multiplier
        mu = min(mu*rejection_multiplier, LM_MAX_DAMPING)
        rejection_multiplier = min(2._dp*rejection_multiplier, LM_MAX_DAMPING)
    end subroutine increase_lm_damping

    !> One scaled, damped and independently bounded pose proposal.
    pure subroutine build_pose_lm_system( gradient, hessian, rotation_scale, mu, active, &
        &scaled_gradient, scaled_hessian, damping_diagonal, solve_matrix, scaled_step, &
        &physical_step, identifiable, stationary, reliable, bounded, shift_step_bound )
        real(dp), intent(in)  :: gradient(5), hessian(5,5), rotation_scale, mu
        logical,  intent(in)  :: active(5)
        real(dp), intent(out) :: scaled_gradient(5), scaled_hessian(5,5)
        real(dp), intent(out) :: damping_diagonal(5), solve_matrix(5,5)
        real(dp), intent(out) :: scaled_step(5), physical_step(5)
        logical,  intent(out) :: identifiable, stationary, reliable, bounded
        real(dp), intent(in)  :: shift_step_bound
        real(dp) :: coordinate_scale(5), ignored_step(5), hessian_scale
        real(dp) :: rotation_norm, shift_norm
        integer  :: axis, jaxis
        coordinate_scale = [rotation_scale, rotation_scale, rotation_scale, 1._dp, 1._dp]
        do axis = 1, 5
            scaled_gradient(axis) = coordinate_scale(axis)*gradient(axis)
            do jaxis = 1, 5
                scaled_hessian(axis,jaxis) = coordinate_scale(axis)*hessian(axis,jaxis)*coordinate_scale(jaxis)
            end do
        end do
        call apply_pose_parameter_mask(scaled_gradient, scaled_hessian, active)
        call solve_pose_cholesky(scaled_hessian, -scaled_gradient, ignored_step, identifiable)
        stationary       = sqrt(dot_product(scaled_gradient, scaled_gradient)) < 1.e-8_dp
        damping_diagonal = 0._dp
        solve_matrix     = scaled_hessian
        scaled_step      = 0._dp
        physical_step    = 0._dp
        reliable         = .false.
        bounded          = .false.
        if( .not. identifiable .or. stationary ) return
        hessian_scale = max(maxval(abs(scaled_hessian)), NUMERIC_FLOOR)
        do axis = 1, 5
            damping_diagonal(axis) = max(scaled_hessian(axis,axis), sqrt(epsilon(1._dp))*hessian_scale, NUMERIC_FLOOR)
            solve_matrix(axis,axis) = solve_matrix(axis,axis) + mu*damping_diagonal(axis)
        end do
        call solve_pose_cholesky(solve_matrix, -scaled_gradient, scaled_step, reliable)
        if( .not. reliable ) return
        physical_step = coordinate_scale*scaled_step
        rotation_norm = sqrt(dot_product(physical_step(1:3), physical_step(1:3)))
        if( rotation_norm > rotation_scale )then
            physical_step(1:3) = physical_step(1:3)*(rotation_scale/rotation_norm)
            bounded = .true.
        endif
        shift_norm = sqrt(dot_product(physical_step(4:5), physical_step(4:5)))
        if( shift_norm > shift_step_bound )then
            physical_step(4:5) = physical_step(4:5)*(shift_step_bound/shift_norm)
            bounded = .true.
        endif
    end subroutine build_pose_lm_system

    !> Freeze inactive pose coordinates while retaining one five-vector LM path.
    pure subroutine apply_pose_parameter_mask( gradient, hessian, active )
        real(dp), intent(inout) :: gradient(5), hessian(5,5)
        logical,  intent(in)    :: active(5)
        integer :: axis
        do axis = 1, 5
            if( active(axis) ) cycle
            gradient(axis)     = 0._dp
            hessian(axis,:)    = 0._dp
            hessian(:,axis)    = 0._dp
            hessian(axis,axis) = 1._dp
        end do
    end subroutine apply_pose_parameter_mask

    !> Cholesky solve with a relative pivot test for a 5x5 symmetric positive-definite block.
    pure subroutine solve_pose_cholesky( matrix, rhs, solution, reliable )
        real(dp), intent(in)  :: matrix(5,5), rhs(5)
        real(dp), intent(out) :: solution(5)
        logical,  intent(out) :: reliable
        real(dp) :: lower(5,5), intermediate(5), pivot, pivot_floor, matrix_scale
        integer  :: i, j
        solution     = 0._dp
        intermediate = 0._dp
        lower        = 0._dp
        reliable     = .false.
        if( any(.not. ieee_is_finite(matrix)) .or. any(.not. ieee_is_finite(rhs)) ) return
        matrix_scale = maxval(abs(matrix))
        if( matrix_scale <= NUMERIC_FLOOR ) return
        pivot_floor = sqrt(epsilon(1._dp))*matrix_scale
        do i = 1, 5
            do j = 1, i - 1
                lower(i,j) = (matrix(i,j) - dot_product(lower(i,1:j-1), lower(j,1:j-1)))/lower(j,j)
            end do
            pivot = matrix(i,i) - dot_product(lower(i,1:i-1), lower(i,1:i-1))
            if( .not. ieee_is_finite(pivot) .or. pivot <= pivot_floor ) return
            lower(i,i) = sqrt(pivot)
        end do
        do i = 1, 5
            intermediate(i) = (rhs(i) - dot_product(lower(i,1:i-1), intermediate(1:i-1)))/lower(i,i)
        end do
        do i = 5, 1, -1
            solution(i) = (intermediate(i) - dot_product(lower(i+1:5,i), solution(i+1:5)))/lower(i,i)
        end do
        reliable = all(ieee_is_finite(solution))
    end subroutine solve_pose_cholesky

end module simple_cartft_pose_opt
