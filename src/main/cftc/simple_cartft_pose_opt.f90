!@descr: bounded Levenberg-Marquardt optimizer of one five-parameter Cartesian particle pose (refine=cont)
! One bounded Levenberg-Marquardt transaction per particle (refine): evaluate the seed, run the
! joint five-parameter solve (lm_solve on a cartft_pose_problem) within trs and athres_cont of the
! seed, and commit only a finite, lower objective; any other outcome returns the seed bit for bit
! with a typed status (CARTFT_*). Rotations update on the right, R exp([omega]x); shifts are
! increments from the stored shift, committed in native-box pixels.
module simple_cartft_pose_opt
use, intrinsic :: ieee_arithmetic, only: ieee_is_finite, ieee_quiet_nan, ieee_value
use, intrinsic :: iso_fortran_env, only: int64
use simple_core_module_api, only: dp, PI, ori, dm2euler, simple_exception
use simple_cartft_calc,     only: cartft_calc
use simple_type_defs,       only: OBJFUN_EUCLID
use simple_lm_solver,       only: lm_problem, lm_options, lm_outcome, lm_solve, LM_ACCEPTED, LM_NO_IMPROVEMENT, &
    &LM_NO_RELIABLE_UPDATE, LM_BOUND_REJECTED, LM_INVALID_NUMERICS, LM_ITERATION_LIMIT
implicit none
private
public :: cartft_pose_opt, cartft_pose_work, right_increment_rotation
! typed per-particle outcomes of a solve or a stage
public :: CARTFT_NOT_ATTEMPTED, CARTFT_INVALID_PREPARATION, CARTFT_ACCEPTED, CARTFT_NO_IMPROVEMENT
public :: CARTFT_NO_RELIABLE_UPDATE, CARTFT_BOUND_REJECTED, CARTFT_INVALID_NUMERICS, CARTFT_ITERATION_LIMIT
public :: run_all_cartft_pose_opt_tests
! stage selectors of get_stage

#include "simple_local_flags.inc"

! the engine's outcomes (simple_lm_solver), one to one, and the transaction's own two
integer, parameter :: CARTFT_NOT_ATTEMPTED       =  0                     !< stage or solve not run
integer, parameter :: CARTFT_INVALID_PREPARATION = -1                     !< the particle slot holds no valid preparation
integer, parameter :: CARTFT_ACCEPTED            = LM_ACCEPTED            !< an objective-reducing pose was committed
integer, parameter :: CARTFT_NO_IMPROVEMENT      = LM_NO_IMPROVEMENT      !< finite, no reduction; seed kept
integer, parameter :: CARTFT_NO_RELIABLE_UPDATE  = LM_NO_RELIABLE_UPDATE  !< unobservable data system; seed kept
integer, parameter :: CARTFT_BOUND_REJECTED      = LM_BOUND_REJECTED      !< every proposal left the step or cumulative bound; seed kept
integer, parameter :: CARTFT_INVALID_NUMERICS    = LM_INVALID_NUMERICS    !< non-finite objective or step; seed kept
integer, parameter :: CARTFT_ITERATION_LIMIT     = LM_ITERATION_LIMIT     !< iteration cap reached without a decision
integer, parameter :: CARTFT_STAGE_SHIFT = 1, CARTFT_STAGE_JOINT = 2

real(dp), parameter :: SHIFT_RESIDUAL_FLOOR = real(epsilon(1.), dp)**2
real(dp), parameter :: NEGLIGIBLE_STEP      = 1.e-8_dp !< radians or cropped pixels: an accepted step this small ends a solve
! default policy
integer,  parameter :: DEFAULT_MAXITS             = 40
real(dp), parameter :: DEFAULT_ROTATION_SCALE     = 0.1_dp
real(dp), parameter :: DEFAULT_SHIFT_STEP_NATIVE  = 1._dp  !< native pixels per proposal
! bench histogram of actual over predicted reduction at accepted steps (>= 0.25 by the acceptance rule)
integer,  parameter :: NRATIO_BINS = 5
real(dp), parameter :: RATIO_EDGES(NRATIO_BINS-1) = [0.75_dp, 1.25_dp, 1.75_dp, 2.5_dp]
character(len=*), parameter :: RATIO_LABELS(NRATIO_BINS) = [character(len=9) :: '0.25-0.75', '0.75-1.25', &
    &'1.25-1.75', '1.75-2.5', '>= 2.5']

!> Record of one LM stage; objectives include the configured prior penalty.
type :: cartft_stage
    integer  :: status         = CARTFT_NOT_ATTEMPTED
    integer  :: niter          = 0
    integer  :: nattempted     = 0     !< trial poses proposed, bound-rejected ones included
    integer  :: ntrials        = 0     !< trial poses evaluated
    integer  :: naccepted      = 0     !< trial poses accepted
    integer  :: nbound_hits    = 0     !< proposals capped or rejected by a bound
    integer  :: nratio(NRATIO_BINS) = 0 !< accepted steps per bin of actual over predicted reduction
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
    logical  :: refine_shifts         = .true. !< false when trs = 0: shifts frozen
    real(dp) :: shift_prior_weight    = 0._dp !< zero disables the Cartesian shift prior
    real(dp) :: rotation_prior_weight = 0._dp !< fixed-width Matrix-Fisher prior
    integer  :: maxits   = DEFAULT_MAXITS  !< iteration cap of each stage
    real(dp) :: rotation_scale     = 0._dp !< radians per joint proposal (scaling and cap)
    real(dp) :: max_total_rotation = 0._dp !< radians from the seed (athres_cont)
    real(dp) :: shift_step_bound   = 0._dp !< cropped pixels per proposal
    real(dp) :: max_total_shift    = 0._dp !< cropped pixels from the seed (trs)
    ! record of the last solve
    integer  :: status = CARTFT_NOT_ATTEMPTED
    real(dp) :: objective_before = -1._dp, objective_after = -1._dp
    real(dp) :: cumulative_rotation = 0._dp, cumulative_shift = 0._dp
    integer  :: npasses = 0                !< calculator evaluations since the last solve began
    type(cartft_stage) :: stages(2)
    real(dp), allocatable :: shell_residual(:) !< whitened |r|^2 per shell at the returned pose (euclid)
    logical  :: exists = .false.
contains
    ! lifecycle
    procedure          :: new
    procedure          :: kill
    procedure          :: set_rotation_prior
    procedure          :: set_shift_prior
    procedure          :: clear_priors
    ! solves
    procedure          :: refine
    procedure          :: refine_pose
    procedure          :: refine_shift
    procedure, private :: refine_joint
    ! record of the last solve
    procedure          :: get_status
    procedure          :: get_objectives
    procedure          :: get_pose_update
    procedure          :: get_stage
    procedure          :: get_shell_residual
    ! private helpers
    procedure, private :: reset_record
end type cartft_pose_opt

!> Work of Cartesian solves summed for a bench record: solves by status, LM counts, calculator
!! passes and pixel-evaluations (passes times the pixels of the particle's shells in the full disk).
type :: cartft_pose_work
    private
    integer(int64) :: nstatus(CARTFT_INVALID_PREPARATION:CARTFT_ITERATION_LIMIT) = 0_int64
    integer(int64) :: nsolved = 0_int64, npasses = 0_int64, npix_evals = 0_int64
    integer(int64) :: niter = 0_int64, ntrials = 0_int64, naccepted = 0_int64, nbound_hits = 0_int64
    integer(int64) :: ncapped = 0_int64                  !< solves whose joint stage used every iteration
    integer(int64) :: nratio(NRATIO_BINS) = 0_int64
contains
    procedure :: add_solve
    procedure :: add
    procedure :: write_bench
end type cartft_pose_work

!> refine_joint's solve for lm_solve: the current and trial pose, the anchor of the cumulative guard,
!! the priors, and the stage statistics the engine does not keep (calc is valid during the solve only)
type, extends(lm_problem) :: cartft_pose_problem
    class(cartft_calc), pointer :: calc => null()
    integer  :: state = 0, iptcl = 0
    logical  :: iseven = .true.
    real(dp) :: rotmat(3,3) = 0._dp, shift(2) = 0._dp, trial_rotmat(3,3) = 0._dp, trial_shift(2) = 0._dp
    logical  :: cumulative_guard = .false.
    real(dp) :: anchor_rotmat(3,3) = 0._dp, anchor_shift(2) = 0._dp
    real(dp) :: shift_prior_weight = 0._dp, rotation_prior_weight = 0._dp
    real(dp) :: rotation_scale = 0._dp, shift_step_bound = 0._dp
    real(dp) :: max_total_rotation = 0._dp, max_total_shift = 0._dp
    real(dp), allocatable :: residual(:), trial_residual(:) !< whitened |r|^2 per shell (with residual_out)
    integer  :: npasses = 0
    integer  :: nratio(NRATIO_BINS) = 0
    real(dp) :: max_rotation_step = 0._dp, max_shift_step = 0._dp
contains
    procedure :: eval            => pose_problem_eval
    procedure :: limit_step      => pose_problem_limit_step
    procedure :: propose         => pose_problem_propose
    procedure :: accept          => pose_problem_accept
    procedure :: penalty         => pose_problem_penalty
    procedure :: observable      => pose_problem_observable
    procedure :: negligible_step => pose_problem_negligible_step
end type cartft_pose_problem

!> refine_shift's solve for lm_solve: the shift at a fixed rotation, bounded from its start
type, extends(lm_problem) :: cartft_shift_problem
    class(cartft_calc), pointer :: calc => null()
    integer  :: state = 0, iptcl = 0
    logical  :: iseven = .true.
    real(dp) :: rotmat(3,3) = 0._dp, shift(2) = 0._dp, trial_shift(2) = 0._dp, anchor_shift(2) = 0._dp
    real(dp) :: shift_step_bound = 0._dp, max_total_shift = 0._dp
    integer  :: npasses = 0
    integer  :: nratio(NRATIO_BINS) = 0
    real(dp) :: max_shift_step = 0._dp
contains
    procedure :: eval            => shift_problem_eval
    procedure :: limit_step      => shift_problem_limit_step
    procedure :: propose         => shift_problem_propose
    procedure :: accept          => shift_problem_accept
    procedure :: negligible_step => shift_problem_negligible_step
end type cartft_shift_problem

interface
    pure module subroutine pose_prior_terms( shift_weight, rotation_weight, rotation, shift, seed_rotation, seed_shift, &
        &penalty, gradient, hessian )
        real(dp), intent(in)  :: shift_weight, rotation_weight, rotation(3,3), shift(2), seed_rotation(3,3), seed_shift(2)
        real(dp), intent(out) :: penalty, gradient(5), hessian(5,5)
    end subroutine pose_prior_terms

    pure module subroutine rotation_prior_terms( weight, rotation, seed, penalty, gradient, hessian )
        real(dp), intent(in)  :: weight, rotation(3,3), seed(3,3)
        real(dp), intent(out) :: penalty, gradient(3), hessian(3,3)
    end subroutine rotation_prior_terms

    pure module subroutine shift_prior_terms( weight, shift, seed, penalty, gradient, hessian )
        real(dp), intent(in)  :: weight, shift(2), seed(2)
        real(dp), intent(out) :: penalty, gradient(2), hessian(2,2)
    end subroutine shift_prior_terms

    !> The data-only observability test of the joint solve: the scaled, masked system factors
    pure module function pose_data_observable( scaled_gradient, scaled_hessian ) result( observable )
        real(dp), intent(in) :: scaled_gradient(:), scaled_hessian(:,:)
        logical :: observable
    end function pose_data_observable

    module subroutine pose_problem_eval( self, f, g, curv )
        class(cartft_pose_problem), intent(inout) :: self
        real(dp),                   intent(out)   :: f, g(:), curv(:,:)
    end subroutine pose_problem_eval

    module subroutine pose_problem_limit_step( self, step, bounded )
        class(cartft_pose_problem), intent(inout) :: self
        real(dp),                   intent(inout) :: step(:)
        logical,                    intent(out)   :: bounded
    end subroutine pose_problem_limit_step

    module function pose_problem_propose( self, step ) result( feasible )
        class(cartft_pose_problem), intent(inout) :: self
        real(dp),                   intent(in)    :: step(:)
        logical :: feasible
    end function pose_problem_propose

    module function pose_problem_accept( self, ratio ) result( changed )
        class(cartft_pose_problem), intent(inout) :: self
        real(dp),                   intent(in)    :: ratio
        logical :: changed
    end function pose_problem_accept

    module function pose_problem_penalty( self, at_trial, penalty, g, curv ) result( has_penalty )
        class(cartft_pose_problem), intent(inout) :: self
        logical,                    intent(in)    :: at_trial
        real(dp),                   intent(out)   :: penalty, g(:), curv(:,:)
        logical :: has_penalty
    end function pose_problem_penalty

    module function pose_problem_observable( self, scaled_g, scaled_curv ) result( observable )
        class(cartft_pose_problem), intent(inout) :: self
        real(dp),                   intent(in)    :: scaled_g(:), scaled_curv(:,:)
        logical :: observable
    end function pose_problem_observable

    module function pose_problem_negligible_step( self, step ) result( negligible )
        class(cartft_pose_problem), intent(inout) :: self
        real(dp),                   intent(in)    :: step(:)
        logical :: negligible
    end function pose_problem_negligible_step

    module subroutine shift_problem_eval( self, f, g, curv )
        class(cartft_shift_problem), intent(inout) :: self
        real(dp),                    intent(out)   :: f, g(:), curv(:,:)
    end subroutine shift_problem_eval

    module subroutine shift_problem_limit_step( self, step, bounded )
        class(cartft_shift_problem), intent(inout) :: self
        real(dp),                    intent(inout) :: step(:)
        logical,                     intent(out)   :: bounded
    end subroutine shift_problem_limit_step

    module function shift_problem_propose( self, step ) result( feasible )
        class(cartft_shift_problem), intent(inout) :: self
        real(dp),                    intent(in)    :: step(:)
        logical :: feasible
    end function shift_problem_propose

    module function shift_problem_accept( self, ratio ) result( changed )
        class(cartft_shift_problem), intent(inout) :: self
        real(dp),                    intent(in)    :: ratio
        logical :: changed
    end function shift_problem_accept

    module function shift_problem_negligible_step( self, step ) result( negligible )
        class(cartft_shift_problem), intent(inout) :: self
        real(dp),                    intent(in)    :: step(:)
        logical :: negligible
    end function shift_problem_negligible_step

    !> white-box tests (submodule simple_cartft_pose_opt_tester, built only with BUILD_TESTS=ON)
    module subroutine run_all_cartft_pose_opt_tests()
    end subroutine run_all_cartft_pose_opt_tests
end interface

contains

    ! LIFECYCLE

    !> Geometry (box of the project shifts, box_crop of the calculator) and bounds (trs the total
    !! shift in native pixels, 0 freezes shifts; rotation_bound the total rotation in degrees).
    subroutine new( self, box, box_crop, trs, rotation_bound, maxits, rotation_scale, shift_step_bound )
        class(cartft_pose_opt), intent(inout) :: self
        integer,                intent(in)    :: box, box_crop
        real,                   intent(in)    :: trs, rotation_bound
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
        self%maxits           = DEFAULT_MAXITS
        self%rotation_scale   = DEFAULT_ROTATION_SCALE
        self%shift_step_bound = DEFAULT_SHIFT_STEP_NATIVE*crop_scale
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
        self%box                   = 0
        self%box_crop              = 0
        self%refine_shifts         = .true.
        self%shift_prior_weight    = 0._dp
        self%rotation_prior_weight = 0._dp
        self%maxits                = DEFAULT_MAXITS
        self%rotation_scale        = 0._dp
        self%max_total_rotation    = 0._dp
        self%shift_step_bound      = 0._dp
        self%max_total_shift       = 0._dp
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

    !> Configure the seed-centred shift prior using native-pixel variance and whitened power.
    !! Gaussian negative-log penalty; SIMPLE crop/power scaling is derived in the plan Numerical contract.
    subroutine set_shift_prior( self, variance, power )
        class(cartft_pose_opt), intent(inout) :: self
        real(dp),              intent(in)    :: variance, power
        real(dp) :: crop_scale, weight
        if( .not. self%exists ) THROW_HARD('cartft_pose_opt used before new')
        if( .not. positive_finite(variance) .or. .not. positive_finite(power) ) &
            &THROW_HARD('cartft_pose_opt shift prior requires positive finite variance and power')
        crop_scale = real(self%box_crop,dp)/real(self%box,dp)
        weight = ((2._dp/variance)/crop_scale)/crop_scale/power
        if( .not. positive_finite(weight) ) THROW_HARD('cartft_pose_opt shift prior weight must be positive and finite')
        self%shift_prior_weight = weight
    end subroutine set_shift_prior

    !> Configure a fixed per-axis rotation width in degrees, converted to physical radians.
    !! Isotropic Matrix-Fisher model: Lee (2018), doi:10.1109/TAC.2018.2797162; local scaling in the plan.
    subroutine set_rotation_prior( self, sigma_degrees, power )
        class(cartft_pose_opt), intent(inout) :: self
        real(dp),              intent(in)    :: sigma_degrees, power
        real(dp) :: sigma, weight
        if( .not. self%exists ) THROW_HARD('cartft_pose_opt used before new')
        if( .not. positive_finite(sigma_degrees) .or. .not. positive_finite(power) ) &
            &THROW_HARD('Cartesian rotation prior requires positive finite width and power')
        sigma = sigma_degrees*(acos(-1._dp)/180._dp)
        weight = (2._dp/sigma)/sigma/power
        if( .not. positive_finite(weight) ) THROW_HARD('Cartesian rotation prior weight must be positive and finite')
        self%rotation_prior_weight = weight
    end subroutine set_rotation_prior

    !> Disable the prior when the optimizer is reused for a different particle.
    subroutine clear_priors( self )
        class(cartft_pose_opt), intent(inout) :: self
        self%shift_prior_weight    = 0._dp
        self%rotation_prior_weight = 0._dp
    end subroutine clear_priors

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
        real(dp) :: staged_rotmat(3,3), staged_shift(2), gradient(5), hessian(5,5)
        real(dp) :: joint_objective_after, sine_half, prior_penalty
        real(dp) :: prior_gradient(5), prior_hessian(5,5)
        real(dp), allocatable :: staged_residual(:), joint_residual(:)
        integer  :: kfromto(2)
        if( .not. self%exists ) THROW_HARD('cartft_pose_opt used before new')
        ! every early exit leaves rotmat and shift at the seed (rollback)
        call self%reset_record
        if( .not. calc%ptcl_is_valid(iptcl) )then
            self%status = CARTFT_INVALID_PREPARATION
            return
        endif
        if( (self%shift_prior_weight > 0._dp .or. self%rotation_prior_weight > 0._dp) .and. &
            &calc%get_ptcl_objfun(iptcl) /= OBJFUN_EUCLID ) &
            &THROW_HARD('cartft_pose_opt shift prior requires Euclidean loss')
        ! the seed objective is the single acceptance baseline of the transaction; its terms also
        ! seed the joint stage, and its shell residual is the returned pose's until a staged pose
        ! is committed
        kfromto = calc%get_ptcl_kfromto(iptcl)
        allocate(self%shell_residual(kfromto(2) - kfromto(1) + 1), joint_residual(kfromto(2) - kfromto(1) + 1))
        call calc%pose_normal_terms(state, iseven, iptcl, rotmat, shift, self%objective_before, gradient, hessian, &
            &self%shell_residual)
        self%npasses = self%npasses + 1
        if( .not. ieee_is_finite(self%objective_before) )then
            self%status = CARTFT_INVALID_NUMERICS
            return
        endif
        staged_rotmat   = rotmat
        staged_shift    = shift
        staged_residual = self%shell_residual
        ! the joint stage from the seed, seeded with its terms and returning the objective at its
        ! endpoint; both cumulative guards are anchored at the seed
        if( self%refine_shifts )then
            call self%refine_joint(calc, state, iseven, iptcl, staged_rotmat, staged_shift, &
                &anchor_rotmat=rotmat, anchor_shift=shift, seed_objective=self%objective_before, &
                &seed_gradient=gradient, seed_hessian=hessian, objective_out=joint_objective_after, &
                &seed_residual=staged_residual, residual_out=joint_residual)
        else
            call self%refine_joint(calc, state, iseven, iptcl, staged_rotmat, staged_shift, &
                &active=[.true., .true., .true., .false., .false.], anchor_rotmat=rotmat, anchor_shift=shift, &
                &seed_objective=self%objective_before, seed_gradient=gradient, seed_hessian=hessian, &
                &objective_out=joint_objective_after, seed_residual=staged_residual, residual_out=joint_residual)
        endif
        associate( st => self%stages(CARTFT_STAGE_JOINT) )
            st%objective_before = self%objective_before
            st%objective_after  = joint_objective_after
            if( aborts_route(st%status) )then
                self%status          = st%status
                self%objective_after = self%objective_before
                return
            endif
        end associate
        ! The transaction seed has zero penalty; final acceptance uses the penalized loss.
        if( (self%shift_prior_weight > 0._dp .or. self%rotation_prior_weight > 0._dp) )then
            call pose_prior_terms(self%shift_prior_weight, self%rotation_prior_weight, staged_rotmat, staged_shift, &
                &rotmat, shift, prior_penalty, prior_gradient, prior_hessian)
            self%stages(CARTFT_STAGE_JOINT)%objective_after = joint_objective_after + prior_penalty
            if( .not. ieee_is_finite(self%stages(CARTFT_STAGE_JOINT)%objective_after) .or. &
                &self%stages(CARTFT_STAGE_JOINT)%objective_after >= self%objective_before )then
                self%status = CARTFT_NO_IMPROVEMENT
                self%objective_after = self%objective_before
                return
            endif
        endif
        ! commit the staged pose only if its final objective beats the seed
        self%objective_after = joint_objective_after
        if( .not. ieee_is_finite(self%objective_before) .or. .not. ieee_is_finite(self%objective_after) )then
            self%status          = CARTFT_INVALID_NUMERICS
            self%objective_after = self%objective_before
            return
        endif
        if( .not.(self%objective_after < self%objective_before) )then
            self%status          = CARTFT_NO_IMPROVEMENT
            self%objective_after = self%objective_before
            return
        endif
        self%status = CARTFT_ACCEPTED
        ! cumulative pose update of the committed transaction
        sine_half = sqrt(sum((staged_rotmat - rotmat)**2))/(2._dp*sqrt(2._dp))
        self%cumulative_rotation = 2._dp*asin(max(0._dp, min(1._dp, sine_half)))
        self%cumulative_shift    = sqrt(sum((staged_shift - shift)**2))
        rotmat = staged_rotmat
        shift  = staged_shift
        self%shell_residual = joint_residual
    end subroutine refine_pose

    !> The shift-only solve at a fixed rotation: lm_solve on a cartft_shift_problem, steps capped at
    !! shift_step_bound and trials within trs of the input shift. shift is updated through accepted steps
    !! only. Overwrites the CARTFT_STAGE_SHIFT record; an invalid particle slot gives CARTFT_INVALID_NUMERICS.
    subroutine refine_shift( self, calc, state, iseven, iptcl, rotmat, shift )
        class(cartft_pose_opt),     intent(inout) :: self
        class(cartft_calc), target, intent(in)    :: calc
        integer,                    intent(in)    :: state, iptcl
        logical,                    intent(in)    :: iseven
        real(dp),                   intent(in)    :: rotmat(3,3)
        real(dp),                   intent(inout) :: shift(2)
        type(cartft_shift_problem) :: problem
        type(lm_options)           :: options
        type(lm_outcome)           :: outcome
        real(dp) :: objective, gradient(2), hessian(2,2)
        if( .not. self%exists ) THROW_HARD('cartft_pose_opt used before new')
        self%stages(CARTFT_STAGE_SHIFT) = cartft_stage(status=CARTFT_ITERATION_LIMIT)
        associate( st => self%stages(CARTFT_STAGE_SHIFT) )
            if( .not. calc%ptcl_is_valid(iptcl) )then
                st%status = CARTFT_INVALID_NUMERICS
                return
            endif
            call calc%shift_normal_terms(state, iseven, iptcl, rotmat, shift, objective, gradient, hessian)
            self%npasses = self%npasses + 1
            if( .not. ieee_is_finite(objective) .or. any(.not. ieee_is_finite(gradient)) .or. &
                &any(.not. ieee_is_finite(hessian)) )then
                st%status = CARTFT_INVALID_NUMERICS
                return
            endif
            ! an exact shift leaves nothing a single-precision observation resolves
            if( objective <= SHIFT_RESIDUAL_FLOOR )then
                st%status = CARTFT_NO_IMPROVEMENT
                return
            endif
            problem%calc             => calc
            problem%state            = state
            problem%iseven           = iseven
            problem%iptcl            = iptcl
            problem%rotmat           = rotmat
            problem%shift            = shift
            problem%trial_shift      = shift
            problem%anchor_shift     = shift
            problem%shift_step_bound = self%shift_step_bound
            problem%max_total_shift  = self%max_total_shift
            options%maxits = self%maxits
            options%scale  = [1._dp, 1._dp]
            options%active = [.true., .true.]
            call lm_solve(problem, options, outcome, objective, gradient, hessian)
            shift = problem%shift
            call copy_outcome(outcome, st)
            st%nratio         = problem%nratio
            st%max_shift_step = problem%max_shift_step
            self%npasses      = self%npasses + problem%npasses
            problem%calc => null()
        end associate
    end subroutine refine_shift

    !> The joint five-parameter stage alone: lm_solve on a cartft_pose_problem, the rotation capped at
    !! rotation_scale and the shift at shift_step_bound per step. active masks parameters (1:3 rotation,
    !! 4:5 shift; default all). With anchor_rotmat and anchor_shift (both or neither) every trial pose must
    !! stay within max_total_rotation and max_total_shift of the anchor. rotmat and shift are updated through
    !! accepted steps only. The seed_* terms (all or none) are the caller's pose_normal_terms at the seed,
    !! which then is not evaluated again; objective_out is the objective at the returned pose, and
    !! residual_out its shell residual (seed_residual, the seed's, with the seed terms). Overwrites the
    !! CARTFT_STAGE_JOINT record; an invalid particle slot gives CARTFT_INVALID_NUMERICS. THROW_HARD for no
    !! active parameter, one anchor without the other, a non-finite anchor or a partial seed.
    subroutine refine_joint( self, calc, state, iseven, iptcl, rotmat, shift, active, anchor_rotmat, anchor_shift, &
        &seed_objective, seed_gradient, seed_hessian, objective_out, seed_residual, residual_out )
        class(cartft_pose_opt),     intent(inout) :: self
        class(cartft_calc), target, intent(in)    :: calc
        integer,                    intent(in)    :: state, iptcl
        logical,                    intent(in)    :: iseven
        real(dp),                   intent(inout) :: rotmat(3,3), shift(2)
        logical,  optional,         intent(in)    :: active(5)
        real(dp), optional,         intent(in)    :: anchor_rotmat(3,3), anchor_shift(2)
        real(dp), optional,         intent(in)    :: seed_objective, seed_gradient(5), seed_hessian(5,5)
        real(dp), optional,         intent(out)   :: objective_out
        real(dp), optional,         intent(in)    :: seed_residual(:)
        real(dp), optional,         intent(out)   :: residual_out(:)
        type(cartft_pose_problem) :: problem
        type(lm_options)          :: options
        type(lm_outcome)          :: outcome
        real(dp) :: objective, gradient(5), hessian(5,5)
        logical  :: l_active(5)
        if( .not. self%exists ) THROW_HARD('cartft_pose_opt used before new')
        l_active = .true.
        if( present(active) ) l_active = active
        if( .not. any(l_active) ) THROW_HARD('cartft_pose_opt joint stage requires one active parameter')
        if( present(anchor_rotmat) .neqv. present(anchor_shift) ) &
            &THROW_HARD('cartft_pose_opt joint stage requires both anchors or none')
        if( (self%shift_prior_weight > 0._dp .or. self%rotation_prior_weight > 0._dp) .and. .not. present(anchor_shift) ) &
            &THROW_HARD('cartft_pose_opt shift prior requires a transaction anchor')
        problem%cumulative_guard = present(anchor_rotmat)
        if( problem%cumulative_guard )then
            if( any(.not. ieee_is_finite(anchor_rotmat)) .or. any(.not. ieee_is_finite(anchor_shift)) ) &
                &THROW_HARD('cartft_pose_opt cumulative anchor must be finite')
            problem%anchor_rotmat = anchor_rotmat
            problem%anchor_shift  = anchor_shift
        endif
        if( (present(seed_objective) .neqv. present(seed_gradient)) .or. (present(seed_objective) .neqv. present(seed_hessian)) ) &
            &THROW_HARD('cartft_pose_opt joint stage requires all seed terms or none')
        if( present(objective_out) ) objective_out = ieee_value(1._dp, ieee_quiet_nan)
        if( present(residual_out) )then
            if( present(seed_objective) .neqv. present(seed_residual) ) &
                &THROW_HARD('cartft_pose_opt joint stage requires the seed residual with the seed terms')
            residual_out = 0._dp
            allocate(problem%residual(size(residual_out)), problem%trial_residual(size(residual_out)), source=0._dp)
        endif
        self%stages(CARTFT_STAGE_JOINT) = cartft_stage(status=CARTFT_ITERATION_LIMIT)
        associate( st => self%stages(CARTFT_STAGE_JOINT) )
            if( .not. calc%ptcl_is_valid(iptcl) )then
                st%status = CARTFT_INVALID_NUMERICS
                return
            endif
            ! the seed's objective, five-vector gradient and 5x5 Gauss-Newton matrix, evaluated once
            if( present(seed_objective) )then
                objective = seed_objective
                gradient  = seed_gradient
                hessian   = seed_hessian
                if( allocated(problem%residual) ) problem%residual = seed_residual
            else
                call calc%pose_normal_terms(state, iseven, iptcl, rotmat, shift, objective, gradient, hessian, &
                    &problem%residual)
                self%npasses = self%npasses + 1
            endif
            if( present(objective_out) ) objective_out = objective
            if( present(residual_out) ) residual_out = problem%residual
            problem%calc                  => calc
            problem%state                 = state
            problem%iseven                = iseven
            problem%iptcl                 = iptcl
            problem%rotmat                = rotmat
            problem%shift                 = shift
            problem%trial_rotmat          = rotmat
            problem%trial_shift           = shift
            problem%shift_prior_weight    = self%shift_prior_weight
            problem%rotation_prior_weight = self%rotation_prior_weight
            problem%rotation_scale        = self%rotation_scale
            problem%shift_step_bound      = self%shift_step_bound
            problem%max_total_rotation    = self%max_total_rotation
            problem%max_total_shift       = self%max_total_shift
            options%maxits = self%maxits
            options%scale  = [self%rotation_scale, self%rotation_scale, self%rotation_scale, 1._dp, 1._dp]
            options%active = l_active
            call lm_solve(problem, options, outcome, objective, gradient, hessian)
            problem%calc => null()
            call copy_outcome(outcome, st)
            st%nratio            = problem%nratio
            st%max_rotation_step = problem%max_rotation_step
            st%max_shift_step    = problem%max_shift_step
            self%npasses         = self%npasses + problem%npasses
            rotmat = problem%rotmat
            shift  = problem%shift
        end associate
        if( present(objective_out) ) objective_out = outcome%f
        if( present(residual_out) ) residual_out = problem%residual
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

    !> Seed-to-result pose update of the last transaction: geodesic rotation (radians) and shift
    !! (cropped pixels); zero unless a pose was committed.
    pure subroutine get_pose_update( self, rotation, shift )
        class(cartft_pose_opt), intent(in)  :: self
        real(dp),               intent(out) :: rotation, shift
        rotation = self%cumulative_rotation
        shift    = self%cumulative_shift
    end subroutine get_pose_update

    !> Record of one stage (CARTFT_STAGE_SHIFT or CARTFT_STAGE_JOINT) of the last solve.
    !! Objectives include the prior penalty; transaction getters retain data-only loss.
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

    !> The whitened |r|^2 per shell of the slot's band at the pose the last transaction returned,
    !! from the evaluation that chose it (calc%sigma_from_shell_residual turns it into the sigma2
    !! contribution); unallocated when the transaction evaluated nothing.
    subroutine get_shell_residual( self, shell_residual )
        class(cartft_pose_opt), intent(in)  :: self
        real(dp), allocatable,  intent(out) :: shell_residual(:)
        if( allocated(self%shell_residual) ) shell_residual = self%shell_residual
    end subroutine get_shell_residual

    subroutine reset_record( self )
        class(cartft_pose_opt), intent(inout) :: self
        self%status              = CARTFT_NOT_ATTEMPTED
        self%objective_before    = -1._dp
        self%objective_after     = -1._dp
        self%cumulative_rotation = 0._dp
        self%cumulative_shift    = 0._dp
        self%npasses             = 0
        self%stages              = cartft_stage()
        if( allocated(self%shell_residual) ) deallocate(self%shell_residual)
    end subroutine reset_record

    ! WORK RECORD

    !> Add the last solve of opt: npix pixels per evaluation.
    subroutine add_solve( self, opt, npix )
        class(cartft_pose_work), intent(inout) :: self
        class(cartft_pose_opt),  intent(in)    :: opt
        integer,                 intent(in)    :: npix
        integer :: istage, npasses
        if( opt%status < lbound(self%nstatus,1) .or. opt%status > ubound(self%nstatus,1) ) &
            &THROW_HARD('cartft_pose_work received an unknown solve status')
        self%nstatus(opt%status) = self%nstatus(opt%status) + 1_int64
        if( opt%npasses > 0 ) self%nsolved = self%nsolved + 1_int64
        npasses         = opt%npasses
        self%npasses    = self%npasses    + int(npasses, int64)
        self%npix_evals = self%npix_evals + int(npasses, int64)*int(npix, int64)
        do istage = CARTFT_STAGE_SHIFT, CARTFT_STAGE_JOINT
            associate( st => opt%stages(istage) )
                self%niter       = self%niter       + int(st%niter, int64)
                self%ntrials     = self%ntrials     + int(st%ntrials, int64)
                self%naccepted   = self%naccepted   + int(st%naccepted, int64)
                self%nbound_hits = self%nbound_hits + int(st%nbound_hits, int64)
                self%nratio      = self%nratio      + int(st%nratio, int64)
            end associate
        end do
        if( opt%stages(CARTFT_STAGE_JOINT)%niter >= opt%maxits ) self%ncapped = self%ncapped + 1_int64
    end subroutine add_solve

    !> Add another record, as the threads' records are merged.
    subroutine add( self, other )
        class(cartft_pose_work), intent(inout) :: self
        class(cartft_pose_work), intent(in)    :: other
        self%nstatus     = self%nstatus     + other%nstatus
        self%nsolved     = self%nsolved     + other%nsolved
        self%npasses     = self%npasses     + other%npasses
        self%npix_evals  = self%npix_evals  + other%npix_evals
        self%niter       = self%niter       + other%niter
        self%ntrials     = self%ntrials     + other%ntrials
        self%naccepted   = self%naccepted   + other%naccepted
        self%nbound_hits = self%nbound_hits + other%nbound_hits
        self%ncapped     = self%ncapped     + other%ncapped
        self%nratio      = self%nratio      + other%nratio
    end subroutine add

    !> The record as a bench section, keys prefixed by label. ts_search is the search's thread-seconds,
    !! so the time per pixel-evaluation covers the whole particle search, sigma2 residual included.
    subroutine write_bench( self, fnr, label, ts_search )
        class(cartft_pose_work), intent(in) :: self
        integer,                 intent(in) :: fnr
        character(len=*),        intent(in) :: label
        real(dp),                intent(in) :: ts_search
        character(len=40) :: key
        integer :: ibin
        write(fnr,'(a)') '*** CONTINUOUS POSE WORK ***'
        call put_count('solves',              sum(self%nstatus))
        call put_count('accepted',            self%nstatus(CARTFT_ACCEPTED))
        call put_count('no improvement',      self%nstatus(CARTFT_NO_IMPROVEMENT))
        call put_count('no reliable update',  self%nstatus(CARTFT_NO_RELIABLE_UPDATE))
        call put_count('bound rejected',      self%nstatus(CARTFT_BOUND_REJECTED))
        call put_count('invalid numerics',    self%nstatus(CARTFT_INVALID_NUMERICS))
        call put_count('iteration limit',     self%nstatus(CARTFT_ITERATION_LIMIT))
        call put_count('invalid preparation', self%nstatus(CARTFT_INVALID_PREPARATION))
        call put_count('not attempted',       self%nstatus(CARTFT_NOT_ATTEMPTED))
        call put_count('LM iterations',       self%niter)
        call put_count('trial evaluations',   self%ntrials)
        call put_count('accepted steps',      self%naccepted)
        call put_count('bound hits',          self%nbound_hits)
        call put_count('calculator passes',   self%npasses)
        call put_count('pixel-evaluations',   self%npix_evals)
        if( self%nsolved > 0_int64 )then
            call put_real('passes per solve',        real(self%npasses, dp)/real(self%nsolved, dp))
            call put_real('LM iterations per solve', real(self%niter, dp)/real(self%nsolved, dp))
        endif
        call put_count('solves at the iteration cap', self%ncapped)
        do ibin = 1, NRATIO_BINS
            call put_count('actual/predicted '//trim(RATIO_LABELS(ibin)), self%nratio(ibin))
        end do
        if( self%npix_evals > 0_int64 ) call put_real('ns per pixel-evaluation', 1.e9_dp*ts_search/real(self%npix_evals, dp))

    contains

        subroutine put_count( name, n )
            character(len=*), intent(in) :: name
            integer(int64),   intent(in) :: n
            key = label//' cont '//name
            write(fnr,'(a,a,1x,i0)') key, ':', n
        end subroutine put_count

        subroutine put_real( name, val )
            character(len=*), intent(in) :: name
            real(dp),         intent(in) :: val
            key = label//' cont '//name
            write(fnr,'(a,a,1x,f0.4)') key, ':', val
        end subroutine put_real

    end subroutine write_bench

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

    !> The bench histogram bin of an accepted step's actual over predicted reduction.
    pure integer function ratio_bin( ratio )
        real(dp), intent(in) :: ratio
        ratio_bin = 1 + count(ratio >= RATIO_EDGES)
    end function ratio_bin

    pure logical function positive_finite( val )
        real(dp), intent(in) :: val
        positive_finite = ieee_is_finite(val) .and. val > 0._dp
    end function positive_finite

    !> A stage outcome other than an accepted or finite non-improving result ends the transaction.
    pure logical function aborts_route( status )
        integer, intent(in) :: status
        select case(status)
            case(CARTFT_ACCEPTED, CARTFT_NO_IMPROVEMENT)
                aborts_route = .false.
            case DEFAULT
                aborts_route = .true.
        end select
    end function aborts_route

    !> The engine's status and counters into a stage record.
    pure subroutine copy_outcome( outcome, st )
        type(lm_outcome),   intent(in)    :: outcome
        type(cartft_stage), intent(inout) :: st
        st%status      = outcome%status
        st%niter       = outcome%niter
        st%nattempted  = outcome%nattempted
        st%ntrials     = outcome%ntrials
        st%naccepted   = outcome%naccepted
        st%nbound_hits = outcome%nbound_hits
    end subroutine copy_outcome


end module simple_cartft_pose_opt
