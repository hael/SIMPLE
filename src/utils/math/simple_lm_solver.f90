!@descr: the library's Levenberg-Marquardt engine: damped Gauss-Newton steps on a client-supplied problem
! lm_solve owns the algorithm only (scaling and masking, the damped system, damping, acceptance, stops);
! the lm_problem client owns its point, its terms, its step limits, its trial point and its constraints.
! The contract: doc/implementation_notes/planned/lm_optimizer.md, decisions 1, 3 and 6, rulings R3 and R4.
module simple_lm_solver
use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
use simple_defs,   only: dp
use simple_error,  only: simple_exception
use simple_linalg, only: spd_solve_rel
implicit none
private
#include "simple_local_flags.inc"

public :: lm_problem, lm_options, lm_outcome, lm_solve, lm_scale_system
public :: LM_ACCEPTED, LM_NO_IMPROVEMENT, LM_NO_RELIABLE_UPDATE, LM_BOUND_REJECTED, LM_INVALID_NUMERICS
public :: LM_ITERATION_LIMIT

! outcomes
integer, parameter :: LM_ACCEPTED           = 1 !< an objective-reducing point was accepted
integer, parameter :: LM_NO_IMPROVEMENT     = 2 !< finite, no reduction; start kept
integer, parameter :: LM_NO_RELIABLE_UPDATE = 3 !< the client's observability test vetoed the solve
integer, parameter :: LM_BOUND_REJECTED     = 4 !< every proposal was capped or infeasible; start kept
integer, parameter :: LM_INVALID_NUMERICS   = 5 !< non-finite terms or step
integer, parameter :: LM_ITERATION_LIMIT    = 6 !< maxits reached without a decision

! step control
real(dp), parameter :: LM_NUMERIC_FLOOR              = epsilon(1._dp)**2
real(dp), parameter :: LM_INITIAL_DAMPING            = 1.e-3_dp
real(dp), parameter :: LM_INITIAL_REJECTION_MULT     = 4._dp
real(dp), parameter :: LM_MAX_DAMPING                = 1._dp/epsilon(1._dp)
real(dp), parameter :: LM_ACCEPTED_RELATIVE_TOL      = 1.e-6_dp !< the default stops (lm_options)
real(dp), parameter :: LM_PREDICTED_RELATIVE_TOL     = 1.e-6_dp
real(dp), parameter :: LM_STATIONARY_TOL             = 1.e-8_dp
real(dp), parameter :: LM_MIN_ACCEPTANCE_RATIO       = 0.25_dp
real(dp), parameter :: LM_DAMPING_DECREASE_RATIO     = 0.75_dp
integer,  parameter :: LM_MAX_CONSECUTIVE_REJECTIONS = 8

!> The client of lm_solve. It holds the current point and a trial point (the start, until the first
!! propose); eval evaluates the trial point, propose forms it from the current one, accept adopts it
type, abstract :: lm_problem
contains
    procedure(lm_eval_i),       deferred :: eval
    procedure(lm_limit_step_i), deferred :: limit_step
    procedure(lm_propose_i),    deferred :: propose
    procedure(lm_accept_i),     deferred :: accept
    procedure :: penalty         => lm_no_penalty
    procedure :: observable      => lm_always_observable
    procedure :: free            => lm_all_free
    procedure :: negligible_step => lm_no_negligible_step
end type lm_problem

!> maxits, the stops, and per coordinate the scale of the damped system and whether the solve may move it
type :: lm_options
    integer               :: maxits         = 40
    real(dp)              :: predicted_tol  = LM_PREDICTED_RELATIVE_TOL !< relative to |f|: a smaller prediction ends it unevaluated
    real(dp)              :: accepted_tol   = LM_ACCEPTED_RELATIVE_TOL  !< relative to |f|: a smaller accepted reduction ends it
    real(dp)              :: stationary_tol = LM_STATIONARY_TOL         !< on the masked, scaled gradient
    real(dp), allocatable :: scale(:)
    logical,  allocatable :: active(:)
end type lm_options

!> The status (LM_*), the counters and the objective at the returned current point (prior excluded)
type :: lm_outcome
    integer  :: status      = LM_ITERATION_LIMIT
    integer  :: niter       = 0 !< iterations begun
    integer  :: nevals      = 0 !< evaluations by eval (seed terms not counted)
    integer  :: nattempted  = 0 !< trial points proposed, infeasible ones included
    integer  :: ntrials     = 0 !< trial points evaluated
    integer  :: naccepted   = 0
    integer  :: nbound_hits = 0 !< steps capped by limit_step or rejected by propose
    real(dp) :: f           = 0._dp
end type lm_outcome

abstract interface
    !> f, gradient and positive semi-definite curvature at the trial point
    subroutine lm_eval_i( self, f, g, curv )
        import :: dp, lm_problem
        class(lm_problem), intent(inout) :: self
        real(dp),          intent(out)   :: f, g(:), curv(:,:)
    end subroutine lm_eval_i

    !> may shorten the step (physical coordinates); bounded reports whether it did
    subroutine lm_limit_step_i( self, step, bounded )
        import :: dp, lm_problem
        class(lm_problem), intent(inout) :: self
        real(dp),          intent(inout) :: step(:)
        logical,           intent(out)   :: bounded
    end subroutine lm_limit_step_i

    !> forms the trial point from the current one and the step; false when it is infeasible
    logical function lm_propose_i( self, step ) result( feasible )
        import :: dp, lm_problem
        class(lm_problem), intent(inout) :: self
        real(dp),          intent(in)    :: step(:)
    end function lm_propose_i

    !> makes the trial point current; ratio is its actual over predicted reduction; true when the
    !! objective at the new current point changed (the engine then evaluates it again)
    logical function lm_accept_i( self, ratio ) result( changed )
        import :: dp, lm_problem
        class(lm_problem), intent(inout) :: self
        real(dp),          intent(in)    :: ratio
    end function lm_accept_i
end interface

contains

    ! DEFAULT HOOKS

    !> no penalty: false, with zero terms
    logical function lm_no_penalty( self, at_trial, penalty, g, curv ) result( has_penalty )
        class(lm_problem), intent(inout) :: self
        logical,           intent(in)    :: at_trial
        real(dp),          intent(out)   :: penalty, g(:), curv(:,:)
        penalty     = 0._dp
        g           = 0._dp
        curv        = 0._dp
        has_penalty = .false.
    end function lm_no_penalty

    !> every scaled, masked data system is solved
    logical function lm_always_observable( self, scaled_g, scaled_curv ) result( observable )
        class(lm_problem), intent(inout) :: self
        real(dp),          intent(in)    :: scaled_g(:), scaled_curv(:,:)
        observable = .true.
    end function lm_always_observable

    !> every coordinate may move
    subroutine lm_all_free( self, g, free )
        class(lm_problem), intent(inout) :: self
        real(dp),          intent(in)    :: g(:)
        logical,           intent(out)   :: free(:)
        free = .true.
    end subroutine lm_all_free

    !> no accepted step is too small to continue from
    logical function lm_no_negligible_step( self, step ) result( negligible )
        class(lm_problem), intent(inout) :: self
        real(dp),          intent(in)    :: step(:)
        negligible = .false.
    end function lm_no_negligible_step

    ! ENGINE

    !> Minimise the problem's f (plus its penalty) from its start by scaled, damped Gauss-Newton steps on
    !! the coordinates options%active and free. The seed terms (all or none) are f, g and curv at the start,
    !! which then is not evaluated. THROW_HARD for inconsistent sizes, maxits < 1 or a partial seed.
    subroutine lm_solve( problem, options, outcome, seed_f, seed_g, seed_curv )
        class(lm_problem),  intent(inout) :: problem
        type(lm_options),   intent(in)    :: options
        type(lm_outcome),   intent(out)   :: outcome
        real(dp), optional, intent(in)    :: seed_f, seed_g(:), seed_curv(:,:)
        real(dp), allocatable :: g(:), curv(:,:), trial_g(:), trial_curv(:,:), total_g(:), total_curv(:,:)
        real(dp), allocatable :: penalty_g(:), penalty_curv(:,:), scaled_g(:), scaled_curv(:,:)
        real(dp), allocatable :: solve_matrix(:,:), scaled_step(:), step(:)
        logical,  allocatable :: lfree(:), mask(:)
        real(dp) :: f, trial_f, total_f, trial_total, penalty, mu, rejection_multiplier
        real(dp) :: predicted, actual, ratio, relative_reduction, curv_scale
        integer  :: n, iteration, axis, consecutive_rejections
        logical  :: bounded_trial, bounded_step, reliable
        n = size(options%scale)
        if( n < 1 ) THROW_HARD('lm_solve needs at least one coordinate')
        if( size(options%active) /= n ) THROW_HARD('lm_solve scale and active differ in size')
        if( options%maxits < 1 ) THROW_HARD('lm_solve requires at least one iteration')
        if( .not. (options%predicted_tol >= 0._dp .and. options%accepted_tol >= 0._dp .and. options%stationary_tol >= 0._dp) ) &
            &THROW_HARD('lm_solve tolerances must be >= 0')
        if( (present(seed_f) .neqv. present(seed_g)) .or. (present(seed_f) .neqv. present(seed_curv)) ) &
            &THROW_HARD('lm_solve requires all seed terms or none')
        allocate(g(n), trial_g(n), total_g(n), penalty_g(n), scaled_g(n), scaled_step(n), step(n), source=0._dp)
        allocate(curv(n,n), trial_curv(n,n), total_curv(n,n), penalty_curv(n,n), scaled_curv(n,n), &
            &solve_matrix(n,n), source=0._dp)
        allocate(lfree(n), mask(n))
        outcome%status = LM_ITERATION_LIMIT
        if( present(seed_f) )then
            if( size(seed_g) /= n .or. any(shape(seed_curv) /= [n,n]) ) THROW_HARD('lm_solve seed terms differ in size')
            f    = seed_f
            g    = seed_g
            curv = seed_curv
        else
            call problem%eval(f, g, curv)
            outcome%nevals = outcome%nevals + 1
        endif
        outcome%f = f
        if( .not. ieee_is_finite(f) .or. any(.not. ieee_is_finite(g)) .or. any(.not. ieee_is_finite(curv)) )then
            outcome%status = LM_INVALID_NUMERICS
            return
        endif
        mu                     = LM_INITIAL_DAMPING
        rejection_multiplier   = LM_INITIAL_REJECTION_MULT
        consecutive_rejections = 0
        bounded_trial          = .false.
        do iteration = 1, options%maxits
            outcome%niter = iteration
            call problem%free(g, lfree)
            mask = options%active .and. lfree
            ! the client's observability test sees the data system before any penalty curvature
            call lm_scale_system(g, curv, options%scale, mask, scaled_g, scaled_curv)
            if( .not. problem%observable(scaled_g, scaled_curv) )then
                outcome%status = merge(LM_ACCEPTED, LM_NO_RELIABLE_UPDATE, outcome%naccepted > 0)
                exit
            endif
            total_f    = f
            total_g    = g
            total_curv = curv
            if( problem%penalty(.false., penalty, penalty_g, penalty_curv) )then
                total_f    = f + penalty
                total_g    = g + penalty_g
                total_curv = curv + penalty_curv
                call lm_scale_system(total_g, total_curv, options%scale, mask, scaled_g, scaled_curv)
                if( .not. ieee_is_finite(total_f) .or. any(.not. ieee_is_finite(total_g)) .or. &
                    &any(.not. ieee_is_finite(total_curv)) )then
                    outcome%status = LM_INVALID_NUMERICS
                    exit
                endif
            endif
            ! keep any improvement already accepted when the solve stops
            if( sqrt(dot_product(scaled_g, scaled_g)) < options%stationary_tol )then
                outcome%status = merge(LM_ACCEPTED, LM_NO_IMPROVEMENT, outcome%naccepted > 0)
                exit
            endif
            ! a coordinate the cost does not depend on (zero gradient and curvature) takes a zero step
            do axis = 1, n
                if( scaled_g(axis) /= 0._dp .or. scaled_curv(axis,axis) /= 0._dp ) cycle
                scaled_curv(axis,:)    = 0._dp
                scaled_curv(:,axis)    = 0._dp
                scaled_curv(axis,axis) = 1._dp
            end do
            ! C + mu diag(d), d_i = max(C_ii, sqrt(eps) max|C|, eps**2): definite even where C is singular
            curv_scale   = max(maxval(abs(scaled_curv)), LM_NUMERIC_FLOOR)
            solve_matrix = scaled_curv
            do axis = 1, n
                solve_matrix(axis,axis) = solve_matrix(axis,axis) + &
                    &mu*max(scaled_curv(axis,axis), sqrt(epsilon(1._dp))*curv_scale, LM_NUMERIC_FLOOR)
            end do
            call spd_solve_rel(solve_matrix, -scaled_g, scaled_step, reliable)
            if( .not. reliable )then
                ! an unreliable damped solve counts as a rejection: damping grows and the step shrinks
                call reject
                if( consecutive_rejections >= LM_MAX_CONSECUTIVE_REJECTIONS ) exit
                cycle
            endif
            step = options%scale*scaled_step
            call problem%limit_step(step, bounded_step)
            if( any(.not. ieee_is_finite(step)) )then
                outcome%status = LM_INVALID_NUMERICS
                exit
            endif
            bounded_trial = bounded_trial .or. bounded_step
            if( bounded_step ) outcome%nbound_hits = outcome%nbound_hits + 1
            ! reduction predicted by the local quadratic model: -g^T d - 1/2 d^T C d
            predicted = -dot_product(total_g, step) - 0.5_dp*dot_product(step, matmul(total_curv, step))
            if( .not. ieee_is_finite(predicted) )then
                outcome%status = LM_INVALID_NUMERICS
                exit
            else if( predicted <= 0._dp )then
                call reject
                if( consecutive_rejections >= LM_MAX_CONSECUTIVE_REJECTIONS ) exit
                cycle
            endif
            ! stop unevaluated below what the objective resolves; the undamped step predicts more,
            ! so this also ends a converged solve
            if( predicted < options%predicted_tol*max(abs(total_f), LM_NUMERIC_FLOOR) )then
                outcome%status = rejection_terminal_status(outcome%naccepted, bounded_trial)
                exit
            endif
            outcome%nattempted = outcome%nattempted + 1
            if( .not. problem%propose(step) )then
                ! an infeasible trial: more damping and a smaller step
                if( .not. bounded_step ) outcome%nbound_hits = outcome%nbound_hits + 1
                bounded_trial = .true.
                call reject
                if( consecutive_rejections >= LM_MAX_CONSECUTIVE_REJECTIONS ) exit
                cycle
            endif
            call problem%eval(trial_f, trial_g, trial_curv)
            outcome%nevals  = outcome%nevals  + 1
            outcome%ntrials = outcome%ntrials + 1
            if( .not. ieee_is_finite(trial_f) .or. any(.not. ieee_is_finite(trial_g)) .or. &
                &any(.not. ieee_is_finite(trial_curv)) )then
                consecutive_rejections = consecutive_rejections + 1
                call increase_lm_damping(mu, rejection_multiplier)
                outcome%status = LM_INVALID_NUMERICS
                if( consecutive_rejections >= LM_MAX_CONSECUTIVE_REJECTIONS ) exit
                cycle
            endif
            trial_total = trial_f
            if( problem%penalty(.true., penalty, penalty_g, penalty_curv) ) trial_total = trial_f + penalty
            actual = total_f - trial_total
            ratio  = actual/predicted
            if( actual > 0._dp .and. ratio >= LM_MIN_ACCEPTANCE_RATIO )then
                relative_reduction = actual/max(abs(total_f), LM_NUMERIC_FLOOR)
                f    = trial_f
                g    = trial_g
                curv = trial_curv
                if( problem%accept(ratio) )then
                    call problem%eval(f, g, curv)
                    outcome%nevals = outcome%nevals + 1
                    if( .not. ieee_is_finite(f) .or. any(.not. ieee_is_finite(g)) .or. &
                        &any(.not. ieee_is_finite(curv)) )then
                        outcome%status = LM_INVALID_NUMERICS
                        exit
                    endif
                endif
                outcome%naccepted = outcome%naccepted + 1
                consecutive_rejections = 0
                if( ratio > LM_DAMPING_DECREASE_RATIO ) mu = max(mu/2._dp, epsilon(1._dp))
                rejection_multiplier = LM_INITIAL_REJECTION_MULT
                outcome%status = LM_ACCEPTED
                ! stop after an accepted but negligible step or reduction
                if( problem%negligible_step(step) .or. relative_reduction < options%accepted_tol ) exit
            else
                call reject
                if( consecutive_rejections >= LM_MAX_CONSECUTIVE_REJECTIONS ) exit
            endif
        end do
        if( outcome%status == LM_ITERATION_LIMIT .and. outcome%naccepted == 0 .and. bounded_trial ) &
            &outcome%status = LM_BOUND_REJECTED
        outcome%f = f

    contains

        !> one more consecutive rejection: more damping; at the limit, the terminal status
        subroutine reject
            consecutive_rejections = consecutive_rejections + 1
            call increase_lm_damping(mu, rejection_multiplier)
            if( consecutive_rejections >= LM_MAX_CONSECUTIVE_REJECTIONS ) &
                &outcome%status = rejection_terminal_status(outcome%naccepted, bounded_trial)
        end subroutine reject

    end subroutine lm_solve

    !> The system in the scaled coordinates s_i x_i, masked: an inactive coordinate gets a zero gradient
    !! and a unit diagonal, so it takes a zero step
    pure subroutine lm_scale_system( g, curv, scale, mask, scaled_g, scaled_curv )
        real(dp), intent(in)  :: g(:), curv(:,:), scale(:)
        logical,  intent(in)  :: mask(:)
        real(dp), intent(out) :: scaled_g(:), scaled_curv(:,:)
        integer :: axis, jaxis
        do axis = 1, size(g)
            scaled_g(axis) = scale(axis)*g(axis)
            do jaxis = 1, size(g)
                scaled_curv(axis,jaxis) = scale(axis)*curv(axis,jaxis)*scale(jaxis)
            end do
        end do
        do axis = 1, size(g)
            if( mask(axis) ) cycle
            scaled_g(axis)         = 0._dp
            scaled_curv(axis,:)    = 0._dp
            scaled_curv(:,axis)    = 0._dp
            scaled_curv(axis,axis) = 1._dp
        end do
    end subroutine lm_scale_system

    ! PRIVATE

    !> a finite rejection limit keeps an earlier accepted point
    pure integer function rejection_terminal_status( naccepted, bounded_trial ) result( status )
        integer, intent(in) :: naccepted
        logical, intent(in) :: bounded_trial
        if( naccepted > 0 )then
            status = LM_ACCEPTED
        else if( bounded_trial )then
            status = LM_BOUND_REJECTED
        else
            status = LM_NO_IMPROVEMENT
        endif
    end function rejection_terminal_status

    !> damping grows faster across consecutive rejections; an acceptance resets the multiplier
    pure subroutine increase_lm_damping( mu, rejection_multiplier )
        real(dp), intent(inout) :: mu, rejection_multiplier
        mu = min(mu*rejection_multiplier, LM_MAX_DAMPING)
        rejection_multiplier = min(2._dp*rejection_multiplier, LM_MAX_DAMPING)
    end subroutine increase_lm_damping

end module simple_lm_solver
