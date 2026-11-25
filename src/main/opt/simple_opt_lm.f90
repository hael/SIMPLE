!@descr: Levenberg-Marquardt optimiser of the factory ('lm'): lm_solve on an opt_spec's fdfcurvfun_8 within its limits
! Double precision. Starts from spec%x (projected into the box, as L-BFGS-B) and returns the solution in x and
! x_8, the counts in nevals, ngevals and niter. Bounds: a coordinate on a bound whose gradient points outward is
! not free, and the free step is clipped to the box (lm_optimizer.md, decision 4). opt_callback runs after
! each accepted step and the objective is then evaluated again unless the callback cleared spec%callback_changed
! (decision 5, ruling R10). The stops take L-BFGS-B's meaning of factr and pgtol (ruling R3), and converged its
! meaning too: a stopping test ended the solve, not maxits or the rejection limit (ruling 2026-11-20).
module simple_opt_lm
use simple_core_module_api
use simple_optimizer, only: optimizer
use simple_opt_spec,  only: opt_spec
use simple_lm_solver, only: lm_problem, lm_options, lm_outcome, lm_solve
implicit none

public :: opt_lm
private
#include "simple_local_flags.inc"

type, extends(optimizer) :: opt_lm
    private
    logical :: exists = .false.
contains
    procedure :: new      => new_opt_lm
    procedure :: minimize => lm_minimize
    procedure :: kill     => kill_opt_lm
end type opt_lm

!> the spec's cost over its box; spec and fun_self are valid during minimize only
type, extends(lm_problem) :: spec_problem
    class(opt_spec), pointer     :: spec     => null()
    class(*),        pointer     :: fun_self => null()
    real(dp),        allocatable :: x(:), trial(:), lower(:), upper(:)
    logical                      :: bounded  = .false.
contains
    procedure :: eval       => spec_problem_eval
    procedure :: limit_step => spec_problem_limit_step
    procedure :: propose    => spec_problem_propose
    procedure :: accept     => spec_problem_accept
    procedure :: free       => spec_problem_free
end type spec_problem

contains

    !> 'lm' needs the cost, gradient and curvature callback (set_fdfcurvfun_8)
    subroutine new_opt_lm( self, spec )
        class(opt_lm),   intent(inout) :: self
        class(opt_spec), intent(inout) :: spec
        call self%kill
        if( .not. associated(spec%fdfcurvfun_8) ) THROW_HARD('the lm optimiser requires set_fdfcurvfun_8; new_opt_lm')
        self%exists = .true.
    end subroutine new_opt_lm

    subroutine lm_minimize( self, spec, fun_self, lowest_cost )
        class(opt_lm),   intent(inout) :: self
        class(opt_spec), intent(inout) :: spec
        class(*),        intent(inout) :: fun_self
        real,            intent(out)   :: lowest_cost
        if( .not. self%exists ) THROW_HARD('lm optimiser used before new; lm_minimize')
        call minimize_targets(spec, fun_self, lowest_cost)
    end subroutine lm_minimize

    !> the solve, with spec and fun_self as targets of the problem's pointers
    subroutine minimize_targets( spec, fun_self, lowest_cost )
        class(opt_spec),         target, intent(inout) :: spec
        class(*),                target, intent(inout) :: fun_self
        real,                    intent(out)           :: lowest_cost
        type(spec_problem) :: problem
        type(lm_options)   :: options
        type(lm_outcome)   :: outcome
        integer            :: n
        if( .not. associated(spec%fdfcurvfun_8) ) THROW_HARD('the lm optimiser requires set_fdfcurvfun_8; lm_minimize')
        n = spec%ndim
        problem%spec     => spec
        problem%fun_self => fun_self
        allocate(problem%x(n), source=real(spec%x, dp))
        problem%bounded = allocated(spec%limits)
        if( problem%bounded )then
            problem%lower = real(spec%limits(:,1), dp)
            problem%upper = real(spec%limits(:,2), dp)
            problem%x     = min(max(problem%x, problem%lower), problem%upper)
        endif
        problem%trial = problem%x
        spec%nevals   = 0
        spec%ngevals  = 0
        options%maxits         = spec%maxits
        ! both relative stops at factr machine epsilons; pgtol on the free (projected) gradient
        options%predicted_tol  = spec%factr*epsilon(1._dp)
        options%accepted_tol   = spec%factr*epsilon(1._dp)
        options%stationary_tol = spec%pgtol
        allocate(options%scale(n), source=1._dp)
        allocate(options%active(n), source=.true.)
        call lm_solve(problem, options, outcome)
        spec%x_8      = problem%x
        spec%x        = real(problem%x)
        spec%niter    = outcome%niter
        ! a stopping test ended the solve; the status (an accepted point at maxits included) is not convergence
        spec%converged = outcome%converged
        lowest_cost   = real(outcome%f)
        problem%spec     => null()
        problem%fun_self => null()
    end subroutine minimize_targets

    subroutine kill_opt_lm( self )
        class(opt_lm), intent(inout) :: self
        self%exists = .false.
    end subroutine kill_opt_lm

    ! THE PROBLEM

    subroutine spec_problem_eval( self, f, g, curv )
        class(spec_problem), intent(inout) :: self
        real(dp),            intent(out)   :: f, g(:), curv(:,:)
        call self%spec%eval_fdfcurv_8(self%fun_self, self%trial, f, g, curv)
    end subroutine spec_problem_eval

    !> a coordinate on a bound whose gradient points outward is frozen this iteration
    subroutine spec_problem_free( self, g, free )
        class(spec_problem), intent(inout) :: self
        real(dp),            intent(in)    :: g(:)
        logical,             intent(out)   :: free(:)
        free = .true.
        if( .not. self%bounded ) return
        free = .not. ((self%x <= self%lower .and. g > 0._dp) .or. (self%x >= self%upper .and. g < 0._dp))
    end subroutine spec_problem_free

    !> clip the step to the box
    subroutine spec_problem_limit_step( self, step, bounded )
        class(spec_problem), intent(inout) :: self
        real(dp),            intent(inout) :: step(:)
        logical,             intent(out)   :: bounded
        integer :: i
        bounded = .false.
        if( .not. self%bounded ) return
        do i = 1, size(step)
            if( self%x(i) + step(i) < self%lower(i) )then
                step(i) = self%lower(i) - self%x(i)
                bounded = .true.
            else if( self%x(i) + step(i) > self%upper(i) )then
                step(i) = self%upper(i) - self%x(i)
                bounded = .true.
            endif
        end do
    end subroutine spec_problem_limit_step

    !> the trial point, kept in the box against rounding
    function spec_problem_propose( self, step ) result( feasible )
        class(spec_problem), intent(inout) :: self
        real(dp),            intent(in)    :: step(:)
        logical :: feasible
        self%trial = self%x + step
        if( self%bounded ) self%trial = min(max(self%trial, self%lower), self%upper)
        feasible = .true.
    end function spec_problem_propose

    !> adopt the trial point; a callback may change the objective, which is then evaluated again unless the callback
    !! cleared the spec's callback_changed
    function spec_problem_accept( self, ratio ) result( changed )
        class(spec_problem), intent(inout) :: self
        real(dp),            intent(in)    :: ratio
        logical :: changed
        self%x = self%trial
        changed = associated(self%spec%opt_callback)
        if( .not. changed ) return
        self%spec%x_8 = self%x
        self%spec%x   = real(self%x)
        self%spec%callback_changed = .true.
        call self%spec%opt_callback(self%fun_self)
        changed = self%spec%callback_changed
    end function spec_problem_accept

end module simple_opt_lm
