!@descr: preconditioned nonlinear conjugate-gradient engine (Polak-Ribiere) with a Gauss-Newton step length
! on a client-supplied loss, gradient, curvature, preconditioner, inner product and constraint projection
module simple_nlcg_solver
use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
use simple_defs,  only: dp, logfhandle
use simple_error, only: simple_exception
implicit none
private
#include "simple_local_flags.inc"

public :: nlcg_objective, nlcg_solver_options, nlcg_solver_outcome, nlcg_solve

!> The client of nlcg_solve: its objective (data misfit plus prior) through contiguous rank-1 views, with the
!! deferred size, eval, curvature, line_losses, precond, dot and project (the contract: wave_model_policy.md,
!! section 6.1); the engine owns only its search vectors
type, abstract :: nlcg_objective
  contains
    procedure(nlcg_size_i),        deferred :: size
    procedure(nlcg_eval_i),        deferred :: eval
    procedure(nlcg_curvature_i),   deferred :: curvature
    procedure(nlcg_line_losses_i), deferred :: line_losses
    procedure(nlcg_precond_i),     deferred :: precond
    procedure(nlcg_dot_i),         deferred :: dot
    procedure(nlcg_project_i),     deferred :: project
end type nlcg_objective

type :: nlcg_solver_options
    integer           :: maxits         = 4       !< iterations (each one full gradient pass)
    real              :: gtol           = 0.0     !< stop when |Pg|/|Pg_0| falls below; 0 runs maxits
    logical           :: refine_step    = .true.  !< refine the Gauss-Newton step by a parabola on the subset
    real              :: tmin_frac      = 0.1     !< the refined step stays within [tmin_frac, tmax_frac] t0
    real              :: tmax_frac      = 3.0
    integer           :: max_backtracks = 3       !< halvings of a step that raised the loss
    integer           :: verbose        = 0
    character(len=96) :: tag            = 'NLCG'
end type nlcg_solver_options

type :: nlcg_solver_outcome
    character(len=24)     :: stop_reason     = 'not_started'
    integer               :: iteration_count = 0
    integer               :: restarts        = 0  !< steepest-descent restarts (direction not downhill, or no descent along it)
    integer               :: backtracks      = 0
    real(dp)              :: initial_loss    = 0.0_dp
    real(dp)              :: final_loss      = 0.0_dp
    real(dp), allocatable :: loss_history(:)      !< L after each iteration
    real(dp), allocatable :: step_history(:)      !< the step taken
    real(dp), allocatable :: gn_step_history(:)   !< the Gauss-Newton step t0
    real(dp), allocatable :: grad_history(:)      !< |Pg|/|Pg_0| after each iteration
  contains
    procedure :: kill => nlcg_solver_outcome_kill
end type nlcg_solver_outcome

abstract interface
    integer function nlcg_size_i( self ) result( n )
        import :: nlcg_objective
        class(nlcg_objective), intent(in) :: self
    end function nlcg_size_i

    subroutine nlcg_eval_i( self, x, loss, g )
        import :: dp, nlcg_objective
        class(nlcg_objective),    intent(inout) :: self
        real, contiguous, target, intent(in)    :: x(:)
        real(dp),                 intent(out)   :: loss
        real, contiguous, target, intent(out)   :: g(:)
    end subroutine nlcg_eval_i

    real(dp) function nlcg_curvature_i( self, x, d ) result( curv )
        import :: dp, nlcg_objective
        class(nlcg_objective),    intent(inout) :: self
        real, contiguous, target, intent(in)    :: x(:), d(:)
    end function nlcg_curvature_i

    subroutine nlcg_line_losses_i( self, x, d, t, losses )
        import :: dp, nlcg_objective
        class(nlcg_objective),    intent(inout) :: self
        real, contiguous, target, intent(in)    :: x(:), d(:)
        real(dp),                 intent(in)    :: t(:)
        real(dp),                 intent(out)   :: losses(:)
    end subroutine nlcg_line_losses_i

    subroutine nlcg_precond_i( self, r, z )
        import :: nlcg_objective
        class(nlcg_objective),    intent(inout) :: self
        real, contiguous, target, intent(in)    :: r(:)
        real, contiguous, target, intent(out)   :: z(:)
    end subroutine nlcg_precond_i

    real(dp) function nlcg_dot_i( self, a, b ) result( value )
        import :: dp, nlcg_objective
        class(nlcg_objective),    intent(in) :: self
        real, contiguous, target, intent(in) :: a(:), b(:)
    end function nlcg_dot_i

    subroutine nlcg_project_i( self, v )
        import :: nlcg_objective
        class(nlcg_objective),    intent(inout) :: self
        real, contiguous, target, intent(inout) :: v(:)
    end subroutine nlcg_project_i
end interface

contains

    !> Minimise L from x by preconditioned Polak-Ribiere conjugate gradients with the Gauss-Newton step
    !! -<g,d> / curvature(d), optionally refined by a parabola, halving steps that raise the loss and
    !! restarting from steepest descent (the algorithm: wave_model_policy.md, section 6.1)
    subroutine nlcg_solve( obj, x, options, outcome )
        class(nlcg_objective),     intent(inout) :: obj
        real, contiguous,          intent(inout) :: x(:)
        type(nlcg_solver_options), intent(in)    :: options
        type(nlcg_solver_outcome), intent(out)   :: outcome
        real,     allocatable :: g(:), gnew(:), r(:), rnew(:), z(:), znew(:), d(:), x_old(:)
        real(dp), allocatable :: lhist(:), thist(:), t0hist(:), ghist(:)
        real(dp) :: loss, loss_new, slope, curv, t0, t, rz, rznew, rold_znew, beta, g0norm
        real(dp) :: phis(3), den
        integer  :: n, iter, nb
        logical  :: l_sd
        n = obj%size()
        if( n < 1 .or. size(x) /= n ) THROW_HARD('NLCG vector size differs from the objective contract')
        if( options%maxits < 1 ) THROW_HARD('NLCG maxits must be at least 1')
        if( .not. ieee_is_finite(options%gtol) ) THROW_HARD('NLCG gradient tolerance must be finite')
        if( options%tmin_frac <= 0.0 .or. options%tmax_frac < options%tmin_frac ) &
            &THROW_HARD('NLCG step bounds must satisfy 0 < tmin_frac <= tmax_frac')
        allocate(g(n), gnew(n), r(n), rnew(n), z(n), znew(n), d(n), x_old(n))
        allocate(lhist(options%maxits), thist(options%maxits), t0hist(options%maxits), &
            &ghist(options%maxits), source=0.0_dp)
        outcome%stop_reason = 'maxits'
        call obj%project(x)
        call obj%eval(x, loss, g)
        if( .not. ieee_is_finite(loss) ) THROW_HARD('NLCG start has a non-finite loss')
        outcome%initial_loss = loss
        g0norm = projected_norm(g)
        if( g0norm <= 0.0_dp )then
            outcome%stop_reason = 'zero_gradient'
            outcome%final_loss  = loss
            call finish(0)
            return
        endif
        r = -g
        call obj%precond(r, z)
        call obj%project(z)
        rz = obj%dot(r, z)
        if( rz <= 0.0_dp ) THROW_HARD('non-positive dot(r,z); NLCG preconditioner is not positive definite')
        d    = z
        l_sd = .true.
        do iter = 1, options%maxits
            slope = obj%dot(g, d)
            if( slope >= 0.0_dp )then
                ! not a descent direction: restart from the preconditioned steepest descent
                d     = z
                slope = -rz
                l_sd  = .true.
                outcome%restarts = outcome%restarts + 1
            endif
            curv = obj%curvature(x, d)
            if( .not. ieee_is_finite(curv) .or. curv <= 0.0_dp )then
                outcome%stop_reason = 'indefinite'
                call finish(iter - 1)
                return
            endif
            t0 = -slope / curv
            t  = t0
            if( options%refine_step )then
                call obj%line_losses(x, d, [0.0_dp, t0, 2.0_dp*t0], phis)
                ! the parabola through (0,phi1), (t0,phi2), (2t0,phi3): minimum at
                ! t0 (3 phi1 - 4 phi2 + phi3) / (2 (phi1 - 2 phi2 + phi3)) when convex
                den = phis(1) - 2.0_dp*phis(2) + phis(3)
                if( ieee_is_finite(den) .and. den > 0.0_dp )then
                    t = t0 * (3.0_dp*phis(1) - 4.0_dp*phis(2) + phis(3)) / (2.0_dp*den)
                    t = min(max(t, real(options%tmin_frac,dp)*t0), real(options%tmax_frac,dp)*t0)
                endif
            endif
            x_old = x
            call try_step(t)
            nb = 0
            do while( rises() .and. nb < options%max_backtracks )
                nb = nb + 1
                t  = 0.5_dp * t
                call try_step(t)
            end do
            outcome%backtracks = outcome%backtracks + nb
            t0hist(iter) = t0
            if( rises() )then
                ! no descent along d: the iterate stays
                x = x_old
                if( l_sd )then
                    outcome%stop_reason = 'no_descent'
                    call finish(iter - 1)
                    return
                endif
                d    = z
                l_sd = .true.
                outcome%restarts = outcome%restarts + 1
                lhist(iter) = loss
                thist(iter) = 0.0_dp
                ghist(iter) = projected_norm(g) / g0norm
                cycle
            endif
            loss = loss_new
            g    = gnew
            rnew = -gnew
            call obj%precond(rnew, znew)
            call obj%project(znew)
            rznew     = obj%dot(rnew, znew)
            rold_znew = obj%dot(r, znew)
            ! preconditioned Polak-Ribiere, clipped: <r_new - r_old, z_new> / <r_old, z_old>
            beta = max(0.0_dp, (rznew - rold_znew) / rz)
            d    = znew + real(beta) * d
            l_sd = beta == 0.0_dp
            r    = rnew
            z    = znew
            rz   = rznew
            lhist(iter) = loss
            thist(iter) = t
            ghist(iter) = projected_norm(g) / g0norm
            if( options%verbose > 0 )then
                write(logfhandle,'(A,A,A,I4,A,ES14.6,A,ES10.3,A,ES10.3,A,ES10.3)') '>>> ', trim(options%tag), ' it', iter, &
                    &'  loss=', loss, '  |Pg|/|Pg0|=', ghist(iter), '  t/t0=', t / t0, '  beta=', beta
            endif
            if( options%gtol > 0.0 .and. ghist(iter) <= real(options%gtol,dp) )then
                outcome%stop_reason = 'gtol'
                call finish(iter)
                return
            endif
            if( rz <= 0.0_dp )then
                outcome%stop_reason = 'converged'
                call finish(iter)
                return
            endif
        end do
        call finish(options%maxits)

    contains

        !> x = P(x_old + t d) and its loss and gradient
        subroutine try_step( tstep )
            real(dp), intent(in) :: tstep
            x = x_old + real(tstep) * d
            call obj%project(x)
            call obj%eval(x, loss_new, gnew)
        end subroutine try_step

        logical function rises()
            rises = .not. ieee_is_finite(loss_new)
            if( .not. rises ) rises = loss_new > loss
        end function rises

        real(dp) function projected_norm( v )
            real, intent(in) :: v(:)
            real, allocatable :: pv(:)
            pv = v
            call obj%project(pv)
            projected_norm = sqrt(obj%dot(pv, pv))
        end function projected_norm

        subroutine finish( ndone )
            integer, intent(in) :: ndone
            outcome%iteration_count = ndone
            outcome%final_loss      = loss
            allocate(outcome%loss_history(ndone),    source=lhist(1:ndone))
            allocate(outcome%step_history(ndone),    source=thist(1:ndone))
            allocate(outcome%gn_step_history(ndone), source=t0hist(1:ndone))
            allocate(outcome%grad_history(ndone),    source=ghist(1:ndone))
            deallocate(g, gnew, r, rnew, z, znew, d, x_old, lhist, thist, t0hist, ghist)
        end subroutine finish

    end subroutine nlcg_solve

    subroutine nlcg_solver_outcome_kill( self )
        class(nlcg_solver_outcome), intent(inout) :: self
        if( allocated(self%loss_history)    ) deallocate(self%loss_history)
        if( allocated(self%step_history)    ) deallocate(self%step_history)
        if( allocated(self%gn_step_history) ) deallocate(self%gn_step_history)
        if( allocated(self%grad_history)    ) deallocate(self%grad_history)
        self%stop_reason     = 'not_started'
        self%iteration_count = 0
        self%restarts        = 0
        self%backtracks      = 0
        self%initial_loss    = 0.0_dp
        self%final_loss      = 0.0_dp
    end subroutine nlcg_solver_outcome_kill

end module simple_nlcg_solver
