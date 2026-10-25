!@descr: unit tests of the nonlinear conjugate-gradient engine on analytic objectives
! A quadratic 1/2 x'Ax - b'x (n = 24, A = B'B/n + I, Jacobi preconditioner): without step refinement the
! engine is linear preconditioned CG, so its iterates after six iterations equal pcg_solve's from the same
! start to 1e-4 relative (single-precision vectors, the bound of a few rounding steps amplified by the
! condition number); with the parabola refinement the same (the quadratic's parabola is exact); run to
! convergence it solves Ax = b to a relative residual of 1e-4, stopping on the gradient tolerance; with a
! support projection (a third of the unknowns held at zero) it solves the restricted system to 1e-4 and
! stops on the projected gradient. A nonconvex least squares r_i = x_i^2 + c x_{i+1} - a_i with a known
! zero-residual minimum x* (a saddle at the origin): from a start inside the basin the engine reaches x*
! to 1e-4 relative with a loss that never rises, also when the curvature it is given is ten times too
! small (every step too long), which three halvings absorb; at twenty times they cannot, and the engine
! ends with no_descent without taking the step.
module simple_nlcg_solver_tester
use simple_defs,        only: dp
use simple_nlcg_solver, only: nlcg_objective, nlcg_solver_options, nlcg_solver_outcome, nlcg_solve
use simple_pcg_solver,  only: pcg_operator, pcg_solver_options, pcg_solver_outcome, pcg_solve
use simple_test_utils,  only: assert_true, assert_int, fix_test_seed
implicit none
private
public :: run_all_nlcg_solver_tests

integer, parameter :: NQ = 24          !< size of the quadratic
integer, parameter :: NL = 16          !< size of the nonlinear least squares
real,    parameter :: CPL = 0.3        !< coupling of the nonlinear residuals

!> 1/2 x'Ax - b'x with a Jacobi preconditioner
type, extends(nlcg_objective) :: quad_obj
    real(dp) :: a(NQ,NQ), b(NQ)
    logical  :: supp(NQ) = .true.   !< the support: unknowns outside it are held at zero
  contains
    procedure :: size        => quad_size
    procedure :: eval        => quad_eval
    procedure :: curvature   => quad_curvature
    procedure :: line_losses => quad_line_losses
    procedure :: precond     => quad_precond
    procedure :: dot         => quad_dot
    procedure :: project     => quad_project
end type quad_obj

!> the same A for the linear engine
type, extends(pcg_operator) :: quad_op
    real(dp) :: a(NQ,NQ)
  contains
    procedure :: size    => qop_size
    procedure :: apply   => qop_apply
    procedure :: precond => qop_precond
    procedure :: dot     => qop_dot
end type quad_op

!> 1/2 sum_i r_i(x)^2, r_i = x_i^2 + CPL x_{i+1} - a_i (x_{NL+1} = 0), identity preconditioner
type, extends(nlcg_objective) :: lsq_obj
    real(dp) :: a(NL)
    real(dp) :: curv_scale = 1.0_dp  !< scales the curvature the engine is given
  contains
    procedure :: size        => lsq_size
    procedure :: eval        => lsq_eval
    procedure :: curvature   => lsq_curvature
    procedure :: line_losses => lsq_line_losses
    procedure :: precond     => lsq_precond
    procedure :: dot         => lsq_dot
    procedure :: project     => lsq_project
end type lsq_obj

contains

    subroutine run_all_nlcg_solver_tests()
        write(*,'(A)') '**** running all nonlinear conjugate-gradient tests ****'
        call fix_test_seed()
        call test_quadratic_is_linear_pcg()
        call test_nonconvex_least_squares()
    end subroutine run_all_nlcg_solver_tests

    subroutine test_quadratic_is_linear_pcg()
        integer, parameter :: K = 6
        class(quad_obj), allocatable :: obj
        class(quad_op),  allocatable :: op
        type(nlcg_solver_options) :: nopts
        type(nlcg_solver_outcome) :: nout
        type(pcg_solver_options)  :: lopts
        type(pcg_solver_outcome)  :: lout
        real(dp) :: bmat(NQ,NQ), resid(NQ)
        real     :: x_lin(NQ), x_nl(NQ), x_ref(NQ), rb(NQ), rbm(NQ,NQ)
        integer  :: i
        allocate(obj, op)
        call random_number(rbm)
        bmat = real(rbm - 0.5, dp)
        obj%a = matmul(transpose(bmat), bmat) / real(NQ,dp)
        do i = 1, NQ
            obj%a(i,i) = obj%a(i,i) + 1.0_dp
        end do
        call random_number(rb)
        obj%b = real(rb - 0.5, dp)
        op%a  = obj%a
        ! linear PCG, K iterations from zero
        lopts%maxits = K
        lopts%rtol   = 0.0
        lopts%track_preconditioned_residual = .false.
        x_lin = 0.0
        call pcg_solve(op, real(obj%b), x_lin, lopts, lout)
        call assert_int(K, lout%iteration_count, 'the linear engine ran its fixed iterations')
        ! the nonlinear engine, exact Gauss-Newton steps
        nopts%maxits      = K
        nopts%refine_step = .false.
        x_nl = 0.0
        call nlcg_solve(obj, x_nl, nopts, nout)
        call assert_int(K, nout%iteration_count, 'the nonlinear engine ran its fixed iterations')
        call assert_true(norm2(x_nl - x_lin) <= 1.0e-4 * norm2(x_lin), &
            &'on a quadratic the nonlinear engine is linear preconditioned CG')
        call assert_true(all(nout%loss_history(2:K) <= nout%loss_history(1:K-1)), 'the quadratic loss never rises')
        call assert_int(0, nout%backtracks, 'exact steps on a quadratic need no backtracking')
        call nout%kill
        ! with the parabola refinement, which is exact on a quadratic
        nopts%refine_step = .true.
        x_ref = 0.0
        call nlcg_solve(obj, x_ref, nopts, nout)
        call assert_true(norm2(x_ref - x_lin) <= 1.0e-4 * norm2(x_lin), &
            &'the parabola refinement keeps the exact step on a quadratic')
        call nout%kill
        ! to convergence
        nopts%maxits = 4 * NQ
        nopts%gtol   = 1.0e-6
        x_nl = 0.0
        call nlcg_solve(obj, x_nl, nopts, nout)
        resid = matmul(obj%a, real(x_nl,dp)) - obj%b
        call assert_true(norm2(resid) <= 1.0e-4_dp * norm2(obj%b), 'run to convergence the engine solves Ax = b')
        call assert_true(trim(nout%stop_reason) == 'gtol', 'it stops on the gradient tolerance')
        call nout%kill
        ! on a support: the restricted system A_SS x_S = b_S, x zero outside S
        do i = 1, NQ
            obj%supp(i) = mod(i,3) /= 0
        end do
        x_nl = 0.0
        call nlcg_solve(obj, x_nl, nopts, nout)
        resid = matmul(obj%a, real(x_nl,dp)) - obj%b
        call assert_true(all(merge(x_nl, 0.0, .not. obj%supp) == 0.0), &
            &'the projection holds the unknowns outside the support at zero')
        call assert_true(norm2(merge(resid, 0.0_dp, obj%supp)) <= 1.0e-4_dp * norm2(merge(obj%b, 0.0_dp, obj%supp)), &
            &'on a support the engine solves the restricted system')
        call assert_true(trim(nout%stop_reason) == 'gtol', 'it stops on the projected gradient')
        call nout%kill
        call lout%kill
        deallocate(obj, op)
    end subroutine test_quadratic_is_linear_pcg

    subroutine test_nonconvex_least_squares()
        class(lsq_obj), allocatable :: obj
        type(nlcg_solver_options) :: opts
        type(nlcg_solver_outcome) :: out
        real     :: xstar(NL), x(NL), pert(NL)
        integer  :: i, nit
        allocate(obj)
        do i = 1, NL
            xstar(i) = 1.0 + 0.5 * sin(real(i))
        end do
        do i = 1, NL
            obj%a(i) = real(xstar(i),dp)**2
            if( i < NL ) obj%a(i) = obj%a(i) + real(CPL,dp) * real(xstar(i+1),dp)
        end do
        call random_number(pert)
        x = xstar + 0.6 * (pert - 0.5)
        opts%maxits = 300
        opts%gtol   = 1.0e-6
        call nlcg_solve(obj, x, opts, out)
        nit = out%iteration_count
        call assert_true(norm2(x - xstar) <= 1.0e-4 * norm2(xstar), 'the engine reaches the known minimum')
        call assert_true(trim(out%stop_reason) == 'gtol' .and. nit < opts%maxits, 'it stops on the gradient tolerance')
        call assert_true(out%loss_history(1) <= out%initial_loss, 'the first step lowers the loss')
        call assert_true(all(out%loss_history(2:nit) <= out%loss_history(1:nit-1)), 'the loss never rises')
        call out%kill
        ! every Gauss-Newton step ten times too long: three halvings bring it within twice the minimiser of
        ! the nearly quadratic loss along the line, where the loss no longer rises
        obj%curv_scale = 0.1_dp
        opts%refine_step = .false.
        x = xstar + 0.6 * (pert - 0.5)
        call nlcg_solve(obj, x, opts, out)
        nit = out%iteration_count
        call assert_true(out%backtracks > 0, 'steps that raise the loss are halved')
        call assert_true(out%loss_history(1) <= out%initial_loss .and. &
            &all(out%loss_history(2:nit) <= out%loss_history(1:nit-1)), 'with overlong steps the loss never rises')
        call assert_true(norm2(x - xstar) <= 1.0e-4 * norm2(xstar), 'with overlong steps the engine reaches the minimum')
        call out%kill
        ! twenty times too long is beyond three halvings: the steepest descent cannot descend, the start stays
        obj%curv_scale = 0.05_dp
        x = xstar + 0.6 * (pert - 0.5)
        call nlcg_solve(obj, x, opts, out)
        call assert_true(trim(out%stop_reason) == 'no_descent' .and. out%iteration_count == 0, &
            &'a step no halving brings downhill ends the solve')
        call assert_true(all(x == xstar + 0.6 * (pert - 0.5)), 'a step that raises the loss is not taken')
        call out%kill
        deallocate(obj)
    end subroutine test_nonconvex_least_squares

    ! ---- the quadratic ----

    integer function quad_size( self ) result( n )
        class(quad_obj), intent(in) :: self
        n = NQ
    end function quad_size

    subroutine quad_eval( self, x, loss, g )
        class(quad_obj),          intent(inout) :: self
        real, contiguous, target, intent(in)    :: x(:)
        real(dp),                 intent(out)   :: loss
        real, contiguous, target, intent(out)   :: g(:)
        real(dp) :: ax(NQ)
        ax   = matmul(self%a, real(x,dp))
        loss = 0.5_dp * dot_product(real(x,dp), ax) - dot_product(self%b, real(x,dp))
        g    = real(ax - self%b)
    end subroutine quad_eval

    real(dp) function quad_curvature( self, x, d ) result( curv )
        class(quad_obj),          intent(inout) :: self
        real, contiguous, target, intent(in)    :: x(:), d(:)
        curv = dot_product(real(d,dp), matmul(self%a, real(d,dp)))
    end function quad_curvature

    subroutine quad_line_losses( self, x, d, t, losses )
        class(quad_obj),          intent(inout) :: self
        real, contiguous, target, intent(in)    :: x(:), d(:)
        real(dp),                 intent(in)    :: t(:)
        real(dp),                 intent(out)   :: losses(:)
        real(dp) :: xt(NQ)
        integer  :: j
        do j = 1, size(t)
            xt = real(x,dp) + t(j) * real(d,dp)
            losses(j) = 0.5_dp * dot_product(xt, matmul(self%a, xt)) - dot_product(self%b, xt)
        end do
    end subroutine quad_line_losses

    subroutine quad_precond( self, r, z )
        class(quad_obj),          intent(inout) :: self
        real, contiguous, target, intent(in)    :: r(:)
        real, contiguous, target, intent(out)   :: z(:)
        integer :: i
        do i = 1, NQ
            z(i) = r(i) / real(self%a(i,i))
        end do
    end subroutine quad_precond

    real(dp) function quad_dot( self, a, b ) result( value )
        class(quad_obj),          intent(in) :: self
        real, contiguous, target, intent(in) :: a(:), b(:)
        value = sum(real(a,dp) * real(b,dp))
    end function quad_dot

    subroutine quad_project( self, v )
        class(quad_obj),          intent(inout) :: self
        real, contiguous, target, intent(inout) :: v(:)
        v = merge(v, 0.0, self%supp)
    end subroutine quad_project

    integer function qop_size( self ) result( n )
        class(quad_op), intent(in) :: self
        n = NQ
    end function qop_size

    subroutine qop_apply( self, x, y )
        class(quad_op),           intent(inout) :: self
        real, contiguous, target, intent(in)    :: x(:)
        real, contiguous, target, intent(out)   :: y(:)
        y = real(matmul(self%a, real(x,dp)))
    end subroutine qop_apply

    subroutine qop_precond( self, r, z )
        class(quad_op),           intent(inout) :: self
        real, contiguous, target, intent(in)    :: r(:)
        real, contiguous, target, intent(out)   :: z(:)
        integer :: i
        do i = 1, NQ
            z(i) = r(i) / real(self%a(i,i))
        end do
    end subroutine qop_precond

    real(dp) function qop_dot( self, a, b ) result( value )
        class(quad_op),           intent(in) :: self
        real, contiguous, target, intent(in) :: a(:), b(:)
        value = sum(real(a,dp) * real(b,dp))
    end function qop_dot

    ! ---- the nonlinear least squares ----

    integer function lsq_size( self ) result( n )
        class(lsq_obj), intent(in) :: self
        n = NL
    end function lsq_size

    pure function lsq_resid( a, x ) result( r )
        real(dp), intent(in) :: a(NL), x(NL)
        real(dp) :: r(NL)
        integer  :: i
        do i = 1, NL
            r(i) = x(i)**2 - a(i)
            if( i < NL ) r(i) = r(i) + real(CPL,dp) * x(i+1)
        end do
    end function lsq_resid

    subroutine lsq_eval( self, x, loss, g )
        class(lsq_obj),           intent(inout) :: self
        real, contiguous, target, intent(in)    :: x(:)
        real(dp),                 intent(out)   :: loss
        real, contiguous, target, intent(out)   :: g(:)
        real(dp) :: r(NL), gd(NL)
        integer  :: j
        r    = lsq_resid(self%a, real(x,dp))
        loss = 0.5_dp * sum(r**2)
        ! g = J^T r: dr_j/dx_j = 2 x_j, dr_{j-1}/dx_j = CPL
        do j = 1, NL
            gd(j) = 2.0_dp * real(x(j),dp) * r(j)
            if( j > 1 ) gd(j) = gd(j) + real(CPL,dp) * r(j-1)
        end do
        g = real(gd)
    end subroutine lsq_eval

    real(dp) function lsq_curvature( self, x, d ) result( curv )
        class(lsq_obj),           intent(inout) :: self
        real, contiguous, target, intent(in)    :: x(:), d(:)
        real(dp) :: jd
        integer  :: i
        curv = 0.0_dp
        do i = 1, NL
            jd = 2.0_dp * real(x(i),dp) * real(d(i),dp)
            if( i < NL ) jd = jd + real(CPL,dp) * real(d(i+1),dp)
            curv = curv + jd**2
        end do
        curv = self%curv_scale * curv
    end function lsq_curvature

    subroutine lsq_line_losses( self, x, d, t, losses )
        class(lsq_obj),           intent(inout) :: self
        real, contiguous, target, intent(in)    :: x(:), d(:)
        real(dp),                 intent(in)    :: t(:)
        real(dp),                 intent(out)   :: losses(:)
        integer :: j
        do j = 1, size(t)
            losses(j) = 0.5_dp * sum(lsq_resid(self%a, real(x,dp) + t(j) * real(d,dp))**2)
        end do
    end subroutine lsq_line_losses

    subroutine lsq_precond( self, r, z )
        class(lsq_obj),           intent(inout) :: self
        real, contiguous, target, intent(in)    :: r(:)
        real, contiguous, target, intent(out)   :: z(:)
        z = r
    end subroutine lsq_precond

    real(dp) function lsq_dot( self, a, b ) result( value )
        class(lsq_obj),           intent(in) :: self
        real, contiguous, target, intent(in) :: a(:), b(:)
        value = sum(real(a,dp) * real(b,dp))
    end function lsq_dot

    subroutine lsq_project( self, v )
        class(lsq_obj),           intent(inout) :: self
        real, contiguous, target, intent(inout) :: v(:)
    end subroutine lsq_project

end module simple_nlcg_solver_tester
