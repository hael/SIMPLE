!@descr: unit tests for the optimiser framework and shared rank-1 PCG engine
! The three optimisers production builds through the factory: L-BFGS-B (shift searches, CTF estimation,
! cavg-quality learning) on a bounded quadratic, a direction-finding problem and Rosenbrock, differential
! evolution (CTF estimation) and the restarted simplex (volume symmetry search) on the quadratic, and the
! bookkeeping of the specification object (defaults, limits, evaluation counters, convergence flag).
module simple_opt_tester
use simple_test_utils      ! assertions etc.
use simple_defs            ! sp, dp
use simple_string_utils,   only: int2str
use simple_optimizer,      only: optimizer
use simple_opt_factory,    only: opt_factory
use simple_opt_spec,       only: opt_spec
use simple_pcg_solver,     only: pcg_operator, pcg_solver_options, pcg_solver_outcome, pcg_solve
implicit none
private
public :: run_all_opt_tests

! the 2D quadratic (x-1)^2 + (y+2)^2 and the direction y_norm of the cosine problem
real, parameter :: XMIN(2)   = [1.0, -2.0]
real, parameter :: YDIR(2)   = [0.75, 0.25]
real, parameter :: LIMS2(2,2) = reshape([-5., -5., 5., 5.], [2,2])

type, extends(pcg_operator) :: dense_pcg_operator
    real :: a(3,3) = 0.0
    real :: minv(3) = 1.0
    integer :: apply_count = 0
    logical :: fail_warm_once = .false.
  contains
    procedure :: size    => dense_pcg_size
    procedure :: apply   => dense_pcg_apply
    procedure :: precond => dense_pcg_precond
    procedure :: dot     => dense_pcg_dot
end type dense_pcg_operator

contains

    subroutine run_all_opt_tests()
        write(*,'(A)') '**** running all optimiser tests ****'
        call test_spec_bookkeeping()
        call test_lbfgsb_quadratic()
        call test_lbfgsb_bounded()
        call test_lbfgsb_direction()
        call test_lbfgsb_rosenbrock()
        call test_de_quadratic()
        call test_simplex_quadratic()
        call test_pcg_dense_spd()
        call test_pcg_diagonal_one_step()
        call test_pcg_cold_restart()
    end subroutine run_all_opt_tests

    !---------------- cost functions ----------------

    function quad1d( fun_self, x, d ) result( r )
        class(*), intent(inout) :: fun_self
        integer,  intent(in)    :: d
        real,     intent(in)    :: x(d)
        real                    :: r
        r = (x(1) - 1.)**2
    end function quad1d

    subroutine quad1d_grad( fun_self, x, grad, d )
        class(*), intent(inout) :: fun_self
        integer,  intent(in)    :: d
        real,     intent(inout) :: x(d)
        real,     intent(out)   :: grad(d)
        grad(1) = 2. * (x(1) - 1.)
    end subroutine quad1d_grad

    function quad2d( fun_self, x, d ) result( r )
        class(*), intent(inout) :: fun_self
        integer,  intent(in)    :: d
        real,     intent(in)    :: x(d)
        real                    :: r
        r = sum((x - XMIN)**2)
    end function quad2d

    subroutine quad2d_grad( fun_self, x, grad, d )
        class(*), intent(inout) :: fun_self
        integer,  intent(in)    :: d
        real,     intent(inout) :: x(d)
        real,     intent(out)   :: grad(d)
        grad = 2. * (x - XMIN)
    end subroutine quad2d_grad

    ! 1 - cos(angle between x and YDIR): the smooth form of the old cosine test (acos has an infinite
    ! derivative at the optimum), minimised by any positive multiple of YDIR
    function dircost( fun_self, x, d ) result( r )
        class(*), intent(inout) :: fun_self
        integer,  intent(in)    :: d
        real,     intent(in)    :: x(d)
        real                    :: r
        real :: ynorm(2)
        ynorm = YDIR / sqrt(sum(YDIR**2))
        r = 1. - sum(x * ynorm) / sqrt(sum(x**2))
    end function dircost

    subroutine dircost_grad( fun_self, x, grad, d )
        class(*), intent(inout) :: fun_self
        integer,  intent(in)    :: d
        real,     intent(inout) :: x(d)
        real,     intent(out)   :: grad(d)
        real :: ynorm(2), nrm, c
        ynorm = YDIR / sqrt(sum(YDIR**2))
        nrm   = sqrt(sum(x**2))
        c     = sum(x * ynorm) / nrm
        grad  = -(ynorm - c * x / nrm) / nrm
    end subroutine dircost_grad

    function rosenbrock( fun_self, x, d ) result( r )
        class(*), intent(inout) :: fun_self
        integer,  intent(in)    :: d
        real,     intent(in)    :: x(d)
        real                    :: r
        r = 100. * (x(2) - x(1)**2)**2 + (1. - x(1))**2
    end function rosenbrock

    subroutine rosenbrock_grad( fun_self, x, grad, d )
        class(*), intent(inout) :: fun_self
        integer,  intent(in)    :: d
        real,     intent(inout) :: x(d)
        real,     intent(out)   :: grad(d)
        grad(1) = -400. * x(1) * (x(2) - x(1)**2) - 2. * (1. - x(1))
        grad(2) =  200. * (x(2) - x(1)**2)
    end subroutine rosenbrock_grad

    !---------------- specification ----------------

    subroutine test_spec_bookkeeping()
        type(opt_spec) :: spec
        real :: lims(2,2)
        write(*,'(A)') 'test_spec_bookkeeping'
        call spec%specify('lbfgsb', 2, limits=LIMS2)
        call assert_char('lbfgsb', spec%str_opt, 'specify stores the optimiser name')
        call assert_int(2, spec%ndim, 'specify stores the dimension')
        call assert_int(1, spec%nrestarts, 'one restart by default')
        call assert_int(100, spec%maxits, '100 iterations by default')
        call assert_real(1.e-5, spec%ftol, 1.e-12, 'default cost tolerance')
        call assert_true(allocated(spec%x), 'specify allocates the solution vector')
        if( allocated(spec%x) )then
            call assert_int(2, size(spec%x), 'the solution vector has ndim entries')
            call assert_true(all(spec%x == 0.), 'the solution vector starts at zero')
        endif
        call assert_true(allocated(spec%limits), 'specify stores the limits')
        if( allocated(spec%limits) ) call assert_real(5., spec%limits(2,2), 1.e-6, 'the limits are stored as given')
        call assert_false(associated(spec%costfun), 'no cost function until set')
        call spec%set_costfun(quad2d)
        call spec%set_gcostfun(quad2d_grad)
        call assert_true(associated(spec%costfun),  'set_costfun attaches the cost function')
        call assert_true(associated(spec%gcostfun), 'set_gcostfun attaches the gradient')
        lims(:,1) = -1.
        lims(:,2) =  1.
        call spec%set_limits(lims)
        call assert_real(1., spec%limits(1,2), 1.e-6, 'set_limits replaces the limits')
        call spec%change_opt('simplex')
        call assert_char('simplex', spec%str_opt, 'change_opt renames the optimiser')
        call spec%specify('de', 3, maxits=250, nrestarts=4, ftol=1.e-3, limits=reshape([-1.,-1.,-1.,1.,1.,1.],[3,2]))
        call assert_int(3,   spec%ndim,      're-specify: new dimension')
        call assert_int(250, spec%maxits,    're-specify: maxits as given')
        call assert_int(4,   spec%nrestarts, 're-specify: restarts as given')
        call assert_real(1.e-3, spec%ftol, 1.e-9, 're-specify: ftol as given')
        call assert_int(3, size(spec%x), 're-specify reallocates the solution vector')
        call spec%kill
    end subroutine test_spec_bookkeeping

    !---------------- L-BFGS-B ----------------

    subroutine test_lbfgsb_quadratic()
        class(optimizer), pointer :: opt => null()
        type(opt_factory) :: ofac
        type(opt_spec)    :: spec
        real :: lims(1,2), lowest_cost
        write(*,'(A)') 'test_lbfgsb_quadratic'
        lims(1,:) = [-5., 5.]
        call spec%specify('lbfgsb', 1, limits=lims)
        call spec%set_costfun(quad1d)
        call spec%set_gcostfun(quad1d_grad)
        call ofac%new(spec, opt)
        spec%x = 0.5
        call opt%minimize(spec, opt, lowest_cost)
        call assert_real(1.0, spec%x(1), 1.e-4, 'the quadratic minimum at 1 from x = 0.5')
        call assert_real(0.0, lowest_cost, 1.e-7, 'the minimum cost is zero')
        call assert_true(spec%converged, 'L-BFGS-B reports convergence')
        call assert_true(spec%nevals > 0 .and. spec%nevals < 50, 'a handful of cost evaluations')
        call assert_int(spec%nevals, spec%ngevals, 'one gradient evaluation per cost evaluation (fdf)')
        call opt%kill
        deallocate(opt)
        call spec%kill
    end subroutine test_lbfgsb_quadratic

    ! the bound is active when it excludes the free minimum
    subroutine test_lbfgsb_bounded()
        class(optimizer), pointer :: opt => null()
        type(opt_factory) :: ofac
        type(opt_spec)    :: spec
        real :: lims(1,2), lowest_cost
        write(*,'(A)') 'test_lbfgsb_bounded'
        lims(1,:) = [-5., 0.5]
        call spec%specify('lbfgsb', 1, limits=lims)
        call spec%set_costfun(quad1d)
        call spec%set_gcostfun(quad1d_grad)
        call ofac%new(spec, opt)
        spec%x = -3.
        call opt%minimize(spec, opt, lowest_cost)
        call assert_real(0.5,  spec%x(1),   1.e-5, 'the solution sits on the upper bound')
        call assert_real(0.25, lowest_cost, 1.e-5, 'the cost is the cost at the bound')
        call assert_true(spec%converged, 'a bound-constrained solution still converges')
        call opt%kill
        deallocate(opt)
        call spec%kill
    end subroutine test_lbfgsb_bounded

    ! the old cosine test: from (-5, -7.5), the direction of YDIR
    subroutine test_lbfgsb_direction()
        class(optimizer), pointer :: opt => null()
        type(opt_factory) :: ofac
        type(opt_spec)    :: spec
        real :: lims(2,2), lowest_cost, xnorm(2), ynorm(2)
        write(*,'(A)') 'test_lbfgsb_direction'
        lims(:,1) = -10.
        lims(:,2) =  10.
        call spec%specify('lbfgsb', 2, limits=lims, maxits=500)
        call spec%set_costfun(dircost)
        call spec%set_gcostfun(dircost_grad)
        call ofac%new(spec, opt)
        spec%x = [-5., -7.5]
        call opt%minimize(spec, opt, lowest_cost)
        xnorm = spec%x / sqrt(sum(spec%x**2))
        ynorm = YDIR / sqrt(sum(YDIR**2))
        call assert_real(ynorm(1), xnorm(1), 1.e-3, 'the solution points along YDIR (x)')
        call assert_real(ynorm(2), xnorm(2), 1.e-3, 'the solution points along YDIR (y)')
        call assert_true(lowest_cost < 1.e-5, 'the angle to YDIR vanishes')
        call assert_true(all(spec%x >= lims(:,1)) .and. all(spec%x <= lims(:,2)), 'the solution respects the box')
        call opt%kill
        deallocate(opt)
        call spec%kill
    end subroutine test_lbfgsb_direction

    subroutine test_lbfgsb_rosenbrock()
        class(optimizer), pointer :: opt => null()
        type(opt_factory) :: ofac
        type(opt_spec)    :: spec
        real :: lowest_cost
        write(*,'(A)') 'test_lbfgsb_rosenbrock'
        call spec%specify('lbfgsb', 2, limits=LIMS2, maxits=1000)
        call spec%set_costfun(rosenbrock)
        call spec%set_gcostfun(rosenbrock_grad)
        call ofac%new(spec, opt)
        spec%x = [-1.2, 1.0]
        call opt%minimize(spec, opt, lowest_cost)
        call assert_real(1.0, spec%x(1), 2.e-3, 'Rosenbrock minimum from the classic start (x)')
        call assert_real(1.0, spec%x(2), 4.e-3, 'Rosenbrock minimum from the classic start (y)')
        call assert_true(lowest_cost < 1.e-5, 'Rosenbrock cost at the minimum')
        call assert_true(spec%converged, 'L-BFGS-B converges on Rosenbrock')
        call opt%kill
        deallocate(opt)
        call spec%kill
    end subroutine test_lbfgsb_rosenbrock

    !---------------- differential evolution ----------------

    ! DE is stochastic (ran3): the population is drawn in the box, the best member ends near the minimum
    subroutine test_de_quadratic()
        class(optimizer), pointer :: opt => null()
        type(opt_factory) :: ofac
        type(opt_spec)    :: spec
        real :: lowest_cost
        write(*,'(A)') 'test_de_quadratic'
        call set_fixed_seed(20260927)
        call spec%specify('de', 2, limits=LIMS2, maxits=3000, ftol=1.e-8)
        call spec%set_costfun(quad2d)
        call ofac%new(spec, opt)
        call opt%minimize(spec, opt, lowest_cost)
        call assert_true(lowest_cost < 1.e-2, 'DE brings the quadratic cost below 1e-2')
        call assert_real(XMIN(1), spec%x(1), 0.1, 'DE solution within 0.1 (x)')
        call assert_real(XMIN(2), spec%x(2), 0.1, 'DE solution within 0.1 (y)')
        call assert_true(spec%nevals > spec%npop, 'DE evaluates the population and then the trials')
        call assert_true(spec%nevals <= spec%npop + spec%maxits, 'one trial per generation at most')
        call assert_true(all(spec%x >= LIMS2(:,1)) .and. all(spec%x <= LIMS2(:,2)), 'DE respects the box')
        ! a preset starting point joins the population and is never lost
        spec%x = XMIN
        call opt%minimize(spec, opt, lowest_cost)
        call assert_true(lowest_cost < 1.e-6, 'a preset optimum survives as the best member')
        call opt%kill
        deallocate(opt)
        call spec%kill
    end subroutine test_de_quadratic

    !---------------- simplex ----------------

    subroutine test_simplex_quadratic()
        class(optimizer), pointer :: opt => null()
        type(opt_factory) :: ofac
        type(opt_spec)    :: spec
        real :: lowest_cost
        write(*,'(A)') 'test_simplex_quadratic'
        call set_fixed_seed(20260928)
        call spec%specify('simplex', 2, limits=LIMS2, maxits=1000, nrestarts=3, ftol=1.e-7)
        call spec%set_costfun(quad2d)
        call ofac%new(spec, opt)
        spec%x = [3., -3.]
        call opt%minimize(spec, opt, lowest_cost)
        call assert_real(XMIN(1), spec%x(1), 2.e-3, 'simplex from (3,-3) (x)')
        call assert_real(XMIN(2), spec%x(2), 2.e-3, 'simplex from (3,-3) (y)')
        call assert_true(lowest_cost < 1.e-5, 'simplex cost at the minimum')
        call assert_true(spec%nevals > 3, 'the simplex vertices and the moves are counted')
        call assert_true(spec%niter > 0, 'the average iteration count over the restarts is reported')
        ! a zero start is replaced by a random point in the box
        spec%x = 0.
        call opt%minimize(spec, opt, lowest_cost)
        call assert_real(XMIN(1), spec%x(1), 2.e-3, 'simplex from a random start (x)')
        call assert_real(XMIN(2), spec%x(2), 2.e-3, 'simplex from a random start (y)')
        call opt%kill
        deallocate(opt)
        call spec%kill
    end subroutine test_simplex_quadratic

    !---------------- shared PCG engine ----------------

    subroutine test_pcg_dense_spd()
        type(dense_pcg_operator) :: op
        type(pcg_solver_options) :: options
        type(pcg_solver_outcome) :: outcome
        real :: truth(3), b(3), x(3)
        write(*,'(A)') 'test_pcg_dense_spd'
        call set_spd_operator(op)
        truth = [1.0, -2.0, 0.5]
        b = matmul(op%a, truth)
        x = 0.0
        options%maxits = 3
        options%rtol = 0.0
        options%xtol = 0.0
        call pcg_solve(op, b, x, options, outcome)
        call assert_true(maxval(abs(x-truth)) < 2.e-5, 'PCG matches the independent dense SPD solution')
        call assert_int(3, outcome%iteration_count, 'fixed three-dimensional PCG uses exactly three iterations')
        call assert_char('fixed_iterations', trim(outcome%stop_reason), 'fixed-iteration stop reason')
        call assert_real(1.0, outcome%initial_rel_residual, 1.e-6, 'zero start has unit relative residual')
        call assert_true(outcome%final_rel_residual < 2.e-5, 'dense SPD final residual is small')
        call assert_int(3, size(outcome%rel_residual_history), 'dense SPD residual history covers every iteration')
        call outcome%kill
    end subroutine test_pcg_dense_spd

    subroutine test_pcg_diagonal_one_step()
        type(dense_pcg_operator) :: op
        type(pcg_solver_options) :: options
        type(pcg_solver_outcome) :: outcome
        real :: truth(3), b(3), x(3)
        write(*,'(A)') 'test_pcg_diagonal_one_step'
        op%a = 0.0
        op%a(1,1) = 2.0
        op%a(2,2) = 3.0
        op%a(3,3) = 5.0
        op%minv = [0.5, 1.0/3.0, 0.2]
        truth = [0.25, -1.5, 2.0]
        b = matmul(op%a, truth)
        x = 0.0
        options%maxits = 5
        options%rtol = 1.e-6
        call pcg_solve(op, b, x, options, outcome)
        call assert_int(1, outcome%iteration_count, 'exact diagonal preconditioner converges in one step')
        call assert_char('rtol', trim(outcome%stop_reason), 'one-step diagonal solve stops on residual')
        call assert_true(outcome%converged, 'one-step diagonal solve reports convergence')
        call assert_true(maxval(abs(x-truth)) < 2.e-6, 'one-step diagonal solution matches the oracle')
        call assert_int(1, size(outcome%rel_update_history), 'one-step outcome records one update')
        call outcome%kill
    end subroutine test_pcg_diagonal_one_step

    subroutine test_pcg_cold_restart()
        type(dense_pcg_operator) :: op
        type(pcg_solver_options) :: options
        type(pcg_solver_outcome) :: outcome
        real :: truth(3), b(3), x(3)
        write(*,'(A)') 'test_pcg_cold_restart'
        call set_spd_operator(op)
        op%fail_warm_once = .true.
        truth = [0.5, -0.25, 1.25]
        b = matmul(op%a, truth)
        x = [0.1, 0.1, 0.1]
        options%maxits = 3
        options%rtol = 0.0
        options%cold_restart = .true.
        call pcg_solve(op, b, x, options, outcome)
        call assert_true(outcome%cold_restart_used, 'indefinite warm curvature triggers the one permitted cold restart')
        call assert_int(1, outcome%restart_trigger_iteration, 'cold restart records the triggering iteration')
        call assert_char('fixed_iterations', trim(outcome%stop_reason), 'successful cold retry reaches its fixed budget')
        call assert_true(maxval(abs(x-truth)) < 2.e-5, 'cold retry recovers the dense SPD solution')
        call outcome%kill
    end subroutine test_pcg_cold_restart

    subroutine set_spd_operator( op )
        type(dense_pcg_operator), intent(inout) :: op
        op%a = 0.0
        op%a(1,:) = [4.0, 1.0, 0.0]
        op%a(2,:) = [1.0, 3.0, 1.0]
        op%a(3,:) = [0.0, 1.0, 2.0]
        op%minv = [0.25, 1.0/3.0, 0.5]
        op%apply_count = 0
        op%fail_warm_once = .false.
    end subroutine set_spd_operator

    integer function dense_pcg_size( self ) result(n)
        class(dense_pcg_operator), intent(in) :: self
        n = 3
    end function dense_pcg_size

    subroutine dense_pcg_apply( self, x, y )
        class(dense_pcg_operator), intent(inout) :: self
        real, contiguous, target,  intent(in)    :: x(:)
        real, contiguous, target,  intent(out)   :: y(:)
        self%apply_count = self%apply_count + 1
        if( self%fail_warm_once .and. self%apply_count == 2 )then
            y = -x
        else
            y = matmul(self%a, x)
        endif
    end subroutine dense_pcg_apply

    subroutine dense_pcg_precond( self, r, z )
        class(dense_pcg_operator), intent(inout) :: self
        real, contiguous, target,  intent(in)    :: r(:)
        real, contiguous, target,  intent(out)   :: z(:)
        z = self%minv * r
    end subroutine dense_pcg_precond

    real(dp) function dense_pcg_dot( self, a, b ) result(value)
        class(dense_pcg_operator), intent(in) :: self
        real, contiguous, target,  intent(in) :: a(:), b(:)
        value = sum(real(a,dp) * real(b,dp))
    end function dense_pcg_dot

end module simple_opt_tester
