!@descr: unit test routines for the optimiser framework (opt_spec, opt_factory, L-BFGS-B, differential evolution, simplex)
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
implicit none
private
public :: run_all_opt_tests

! the 2D quadratic (x-1)^2 + (y+2)^2 and the direction y_norm of the cosine problem
real, parameter :: XMIN(2)   = [1.0, -2.0]
real, parameter :: YDIR(2)   = [0.75, 0.25]
real, parameter :: LIMS2(2,2) = reshape([-5., -5., 5., 5.], [2,2])

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
        call spec%specify('de', 2, limits=LIMS2, maxits=3000, ftol=1.e-8)
        call spec%set_costfun(quad2d)
        call ofac%new(spec, opt)
        call opt%minimize(spec, opt, lowest_cost)
        call assert_true(lowest_cost < 1.e-2, 'DE brings the quadratic cost below 1e-2 (population tolerance stop)')
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

end module simple_opt_tester
