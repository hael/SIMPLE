!@descr: unit tests for the optimiser framework and shared rank-1 PCG engine
! The optimisers production builds through the factory: L-BFGS-B (shift searches, CTF estimation,
! cavg-quality learning) on a bounded quadratic, a direction-finding problem and Rosenbrock, differential
! evolution (CTF estimation) and the restarted simplex (volume symmetry search) on the quadratic, the
! Levenberg-Marquardt 'lm' (L1-L6 of lm_optimizer.md), and the bookkeeping of the specification object.
module simple_opt_tester
use simple_test_utils      ! assertions etc.
use simple_defs            ! sp, dp
use simple_string_utils,   only: int2str
use simple_optimizer,      only: optimizer
use simple_opt_factory,    only: opt_factory
use simple_opt_spec,       only: opt_spec
use simple_pcg_solver,     only: pcg_operator, pcg_solver_options, pcg_solver_outcome, pcg_solve
use simple_lm_solver,      only: lm_problem, lm_options, lm_outcome, lm_solve, LM_ACCEPTED
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

! L4: the objective (x1 - c)^2/2 + (x2 + 2)^2/2, whose centre c the callback moves once (a discrete
! choice, as the polar search's angle), and the record of the evaluations and callbacks
type :: lm_recorder
    real(dp) :: c = 1._dp
    real(dp) :: cb_x(2) = huge(1._dp)
    integer  :: nev = 0, ncb = 0, nreeval = 0, nrepeat = 0
    logical  :: pending = .false.
    logical  :: clear_unchanged = .false. !< L4b: after its first call the callback reports the objective unchanged
    type(opt_spec), pointer :: spec => null() !< the callback reads the accepted point here
end type lm_recorder

! the engine on a tester cost with an optional box clip and no free hook: L2's plain-clipping control
! (the L2 cost, clipped) and L6's status (the null-coordinate cost)
type, extends(lm_problem) :: clip_only_problem
    real(dp) :: x(2) = 0._dp, trial(2) = 0._dp, lower(2) = -huge(1._dp), upper(2) = huge(1._dp)
    logical  :: null_cost = .false.
  contains
    procedure :: eval       => clip_only_eval
    procedure :: limit_step => clip_only_limit_step
    procedure :: propose    => clip_only_propose
    procedure :: accept     => clip_only_accept
end type clip_only_problem

! the first point a minimiser evaluates once armed: simplex's start, DE's first member
type :: eval_recorder
    real    :: first(2) = huge(1.)
    logical :: armed    = .false.
end type eval_recorder

contains

    subroutine run_all_opt_tests()
        write(*,'(A)') '**** running all optimiser tests ****'
        call test_spec_bookkeeping()
        call test_lbfgsb_quadratic()
        call test_lbfgsb_bounded()
        call test_lbfgsb_direction()
        call test_lbfgsb_rosenbrock()
        call test_lm_rosenbrock()
        call test_lm_coupled_bound()
        call test_lm_flipped_gradient()
        call test_lm_callback()
        call test_lm_callback_unchanged()
        call test_lm_rank_deficient()
        call test_lm_null_coordinate()
        call test_de_quadratic()
        call test_simplex_quadratic()
        call test_start_at_origin()
        call test_de_npop()
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

    ! Rosenbrock in residual form r = (10(y - x^2), 1 - x), f = |r|^2/2: J^T r and J^T J
    subroutine rosen_res_fdfcurv( fun_self, x, f, grad, curv, d )
        class(*),     intent(inout) :: fun_self
        integer,      intent(in)    :: d
        real(kind=8), intent(inout) :: x(d)
        real(kind=8), intent(out)   :: f, grad(d), curv(d,d)
        real(dp) :: r(2), jac(2,2)
        r      = [10._dp*(x(2) - x(1)**2), 1._dp - x(1)]
        jac(1,:) = [-20._dp*x(1), 10._dp]
        jac(2,:) = [-1._dp, 0._dp]
        call lsq_terms(r, jac, f, grad, curv)
    end subroutine rosen_res_fdfcurv

    ! the same with the gradient's sign flipped (L3)
    subroutine rosen_res_flipped( fun_self, x, f, grad, curv, d )
        class(*),     intent(inout) :: fun_self
        integer,      intent(in)    :: d
        real(kind=8), intent(inout) :: x(d)
        real(kind=8), intent(out)   :: f, grad(d), curv(d,d)
        call rosen_res_fdfcurv(fun_self, x, f, grad, curv, d)
        grad = -grad
    end subroutine rosen_res_flipped

    function rosen_res_cost( fun_self, x, d ) result( f )
        class(*),     intent(inout) :: fun_self
        integer,      intent(in)    :: d
        real(kind=8), intent(in)    :: x(d)
        real(kind=8)                :: f
        f = 0.5_dp*((10._dp*(x(2) - x(1)**2))**2 + (1._dp - x(1))**2)
    end function rosen_res_cost

    subroutine rosen_res_grad( fun_self, x, grad, d )
        class(*),     intent(inout) :: fun_self
        integer,      intent(in)    :: d
        real(kind=8), intent(inout) :: x(d)
        real(kind=8), intent(out)   :: grad(d)
        real(dp) :: f, curv(2,2)
        call rosen_res_fdfcurv(fun_self, x, f, grad, curv, d)
    end subroutine rosen_res_grad

    ! L2: ((x + 2y + 1)^2 + (y - 1.9)^2)/2
    subroutine coupled_fdfcurv( fun_self, x, f, grad, curv, d )
        class(*),     intent(inout) :: fun_self
        integer,      intent(in)    :: d
        real(kind=8), intent(inout) :: x(d)
        real(kind=8), intent(out)   :: f, grad(d), curv(d,d)
        real(dp) :: jac(2,2)
        jac(1,:) = [1._dp, 2._dp]
        jac(2,:) = [0._dp, 1._dp]
        call lsq_terms([x(1) + 2._dp*x(2) + 1._dp, x(2) - 1.9_dp], jac, f, grad, curv)
    end subroutine coupled_fdfcurv

    ! L5: r = x + y - 1, a singular curvature
    subroutine rank1_fdfcurv( fun_self, x, f, grad, curv, d )
        class(*),     intent(inout) :: fun_self
        integer,      intent(in)    :: d
        real(kind=8), intent(inout) :: x(d)
        real(kind=8), intent(out)   :: f, grad(d), curv(d,d)
        real(dp) :: jac(1,2)
        jac(1,:) = [1._dp, 1._dp]
        call lsq_terms([x(1) + x(2) - 1._dp], jac, f, grad, curv)
    end subroutine rank1_fdfcurv

    ! L6: r = x - 1 does not depend on y
    subroutine null_fdfcurv( fun_self, x, f, grad, curv, d )
        class(*),     intent(inout) :: fun_self
        integer,      intent(in)    :: d
        real(kind=8), intent(inout) :: x(d)
        real(kind=8), intent(out)   :: f, grad(d), curv(d,d)
        real(dp) :: jac(1,2)
        jac(1,:) = [1._dp, 0._dp]
        call lsq_terms([x(1) - 1._dp], jac, f, grad, curv)
    end subroutine null_fdfcurv

    ! L4: the recorder's moving quadratic
    subroutine moving_fdfcurv( fun_self, x, f, grad, curv, d )
        class(*),     intent(inout) :: fun_self
        integer,      intent(in)    :: d
        real(kind=8), intent(inout) :: x(d)
        real(kind=8), intent(out)   :: f, grad(d), curv(d,d)
        real(dp) :: jac(2,2)
        select type(fun_self)
            type is(lm_recorder)
                fun_self%nev = fun_self%nev + 1
                if( fun_self%pending )then
                    if( all(x == fun_self%cb_x) ) fun_self%nreeval = fun_self%nreeval + 1
                    fun_self%pending = .false.
                endif
                jac      = 0._dp
                jac(1,1) = 1._dp
                jac(2,2) = 1._dp
                call lsq_terms([x(1) - fun_self%c, x(2) + 2._dp], jac, f, grad, curv)
            class default
                f = 0._dp; grad = 0._dp; curv = 0._dp
        end select
    end subroutine moving_fdfcurv

    ! L4: the callback moves the centre at its first call and records where it ran
    subroutine moving_callback( fun_self )
        class(*), intent(inout) :: fun_self
        select type(fun_self)
            type is(lm_recorder)
                fun_self%ncb = fun_self%ncb + 1
                fun_self%c   = 1.5_dp
                if( fun_self%clear_unchanged .and. fun_self%ncb > 1 ) fun_self%spec%callback_changed = .false.
                if( all(fun_self%cb_x == fun_self%spec%x_8) ) fun_self%nrepeat = fun_self%nrepeat + 1
                fun_self%cb_x    = fun_self%spec%x_8
                fun_self%pending = .true.
        end select
    end subroutine moving_callback

    pure subroutine lsq_terms( r, jac, f, grad, curv )
        real(dp), intent(in)  :: r(:), jac(:,:)
        real(dp), intent(out) :: f, grad(:), curv(:,:)
        f    = 0.5_dp*sum(r**2)
        grad = matmul(transpose(jac), r)
        curv = matmul(transpose(jac), jac)
    end subroutine lsq_terms

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
        call assert_true(associated(spec%costfun),  'set_costfun attaches the cost function')
        lims(:,1) = -1.
        lims(:,2) =  1.
        call spec%set_limits(lims)
        call assert_real(1., spec%limits(1,2), 1.e-6, 'set_limits replaces the limits')
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
        spec%gcostfun => quad1d_grad
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
        spec%gcostfun => quad1d_grad
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
        spec%gcostfun => dircost_grad
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
        real :: lowest_cost, first_cost, first_x(2)
        integer :: first_nevals, i
        write(*,'(A)') 'test_lbfgsb_rosenbrock'
        call spec%specify('lbfgsb', 2, limits=LIMS2, maxits=1000)
        call spec%set_costfun(rosenbrock)
        spec%gcostfun => rosenbrock_grad
        call ofac%new(spec, opt)
        spec%x = [-1.2, 1.0]
        call opt%minimize(spec, opt, lowest_cost)
        call assert_real(1.0, spec%x(1), 2.e-3, 'Rosenbrock minimum from the classic start (x)')
        call assert_real(1.0, spec%x(2), 4.e-3, 'Rosenbrock minimum from the classic start (y)')
        call assert_true(lowest_cost < 1.e-5, 'Rosenbrock cost at the minimum')
        call assert_true(spec%converged, 'L-BFGS-B converges on Rosenbrock')
        first_x = spec%x
        first_cost = lowest_cost
        first_nevals = spec%nevals
        do i = 1, 3
            ! Leave different Hessian history in a reused optimizer, then
            ! replay the same problem: START must be independent of history.
            spec%x = [4.0, -4.0]
            call opt%minimize(spec, opt, lowest_cost)
            spec%x = [-1.2, 1.0]
            call opt%minimize(spec, opt, lowest_cost)
            call assert_true(all(spec%x == first_x), 'reused L-BFGS-B replays the solution exactly')
            call assert_true(lowest_cost == first_cost, 'reused L-BFGS-B replays the cost exactly')
            call assert_int(first_nevals, spec%nevals, 'reused L-BFGS-B replays the evaluation count')
        enddo
        call opt%kill
        deallocate(opt)
        call spec%kill
    end subroutine test_lbfgsb_rosenbrock

    !---------------- Levenberg-Marquardt ('lm', lm_optimizer.md section 4) ----------------

    ! L1: Rosenbrock in residual form from the classic start reaches (1, 1) within 1e-6 (guessed); the
    ! evaluations of L-BFGS-B on the same cost are reported, not gated. factr 1e7, pgtol 1e-8 (ruling R3)
    subroutine test_lm_rosenbrock()
        class(optimizer), pointer :: opt => null()
        type(opt_factory) :: ofac
        type(opt_spec)    :: spec
        real    :: lowest_cost
        integer :: nevals_lm
        write(*,'(A)') 'test_lm_rosenbrock'
        call spec%specify('lm', 2, limits=LIMS2, maxits=1000, factr=1.d7, pgtol=1.d-8)
        call spec%set_fdfcurvfun_8(rosen_res_fdfcurv)
        call ofac%new(spec, opt)
        spec%x = [-1.2, 1.0]
        call opt%minimize(spec, opt, lowest_cost)
        call assert_true(abs(spec%x_8(1) - 1._dp) < 1.e-6_dp .and. abs(spec%x_8(2) - 1._dp) < 1.e-6_dp, &
            &'L1: lm reaches the Rosenbrock minimum within 1e-6')
        call assert_true(spec%converged, 'L1: lm reports convergence')
        call assert_int(spec%nevals, spec%ngevals, 'L1: one gradient (and curvature) per cost evaluation')
        nevals_lm = spec%nevals
        call opt%kill
        deallocate(opt)
        call spec%specify('lbfgsb', 2, limits=LIMS2, maxits=1000)
        call spec%set_costfun_8(rosen_res_cost)
        call spec%set_gcostfun_8(rosen_res_grad)
        call ofac%new(spec, opt)
        spec%x = [-1.2, 1.0]
        call opt%minimize(spec, opt, lowest_cost)
        write(*,'(A,I0,A,I0,A,2ES12.4)') '  L1 evaluations: lm ', nevals_lm, ', lbfgsb ', spec%nevals, &
            &'; lbfgsb error ', abs(spec%x_8 - 1._dp)
        call opt%kill
        deallocate(opt)
        call spec%kill
    end subroutine test_lm_rosenbrock

    ! L2 (review 1): the coupled quadratic on [0,10] x [-10,10] from (0,0) on the bound reaches the constrained
    ! minimum (0,-0.02) within 1e-6 (guessed) with a projected gradient below 1e-8; plain clipping (the
    ! engine with the box clip and no free hook) must fail it: the negative control. factr 0 (no relative
    ! stop, as in L-BFGS-B), pgtol 1e-10 (ruling R3)
    subroutine test_lm_coupled_bound()
        class(optimizer), pointer :: opt => null()
        type(opt_factory)       :: ofac
        type(opt_spec)          :: spec
        type(clip_only_problem) :: clip_only
        type(lm_options)        :: options
        type(lm_outcome)        :: outcome
        real(dp) :: f, grad(2), curv(2,2), pgrad(2)
        real     :: lims(2,2), lowest_cost
        write(*,'(A)') 'test_lm_coupled_bound'
        lims(:,1) = [ 0., -10.]
        lims(:,2) = [10.,  10.]
        call spec%specify('lm', 2, limits=lims, factr=0.d0, pgtol=1.d-10)
        call spec%set_fdfcurvfun_8(coupled_fdfcurv)
        call ofac%new(spec, opt)
        spec%x = 0.
        call opt%minimize(spec, opt, lowest_cost)
        call coupled_fdfcurv(opt, spec%x_8, f, grad, curv, 2)
        pgrad = grad
        if( spec%x_8(1) <= 0._dp .and. grad(1) > 0._dp ) pgrad(1) = 0._dp
        write(*,'(A,2ES14.6,A,2ES12.4,A,I0)') '  L2 lm solution ', spec%x_8, ', projected gradient ', pgrad, &
            &', evaluations ', spec%nevals
        call assert_true(abs(spec%x_8(1)) < 1.e-6_dp .and. abs(spec%x_8(2) + 0.02_dp) < 1.e-6_dp, &
            &'L2: lm reaches the constrained minimum (0,-0.02) within 1e-6')
        call assert_true(sqrt(sum(pgrad**2)) < 1.e-8_dp, 'L2: the projected gradient falls below 1e-8')
        call opt%kill
        deallocate(opt)
        call spec%kill
        ! the negative control, at L2's tolerances: it differs from L2 only by the free hook
        clip_only%lower = real(lims(:,1), dp)
        clip_only%upper = real(lims(:,2), dp)
        options%maxits         = 100
        options%predicted_tol  = 0._dp
        options%accepted_tol   = 0._dp
        options%stationary_tol = 1.e-10_dp
        options%scale  = [1._dp, 1._dp]
        options%active = [.true., .true.]
        call lm_solve(clip_only, options, outcome)
        write(*,'(A,2ES14.6,A,I0,A,I0)') '  L2 plain clipping solution ', clip_only%x, ', status ', outcome%status, &
            &', accepted ', outcome%naccepted
        call assert_false(abs(clip_only%x(1)) < 1.e-6_dp .and. abs(clip_only%x(2) + 0.02_dp) < 1.e-6_dp, &
            &'L2 negative control: plain clipping does not reach the constrained minimum')
    end subroutine test_lm_coupled_bound

    ! L3, negative control: with the gradient's sign flipped, lm accepts no step and returns the start after
    ! 8 rejections (1 + 8 evaluations). factr 1e7, pgtol 1e-5, L-BFGS-B's defaults (ruling R3)
    subroutine test_lm_flipped_gradient()
        class(optimizer), pointer :: opt => null()
        type(opt_factory) :: ofac
        type(opt_spec)    :: spec
        real :: lowest_cost
        write(*,'(A)') 'test_lm_flipped_gradient'
        call spec%specify('lm', 2, limits=LIMS2, maxits=1000, factr=1.d7, pgtol=1.d-5)
        call spec%set_fdfcurvfun_8(rosen_res_flipped)
        call ofac%new(spec, opt)
        spec%x = [-1.2, 1.0]
        call opt%minimize(spec, opt, lowest_cost)
        write(*,'(A,2ES14.6,A,I0,A,I0)') '  L3 solution ', spec%x_8, ', evaluations ', spec%nevals, &
            &', iterations ', spec%niter
        call assert_true(all(spec%x_8 == real([-1.2, 1.0], dp)), 'L3: a flipped gradient returns the start')
        call assert_real(real(rosen_res_cost(opt, real([-1.2, 1.0], dp), 2)), lowest_cost, 0., &
            &'L3: the cost is the start''s')
        call assert_int(9, spec%nevals, 'L3: the start and 8 rejected trials are evaluated')
        call opt%kill
        deallocate(opt)
        call spec%kill
    end subroutine test_lm_flipped_gradient

    ! L4: opt_callback runs once per accepted step, at the accepted point, and the objective it changes is
    ! evaluated again there and counted in nevals; the solve continues on the changed objective
    subroutine test_lm_callback()
        class(optimizer), pointer :: opt => null()
        type(opt_factory)       :: ofac
        type(opt_spec), target  :: spec
        type(lm_recorder)       :: rec
        real :: lowest_cost
        write(*,'(A)') 'test_lm_callback'
        call spec%specify('lm', 2, limits=LIMS2, factr=1.d7, pgtol=1.d-8)
        call spec%set_fdfcurvfun_8(moving_fdfcurv)
        spec%opt_callback => moving_callback
        call ofac%new(spec, opt)
        rec%spec => spec
        spec%x = [-3., 3.]
        call opt%minimize(spec, rec, lowest_cost)
        write(*,'(A,I0,A,I0,A,I0)') '  L4 callbacks ', rec%ncb, ', evaluations ', spec%nevals, ', re-evaluations ', rec%nreeval
        call assert_true(rec%ncb > 0, 'L4: the callback runs after an accepted step')
        call assert_int(0, rec%nrepeat, 'L4: the callback runs once per accepted point')
        call assert_int(rec%ncb, rec%nreeval, 'L4: every callback is followed by an evaluation at its point')
        call assert_int(rec%nev, spec%nevals, 'L4: nevals counts every evaluation, the re-evaluations included')
        call assert_true(all(spec%x_8 == rec%cb_x), 'L4: the solution is the last accepted point')
        call assert_true(abs(spec%x_8(1) - 1.5_dp) < 1.e-6_dp .and. abs(spec%x_8(2) + 2._dp) < 1.e-6_dp, &
            &'L4: the solve continues on the objective the callback changed')
        rec%spec => null()
        call opt%kill
        deallocate(opt)
        call spec%kill
    end subroutine test_lm_callback

    ! L4b (R10): the same solve with a callback that clears callback_changed after its first call (the centre no
    ! longer moves): those callbacks add no evaluation, the path and the minimum are L4's
    subroutine test_lm_callback_unchanged()
        integer  :: ncb(2), nevals(2), nreeval(2)
        real(dp) :: x(2,2)
        integer  :: irun
        write(*,'(A)') 'test_lm_callback_unchanged'
        do irun = 1, 2
            call run_moving(irun == 2, ncb(irun), nevals(irun), nreeval(irun), x(:,irun))
        enddo
        write(*,'(A,I0,A,I0,A,I0)') '  L4b callbacks ', ncb(2), ', evaluations ', nevals(2), ', re-evaluations ', nreeval(2)
        call assert_true(ncb(2) > 1, 'L4b: the callback runs at more than one accepted step')
        call assert_int(ncb(1), ncb(2), 'L4b: the same accepted steps as L4')
        call assert_int(1, nreeval(2), 'L4b: only the callback that changed the objective is followed by an evaluation')
        call assert_int(nevals(1) - (ncb(2) - 1), nevals(2), 'L4b: a callback that clears the flag adds no evaluation')
        call assert_true(all(x(:,1) == x(:,2)), 'L4b: the solve reaches L4''s minimum')
    end subroutine test_lm_callback_unchanged

    subroutine run_moving( clear_unchanged, ncb, nevals, nreeval, x )
        logical,  intent(in)  :: clear_unchanged
        integer,  intent(out) :: ncb, nevals, nreeval
        real(dp), intent(out) :: x(2)
        class(optimizer), pointer :: opt => null()
        type(opt_factory)       :: ofac
        type(opt_spec), target  :: spec
        type(lm_recorder)       :: rec
        real :: lowest_cost
        call spec%specify('lm', 2, limits=LIMS2, factr=1.d7, pgtol=1.d-8)
        call spec%set_fdfcurvfun_8(moving_fdfcurv)
        spec%opt_callback => moving_callback
        call ofac%new(spec, opt)
        rec%spec            => spec
        rec%clear_unchanged = clear_unchanged
        spec%x = [-3., 3.]
        call opt%minimize(spec, rec, lowest_cost)
        ncb     = rec%ncb
        nevals  = spec%nevals
        nreeval = rec%nreeval
        x       = spec%x_8
        rec%spec => null()
        call opt%kill
        deallocate(opt)
        call spec%kill
    end subroutine run_moving

    ! L5 (review 2): r = x + y - 1 has a singular curvature; lm reaches f < 1e-12 and x - y keeps its start
    ! value (to rounding, 1e-12, guessed): the null direction gets no step. factr 1e7, pgtol 1e-8
    subroutine test_lm_rank_deficient()
        class(optimizer), pointer :: opt => null()
        type(opt_factory) :: ofac
        type(opt_spec)    :: spec
        real :: lowest_cost
        write(*,'(A)') 'test_lm_rank_deficient'
        call spec%specify('lm', 2, factr=1.d7, pgtol=1.d-8)
        call spec%set_fdfcurvfun_8(rank1_fdfcurv)
        call ofac%new(spec, opt)
        spec%x = [0.25, -0.5]
        call opt%minimize(spec, opt, lowest_cost)
        write(*,'(A,2ES24.16,A,I0)') '  L5 solution ', spec%x_8, ', evaluations ', spec%nevals
        call assert_true(0.5_dp*(spec%x_8(1) + spec%x_8(2) - 1._dp)**2 < 1.e-12_dp, 'L5: lm reaches f < 1e-12')
        call assert_true(abs((spec%x_8(1) - spec%x_8(2)) - 0.75_dp) < 1.e-12_dp, 'L5: the null direction x - y keeps its value')
        call opt%kill
        deallocate(opt)
        call spec%kill
    end subroutine test_lm_rank_deficient

    ! L6 (ruling R4): r = x - 1 in two variables leaves y without gradient or curvature; lm reaches f < 1e-12,
    ! y keeps its start value exactly, and the engine reports an accepted solve, not "no improvement".
    ! factr 1e7, pgtol 1e-8
    subroutine test_lm_null_coordinate()
        class(optimizer), pointer :: opt => null()
        type(opt_factory)       :: ofac
        type(opt_spec)          :: spec
        type(clip_only_problem) :: null_problem
        type(lm_options)        :: options
        type(lm_outcome)        :: outcome
        real :: lowest_cost
        write(*,'(A)') 'test_lm_null_coordinate'
        call spec%specify('lm', 2, factr=1.d7, pgtol=1.d-8)
        call spec%set_fdfcurvfun_8(null_fdfcurv)
        call ofac%new(spec, opt)
        spec%x = [0.25, -0.5]
        call opt%minimize(spec, opt, lowest_cost)
        write(*,'(A,2ES24.16,A,I0)') '  L6 solution ', spec%x_8, ', evaluations ', spec%nevals
        call assert_true(0.5_dp*(spec%x_8(1) - 1._dp)**2 < 1.e-12_dp, 'L6: lm reaches f < 1e-12')
        call assert_true(spec%x_8(2) == -0.5_dp, 'L6: the null coordinate keeps its start value')
        call opt%kill
        deallocate(opt)
        call spec%kill
        ! the engine's status on the same cost
        null_problem%null_cost = .true.
        null_problem%x         = [0.25_dp, -0.5_dp]
        null_problem%trial     = null_problem%x
        options%stationary_tol = 1.e-8_dp
        options%scale  = [1._dp, 1._dp]
        options%active = [.true., .true.]
        call lm_solve(null_problem, options, outcome)
        call assert_int(LM_ACCEPTED, outcome%status, 'L6: the engine accepts the solve (not no improvement)')
        call assert_true(null_problem%x(2) == -0.5_dp, 'L6: the engine leaves the null coordinate')
    end subroutine test_lm_null_coordinate

    subroutine clip_only_eval( self, f, g, curv )
        class(clip_only_problem), intent(inout) :: self
        real(dp),                 intent(out)   :: f, g(:), curv(:,:)
        if( self%null_cost )then
            call null_fdfcurv(self, self%trial, f, g, curv, 2)
        else
            call coupled_fdfcurv(self, self%trial, f, g, curv, 2)
        endif
    end subroutine clip_only_eval

    subroutine clip_only_limit_step( self, step, bounded )
        class(clip_only_problem), intent(inout) :: self
        real(dp),                 intent(inout) :: step(:)
        logical,                  intent(out)   :: bounded
        real(dp) :: clipped(2)
        clipped = min(max(self%x + step, self%lower), self%upper) - self%x
        bounded = any(clipped /= step)
        step    = clipped
    end subroutine clip_only_limit_step

    logical function clip_only_propose( self, step ) result( feasible )
        class(clip_only_problem), intent(inout) :: self
        real(dp),                 intent(in)    :: step(:)
        self%trial = self%x + step
        feasible   = .true.
    end function clip_only_propose

    logical function clip_only_accept( self, ratio ) result( changed )
        class(clip_only_problem), intent(inout) :: self
        real(dp),                 intent(in)    :: ratio
        self%x  = self%trial
        changed = .false.
    end function clip_only_accept

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
        ! no start: a random point in the box
        spec%x       = 0.
        spec%nostart = .true.
        call opt%minimize(spec, opt, lowest_cost)
        call assert_real(XMIN(1), spec%x(1), 2.e-3, 'simplex from a random start (x)')
        call assert_real(XMIN(2), spec%x(2), 2.e-3, 'simplex from a random start (y)')
        call opt%kill
        deallocate(opt)
        call spec%kill
    end subroutine test_simplex_quadratic

    ! quad2d, recording the first point evaluated after arming
    function quad2d_rec( fun_self, x, d ) result( r )
        class(*), intent(inout) :: fun_self
        integer,  intent(in)    :: d
        real,     intent(in)    :: x(d)
        real                    :: r
        select type(fun_self)
            type is(eval_recorder)
                if( fun_self%armed )then
                    fun_self%first = x(1:2)
                    fun_self%armed = .false.
                endif
        end select
        r = sum((x - XMIN)**2)
    end function quad2d_rec

    ! an all-zero x is a start (the origin); nostart says there is none
    subroutine test_start_at_origin()
        class(optimizer), pointer :: opt => null()
        type(opt_factory)   :: ofac
        type(opt_spec)      :: spec
        type(eval_recorder) :: rec
        character(len=7), parameter :: OPTS(2) = ['simplex', 'de     ']
        real    :: lowest_cost
        integer :: iopt
        write(*,'(A)') 'test_start_at_origin'
        call set_fixed_seed(20261105)
        do iopt = 1,2
            call spec%specify(trim(OPTS(iopt)), 2, limits=LIMS2, maxits=200, ftol=1.e-7)
            call spec%set_costfun(quad2d_rec)
            call ofac%new(spec, opt)
            spec%x    = 0.
            rec%armed = .true.
            call opt%minimize(spec, rec, lowest_cost)
            call assert_true(all(rec%first == 0.), trim(OPTS(iopt))//' starts at the origin when x is zero')
            spec%nostart = .true.
            rec%armed    = .true.
            call opt%minimize(spec, rec, lowest_cost)
            call assert_true(any(rec%first /= 0.), trim(OPTS(iopt))//' draws its start under nostart')
            call opt%kill
            deallocate(opt)
            call spec%kill
        end do
    end subroutine test_start_at_origin

    ! DE takes the npop given to specify
    subroutine test_de_npop()
        class(optimizer), pointer :: opt => null()
        type(opt_factory)   :: ofac
        type(opt_spec)      :: spec
        type(eval_recorder) :: rec
        real :: lowest_cost
        write(*,'(A)') 'test_de_npop'
        call set_fixed_seed(20261106)
        call spec%specify('de', 2, limits=LIMS2, maxits=1, npop=20)
        call spec%set_costfun(quad2d_rec)
        call ofac%new(spec, opt)
        call assert_int(20, spec%npop, 'DE keeps the npop given to specify')
        spec%nostart = .true.
        call opt%minimize(spec, rec, lowest_cost)
        call assert_int(21, spec%nevals, 'twenty members and one trial')
        call opt%kill
        deallocate(opt)
        call spec%kill
    end subroutine test_de_npop

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
