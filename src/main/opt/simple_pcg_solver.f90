!@descr: rank-1 preconditioned conjugate-gradient engine with client-supplied operator, preconditioner and inner product
module simple_pcg_solver
use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
use simple_core_module_api, only: dp, logfhandle, timer_int_kind, simple_exception
implicit none
private
#include "simple_local_flags.inc"

public :: pcg_operator, pcg_solver_options, pcg_solver_outcome, pcg_solve
public :: PCG_STOP_INDEFINITE, PCG_XTOL, PCG_RESID_REPLACE, PCG_RHO_FLOOR_FRAC

character(len=*), parameter :: PCG_STOP_INDEFINITE = 'indefinite'
real,             parameter :: PCG_XTOL            = 1.5e-2
integer,          parameter :: PCG_RESID_REPLACE   = 25
real,             parameter :: PCG_RHO_FLOOR_FRAC  = 1.0e-2

!> A client exposes its natural operator through contiguous rank-1 views. The
!! engine owns only Krylov vectors; it never copies the caller's solution.
type, abstract :: pcg_operator
  contains
    procedure(pcg_size_i),    deferred :: size
    procedure(pcg_apply_i),   deferred :: apply
    procedure(pcg_precond_i), deferred :: precond
    procedure(pcg_dot_i),     deferred :: dot
end type pcg_operator

type :: pcg_solver_options
    integer           :: maxits = 50
    real              :: rtol = 1.0e-4
    real              :: xtol = PCG_XTOL
    integer           :: residual_replace = PCG_RESID_REPLACE
    real              :: start_max_rel_residual = -1.0
    logical           :: rescale_start = .false.
    logical           :: cold_restart = .false.
    logical           :: track_preconditioned_residual = .true.
    logical           :: record_history = .true.
    integer           :: verbose = 0
    character(len=96) :: tag = 'PCG'
end type pcg_solver_options

!> Solver facts shared by reconstruction and coupled FLEX solves. The
!! closed-form fields are populated by reconstruction after this engine returns.
type :: pcg_solver_outcome
    character(len=24) :: stop_reason          = 'not_started'
    integer           :: iteration_count      = 0
    integer           :: requested_maxits     = 0
    real              :: initial_rel_residual = 0.0
    real              :: final_rel_residual   = 0.0
    real              :: final_rel_residual_m = -1.0
    real              :: final_rel_update     = 0.0
    real(dp)          :: failure_curvature    = 0.0_dp
    integer           :: failure_iteration    = 0
    real(dp)          :: restart_trigger_curvature = 0.0_dp
    integer           :: restart_trigger_iteration = 0
    logical           :: cold_restart_used    = .false.
    logical           :: start_rejected       = .false.
    real              :: rejected_start_initial = 0.0
    logical           :: converged            = .false.
    real              :: rhs_norm             = 0.0
    real              :: start_norm           = 0.0
    real              :: start_corr           = 0.0
    real              :: start_scale          = 0.0
    integer           :: iters_to_1e2          = 0
    real              :: seconds               = 0.0
    real              :: closed_form_rel_residual   = -1.0
    real              :: closed_form_rel_residual_m = -1.0
    real              :: closed_form_fsc05_res      = 0.0
    real              :: closed_form_fsc0143_res    = 0.0
    real              :: closed_form_min_fsc_inband = 0.0
    integer           :: closed_form_band_shell     = 0
    real, allocatable :: rel_residual_history(:)
    real, allocatable :: rel_update_history(:)
    real, allocatable :: preconditioned_residual_history(:)
    real, allocatable :: iteration_seconds(:)
  contains
    procedure :: kill => pcg_solver_outcome_kill
end type pcg_solver_outcome

abstract interface
    integer function pcg_size_i( self ) result(n)
        import :: pcg_operator
        class(pcg_operator), intent(in) :: self
    end function pcg_size_i

    subroutine pcg_apply_i( self, x, y )
        import :: pcg_operator
        class(pcg_operator), intent(inout) :: self
        real, contiguous, target, intent(in)  :: x(:)
        real, contiguous, target, intent(out) :: y(:)
    end subroutine pcg_apply_i

    subroutine pcg_precond_i( self, r, z )
        import :: pcg_operator
        class(pcg_operator), intent(inout) :: self
        real, contiguous, target, intent(in)  :: r(:)
        real, contiguous, target, intent(out) :: z(:)
    end subroutine pcg_precond_i

    real(dp) function pcg_dot_i( self, a, b ) result(value)
        import :: dp, pcg_operator
        class(pcg_operator), intent(in) :: self
        real, contiguous, target, intent(in) :: a(:), b(:)
    end function pcg_dot_i
end interface

contains

    !> Solve H x = b. Client arrays are contiguous rank-1 views; only r, p, z
    !! and Hp are allocated here. A cold retry, when requested, is permitted
    !! only for an indefinite warm attempt.
    subroutine pcg_solve( op, b, x, options, outcome )
        class(pcg_operator),       intent(inout) :: op
        real, contiguous,          intent(in)    :: b(:)
        real, contiguous,          intent(inout) :: x(:)
        type(pcg_solver_options),  intent(in)    :: options
        type(pcg_solver_outcome),  intent(out)   :: outcome
        real, allocatable :: r(:), p(:), hp(:), z(:)
        real, allocatable :: rhist(:), xhist(:), mhist(:), thist(:)
        real(dp) :: rho, rho_new, rho0, alpha, beta, pHp
        real(dp) :: bnorm, rnorm, xnorm, dxnorm, mnorm, dxx, hnorm, bh
        integer :: n, iter, n_done, attempt, max_attempts
        logical :: l_warm, retry, stop_rtol, stop_xtol
        integer(timer_int_kind) :: t_it

        n = op%size()
        if( n < 1 .or. size(b) /= n .or. size(x) /= n ) THROW_HARD('PCG vector size differs from the operator contract')
        if( options%maxits < 1 ) THROW_HARD('PCG maxits must be at least 1')
        if( .not. ieee_is_finite(options%rtol) .or. .not. ieee_is_finite(options%xtol) ) &
            &THROW_HARD('PCG stopping tolerances must be finite')
        if( options%residual_replace < 1 ) THROW_HARD('PCG residual replacement interval must be positive')

        outcome%requested_maxits = options%maxits
        outcome%stop_reason = 'maxits'
        if( options%rtol <= 0.0 ) outcome%stop_reason = 'fixed_iterations'
        allocate(r(n), p(n), hp(n), z(n))
        allocate(rhist(options%maxits), xhist(options%maxits), thist(options%maxits), source=0.0)
        allocate(mhist(options%maxits), source=-1.0)

        bnorm = sqrt(op%dot(b,b))
        if( bnorm <= 0.0_dp ) THROW_HARD('zero right-hand side; nothing to solve with PCG')
        outcome%rhs_norm   = real(bnorm)
        outcome%start_norm = real(sqrt(op%dot(x,x)))
        l_warm = outcome%start_norm > 0.0
        max_attempts = merge(2, 1, options%cold_restart)
        n_done = 0
        dxx = 0.0_dp
        rnorm = bnorm

        do attempt = 1, max_attempts
            retry = .false.
            if( l_warm )then
                call op%apply(x, hp)
                r = b - hp
            else
                hp = 0.0
                r  = b
            endif
            rnorm = sqrt(op%dot(r,r))

            if( attempt == 1 )then
                outcome%initial_rel_residual = real(rnorm / bnorm)
                hnorm = sqrt(op%dot(hp,hp))
                bh = op%dot(b,hp)
                outcome%start_corr  = real(bh / max(bnorm*hnorm, 1.0e-30_dp))
                outcome%start_scale = real(bh / max(hnorm*hnorm, 1.0e-30_dp))
                if( l_warm .and. options%rescale_start .and. &
                    &ieee_is_finite(real(outcome%start_scale,dp)) .and. outcome%start_scale > 0.0 )then
                    x  = outcome%start_scale * x
                    hp = outcome%start_scale * hp
                    r  = b - hp
                    rnorm = sqrt(op%dot(r,r))
                    outcome%start_norm = real(sqrt(op%dot(x,x)))
                    outcome%initial_rel_residual = real(rnorm / bnorm)
                endif
                if( l_warm .and. options%start_max_rel_residual >= 0.0 .and. &
                    &outcome%initial_rel_residual > options%start_max_rel_residual )then
                    outcome%start_rejected = .true.
                    outcome%rejected_start_initial = outcome%initial_rel_residual
                    x = 0.0
                    r = b
                    rnorm = bnorm
                    l_warm = .false.
                    outcome%start_norm = 0.0
                    outcome%initial_rel_residual = 1.0
                endif
                if( options%verbose > 0 )then
                    write(logfhandle,'(A,A,A,ES10.3,A,F8.4,A,ES10.3)') '>>> ', trim(options%tag), &
                        &' start: rel resid=', outcome%initial_rel_residual, '  corr(b,Hx0)=', &
                        &outcome%start_corr, '  scale=', outcome%start_scale
                endif
            endif

            call op%precond(r, z)
            p    = z
            rho  = op%dot(r,z)
            rho0 = rho
            if( rho0 <= 0.0_dp ) THROW_HARD('non-positive initial dot(r,z); PCG preconditioner is not positive definite')
            n_done = 0
            dxx = 0.0_dp

            do iter = 1, options%maxits
                if( options%record_history ) t_it = pcg_tic()
                call op%apply(p, hp)
                pHp = op%dot(p,hp)
                if( .not. ieee_is_finite(pHp) .or. pHp <= 0.0_dp )then
                    outcome%failure_curvature = pHp
                    outcome%failure_iteration = iter
                    if( l_warm .and. options%cold_restart .and. attempt == 1 )then
                        outcome%restart_trigger_curvature = pHp
                        outcome%restart_trigger_iteration = iter
                        outcome%cold_restart_used = .true.
                        x = 0.0
                        l_warm = .false.
                        retry = .true.
                        write(logfhandle,'(A,A,ES11.3,A,I0,A)') '>>> ', trim(options%tag), real(pHp), &
                            &' curvature at iteration ', iter, ': warm start discarded, cold restart'
                    else
                        outcome%stop_reason = PCG_STOP_INDEFINITE
                        outcome%converged = .false.
                    endif
                    exit
                endif

                alpha = rho / pHp
                x = x + real(alpha) * p
                r = r - real(alpha) * hp
                if( mod(iter, options%residual_replace) == 0 )then
                    call op%apply(x, hp)
                    r = b - hp
                endif
                n_done = iter
                rnorm  = sqrt(op%dot(r,r))
                xnorm  = sqrt(op%dot(x,x))
                dxnorm = abs(alpha) * sqrt(op%dot(p,p))
                dxx    = dxnorm / max(xnorm, epsilon(1.0_dp))
                rhist(iter) = real(rnorm / bnorm)
                xhist(iter) = real(dxx)
                if( outcome%iters_to_1e2 == 0 .and. rnorm / bnorm <= 1.0e-2_dp ) outcome%iters_to_1e2 = iter
                stop_rtol = options%rtol > 0.0 .and. rnorm / bnorm <= real(options%rtol,dp)
                stop_xtol = options%rtol > 0.0 .and. dxx <= real(options%xtol,dp)
                if( options%verbose > 0 )then
                    write(logfhandle,'(A,A,A,I4,A,ES10.3,A,ES10.3,A,ES10.3)') '>>> ', trim(options%tag), &
                        &' it', iter, '  rel resid=', real(rnorm/bnorm), '  dx/x=', real(dxx), &
                        &'  alpha=', real(alpha)
                endif

                if( stop_rtol .or. stop_xtol .or. iter == options%maxits )then
                    if( options%track_preconditioned_residual )then
                        call op%precond(r, z)
                        rho_new = op%dot(r,z)
                        mnorm = sqrt(abs(rho_new) / rho0)
                        mhist(iter) = real(mnorm)
                    endif
                    if( options%record_history ) thist(iter) = real(pcg_toc(t_it))
                    if( stop_rtol )then
                        outcome%stop_reason = 'rtol'
                        outcome%converged = .true.
                    else if( stop_xtol )then
                        outcome%stop_reason = 'xtol'
                        outcome%converged = .true.
                    endif
                    exit
                endif

                call op%precond(r, z)
                rho_new = op%dot(r,z)
                if( options%track_preconditioned_residual )then
                    mnorm = sqrt(abs(rho_new) / rho0)
                    mhist(iter) = real(mnorm)
                endif
                if( options%record_history ) thist(iter) = real(pcg_toc(t_it))
                beta = rho_new / rho
                p = z + real(beta) * p
                rho = rho_new
            end do

            if( retry ) cycle
            exit
        end do

        outcome%iteration_count = n_done
        outcome%final_rel_residual = real(rnorm / bnorm)
        outcome%final_rel_update = real(dxx)
        if( n_done > 0 .and. options%track_preconditioned_residual ) &
            &outcome%final_rel_residual_m = mhist(n_done)
        if( options%record_history )then
            allocate(outcome%rel_residual_history(n_done), source=rhist(1:n_done))
            allocate(outcome%rel_update_history(n_done), source=xhist(1:n_done))
            allocate(outcome%preconditioned_residual_history(n_done), source=mhist(1:n_done))
            allocate(outcome%iteration_seconds(n_done), source=thist(1:n_done))
        endif
        deallocate(r, p, hp, z, rhist, xhist, mhist, thist)
    end subroutine pcg_solve

    subroutine pcg_solver_outcome_kill( self )
        class(pcg_solver_outcome), intent(inout) :: self
        if( allocated(self%rel_residual_history) ) deallocate(self%rel_residual_history)
        if( allocated(self%rel_update_history) ) deallocate(self%rel_update_history)
        if( allocated(self%preconditioned_residual_history) ) deallocate(self%preconditioned_residual_history)
        if( allocated(self%iteration_seconds) ) deallocate(self%iteration_seconds)
        self%stop_reason = 'not_started'
        self%iteration_count = 0
        self%requested_maxits = 0
        self%initial_rel_residual = 0.0
        self%final_rel_residual = 0.0
        self%final_rel_residual_m = -1.0
        self%final_rel_update = 0.0
        self%failure_curvature = 0.0_dp
        self%failure_iteration = 0
        self%restart_trigger_curvature = 0.0_dp
        self%restart_trigger_iteration = 0
        self%cold_restart_used = .false.
        self%start_rejected = .false.
        self%rejected_start_initial = 0.0
        self%converged = .false.
        self%rhs_norm = 0.0
        self%start_norm = 0.0
        self%start_corr = 0.0
        self%start_scale = 0.0
        self%iters_to_1e2 = 0
        self%seconds = 0.0
        self%closed_form_rel_residual = -1.0
        self%closed_form_rel_residual_m = -1.0
        self%closed_form_fsc05_res = 0.0
        self%closed_form_fsc0143_res = 0.0
        self%closed_form_min_fsc_inband = 0.0
        self%closed_form_band_shell = 0
    end subroutine pcg_solver_outcome_kill

    integer(timer_int_kind) function pcg_tic() result(tstart)
        call system_clock(count=tstart)
    end function pcg_tic

    real(dp) function pcg_toc( tstart ) result(seconds)
        integer(timer_int_kind), intent(in) :: tstart
        integer(timer_int_kind) :: tend, rate
        call system_clock(count=tend, count_rate=rate)
        seconds = real(tend-tstart,dp) / real(rate,dp)
    end function pcg_toc

end module simple_pcg_solver
