! function minimization by BFGS algorithm, translated from gsl 2.4 (vector_bfgs2.c)

module simple_bfgs2_opt
#include "simple_lib.f08"

use simple_optimizer, only: optimizer
use simple_opt_helpers
implicit none

public :: bfgs2_opt
private

type :: bfgs2_wrapper
    ! cached values, for x(alpha) = x + alpha * p
    real(dp) :: f_alpha, df_alpha
    real(dp), allocatable :: x_alpha(:), g_alpha(:)
    ! cache "keys"
    real(dp) :: f_cache_key, df_cache_key, x_cache_key, g_cache_key
end type bfgs2_wrapper

type, extends(optimizer) :: bfgs2_opt
    private
    integer :: bfgs2_iter       !< internal iter for bfgs2 algorithm
    real(dp) :: step, g0norm, pnorm, delta_f, f
    real(dp) :: fp0 ! f'(0) for f(x-alpha*p)
    real(dp), allocatable :: x0(:), g0(:), p(:), gradient(:)
    ! work space
    real(dp), allocatable :: dx0(:), dg0(:), dx(:)
    ! wrapper function
    type(bfgs2_wrapper) :: wrapper
    ! minimization parameters
    real(dp) :: rho, sigma, tau1, tau2, tau3
    integer :: order
    logical :: exists=.false.
contains
    procedure :: new      => new_bfgs2_opt
    procedure :: minimize => bfgs2_minimize
    procedure :: kill     => kill_bfgs2_opt
end type

contains

    !> \brief  is a constructor
    subroutine new_bfgs2_opt( self, spec )
        use simple_opt_spec, only: opt_spec
        use simple_syslib,   only: alloc_errchk
        class(bfgs2_opt), intent(inout) :: self !< instance
        class(opt_spec), intent(inout)  :: spec !< specification
        call self%kill
        allocate(self%x0(spec%ndim),self%p(spec%ndim),self%g0(spec%ndim),&
            & self%gradient(spec%ndim), self%dx0(spec%ndim), self%dg0(spec%ndim),&
            & self%dx(spec%ndim), self%wrapper%x_alpha(spec%ndim), &
            & self%wrapper%g_alpha(spec%ndim), stat=alloc_stat)
        allocchk('In: new_bfgs2_opt; simple_bfgs2_opt')
        self%exists = .true.
    end subroutine new_bfgs2_opt

    !>  \brief  nonlinear conjugate gradient minimizer
    subroutine bfgs2_minimize( self, spec, fun_self, lowest_cost )
        use simple_opt_spec, only: opt_spec
        use simple_opt_subs, only: lnsrch
        class(bfgs2_opt), intent(inout) :: self        !< instance
        class(opt_spec),  intent(inout) :: spec        !< specification
        class(*),         intent(inout) :: fun_self    !< self-pointer for cost function
        real, intent(out)               :: lowest_cost !< minimum function value
        if( .not. associated(spec%costfun) )then
            stop 'cost function not associated in opt_spec; bfgs2_minimize; simple_bfgs2_opt'
        endif
        if( .not. associated(spec%gcostfun) )then
            stop 'gradient of cost function not associated in opt_spec; bfgs2_minimize; simple_bfgs2_opt'
        endif
        spec%x_8 = spec%x
        ! initialise nevals counters
        spec%nevals  = 0
        spec%ngevals = 0
        call bfgs2min
        
    contains

        !>  \brief  nonlinear conjugate gradient minimizer
        subroutine bfgs2min
            integer :: status
            integer :: iter
            iter        = 0
            self%step   = spec%max_step
            self%x0     = 0.0_8
            self%g0     = 0.0_8
            self%p      = 0.0_8
            call bfgs2_set
            do
                iter = iter + 1
                status = bfgs2_iterate()
                if (status == OPT_STATUS_ERROR) then
                    write (*,*) 'simple_bfgs2_opt: error in minimizer routine'
                    return
                end if
                status = test_gradient(self%gradient, real(spec%gtol, dp))
                if ((global_debug).and.(global_verbose)) then
                    if (status == OPT_STATUS_SUCCESS) then
                        write (*,*) 'Minimum found at:'
                    end if
                    write (*,*) iter, 'x = ', spec%x_8, 'f = ', self%f
                end if
                if ((status .ne. OPT_STATUS_CONTINUE) .or. (iter > spec%maxits)) then
                    exit
                end if
            end do
            if (status == OPT_STATUS_SUCCESS) then
                spec%converged = .true.
            else
                spec%converged = .false.
            end if
            lowest_cost = real(self%f, kind=4)
        end subroutine bfgs2min

        subroutine bfgs2_set
            self%bfgs2_iter = 0
            self%step       = spec%max_step
            self%delta_f    = 0.0_8
            ! Use the gradient as the initial direction
            call spec%eval_fdf(fun_self, spec%x_8, self%f, self%gradient)
            self%x0         = spec%x_8
            self%g0         = self%gradient
            self%g0norm     = norm2(self%g0)
            self%p          = -self%gradient / self%g0norm
            self%pnorm      = norm2(self%p)
            self%fp0        = -self%g0norm
            ! Prepare the wrapper
            call prepare_wrapper
            ! Prepare 1d minimisation parameters
            self%rho        = 0.01_8
            self%sigma      = spec%ftol
            self%tau1       = 9.0_8
            self%tau2       = 0.05_8
            self%tau3       = 0.5_8
            self%order      = 3  ! use cubic interpolation where possible
        end subroutine bfgs2_set

        subroutine prepare_wrapper
            self%wrapper%x_alpha      = self%x0
            self%wrapper%x_cache_key  = 0.0_8
            self%wrapper%f_alpha      = self%f
            self%wrapper%f_cache_key  = 0.0_8
            self%wrapper%g_alpha      = self%g0
            self%wrapper%g_cache_key  = 0.0_8
            self%wrapper%df_alpha     = slope()
            self%wrapper%df_cache_key = 0.0_8
        end subroutine prepare_wrapper

        function slope() result(res)
            real(dp) :: res
            res = dot_product(self%wrapper%g_alpha, self%p)
        end function slope

        function bfgs2_iterate() result(status)
            real(dp) :: alpha, alpha1, pg, dir, f0, del
            real(dp) :: dxg, dgg, dxdg, dgnorm, A, B;
            integer :: status
            alpha = 0.0_8
            f0 = self%f
            if ((self%pnorm == 0.0_8) .or. (self%g0norm == 0.0_8) .or. (self%fp0 == 0.0_8)) then
                status = OPT_STATUS_ERROR
                return
            end if
            if (self%delta_f < 0.0_8) then
                del = max(-self%delta_f, EPSILON(alpha) * abs(f0))
                alpha1 = min(1.0_8, 2.0_8 * del / (-self%fp0))
            else
                alpha1 = abs(self%step)
            end if
            ! line minimisation, with cubic interpolation (order = 3)
            status = linear_minimize(self%rho, self%sigma, self%tau1, self%tau2, self%tau3, &
                & self%order, alpha1,  alpha)
            if (status .ne. OPT_STATUS_SUCCESS) then
                return
            end if
            call update_position(alpha, spec%x_8, self%f, self%gradient)
            spec%x = real(spec%x_8, kind=4)
            self%delta_f = self%f - f0
            ! Choose a new direction for the next step
            ! This is the BFGS update:
            ! p' = g1 - A dx - B dg
            ! A = - (1+ dg.dg/dx.dg) B + dg.g/dx.dg
            ! B = dx.g/dx.dg

            ! dx0 = x - x0
            self%dx0 = spec%x_8 - self%x0
            self%dx  = self%dx0 ! keep a copy
            ! dg0 = g - g0
            self%dg0 = self%gradient - self%g0
            dxg    = dot_product(self%dx0, self%gradient)
            dgg    = dot_product(self%dg0, self%gradient)
            dxdg   = dot_product(self%dx0, self%dg0)
            dgnorm = norm2(self%dg0)
            if (dxdg .ne. 0.0_8) then
                B = dxg / dxdg
                A = -(1.0_8 + dgnorm * dgnorm / dxdg) * B + dgg / dxdg
            else
                B = 0.0_8
                A = 0.0_8
            end if
            self%p      = self%gradient - A*self%dx0 - B*self%dg0
            self%g0     = self%gradient
            self%x0     = spec%x_8
            self%g0norm = norm2(self%g0)
            self%pnorm  = norm2(self%p)
            ! update direction and fp0
            pg = dot_product(self%p, self%gradient)
            if (pg >= 0.) then
                dir = -1.0_8
            else
                dir = 1.0_8
            end if
            self%p     = self%p * dir / self%pnorm
            self%pnorm = norm2(self%p)
            self%fp0   = dot_product(self%p, self%g0)
            call change_direction
            status = OPT_STATUS_SUCCESS
        end function bfgs2_iterate

        function linear_minimize(rho, sigma, tau1, tau2, tau3, order, alpha1, alpha_new) result(status)
            ! recommended values from Fletcher are
            ! rho = 0.01, sigma = 0.1, tau1 = 9, tau2 = 0.05, tau3 = 0.5
            real(dp), intent(in) :: rho, sigma, tau1, tau2, tau3, alpha1
            integer, intent(in) :: order
            real(dp), intent(out) :: alpha_new
            integer :: status
            real(dp), volatile :: f0, fp0, falpha, falpha_prev, fpalpha, fpalpha_prev, delta, alpha_next, alpha, alpha_prev
            real(dp), volatile :: a, b, fa, fb, fpa, fpb
            real(dp), volatile :: lower, upper
            logical            :: fpb_nan
            integer, parameter :: bracket_iters = 100, section_iters = 100
            integer :: i
            alpha = alpha1
            alpha_prev = 0.0_8
            i = 0
            call  wrap_fdf(0.0_8, f0, fp0)
            falpha_prev = f0
            fpalpha_prev = fp0
            ! Avoid uninitialized variables morning
            a = 0.0_8
            b = alpha
            fa = f0
            fb = 0.0_8
            fpa = fp0
            fpb = 0.0_8
            fpb_nan = .false.
            ! Begin bracketing
            do while (i < bracket_iters-1)
                i = i + 1
                falpha = wrap_f(alpha)
                ! Fletcher's rho test
                if ((falpha > f0 + alpha * rho * fp0).or.(falpha >= falpha_prev)) then
                    a = alpha_prev
                    fa = falpha_prev
                    fpa = fpalpha_prev
                    b = alpha
                    fb = falpha
                    fpb_nan = .true.
                    exit ! goto sectioning
                end if
                fpalpha = wrap_df(alpha)
                ! Fletcher's sigma test
                if (abs(fpalpha) <= -sigma * fp0) then
                    alpha_new = alpha
                    status = OPT_STATUS_SUCCESS
                    return
                end if
                if (fpalpha >= 0.0_8) then
                    a = alpha
                    fa = falpha
                    fpa = fpalpha
                    b = alpha_prev
                    fb = falpha_prev
                    fpb = fpalpha_prev
                    exit ! goto sectioning
                end if
                delta = alpha - alpha_prev
                lower = alpha + delta
                upper = alpha + tau1 * delta
                alpha_next = interpolate(alpha_prev, falpha_prev, fpalpha_prev, &
                    & alpha, falpha, fpalpha, .true., lower, upper, order)
                alpha_prev = alpha
                falpha_prev = falpha
                fpalpha_prev = fpalpha
                alpha = alpha_next
            end do
            ! Sectioning of bracket [a,b]
            do while (i < section_iters-1)
                i = i + 1
                delta = b - a
                lower = a + tau2 * delta
                upper = b - tau3 * delta
                alpha = interpolate(a, fa, fpa, b, fb, fpb, fpb_nan, lower, upper, order)
                falpha = wrap_f(alpha)
                if ((a-alpha)*fpa <= EPSILON(rho)) then
                    ! roundoff prevents progress
                    write (*,*) 'simple_bfgs2_opt: roundoff prevents progress'
                    status = OPT_STATUS_ERROR
                    return
                end if
                if ((falpha > f0 + rho * alpha * fp0).or.(falpha >= fa)) then
                    ! a_next = a
                    b = alpha
                    fb = falpha
                    fpb_nan = .true.
                else
                    fpalpha = wrap_df(alpha)
                    if (abs(fpalpha) <= -sigma * fp0) then
                        alpha_new = alpha
                        status = OPT_STATUS_SUCCESS
                        return ! terminate
                    end if
                    if ( ((b-a >= 0.0_8) .and. (fpalpha >= 0.0_8)) &
                        & .or. ((b-a <= 0.0_8) .and. (fpalpha <= 0.0_8))) then
                        b = a
                        fb = fa
                        fpb = fpa
                        a = alpha
                        fa = falpha
                        fpa = fpalpha
                    else
                        a = alpha
                        fa = falpha
                        fpa = fpalpha
                    end if
                end if
            end do
            status = OPT_STATUS_SUCCESS
        end function linear_minimize

        function wrap_f(alpha) result(res)
            real(dp), intent(in) :: alpha
            real(dp) :: res
            if (alpha == self%wrapper%f_cache_key) then ! using previously cached f(alpha)
                res = self%wrapper%f_alpha
                return
            end if
            call moveto(alpha)
            self%wrapper%f_alpha     = spec%eval_f(fun_self, self%wrapper%x_alpha)
            self%wrapper%f_cache_key = alpha
            res = self%wrapper%f_alpha
        end function wrap_f

        function wrap_df(alpha) result(res)
            real(dp), intent(in) :: alpha
            real(dp) :: res
            if (alpha == self%wrapper%df_cache_key) then ! using previously cached df(alpha)
                res = self%wrapper%df_alpha
                return
            end if
            call moveto(alpha)
            if (alpha .ne. self%wrapper%g_cache_key) then
                call spec%eval_df(fun_self, self%wrapper%x_alpha, self%wrapper%g_alpha)
                self%wrapper%g_cache_key = alpha
            end if
            self%wrapper%df_alpha = slope()
            self%wrapper%df_cache_key = alpha
            res = self%wrapper%df_alpha
        end function wrap_df

        subroutine wrap_fdf(alpha, f, df)
            real(dp), intent(in) :: alpha
            real(dp), intent(out) :: f, df
            ! Check for previously cached values
            if ((alpha == self%wrapper%f_cache_key).and.(alpha == self%wrapper%df_cache_key)) then
                f  = self%wrapper%f_alpha
                df = self%wrapper%df_alpha
                return
            end if
            if ((alpha == self%wrapper%f_cache_key).or.(alpha == self%wrapper%df_cache_key)) then
                f  = wrap_f(alpha)
                df = wrap_df(alpha)
                return
            end if
            call moveto(alpha)
            call spec%eval_fdf(fun_self, self%wrapper%x_alpha, self%wrapper%f_alpha, self%wrapper%g_alpha)
            self%wrapper%f_cache_key  = alpha
            self%wrapper%g_cache_key  = alpha
            self%wrapper%df_alpha     = slope()
            self%wrapper%df_cache_key = alpha
            f  = self%wrapper%f_alpha
            df = self%wrapper%df_alpha
        end subroutine wrap_fdf

        subroutine update_position(alpha, x, f, g)
            real(dp), intent(in) :: alpha
            real(dp), intent(out) :: x(:), f, g(:)
            real(dp) :: f_alpha, df_alpha
            ! ensure that everything is fully cached
            call wrap_fdf(alpha, f_alpha, df_alpha)
            f = self%wrapper%f_alpha
            x = self%wrapper%x_alpha
            g = self%wrapper%g_alpha
        end subroutine update_position

        subroutine moveto(alpha)
            real(dp), intent(in) :: alpha
            if (alpha == self%wrapper%x_cache_key) then ! using previously cached position
                return
            end if
            ! set x_alpha = x + alpha * p
            self%wrapper%x_alpha = self%x0
            self%wrapper%x_alpha = self%wrapper%x_alpha + alpha*self%p
            self%wrapper%x_cache_key = alpha
        end subroutine moveto

        function interpolate(a, fa, fpa, b, fb, fpb, fpb_nan, xmin, xmax, order) result(alpha)
            !use, intrinsic :: IEEE_ARITHMETIC, only: IEEE_IS_FINITE
            real(dp), intent(in) :: a, fa, fpa, b, fb, fpb, xmin, xmax
            integer, intent(in) :: order
            logical, intent(in) :: fpb_nan
            real(dp) :: z, alpha, zmin, zmax, tmp
            ! Map [a,b] to [0,1]
            zmin = (xmin - a) / (b - a)
            zmax = (xmax - a) / (b - a)
            if (zmin > zmax) then
                tmp = zmin
                zmin = zmax
                zmax = tmp
            end if
            if ((order > 2) .and. (.not. fpb_nan)) then
                z = interp_cubic(fa, fpa * (b - a), fb, fpb * (b - a), zmin, zmax)
            else
                z = interp_quad(fa, fpa * (b - a), fb, zmin, zmax)
            end if
            alpha = a + z * (b - a)
        end function interpolate

        function cubic(c0, c1, c2, c3, z) result(res)
            ! Find a minimum in x=[0,1] of the interpolating cubic through
            ! (0,f0) (1,f1) with derivatives fp0 at x=0 and fp1 at x=1.
            !
            ! The interpolating polynomial is:
            !
            ! c(x) = f0 + fp0 * z + eta * z^2 + xi * z^3
            !
            ! where eta=3*(f1-f0)-2*fp0-fp1, xi=fp0+fp1-2*(f1-f0).
            real(dp), intent(in) :: c0, c1, c2, c3, z
            real(dp) :: res
            res = c0 + z * (c1 + z * (c2 + z * c3))
        end function cubic

        subroutine check_extremum(c0, c1, c2, c3, z, zmin, fmin)
            ! could make an early return by testing curvature >0 for minimum
            real(dp), intent(in) :: c0, c1, c2, c3, z
            real(dp), intent(inout) :: fmin
            real(dp), intent(out) :: zmin
            real(dp) :: y
            y = cubic(c0, c1, c2, c3, z)
            if (y < fmin) then
                zmin = z  ! accepted new point
                fmin = y
            end if
        end subroutine check_extremum

        function interp_cubic(f0, fp0, f1, fp1, zl, zh) result(zmin)
            real(dp), intent(in) :: f0, fp0, f1, fp1, zl, zh
            real(dp), volatile :: zmin
            real(dp), volatile :: eta, xi, c0, c1, c2, c3, fmin, z0, z1
            integer :: n
            eta = 3.0_8 * (f1 - f0) - 2.0_8 * fp0 - fp1;
            xi = fp0 + fp1 - 2.0_8 * (f1 - f0);
            c0 = f0
            c1 = fp0
            c2 = eta
            c3 = xi
            zmin = zl
            fmin = cubic(c0, c1, c2, c3, zl)
            call check_extremum(c0, c1, c2, c3, zh, zmin, fmin)
            n = solve_quadratic(3.0_8 * c3, 2.0_8 * c2, c1, z0, z1)
            if (n == 2) then  ! found 2 roots
                if ((z0 > zl) .and. (z0 < zh)) then
                    call check_extremum(c0, c1, c2, c3, z0, zmin, fmin)
                end if
                if ((z1 > zl) .and. (z1 < zh)) then
                    call check_extremum(c0, c1, c2, c3, z1, zmin, fmin)
                end if
            else if (n == 1) then  ! found 1 root
                if ((z0 > zl) .and. (z0 < zh)) then
                    call check_extremum(c0, c1, c2, c3, z0, zmin, fmin);
                end if
            end if
        end function interp_cubic

        function interp_quad(f0, fp0, f1, zl, zh) result(zmin)
            ! Find a minimum in x=[0,1] of the interpolating quadratic through
            ! (0,f0) (1,f1) with derivative fp0 at x=0.  The interpolating
            ! polynomial is q(x) = f0 + fp0 * z + (f1-f0-fp0) * z^2
            real(dp), intent(in) :: f0, fp0, f1, zl, zh
            real(dp) :: fl, fh, c, zmin, fmin, z, f
            fl = f0 + zl*(fp0 + zl*(f1 - f0 -fp0))
            fh = f0 + zh*(fp0 + zh*(f1 - f0 -fp0))
            c = 2.0_8 * (f1 - f0 - fp0) ! curvature
            zmin = zl
            fmin = fl
            if (fh < fmin) then
                zmin = zh
                fmin = fh
            end if
            if (c > 0.0_8) then ! positive curvature required for a minimum
                z = -fp0 / c      ! location of minimum
                if ((z > zl) .and. (z < zh))  then
                    f = f0 + z*(fp0 + z*(f1 - f0 -fp0))
                    if (f < fmin) then
                        zmin = z
                        fmin = f
                    end if
                end if
            end if
        end function interp_quad

        function solve_quadratic(a, b, c, x0, x1) result(nsol)
            real(dp), intent(in) :: a, b, c
            real(dp), intent(out) :: x0, x1
            integer :: nsol  ! number of solutions
            real(dp) :: disc, r, sgnb, temp, r1, r2
            if (a == 0.0_8) then ! Handle linear case
                if (b == 0.0_8) then
                    nsol = 0
                    return
                else
                    x0 = -c / b
                    nsol = 1
                    return
                end if
            end if
            disc = b * b - 4.0_8 * a * c;
            if (disc > 0.0_8) then
                if (b == 0.0_8) then
                    r = sqrt(-c / a)
                    x0 = -r
                    x1 =  r
                else
                    if (b > 0.0_8) then
                        sgnb = 1.0_8
                    else
                        sgnb = -1.0_8
                    end if
                    temp = -0.5_8 * (b + sgnb * sqrt(disc))
                    r1 = temp / a
                    r2 = c / temp
                    if (r1 < r2) then
                        x0 = r1
                        x1 = r2
                    else
                        x0 = r2
                        x1 = r1
                    end if
                end if
                nsol = 2
                return
            else if (disc == 0.0_8) then
                x0 = -0.5_8 * b / a
                x1 = -0.5_8 * b / a
                nsol = 2
                return
            else
                nsol = 0
                return
            end if
        end function solve_quadratic

        subroutine change_direction
            ! Convert the cache values from the end of the current minimisation
            ! to those needed for the start of the next minimisation, alpha=0
            ! The new x_alpha for alpha=0 is the current position
            self%wrapper%x_alpha = spec%x_8
            self%wrapper%x_cache_key = 0.0_8
            ! The function value does not change
            self%wrapper%f_cache_key = 0.0_8
            ! The new g_alpha for alpha=0 is the current gradient at the endpoint
            self%wrapper%g_alpha = self%gradient
            self%wrapper%g_cache_key = 0.0_8
            ! Calculate the slope along the new direction vector, p
            self%wrapper%df_alpha = slope()
            self%wrapper%df_cache_key = 0.0_8
        end subroutine change_direction

    end subroutine bfgs2_minimize

    !> \brief  is a destructor
    subroutine kill_bfgs2_opt( self )
        class(bfgs2_opt), intent(inout) :: self !< instance
        if( self%exists )then
            deallocate(self%x0, self%g0, self%p, self%gradient, &
                & self%dx0, self%dg0, self%dx, self%wrapper%x_alpha, self%wrapper%g_alpha)
            self%exists = .false.
        endif
    end subroutine kill_bfgs2_opt

end module simple_bfgs2_opt
