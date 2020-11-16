! helper routines and definitions for certain optimizer routines (those taken from the gsl library)
module simple_opt_helpers
include 'simple_lib.f08'
implicit none

!< constants for return values
integer, parameter  :: OPT_STATUS_SUCCESS = 1, OPT_STATUS_CONTINUE = 0, OPT_STATUS_ERROR = -1

contains

    subroutine take_step (x, p, step, lambda, x1)
        real(dp), intent(in) :: x(:), p(:), step, lambda
        real(dp), intent(out) :: x1(:)
        x1 = x - step*lambda*p
    end subroutine take_step

    subroutine intermediate_point(x, p, lambda,  pg, stepa, stepc, fa, fc, &
        & x1, gradient, step, f, spec, fun_self)
        use simple_opt_spec, only: opt_spec
        real(dp), intent(in)       :: x(:), p(:), lambda, pg, stepa,  fa
        real(dp), intent(out)      :: x1(:), gradient(:), step, f, fc, stepc
        class(opt_spec), intent(inout) :: spec
        class(*),        intent(inout) :: fun_self
        real(dp) :: stepb, fb, u
10      u = abs(pg * lambda * stepc)  ! label: trial
        stepb = 0.5_8 * stepc * u / ((fc - fa) + u)
        call take_step(x, p , stepb, lambda, x1)
        if (vect_equal(x, x1)) then
            ! Take fast exit if trial point does not move from initial point
            if (global_debug .and. global_verbose) then
                write (*,*) 'fast exit x == x1 for stepb = ', stepb
            end if
            step = 0
            f = fa
            call spec%eval_df(fun_self, x1, gradient)
            return
        end if
        fb = spec%eval_f(fun_self,x1)
        if (global_debug .and. global_verbose) then
            write (*,*) 'trying stepb = ', stepb, 'fb = ', fb
        end if
        if ((fb >= fa) .and. (is_gt_zero(stepb))) then
            ! downhill step failed, reduce step-size and try again
            fc = fb
            stepc = stepb
            go to 10  ! goto trial
        end if
        if (global_debug .and. global_verbose) then
            write (*,*) 'ok!'
        end if
        step = stepb
        f = fb
        call spec%eval_df(fun_self, x1, gradient)
    end subroutine intermediate_point

    subroutine minimize(x, p, lambda, stepa, stepb, stepc, fa, fb, fc, &
        & x1, x2, gradient, step, f, gnorm, spec, fun_self)
        use simple_opt_spec, only: opt_spec
        real(dp), intent(in) :: x(:), p(:), lambda
        real(dp), intent(out) :: x1(:), x2(:), gradient(:), step, stepa, stepb, stepc, fa, fb, fc, f, gnorm
        class(opt_spec), intent(inout) :: spec
        class(*),        intent(inout) :: fun_self
        real(dp) :: u, v, w, dw, dv, du, fu, fv, fw, old1, old2, stepm, fm, pg, gnorm1, iter, e1, e2
        u = stepb
        v = stepa
        w = stepc
        fu = fb
        fv = fa
        fw = fc
        old2 = abs(w - v)
        old1 = abs(v - u)
        iter = 0
        x2 = x1
        f = fb
        step = stepb
        gnorm = norm_2(gradient)
10      iter = iter + 1 ! label: mid_trial
        if (iter > 10) then
            return ! MAX ITERATIONS
        end if
        dw = w - u
        dv = v - u
        du = 0.0_8
        e1 = ((fv - fu) * dw * dw + (fu - fw) * dv * dv)
        e2 = 2.0_8 * ((fv - fu) * dw + (fu - fw) * dv)
        if (.not. is_zero(e2)) then
            du = e1 / e2
        end if
        if ((du > 0.0_8) .and. (du < (stepc - stepb)) .and. (abs(du) < 0.5_8 * old2)) then
            stepm = u + du
        else if ((du < 0.0_8) .and. (du > (stepa - stepb)) .and. (abs(du) < 0.5_8 * old2)) then
            stepm = u + du
        else if ((stepc - stepb) > (stepb - stepa)) then
            stepm = 0.38_8 * (stepc - stepb) + stepb
        else
            stepm = stepb - 0.38_8 * (stepb - stepa)
        end if
        call take_step(x, p, stepm, lambda, x1)
        fm = spec%eval_f(fun_self,x1)
        if (global_debug .and. global_verbose) then
            write (*,*) 'trying stepm = ', stepm, ' fm = ', fm
        end if
        if (fm > fb) then
            if (fm < fv) then
                w = v
                v = stepm
                fw = fv
                fv = fm
            else if (fm < fw) then
                w = stepm
                fw = fm
            end if
            if (stepm < stepb) then
                stepa = stepm
                fa = fm
            else
                stepc = stepm
                fc = fm
            end if
            go to 10    ! goto mid_trial
        else if (fm <= fb) then
            old2 = old1
            old1 = abs(u - stepm)
            w = v
            v = u
            u = stepm
            fw = fv
            fv = fu
            fu = fm
            x2 = x1
            call spec%eval_df(fun_self, x1, gradient)
            pg = dot_product(p, gradient)
            gnorm1 = norm_2(gradient)
            if (global_debug .and. global_verbose) then
                write (*,*) 'p: ', p
                write (*,*) 'g: ', gradient
                write (*,*) 'gnorm: ', gnorm1
                write (*,*) 'pg: ', pg
                write (*,*) 'orth: ', abs(pg * lambda / gnorm1)
            end if
            f = fm
            step = stepm
            gnorm = gnorm1
            if (abs(pg * lambda / gnorm1) < spec%ftol) then
                if (global_debug .and. global_verbose) then
                    write (*,*) 'ok!'
                end if
                return  ! SUCCESS
            end if
            if (stepm < stepb) then
                stepc = stepb
                fc = fb
                stepb = stepm
                fb = fm
            else
                stepa = stepb
                fa = fb
                stepb = stepm
                fb = fm
            end if
            go to 10    ! goto mid_trial
        end if
    end subroutine minimize

    function test_gradient(g, epsabs) result(status)
        real(dp), intent(in) :: g(:), epsabs
        integer                  :: status
        if (norm_2(g) < epsabs) then
            status = OPT_STATUS_SUCCESS
        else
            status = OPT_STATUS_CONTINUE
        end if
    end function test_gradient

    logical function vect_equal_4( array1, array2 )
        real, intent(in) :: array1(:), array2(:)
        integer :: i
        vect_equal_4 = size(array1) == size(array2)
        if ( vect_equal_4 ) then
            do i = 1,size(array1)
                vect_equal_4 = is_equal(array1(i), array2(i))
                if ( .not. vect_equal_4 )exit
            enddo
        endif
    end function vect_equal_4

    logical function vect_equal( array1, array2 )
        real(dp), intent(in) :: array1(:), array2(:)
        integer :: i
        vect_equal = size(array1) == size(array2)
        if ( vect_equal ) then
            do i = 1,size(array1)
                vect_equal = is_equal(array1(i), array2(i))
                if ( .not. vect_equal )exit
            enddo
        endif
    end function vect_equal

end module simple_opt_helpers
