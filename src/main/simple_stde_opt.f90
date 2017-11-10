! function minimization by steepest descent algorithm, translated from gsl 2.4 (steepest_descent.c)

module simple_stde_opt
#include "simple_lib.f08"

use simple_optimizer, only: optimizer
use simple_opt_helpers
implicit none

public :: stde_opt
private

type, extends(optimizer) :: stde_opt
    private
    real(dp) :: step, f
    real(dp), allocatable :: x1(:), g1(:), gradient(:)
    logical :: exists=.false.
contains
    procedure :: new          => new_stde_opt
    procedure :: minimize     => stde_minimize
    procedure :: kill         => kill_stde_opt
end type

contains

    !> \brief  is a constructor
    subroutine new_stde_opt( self, spec )
        use simple_opt_spec, only: opt_spec
        use simple_syslib,   only: alloc_errchk
        class(stde_opt), intent(inout)  :: self !< instance
        class(opt_spec), intent(inout)  :: spec !< specification
        call self%kill
        allocate(self%x1(spec%ndim),self%g1(spec%ndim), self%gradient(spec%ndim), &
            & stat=alloc_stat)
        allocchk('In: new_stde_opt; simple_stde_opt')
        self%exists = .true.
    end subroutine new_stde_opt

    !>  \brief  nonlinear conjugate gradient minimizer
    subroutine stde_minimize( self, spec, fun_self, lowest_cost )
        use simple_opt_spec, only: opt_spec
        use simple_opt_subs, only: lnsrch
        class(stde_opt), intent(inout)  :: self        !< instance
        class(opt_spec), intent(inout)  :: spec        !< specification
        class(*),        intent(inout)  :: fun_self    !< self-pointer for cost function
        real, intent(out)               :: lowest_cost !< minimum function value
        if( .not. associated(spec%costfun) )then
            stop 'cost function not associated in opt_spec; stde_minimize; simple_stde_opt'
        endif
        if( .not. associated(spec%gcostfun) )then
            stop 'gradient of cost function not associated in opt_spec; stde_minimize; simple_stde_opt'
        endif
        spec%x_8     = spec%x
        ! initialise nevals counters
        spec%nevals  = 0
        spec%ngevals = 0
        call stdemin

    contains

        !>  \brief  nonlinear conjugate gradient minimizer
        subroutine stdemin
            integer :: status
            integer :: iter
            iter        = 0
            self%step   = spec%max_step
            self%x1     = 0.0_8
            self%g1     = 0.0_8
            call stde_set
            do
                iter = iter + 1
                status = stde_iterate()
                if (status == OPT_STATUS_ERROR) then
                    write (*,*) 'simple_stde_opt: error in minimizer routine'
                    return
                end if
                status = test_gradient(self%gradient, real(spec%gtol, kind=8))
                if ((global_debug).and.(global_verbose)) then
                    if (status == OPT_STATUS_SUCCESS) then
                        write (*,*) 'Minimum found at:'
                    end if
                    write (*,*) iter, 'x = ', spec%x_8, 'f = ', self%f
                end if
                if ((status .ne. OPT_STATUS_CONTINUE) .or. (iter > spec%maxits)) then
                    exit
                end if
                if (associated(spec%opt_callback)) then
                    call spec%opt_callback(fun_self)
                end if
            end do
            if (status == OPT_STATUS_SUCCESS) then
                spec%converged = .true.
            else
                spec%converged = .false.
            end if
            lowest_cost = real(self%f, kind=4)
        end subroutine stdemin

        subroutine stde_set
            self%step    = spec%max_step
            call spec%eval_fdf_8(fun_self, spec%x_8, self%f, self%gradient)
        end subroutine stde_set

        function stde_iterate() result(status)
            integer      :: status
            real(dp) :: f0, f1, gnorm
            logical      :: failed
            f0 = self%f
            failed = .false.
            ! compute new trial point at x1= x - step * dir, where dir is the
            ! normalized gradient
            gnorm = norm2(self%gradient)
            if (gnorm == 0.0_8) then
                status = OPT_STATUS_ERROR
                return
            end if
10          self%x1 = spec%x_8 - self%step / gnorm * self%gradient ! label: trial
            if (vect_equal(spec%x_8, self%x1)) then
                status = OPT_STATUS_ERROR
                return
            end if
            ! evaluate function and gradient at new point x1
            call spec%eval_fdf(fun_self, self%x1, f1, self%g1)
            if (f1 > f0) then
                ! downhill step failed, reduce step-size and try again
                failed = .true.
                self%step = self%step * spec%ftol
                go to 10  ! goto label trial
            end if
            if (failed) then
                self%step = self%step * spec%ftol
            else
                self%step = self%step * 2.0_8
            end if
            !write (*,*) 'step = ', self%step
            spec%x_8      = self%x1
            spec%x        = real(spec%x_8, kind=4)
            self%gradient = self%g1
            self%f = f1
            status = OPT_STATUS_SUCCESS
        end function stde_iterate

    end subroutine stde_minimize

        !> \brief  is a destructor
    subroutine kill_stde_opt( self )
        class(stde_opt), intent(inout) :: self !< instance
        if( self%exists )then
            deallocate(self%x1,self%g1, self%gradient)
            self%exists = .false.
        endif
    end subroutine kill_stde_opt
end module simple_stde_opt
