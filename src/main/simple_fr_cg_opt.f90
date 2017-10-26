! function minimization by Fletcher-Reeves conjugate gradient algorithm, translated from gsl 2.4

module simple_fr_cg_opt
#include "simple_lib.f08"

use simple_optimizer, only: optimizer
use simple_opt_helpers
implicit none

public :: fr_cg_opt
private

type, extends(optimizer) :: fr_cg_opt
    private
    real, allocatable :: x1(:), x2(:), p(:), g0(:), gradient(:)
    real :: step, pnorm, g0norm, f
    integer :: fr_cg_iter       !< internal iter for fr_cg algorithm
    logical :: exists=.false.
contains
    procedure :: new          => new_fr_cg_opt
    procedure :: minimize     => fr_cg_minimize
    procedure :: kill         => kill_fr_cg_opt
end type

contains

    !> \brief  is a constructor
    subroutine new_fr_cg_opt( self, spec )
        use simple_opt_spec, only: opt_spec
        use simple_syslib,   only: alloc_errchk
        class(fr_cg_opt), intent(inout) :: self !< instance
        class(opt_spec), intent(inout)  :: spec !< specification
        call self%kill
        allocate(self%x1(spec%ndim),self%x2(spec%ndim),self%p(spec%ndim),self%g0(spec%ndim),&
            & self%gradient(spec%ndim),stat=alloc_stat)
        allocchk('In: new_fr_cg_opt; simple_fr_cg_opt')
        self%exists = .true.
    end subroutine

    !>  \brief  nonlinear conjugate gradient minimizer
    subroutine fr_cg_minimize( self, spec, lowest_cost )
        use simple_opt_spec, only: opt_spec
        use simple_opt_subs, only: lnsrch
        class(fr_cg_opt), intent(inout) :: self        !< instance
        class(opt_spec), intent(inout)  :: spec        !< specification
        real, intent(out)               :: lowest_cost !< minimum function value
        integer                         :: avgniter,i
        if( .not. associated(spec%costfun) )then
            stop 'cost function not associated in opt_spec; fr_cg_minimize; simple_fr_cg_opt'
        endif
        if( .not. associated(spec%gcostfun) )then
            stop 'gradient of cost function not associated in opt_spec; fr_cg_minimize; simple_fr_cg_opt'
        endif
        ! run nrestarts restarts
        avgniter = 0
        spec%nevals = 0
        do i=1,spec%nrestarts
            call fr_cgmin
            avgniter = avgniter+spec%niter
        end do
        spec%niter  = avgniter/spec%nrestarts
        spec%nevals = spec%nevals/spec%nrestarts
        !spec%x      = self%p

        contains

            !>  \brief  nonlinear conjugate gradient minimizer
            subroutine fr_cgmin
                integer :: status
                integer :: iter
                iter        = 0
                self%step   = spec%max_step
                self%x1     = 0.
                self%x2     = 0.
                self%p      = 0.
                self%g0     = 0.
                call fr_cg_set
                do
                    iter = iter + 1
                    status = fr_cg_iterate()
                    if (status == OPT_STATUS_ERROR) then
                        write (*,*) 'simple_fr_cg_opt: error in minimizer routine'
                        return
                    end if
                    status = test_gradient(self%gradient, spec%gtol)                    
                    if ((global_debug).and.(global_verbose)) then
                        if (status == OPT_STATUS_SUCCESS) then
                            write (*,*) 'Minimum found at:'
                        end if
                        write (*,*) iter, 'x = ', spec%x, 'f = ', self%f
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
                lowest_cost = self%f
            end subroutine fr_cgmin

            subroutine fr_cg_set
                real :: gnorm
                self%fr_cg_iter = 0
                self%step = spec%max_step
                ! Use the gradient as the initial direction
                call spec%eval_fdf(spec%x, self%f, self%gradient)
                spec%nevals  = spec%nevals  + 1
                spec%ngevals = spec%ngevals + 1
                self%p       = self%gradient
                self%g0      = self%gradient
                gnorm        = norm2(self%gradient)
                self%pnorm   = gnorm
                self%g0norm  = gnorm
            end subroutine fr_cg_set

            function fr_cg_iterate() result(status)
                integer :: status
                real :: fa, fb, fc, dir, stepa, stepb, stepc, g1norm, pg, beta
                fa = self%f
                stepa = 0.
                stepc = self%step
                if ((self%pnorm == 0.).or.(self%g0norm == 0.)) then
                    status = OPT_STATUS_ERROR
                    return
                end if
                ! Determine which direction is downhill, +p or -p
                pg = dot_product(self%p, self%gradient)
                if (pg >= 0.) then
                    dir = 1.
                else
                    dir = -1.
                end if
                call take_step(spec%x, self%p, stepc, dir / self%pnorm, self%x1)
                fc = spec%costfun(self%x1,spec%ndim)
                spec%nevals = spec%nevals  + 1                
                if (fc < fa) then
                    ! Success, reduced the function value
                    self%step = stepc * 2.0
                    self%f = fc
                    spec%x = self%x1
                    self%gradient = spec%gcostfun(self%x1,spec%ndim)
                    spec%ngevals = spec%ngevals + 1
                    status = OPT_STATUS_CONTINUE
                    return
                end if
                if (global_debug .and. global_verbose) then
                    write (*,*) 'got stepc = ', stepc, 'fc = ', fc
                end if
                ! Do a line minimisation in the region (xa,fa) (xc,fc) to find an
                ! intermediate (xb,fb) satisifying fa > fb < fc.  Choose an initial
                ! xb based on parabolic interpolation
                call intermediate_point (spec%x, self%p, dir / self%pnorm, pg, &
                    & stepa, stepc, fa, fc, self%x1, self%gradient, stepb, fb, spec)
                if (stepb == 0.0) then
                    status = OPT_STATUS_ERROR
                    return
                end if
                call minimize (spec%x, self%p, dir / self%pnorm, &
                    & stepa, stepb, stepc, fa, fb, fc, &
                    & self%x1, self%x2, self%gradient, self%step, self%f, g1norm, spec)
                spec%x = self%x2
                ! Choose a new conjugate direction for the next step
                self%fr_cg_iter = modulo(self%fr_cg_iter + 1, spec%ndim)
                if (self%fr_cg_iter == 0) then
                    self%p = self%gradient
                    self%pnorm = g1norm
                else
                    ! p' = g1 - beta * p
                    beta = - (g1norm / self%g0norm)**2
                    self%p = -beta * self%p
                    self%p = self%p + self%gradient
                    self%pnorm = norm2(self%p);
                end if
                self%g0norm = g1norm
                self%g0 = self%gradient
                if (global_debug .and. global_verbose) then
                    write (*,*) 'updated conjugate directions'
                    write (*,*) 'p: ', self%p
                    write (*,*) 'g: ', self%gradient
                end if
                status = OPT_STATUS_CONTINUE
            end function fr_cg_iterate            
    end subroutine

    !> \brief  is a destructor
    subroutine kill_fr_cg_opt( self )
        class(fr_cg_opt), intent(inout) :: self !< instance
        if( self%exists )then
            deallocate(self%x1, self%x2, self%p, self%g0, self%gradient)
            self%exists = .false.
        endif
    end subroutine

end module simple_fr_cg_opt
