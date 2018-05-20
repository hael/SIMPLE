! function minimization by Fletcher-Reeves conjugate gradient algorithm, translated from gsl 2.4
module simple_opt_fr_cg
include 'simple_lib.f08'
use simple_optimizer, only: optimizer
use simple_opt_spec, only: opt_spec
use simple_opt_helpers
implicit none

public :: opt_fr_cg
private

type, extends(optimizer) :: opt_fr_cg
    private
    real(kind=8), allocatable :: x1(:), x2(:), p(:), g0(:), gradient(:)
    real(kind=8) :: step, pnorm, g0norm, f
    integer :: fr_cg_iter       !< internal iter for fr_cg algorithm
    logical :: exists=.false.
contains
    procedure :: new          => new_opt_fr_cg
    procedure :: minimize     => fr_cg_minimize
    procedure :: kill         => kill_opt_fr_cg
end type

contains

    !> \brief  is a constructor
    subroutine new_opt_fr_cg( self, spec )
        use simple_opt_spec, only: opt_spec
        class(opt_fr_cg), intent(inout) :: self        !< instance
        class(opt_spec),  intent(inout) :: spec        !< specification
        call self%kill
        allocate(self%x1(spec%ndim),self%x2(spec%ndim),self%p(spec%ndim),self%g0(spec%ndim),&
            & self%gradient(spec%ndim),stat=alloc_stat)
        if(alloc_stat/=0)call allocchk('In: new_fr_cg_opt; simple_opt_fr_cg')
        self%exists = .true.
    end subroutine

    !>  \brief  nonlinear conjugate gradient minimizer
    subroutine fr_cg_minimize( self, spec, fun_self, lowest_cost )
        use simple_opt_subs, only: lnsrch
        class(opt_fr_cg), intent(inout) :: self        !< instance
        class(opt_spec),  intent(inout) :: spec        !< specification
        class(*),         intent(inout) :: fun_self    !< self-pointer for cost function
        real, intent(out)               :: lowest_cost !< minimum function value
        if ( (.not. associated(spec%costfun) ).and.( .not. associated(spec%costfun_8)) ) then
            stop 'cost function not associated in opt_spec; fr_cg_minimize; simple_opt_fr_cg'
        endif
        if ( (.not. associated(spec%gcostfun) ).and.(.not. associated(spec%gcostfun_8)) ) then
            stop 'gradient of cost function not associated in opt_spec; fr_cg_minimize; simple_opt_fr_cg'
        endif
        spec%x_8 = spec%x
        ! initialise nevals counters
        spec%nevals  = 0
        spec%ngevals = 0
        call fr_cgmin

        contains

            !>  \brief  nonlinear conjugate gradient minimizer
            subroutine fr_cgmin
                integer :: status
                integer :: iter
                iter        = 0
                self%step   = spec%max_step
                self%x1     = 0.0_8
                self%x2     = 0.0_8
                self%p      = 0.0_8
                self%g0     = 0.0_8
                call fr_cg_set
                do
                    iter = iter + 1
                    status = fr_cg_iterate()
                    if (status == OPT_STATUS_ERROR) then
                        write (*,*) 'simple_opt_fr_cg: error in minimizer routine'
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
            end subroutine fr_cgmin

            subroutine fr_cg_set
                real(kind=8) :: gnorm
                self%fr_cg_iter = 0
                self%step = spec%max_step
                ! Use the gradient as the initial direction
                call spec%eval_fdf(fun_self, spec%x_8, self%f, self%gradient)
                self%p       = self%gradient
                self%g0      = self%gradient
                gnorm        = norm_2(self%gradient)
                self%pnorm   = gnorm
                self%g0norm  = gnorm
            end subroutine fr_cg_set

            function fr_cg_iterate() result(status)
                integer :: status
                real(kind=8) :: fa, fb, fc, dir, stepa, stepb, stepc, g1norm, pg, beta
                fa = self%f
                stepa = 0.0_8
                stepc = self%step
                if ((self%pnorm == 0.0_8).or.(self%g0norm == 0.0_8)) then
                    status = OPT_STATUS_ERROR
                    return
                end if
                ! Determine which direction is downhill, +p or -p
                pg = dot_product(self%p, self%gradient)
                if (pg >= 0.0_8) then
                    dir = 1.0_8
                else
                    dir = -1.0_8
                end if
                call take_step(spec%x_8, self%p, stepc, dir / self%pnorm, self%x1)
                fc = spec%eval_f(fun_self,self%x1)
                if (fc < fa) then
                    ! Success, reduced the function value
                    self%step = stepc * 2.0_8
                    self%f    = fc
                    spec%x_8  = self%x1
                    spec%x    = real(spec%x_8, kind=4)
                    call spec%eval_df(fun_self, self%x1, self%gradient)
                    status = OPT_STATUS_CONTINUE
                    return
                end if
                if (global_debug .and. global_verbose) then
                    write (*,*) 'got stepc = ', stepc, 'fc = ', fc
                end if
                ! Do a line minimisation in the region (xa,fa) (xc,fc) to find an
                ! intermediate (xb,fb) satisifying fa > fb < fc.  Choose an initial
                ! xb based on parabolic interpolation
                call intermediate_point(spec%x_8, self%p, dir / self%pnorm, pg, &
                    & stepa, stepc, fa, fc, self%x1, self%gradient, stepb, fb, spec, fun_self)
                if (stepb == 0.0_8) then
                    status = OPT_STATUS_ERROR
                    return
                end if
                call minimize(spec%x_8, self%p, dir / self%pnorm, &
                    & stepa, stepb, stepc, fa, fb, fc, &
                    & self%x1, self%x2, self%gradient, self%step, self%f, g1norm, spec, fun_self)
                spec%x_8 = self%x2
                spec%x   = real(spec%x_8, kind=4)
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
                    self%pnorm = norm_2(self%p)
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
    subroutine kill_opt_fr_cg( self )
        class(opt_fr_cg), intent(inout) :: self !< instance
        if( self%exists )then
            deallocate(self%x1, self%x2, self%p, self%g0, self%gradient)
            self%exists = .false.
        endif
    end subroutine

end module simple_opt_fr_cg
