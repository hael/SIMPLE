!==Class simple_powell_opt
!
! Minimization of an externally defined function by Powell's direction set method
! The code is distributed with the hope that it will be useful, but _WITHOUT_ _ANY_ _WARRANTY_. 
! Redistribution or modification is regulated by the GNU General Public License. 
! *Author:* Hans Elmlund, 2013-10-15
!
module simple_powell_opt
use simple_optimizer, only: optimizer
use simple_opt_spec,  only: opt_spec
use simple_jiffys,    only: alloc_err
implicit none

public :: powell_opt
private
  
type, extends(optimizer) :: powell_opt
    private
    type(opt_spec)    :: spec_linmin        !< specification of linmin optimizer
    real, allocatable :: direction_set(:,:) !< set of search directions
    real              :: yb=0.              !< best cost function value
    logical           :: exists=.false.     !< to indicate existence
  contains
    procedure :: new      => new_powell_opt
    procedure :: minimize => powell_minimize
    procedure :: kill     => kill_powell_opt
end type

logical :: warn=.false.

contains
    
    !> \brief  is a constructor
    subroutine new_powell_opt( self, spec )
        class(powell_opt), intent(inout) :: self !< instance
        class(opt_spec), intent(inout)   :: spec !< specification
        integer                          :: alloc_stat
        real                             :: x
        call self%kill
        allocate(self%direction_set(spec%ndim,spec%ndim), stat=alloc_stat)
        call alloc_err("In: new; simple_powell_opt", alloc_stat)
        self%direction_set = 0.
        self%yb            = huge(x) ! initialize best cost to huge number
        ! make line minimizer
        call self%spec_linmin%specify('linmin', spec%ndim, maxits=spec%maxits, ftol=spec%ftol)
        call self%spec_linmin%set_costfun(spec%costfun)
        self%exists = .true.
    end subroutine
    
    !>  \brief  the high-level minimization routine
    subroutine powell_minimize(self,spec,lowest_cost)
        use simple_rnd, only: ran3
        class(powell_opt), intent(inout) :: self        !< instance
        class(opt_spec), intent(inout)   :: spec        !< specification
        real, intent(out)                :: lowest_cost !< lowest cost
        logical :: found_better
        real    :: cost
        integer :: i, j
        logical :: arezero(spec%ndim)
        if( .not. associated(spec%costfun) )then 
            stop 'cost function not associated in opt_spec; powell_minimize; simple_powell_opt'
        endif
        ! initialize the direction set
        self%direction_set = 0.
        do i=1,spec%ndim
            self%direction_set(i,i) = 1.
        end do
        ! test if best point in spec is set
        arezero = .false.
        do i=1,spec%ndim
            if( spec%x(i) == 0. ) arezero(i) = .true.
        end do
        ! generate initial vector
        if( all(arezero) )then
            do i=1,spec%ndim
                ! initialize each variable by randomized bounds
                spec%x(i) = spec%limits(i,1)+ran3()*(spec%limits(i,2)-spec%limits(i,1))
            end do
        endif
        ! set best point in line minimizer to point in spec
        self%spec_linmin%x = spec%x
        ! set best cost
        spec%nevals  = 0
        self%yb      = spec%costfun(spec%x, spec%ndim)
        spec%nevals  = spec%nevals+1
        ! run nrestarts
        found_better = .false.
        do i=1,spec%nrestarts
            call powell(cost)
            if( cost <= self%yb )then
                self%yb = cost
                spec%x = self%spec_linmin%x
                found_better = .true.
            endif
            if( found_better )then
                ! re-initialize the direction set
                self%direction_set = 0.
                do j=1,spec%ndim
                    self%direction_set(j,j) = 1.
                end do
            else
                exit
            endif
        end do
        spec%nevals = spec%nevals/spec%nrestarts
        lowest_cost = self%yb
        
        contains
        
            !>  \brief  minimization of a function spec%costfun of spec%ndim variables. 
            !!          Input consists of an initial starting point self%spec_linmin%x 
            !!          that is a vector of length spec%ndim; an initial matrix 
            !!          self%direction_set whose logical dimensions are spec%ndim by 
            !!          spec%ndim, physical dimensions np by np, and whose columns
            !!          contain the initial set of directions (usually the n unit
            !!          vectors); and spec%ftol, the fractional tolerance in the func-
            !!          tion value such that failure to decrease by more than this
            !!          amount on one iteration signals doneness. On output, 
            !!          self%spec_linmin%x is
            !!          set to the best point found, xi is the then-current direc-
            !!          tion set,  cost is the returned function value at p,  and
            !!          iter is the number of iterations taken. the routine linmin
            !!          is used
            subroutine powell(cost)
                use simple_opt_subs, only: linmin
                real, intent(out) :: cost
                real, allocatable :: pt(:),ptt(:)
                integer :: i,ibig,j,iter,alloc_stat
                real    :: del,fp,fptt,t
                allocate(pt(spec%ndim),ptt(spec%ndim), stat=alloc_stat)
                call alloc_err("In: powell; simple_powell_opt", alloc_stat)
                cost=spec%costfun(self%spec_linmin%x,spec%ndim) ! set initial costfun val
                spec%nevals = spec%nevals+1
                do j=1,spec%ndim                
                    pt(j)=self%spec_linmin%x(j) ! save initial pont
                end do
                iter=0
                do
                    iter=iter+1
                    fp=cost
                    ibig=1
                    del=0.
                    do i=1,spec%ndim     ! in each iteration, loop over all directions in the set
                        do j=1,spec%ndim ! copy the direction
                            self%spec_linmin%xi(j)=self%direction_set(j,i)
                        end do
                        fptt=cost
                        self%spec_linmin%nevals = 0
                        call linmin(self%spec_linmin,cost) ! minimize along it
                        spec%nevals = spec%nevals+self%spec_linmin%nevals
                        if (abs(fptt-cost).gt.del) then
                            del=abs(fptt-cost)
                            ibig=i
                        end if
                    end do
                    if(2.*abs(fp-cost).le.spec%ftol*(abs(fp)+abs(cost))) return ! termination criterion
                    if(iter.eq.spec%maxits)then
                        if( warn ) write(*,'(a)') 'powell exceeding maximum iterations; simple_powell_opt'
                        return
                    end if         
                    do j=1,spec%ndim
                        ! construct the extrapolated point and the average
                        ptt(j)=2.*self%spec_linmin%x(j)-pt(j)
                        ! direction moved. Save the old starting point              
                        self%spec_linmin%xi(j)=self%spec_linmin%x(j)-pt(j) 
                        pt(j)=self%spec_linmin%x(j)
                    end do
                    fptt=spec%costfun(ptt,spec%ndim)   ! function value at extrapolated point
                    spec%nevals = spec%nevals+1        ! increment nr ov cost fun evals counter
                    if (fptt.ge.fp) cycle              ! one reason not to use new direction 
                    t=2.*(fp-2.*cost+fptt)*(fp-cost-del)**2-del*(fp-fptt)**2
                    if (t.ge.0.) cycle                 ! another reason not to use new direction
                    self%spec_linmin%nevals = 0
                    call linmin(self%spec_linmin,cost) ! move to the minimum of the new direction
                    spec%nevals = spec%nevals+self%spec_linmin%nevals
                    do j=1,spec%ndim                   ! and save the new direction
                        self%direction_set(j,ibig)=self%direction_set(j,spec%ndim)
                        self%direction_set(j,spec%ndim)=self%spec_linmin%xi(j)
                    end do
                end do
                deallocate(pt,ptt)
            end subroutine

    end subroutine
        
    ! DESTRUCTOR
    
    !> \brief  is a destructor
    subroutine kill_powell_opt( self )
        class(powell_opt), intent(inout) :: self
        if( self%exists )then
            deallocate(self%direction_set)
            call self%spec_linmin%kill
            self%exists = .false.
        endif
    end subroutine
    
end module simple_powell_opt       
