!==Class simple_bfgs_opt
!
! Minimization of an externally defined function by the nonlinear conjugate gradient method
! The code is distributed with the hope that it will be useful, but _WITHOUT_ _ANY_ _WARRANTY_. 
! Redistribution or modification is regulated by the GNU General Public License. 
! *Author:* Hans Elmlund, 2013-10-15
!
module simple_bfgs_opt
use simple_optimizer, only: optimizer
implicit none

public :: bfgs_opt
private

type, extends(optimizer) :: bfgs_opt
    private
    real, allocatable :: p(:),dg(:),g(:),hdg(:),hessin(:,:),pnew(:),xi(:)
    logical :: exists=.false.
  contains
    procedure :: new      => new_bfgs_opt
    procedure :: minimize => bfgs_minimize
    procedure :: kill     => kill_bfgs_opt
end type

contains
    
    !> \brief  is a constructor
    subroutine new_bfgs_opt( self, spec )
        use simple_opt_spec, only: opt_spec
        use simple_jiffys,   only: alloc_err
        class(bfgs_opt), intent(inout) :: self !< instance
        class(opt_spec), intent(inout) :: spec !< specification
        integer :: alloc_stat
        call self%kill
        allocate(self%p(spec%ndim),self%dg(spec%ndim),self%hdg(spec%ndim),&
        self%hessin(spec%ndim,spec%ndim),self%pnew(spec%ndim),self%xi(spec%ndim),stat=alloc_stat)
        call alloc_err('In: new_bfgs_opt; simple_bfgs_opt', alloc_stat)
        self%dg     = 0.
        self%hdg    = 0.
        self%hessin = 0.
        self%pnew   = 0.
        self%xi     = 0.
        self%p      = 0.
        self%exists = .true.
    end subroutine
    
    !>  \brief  nonlinear conjugate gradient minimizer
    subroutine bfgs_minimize( self, spec, lowest_cost )
        use simple_opt_spec, only: opt_spec
        use simple_opt_subs, only: lnsrch
        class(bfgs_opt), intent(inout) :: self        !< instance
        class(opt_spec), intent(inout) :: spec        !< specification
        real, intent(out)              :: lowest_cost !< minimum function value
        integer                        :: avgniter,i
        if( .not. associated(spec%costfun) )then 
            stop 'cost function not associated in opt_spec; bfgs_minimize; simple_bfgs_opt'
        endif
        if( .not. associated(spec%gcostfun) )then
            stop 'gradient of cost function not associated in opt_spec; bfgs_minimize; simple_bfgs_opt'
        endif
        ! init allocated variables
        self%dg     = 0.
        self%hdg    = 0.
        self%hessin = 0.
        self%pnew   = 0.
        self%xi     = 0.
        self%p      = spec%x
        ! run nrestarts restarts
        avgniter = 0
        spec%nevals = 0
        do i=1,spec%nrestarts
            call bfgsmin
            avgniter = avgniter+spec%niter
        end do
        spec%niter  = avgniter/spec%nrestarts
        spec%nevals = spec%nevals/spec%nrestarts
        spec%x      = self%p

        contains
        
            !>  \brief  nonlinear conjugate gradient minimizer
            subroutine bfgsmin
                real, parameter :: STPMX=100.,EPS=3.e-8,TOLX=4.*EPS
                integer         :: i,its,j
                logical         :: check
                real            :: den,fac,fad,fae,fp
                real            :: stpmax,sum,sumdg,sumxi,temp,test,dgg,gnorm
                fp=spec%costfun(self%p,spec%ndim) ! calc start costfun val & gradient
                spec%nevals = spec%nevals+1       ! increment the nevals counter
                self%g=spec%gcostfun(self%p,spec%ndim)
                sum=0.
                do i=1,spec%ndim     ! init Hessian 2 unit matrix
                    do j=1,spec%ndim
                        self%hessin(i,j)=0.
                    end do
                    self%hessin(i,i)=1.
                    self%xi(i)=-self%g(i) ! init line direction, direction is simply negative gradient
                    sum=sum+self%p(i)**2
                end do
                stpmax=STPMX*max(sqrt(sum),real(spec%ndim))
                do its=1,spec%maxits                       ! main loop
                    spec%niter=its
                    call lnsrch(spec%ndim,self%p,fp,self%g,self%xi,&
                    self%pnew,lowest_cost,stpmax,check,spec%costfun,spec%nevals)
                    fp=lowest_cost
                    do i=1,spec%ndim
                        self%xi(i)=self%pnew(i)-self%p(i)  ! update line direction
                        self%p(i)=self%pnew(i)             ! and current point
                    end do
                    test=0.                                ! test for convergence (diff x)
                    do i=1,spec%ndim
                        temp=abs(self%xi(i))/max(abs(self%p(i)),1.)
                        if(temp.gt.test)test=temp
                    end do
                    if(test.lt.TOLX)return
                    do i=1,spec%ndim                       ! save the old gradient
                        self%dg(i)=self%g(i)
                    end do
                    self%g=spec%gcostfun(self%p,spec%ndim) ! and calculate the new gradient
                    test=0.                                ! test for convergence (zero grad)
                    den=max(lowest_cost,1.)
                    do i=1,spec%ndim
                        temp=abs(self%g(i))*max(abs(self%p(i)),1.)/den
                        if(temp.gt.test)test=temp
                    end do
                    if(test.lt.spec%gtol)return
                    dgg=0.
                    gnorm=0.
                    do i=1,spec%ndim                       ! compute difference of gradients
                        dgg=dgg+self%dg(i)*self%g(i)       ! scalarprod(old_grad,grad)
                        gnorm=gnorm+self%g(i)*self%g(i)    ! gradient norm
                        self%dg(i)=self%g(i)-self%dg(i)
                    end do
                    do i=1,spec%ndim
                        self%hdg(i)=0.
                        do j=1,spec%ndim                   ! and difference times current matrix
                            self%hdg(i)=self%hdg(i)+self%hessin(i,j)*self%dg(j)
                        end do
                    end do
                    fac=0. ! calc dot products 4 the denominators
                    fae=0.
                    sumdg=0.
                    sumxi=0.
                    do i=1,spec%ndim
                        fac=fac+self%dg(i)*self%xi(i)
                        fae=fae+self%dg(i)*self%hdg(i)
                        sumdg=sumdg+self%dg(i)**2
                        sumxi=sumxi+self%xi(i)**2
                    end do
                    if(fac.gt.sqrt(EPS*sumdg*sumxi))then ! skip update if fac not sufficiently positive
                        fac=1./fac
                        fad=1./fae
                        do i=1,spec%ndim                 ! the vector that makes BFGS different from DFP
                            self%dg(i)=fac*self%xi(i)-fad*self%hdg(i)
                        end do
                        do i=1,spec%ndim                 ! the BFGS update formula:
                            do j=i,spec%ndim
                                self%hessin(i,j)=self%hessin(i,j)+fac*self%xi(i)*self%xi(j)-&
                                fad*self%hdg(i)*self%hdg(j)+fae*self%dg(i)*self%dg(j)
                                self%hessin(j,i)=self%hessin(i,j)
                            end do
                        end do
                    endif
                    do i=1,spec%ndim ! calc direction to go
                        self%xi(i)=0.
                        do j=1,spec%ndim
                            self%xi(i)=self%xi(i)-self%hessin(i,j)*self%g(j)
                        end do
                    end do
                end do
                if(spec%warn) write(*,'(a)') 'too many iterations in minimize; simple_bfgs_opt'
            end subroutine

    end subroutine
    
    !> \brief  is a destructor
    subroutine kill_bfgs_opt( self )
        class(bfgs_opt), intent(inout) :: self !< instance
        if( self%exists )then
            deallocate(self%p,self%dg,self%hdg,self%hessin,self%pnew,self%xi)
            if( allocated(self%g) ) deallocate(self%g)
            self%exists = .false.
        endif
    end subroutine
 
end module simple_bfgs_opt
