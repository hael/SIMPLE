! function minimization by L-BFGS (Limited memory Broyden–Fletcher–Goldfarb–Shannon optimisation)

module simple_opt_bfgs
include 'simple_lib.f08'
use simple_optimizer, only: optimizer
implicit none

public :: opt_bfgs
private
#include "simple_local_flags.inc"

type, extends(optimizer) :: opt_bfgs
    private
    real, allocatable :: p(:),dg(:),g(:),hdg(:),hessin(:,:),pnew(:),xi(:)
    logical :: exists=.false.
  contains
    procedure :: new          => new_opt_bfgs
    procedure :: minimize     => bfgs_minimize
    procedure :: kill         => kill_opt_bfgs
end type

contains

    !> \brief  is a constructor
    subroutine new_opt_bfgs( self, spec )
        use simple_opt_spec, only: opt_spec
        class(opt_bfgs), intent(inout) :: self !< instance
        class(opt_spec), intent(inout) :: spec !< specification
        call self%kill
        allocate(self%p(spec%ndim),self%dg(spec%ndim),self%hdg(spec%ndim),&
        self%hessin(spec%ndim,spec%ndim),self%pnew(spec%ndim),self%xi(spec%ndim))
        self%dg     = 0.
        self%hdg    = 0.
        self%hessin = 0.
        self%pnew   = 0.
        self%xi     = 0.
        self%p      = 0.
        self%exists = .true.
    end subroutine

    !>  \brief  nonlinear conjugate gradient minimizer
    subroutine bfgs_minimize( self, spec, fun_self, lowest_cost )
        use simple_opt_spec, only: opt_spec
        use simple_opt_subs, only: lnsrch
        class(opt_bfgs), intent(inout) :: self        !< instance
        class(opt_spec), intent(inout) :: spec        !< specification
        class(*),        intent(inout) :: fun_self    !< self-pointer for cost function
        real, intent(out)              :: lowest_cost !< minimum function value
        integer                        :: avgniter,i
        if( .not. associated(spec%costfun) )then
            THROW_HARD('cost function not associated in opt_spec; bfgs_minimize')
        endif
        if( .not. associated(spec%gcostfun) )then
            THROW_HARD('gradient of cost function not associated in opt_spec; bfgs_minimize')
        endif
        ! initialise nevals counters
        spec%nevals  = 0
        spec%ngevals = 0
        avgniter     = 0
        ! init allocated variables
        self%dg     = 0.
        self%hdg    = 0.
        self%hessin = 0.
        self%pnew   = 0.
        self%xi     = 0.
        self%p      = spec%x
        ! run nrestarts restarts
        do i=1,spec%nrestarts
            call bfgsmin
            avgniter = avgniter+spec%niter
        end do
        spec%niter = avgniter/spec%nrestarts
        spec%x     = self%p

        contains

            !>  \brief  nonlinear conjugate gradient minimizer
            subroutine bfgsmin
                real, parameter :: STPMX=100.,EPS=3.e-8,TOLX=4.*EPS
                integer         :: i,its,j
                logical         :: check
                real            :: den,fac,fad,fae,fp
                real            :: stpmax,sum,sumdg,sumxi,temp,test,dgg,gnorm
                fp=spec%costfun(fun_self,self%p,spec%ndim) ! calc start costfun val & gradient
                spec%nevals = spec%nevals+1       ! increment the nevals counter
                call spec%gcostfun(fun_self,self%p,self%g,spec%ndim)
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
                    self%pnew,lowest_cost,stpmax,check,spec%costfun,fun_self,spec%nevals)
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
                    call spec%gcostfun(fun_self,self%p,self%g,spec%ndim) ! and calculate the new gradient
                    test=0.                                     ! test for convergence (zero grad)
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
                if(spec%warn) THROW_WARN('too many iterations in minimize; simple_opt_bfgs')
            end subroutine

    end subroutine

    !> \brief  is a destructor
    subroutine kill_opt_bfgs( self )
        class(opt_bfgs), intent(inout) :: self !< instance
        if( self%exists )then
            deallocate(self%p,self%dg,self%hdg,self%hessin,self%pnew,self%xi)
            if( allocated(self%g) ) deallocate(self%g)
            self%exists = .false.
        endif
    end subroutine

end module simple_opt_bfgs
