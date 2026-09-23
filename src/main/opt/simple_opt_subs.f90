!@descr: the downhill simplex step (amoeba) behind the restarted simplex optimiser
! The line searches, the hill-climbing selector and the limit corrector that lived here served the
! BFGS and steepest-descent optimisers removed on 2026-09-23 (no production caller).
module simple_opt_subs
use simple_core_module_api
use simple_opt_spec, only: costfun
implicit none

public :: amoeba
private
#include "simple_local_flags.inc"

logical :: warn=.false.

contains

    !> \brief multidimensional minimization of the function func(x) (x(1:ndim)
    !>          is a vector in ndim dimensions) by the downhill simplex method
    !>          of Nelder and Mead.
    !!          The matrix p(1:ndim+1,1:ndim) is input/output. Its ndim+1 rows
    !!          are ndim-dimensional vectors which are the vertices of the
    !!          starting simplex. Input is also the vector y(1:ndim+1), whose
    !!          components must be pre-initialized to the values of funk
    !!          evaluated at the ndim+1 vertices (rows) of p. ftol is the
    !!          fractional convergence tolerance to be achieved in the function
    !!          value. On output, p and y will have been reset to ndim+1 new
    !!          points all within ftol of a minimum function value, and iter
    !!          gives the number of function evaluations taken.
    subroutine amoeba(p,y,pb,yb,ftol,func,fun_self,iter,itmax,nevals)
        real,     intent(inout) :: p(:,:)   !< the ndim+1 rows of p are ndim vec:s which are the vertices of the starting simplex
                                         !! the best point is put in slot 1 upon convergence
        real,     intent(inout) :: y(:)     !< must be pre-initialized to the values of func evaluated at the ndim+1 vertices (rows) of p
        real,     intent(inout) :: pb(:)    !< for updating the best point
        real,     intent(inout) :: yb       !< for updating the cost of best point
        real,     intent(in)    :: ftol     !< fractional convergence tolerance to be achieved in the function value (0.0005)
        class(*), intent(inout) :: fun_self !< self-pointer for cost function
        integer,  intent(out)   :: iter     !< number of exectuted iterations
        integer,  intent(in)    :: itmax    !< maximum number of iterations
        integer,  intent(inout) :: nevals   !< number of costfun evals counter
        procedure(costfun), pointer :: func
        real, parameter :: tiny=1.0e-10
        integer         :: ihi,ndim
        real, dimension(size(p,2)) :: psum
        call amoeba_private

        contains

            subroutine amoeba_private
                integer :: i,ilo,inhi,loc(1)
                real :: rtol,ysave,ytry,ytmp
                if(size(p,2) == size(p,1)-1 .and. size(p,1)-1 == size(y)-1 .and. size(y)-1 == size(pb)) then
                    ndim = size(p,2)
                else
                    THROW_HARD('assert eq failed; amoeba private')
                end if
                iter=0
                psum(:)=sum(p(:,:),dim=1)
                do
                    loc=minloc(y) ! determine which point has the highest (worst), next-highest, and lowest (best)
                    ilo=loc(1)
                    loc=maxloc(y)
                    ihi=loc(1)
                    ytmp=y(ihi)
                    y(ihi)=y(ilo)
                    loc=maxloc(y)
                    inhi=loc(1)
                    y(ihi)=ytmp
                    ! Compute the fractional range from highest to lowest and return if satisfactory
                    rtol=2.0*abs(y(ihi)-y(ilo))/(abs(y(ihi))+abs(y(ilo))+tiny) ! relative tolerance
                    if(rtol < ftol)then
                        call swap(y(1),y(ilo))
                        call swap(p(1,:),p(ilo,:))
                        exit
                    end if
                    if(iter >= itmax)then
                        if(warn) THROW_WARN('itmax exceeded in amoeba')
                        exit
                    endif
                    ! Begin a new iteration. First extrapolate by a factor -1 through the face of the
                    ! simplex across from the high point, i.e. reflect the simplex from the high point
                    ytry=amotry(-1.0)
                    iter=iter+1
                    if(ytry <= y(ilo))then
                        ! gives a result better than the best point, so try an additional extrapolation by a factor 2
                        ytry=amotry(2.0)
                        iter=iter+1
                    else if(ytry >= y(inhi))then
                        ! the reflected point is worse than the second-highest, so look for an intermediate
                        ! lower point, i.e. do a one-dimensional contraction
                        ysave=y(ihi)
                        ytry=amotry(0.5)
                        iter=iter+1
                        if(ytry >= ysave)then ! can't seem to get rid of that high point
                            p(:,:)=0.5*(p(:,:)+spread(p(ilo,:),1,size(p,1))) ! better contract around the lowest (best) point
                            do i=1,ndim+1
                                if(i /= ilo)then
                                    y(i)=func(fun_self,p(i,:),ndim)
                                    nevals = nevals+1
                                endif
                            end do
                            iter=iter+ndim ! keep track of function evaluations
                            psum(:)=sum(p(:,:),dim=1) ! recompute psum
                        end if
                    end if
                end do
                ! store best
                if( y(1) <= yb )then
                    pb(:) = p(1,:)
                    yb    = y(1)
                endif
            end subroutine amoeba_private

            !>  \brief  extrapolates by a factor fac through the face of the simplex across from the
            !!          high point, tries it, and replaces the high point if the new point is better
            function amotry(fac)
                real, intent(in) :: fac
                real :: amotry
                real :: fac1,fac2,ytry
                real, dimension(size(p,2)) :: ptry
                fac1=(1.0-fac)/ndim
                fac2=fac1-fac
                ptry(:)=psum(:)*fac1-p(ihi,:)*fac2
                ytry=func(fun_self,ptry,ndim)  ! evaluate the function at the trial point
                nevals = nevals+1
                if(ytry < y(ihi))then ! if it is better than the highest, then replace the highest
                    y(ihi)=ytry
                    psum(:)=psum(:)-p(ihi,:)+ptry(:)
                    p(ihi,:)=ptry(:)
                end if
                amotry=ytry
            end function amotry

    end subroutine amoeba

end module simple_opt_subs
