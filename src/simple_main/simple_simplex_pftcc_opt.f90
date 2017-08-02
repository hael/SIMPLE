!------------------------------------------------------------------------------!
! SIMPLE v2.5         Elmlund & Elmlund Lab          simplecryoem.com          !
!------------------------------------------------------------------------------!
!> simple optimisation module: Simplex optimisation on pft correlations
!
!! Minimization of an externally defined function by the simplex method of Nelder and Mead
! The code is distributed with the hope that it will be useful, but _WITHOUT_ _ANY_ _WARRANTY_.
! Redistribution or modification is regulated by the GNU General Public License.
! *Author:* Hans Elmlund, 2013-10-15
module simple_simplex_pftcc_opt
use simple_pftcc_opt, only: pftcc_opt
implicit none

public :: simplex_pftcc_opt
private

logical, parameter :: WARN = .false.

type :: simplex_pftcc_opt
    private
    real, allocatable :: p(:,:)         !< vertices of the simplex
    real, allocatable :: y(:)           !< cost function vals
    real, allocatable :: pb(:)          !< best point
    real              :: yb=0.          !< best cost function value
    logical           :: exists=.false. !< to indicate existence
  contains
    procedure :: new
    procedure :: minimize
    procedure :: kill
end type simplex_pftcc_opt

contains

    !> \brief  is a constructor
    subroutine new( self, spec )
        use simple_opt_spec, only: opt_spec
        use simple_jiffys,   only: alloc_err
        class(simplex_pftcc_opt), intent(inout) :: self !< instance
        class(opt_spec),          intent(inout) :: spec   !< specification
        integer :: alloc_stat
        real    :: x
        call self%kill
        allocate(self%p(spec%ndim+1,spec%ndim), self%y(spec%ndim+1), self%pb(spec%ndim), stat=alloc_stat)
        call alloc_err("In: new_simplex_opt_c", alloc_stat)
        ! initialize best cost to huge number
        self%yb = huge(x)
        self%exists = .true. ! indicates existence
    end subroutine new

    !> \brief  restarted simplex minimization
    subroutine minimize( self, spec, funcontainer, lowest_cost )
        use simple_opt_spec, only: opt_spec
        use simple_opt_subs, only: amoeba
        use simple_rnd,      only: ran3
        class(simplex_pftcc_opt), intent(inout) :: self         !< instance
        class(opt_spec),          intent(inout) :: spec         !< specification
        class(pftcc_opt),         intent(inout) :: funcontainer !< container for the cost function
        real,                     intent(out)   :: lowest_cost  !< lowest cost
        integer :: i, avgniter
        integer :: niters(spec%nrestarts)
        ! generate initial vector
        if( all(spec%x == 0.) )then
            do i=1,spec%ndim
                ! initialize each variable by randomized bounds
                spec%x(i) = spec%limits(i,1)+ran3()*(spec%limits(i,2)-spec%limits(i,1))
            end do
        endif
        ! set best point to best point in spec
        self%pb = spec%x
        ! set best cost
        spec%nevals = 0
        self%yb = funcontainer%costfun(self%pb, spec%ndim)
        spec%nevals = spec%nevals+1
        ! run nrestarts
        do i=1,spec%nrestarts
            call init
            ! run the amoeba routine
            call amoeba_pftcc_opt(self%p,self%y,self%pb,self%yb,spec%ftol,funcontainer,niters(i),spec%maxits,spec%nevals)
        end do
        avgniter    = sum(niters)/spec%nrestarts
        spec%niter  = avgniter
        spec%nevals = spec%nevals/spec%nrestarts
        spec%x      = self%pb
        lowest_cost = self%yb

        contains

            !> \brief  initializes the simplex using randomized bounds
            subroutine init
                use simple_rnd, only: ran3
                integer :: i, j
                ! first vertex is the best-so-far solution solution
                self%p(1,:) = self%pb
                ! the others are obtained by randomized bounds
                do i=2,spec%ndim+1
                    do j=1,spec%ndim
                        self%p(i,j) = spec%limits(j,1)+ran3()*(spec%limits(j,2)-spec%limits(j,1))
                    end do
                end do
                ! calculate costs
                do i=1,spec%ndim+1
                    self%y(i) = funcontainer%costfun(self%p(i,:), spec%ndim)
                end do
            end subroutine init

    end subroutine minimize

    !> \brief  is a destructor
    subroutine kill( self )
        class(simplex_pftcc_opt), intent(inout) :: self
        if( allocated(self%p) )  deallocate(self%p)
        if( allocated(self%y) )  deallocate(self%y)
        if( allocated(self%pb) ) deallocate(self%pb)
        self%exists = .false.
    end subroutine kill

    !> \brief  is a private optimisation routine
    subroutine amoeba_pftcc_opt(p,y,pb,yb,ftol,funcontainer,iter,itmax,nevals)
        use simple_pftcc_opt, only: pftcc_opt
         use simple_jiffys,   only: assert_eq, swap
        real,    intent(inout)  :: p(:,:)       !< the ndim+1 rows of p are ndim vec:s which are the vertices of the starting simplex
                                                !! the best point is put in slot 1 upon convergence
        real,    intent(inout)  :: y(:)         !< must be pre-initialized to the values of func evaluated at the ndim+1 vertices (rows) of p
        real,    intent(inout)  :: pb(:)        !< for updating the best point
        real,    intent(inout)  :: yb           !< for updating the cost of best point
        real,    intent(in)     :: ftol         !< fractional convergence tolerance to be achieved in the function value (0.0005)
        class(pftcc_opt)        :: funcontainer !< container object for the cost function
        integer, intent(out)    :: iter         !< number of exectuted iterations
        integer, intent(in)     :: itmax        !< maximum number of iterations
        integer, intent(inout)  :: nevals       !< number of costfun evals counter
        real, parameter :: tiny=1.0e-10
        real, dimension(size(p,2)) :: psum
        integer :: i,ilo,inhi,loc(1),ihi,ndim
        real :: rtol,ysave,ytry,ytmp
        ndim=assert_eq(size(p,2),size(p,1)-1,size(y)-1,size(pb),'amoeba; simple_opt_subs')
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
                if(WARN) write (*,*) 'itmax exceeded in amoeba; simple_opt_subs'
                exit
            endif
            ! Begin a new iteration. First extrapolate by a factor −1 through the face of the
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
                if(ytry >= ysave)then ! can’t seem to get rid of that high point
                    p(:,:)=0.5*(p(:,:)+spread(p(ilo,:),1,size(p,1))) ! better contract around the lowest (best) point
                    do i=1,ndim+1
                        if(i /= ilo)then
                            y(i)=funcontainer%costfun(p(i,:),ndim)
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

      contains

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
            ytry=funcontainer%costfun(ptry,ndim) ! evaluate the function at the trial point
            nevals = nevals+1
            if(ytry < y(ihi))then ! if it is better than the highest, then replace the highest
                y(ihi)=ytry
                psum(:)=psum(:)-p(ihi,:)+ptry(:)
                p(ihi,:)=ptry(:)
            end if
            amotry=ytry
        end function amotry

    end subroutine amoeba_pftcc_opt

end module simple_simplex_pftcc_opt
