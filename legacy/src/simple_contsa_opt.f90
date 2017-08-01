!==Class simple_contsa_opt
!
! Minimization of an externally defined function by continuous simulated annealing
! The code is distributed with the hope that it will be useful, but _WITHOUT_ _ANY_ _WARRANTY_. 
! Redistribution or modification is regulated by the GNU General Public License. 
! *Author:* Hans Elmlund, 2013-10-15
module simple_contsa_opt
use simple_optimizer, only: optimizer
implicit none

public :: contsa_opt
private

type, extends(optimizer) :: contsa_opt
    private
    real, allocatable :: p(:,:)         !< vertices of the contsa
    real, allocatable :: y(:)           !< cost function vals
    real, allocatable :: pb(:)          !< best point
    real              :: yb=0.          !< best cost function value
    logical           :: exists=.false. !< to indicate existence
  contains
    procedure :: new      => new_contsa_opt
    procedure :: minimize => contsa_minimize
    procedure :: kill     => kill_contsa_opt
end type

contains

    ! CONSTRUCTORS

    !> \brief  is a constructor
    subroutine new_contsa_opt( self, spec )
        use simple_opt_spec, only: opt_spec
        use simple_jiffys,   only: alloc_err
        class(contsa_opt), intent(inout) :: self
        class(opt_spec), intent(in)      :: spec
        integer                          :: alloc_stat
        real                             :: x
        call self%kill
        allocate(self%p(spec%ndim+1,spec%ndim), self%y(spec%ndim+1), self%pb(spec%ndim), stat=alloc_stat)
        call alloc_err("In: new_contsa_opt", alloc_stat)
        ! initialize best cost to huge number
        self%yb = huge(x)
        self%exists = .true. ! indicates existence
    end subroutine
    
    ! THE STUFF THAT DO REAL WORK
    
    !> \brief  restarted contsa minimization
    subroutine contsa_minimize( self, spec, lowest_cost )
        use simple_opt_spec, only: opt_spec
        use simple_opt_subs, only: amebsa
        class(contsa_opt), intent(inout) :: self        !< instance
        class(opt_spec), intent(inout)   :: spec        !< specification
        real, intent(out)                :: lowest_cost !< lowest cost
        integer                          :: iter
        real                             :: t
        if( .not. associated(spec%costfun) )then
            stop 'cost function not associated in opt_spec; contsa_minimize; simple_contsa_opt'
        endif
        ! set best point to point in spec
        self%pb = spec%x
        ! set best cost
        self%yb = spec%costfun(self%pb, spec%ndim)
        ! make an initial simplex
        call init
        t = spec%tinit 
        do
            ! run the amebsa routine
            iter = spec%maxits
            spec%converged = .false. 
            call amebsa(self%p,self%y,self%pb,self%yb,spec%ftol,spec%costfun,iter,t)
            ! upodate temp
            t = t*spec%tc
            ! report early convergence via spec
            if( iter >= 0 .or. t < spec%tmin ) exit
        end do
        spec%x      = self%pb
        lowest_cost = self%yb
        
        contains
        
            !> \brief  initializes the contsa using randomized bounds
            subroutine init
                use simple_rnd, only: ran3
                integer :: i, j
                ! first vertex is the best-so-far solution solution
                self%p(1,:) = self%pb
                ! randomly permute init solution to obtain the other vertices
                do i=2,spec%ndim+1
                    do j=1,spec%ndim
                        self%p(i,j) = spec%limits(j,1)+ran3()*(spec%limits(j,2)-spec%limits(j,1))
                    end do
                end do
                ! calculate costs
                do i=1,spec%ndim+1
                    self%y(i) = spec%costfun(self%p(i,:), spec%ndim)
                end do
            end subroutine
        
    end subroutine
    
    ! DESTRUCTOR
    
    !> \brief  is a destructor
    subroutine kill_contsa_opt( self )
        class(contsa_opt), intent(inout) :: self
        if( self%exists )then
            deallocate( self%p, self%y )
            self%exists = .false.
        endif
    end subroutine
    
end module simple_contsa_opt
