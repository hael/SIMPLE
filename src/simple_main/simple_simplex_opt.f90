!==Class simple_simplex_opt
!
! Minimization of an externally defined function by the simplex method of Nelder and Mead
! The code is distributed with the hope that it will be useful, but _WITHOUT_ _ANY_ _WARRANTY_. 
! Redistribution or modification is regulated by the GNU General Public License. 
! *Author:* Hans Elmlund, 2013-10-15
module simple_simplex_opt
use simple_optimizer, only: optimizer
implicit none

public :: simplex_opt
private

type, extends(optimizer) :: simplex_opt
    private
    real, allocatable :: p(:,:)         !< vertices of the simplex
    real, allocatable :: y(:)           !< cost function vals
    real, allocatable :: pb(:)          !< best point
    real              :: yb=0.          !< best cost function value
    logical           :: exists=.false. !< to indicate existence
  contains
    procedure :: new      => new_simplex_opt
    procedure :: minimize => simplex_minimize
    procedure :: kill     => kill_simplex_opt
end type simplex_opt

contains

    !> \brief  is a constructor
    subroutine new_simplex_opt( self, spec )
        use simple_opt_spec, only: opt_spec
        use simple_jiffys,   only: alloc_err
        class(simplex_opt), intent(inout) :: self !< instance
        class(opt_spec), intent(inout)    :: spec !< specification
        integer                           :: alloc_stat
        real                              :: x
        call self%kill
        allocate(self%p(spec%ndim+1,spec%ndim), self%y(spec%ndim+1), self%pb(spec%ndim), stat=alloc_stat)
        call alloc_err("In: new_simplex_opt", alloc_stat)
        ! initialize best cost to huge number
        self%yb = huge(x)
        self%exists = .true. ! indicates existence
    end subroutine new_simplex_opt
    
    !> \brief  restarted simplex minimization
    subroutine simplex_minimize( self, spec, lowest_cost )
        !$ use omp_lib
        !$ use omp_lib_kinds
        use simple_opt_spec, only: opt_spec
        use simple_opt_subs, only: amoeba
        use simple_rnd,      only: ran3
        class(simplex_opt), intent(inout) :: self        !< instance
        class(opt_spec), intent(inout)    :: spec        !< specification
        real, intent(out)                 :: lowest_cost !< lowest cost
        integer                           :: i, avgniter
        integer                           :: niters(spec%nrestarts)
        logical                           :: arezero(spec%ndim)
        if( .not. associated(spec%costfun) )then
            stop 'cost function not associated in opt_spec; simplex_minimize; simple_simplex_opt'
        endif
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
        ! set best point to best point in spec
        self%pb = spec%x
        ! set best cost
        spec%nevals = 0
        self%yb = spec%costfun(self%pb, spec%ndim)
        spec%nevals = spec%nevals+1
        ! run nrestarts
        do i=1,spec%nrestarts
            call init
            ! run the amoeba routine
            call amoeba(self%p,self%y,self%pb,self%yb,spec%ftol,spec%costfun,niters(i),spec%maxits,spec%nevals)
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
                    self%y(i) = spec%costfun(self%p(i,:), spec%ndim)
                end do
            end subroutine init
        
    end subroutine simplex_minimize
    
    !> \brief  is a destructor
    subroutine kill_simplex_opt( self )
        class(simplex_opt), intent(inout) :: self
        if( allocated(self%p) )  deallocate(self%p)
        if( allocated(self%y) )  deallocate(self%y)
        if( allocated(self%pb) ) deallocate(self%pb)
        self%exists = .false.
    end subroutine kill_simplex_opt
    
end module simple_simplex_opt
