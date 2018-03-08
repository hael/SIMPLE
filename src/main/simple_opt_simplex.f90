! The Nelder-Mead simplex method for continuous function minimisation
module simple_opt_simplex
include 'simple_lib.f08'
use simple_optimizer, only: optimizer
use simple_opt_spec, only: opt_spec
implicit none

public :: opt_simplex
private

type, extends(optimizer) :: opt_simplex
    private
    real, allocatable :: p(:,:)         !< vertices of the simplex
    real, allocatable :: y(:)           !< cost function vals
    real, allocatable :: pb(:)          !< best point
    real              :: yb=0.          !< best cost function value
    logical           :: exists=.false. !< to indicate existence
  contains
    procedure :: new      => new_opt_simplex
    procedure :: minimize => simplex_minimize
    procedure :: kill     => kill_opt_simplex
end type opt_simplex

contains

    !> \brief  is a constructor
    subroutine new_opt_simplex( self, spec )
        use simple_opt_spec, only: opt_spec
        class(opt_simplex), intent(inout) :: self !< instance
        class(opt_spec),    intent(inout) :: spec !< specification
        real :: x
        call self%kill
        allocate(self%p(spec%ndim+1,spec%ndim), self%y(spec%ndim+1), self%pb(spec%ndim), stat=alloc_stat)
        if(alloc_stat.ne.0)call allocchk("In: new_opt_simplex",alloc_stat)
        ! initialize best cost to huge number
        self%yb = huge(x)
        self%exists = .true. ! indicates existence
    end subroutine new_opt_simplex

    !> \brief  restarted simplex minimization
    subroutine simplex_minimize( self, spec, fun_self, lowest_cost )
        use simple_opt_spec, only: opt_spec
        use simple_opt_subs, only: amoeba
        class(opt_simplex), intent(inout) :: self        !< instance
        class(opt_spec),    intent(inout) :: spec        !< specification
        class(*),           intent(inout) :: fun_self    !< self-pointer for cost function
        real,               intent(out)   :: lowest_cost !< lowest cost
        real, allocatable :: lims_dyn(:,:)
        integer           :: i, avgniter
        integer           :: niters(spec%nrestarts)
        logical           :: arezero(spec%ndim)
        if( .not. associated(spec%costfun) )then
           call simple_stop ('cost function not associated in opt_spec; simplex_minimize; simple_opt_simplex')
        end if
        ! initialise nevals counter
        spec%nevals = 0
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
        self%yb     = spec%costfun(fun_self, self%pb, spec%ndim)
        spec%nevals = spec%nevals + 1
        ! initialise dynamic bounds
        if( allocated(spec%limits_init) )then
            allocate(lims_dyn(spec%ndim,2), source=spec%limits_init)
        endif
        ! run nrestarts
        do i=1,spec%nrestarts
            if( allocated(spec%limits_init) )then
                call init( lims_dyn )
            else
                call init( spec%limits )
            endif
            ! run the amoeba routine
            call amoeba(self%p,self%y,self%pb,self%yb,spec%ftol,spec%costfun,fun_self,niters(i),spec%maxits,spec%nevals)
            if( allocated(spec%limits_init) )then
                ! movie the shift limits to around the best point
                ! left limit
                lims_dyn(:,1) = self%pb - (lims_dyn(:,2) - lims_dyn(:,1)) / 2.0
                ! right limit
                lims_dyn(:,2) = self%pb + (lims_dyn(:,2) - lims_dyn(:,1)) / 2.0
                ! bound according to spec%limits
                where( lims_dyn(:,1) < spec%limits(:,1) ) lims_dyn(:,1) = spec%limits(:,1)
                where( lims_dyn(:,2) > spec%limits(:,2) ) lims_dyn(:,2) = spec%limits(:,2)
            endif
        end do
        avgniter    = sum(niters)/spec%nrestarts
        spec%niter  = avgniter
        spec%x      = self%pb
        lowest_cost = self%yb

        contains

             !> \brief  initializes the simplex using randomized bounds
            subroutine init( limits )
                real, intent(in) :: limits(spec%ndim,2)
                integer :: i, j
                ! first vertex is the best-so-far solution solution
                self%p(1,:) = self%pb
                ! the others are obtained by randomized bounds
                do i=2,spec%ndim + 1
                    do j=1,spec%ndim
                        self%p(i,j) = limits(j,1) + ran3() * (limits(j,2) - limits(j,1))
                    end do
                end do
                ! calculate costs
                do i=1,spec%ndim + 1
                    self%y(i) = spec%costfun(fun_self, self%p(i,:), spec%ndim)
                    spec%nevals = spec%nevals + 1
                end do
            end subroutine init

    end subroutine simplex_minimize

    !> \brief  is a destructor
    subroutine kill_opt_simplex( self )
        class(opt_simplex), intent(inout) :: self
        alloc_stat=0
        if( allocated(self%p) )then
            deallocate(self%p, stat=alloc_stat)
            if(alloc_stat.ne.0)call allocchk("In: kill_simplex_opt p",alloc_stat)
        end if
        if( allocated(self%y) )then
            deallocate(self%y, stat=alloc_stat)
            if(alloc_stat.ne.0)call allocchk("In: kill_simplex_opt y ",alloc_stat)
        end if
        if( allocated(self%pb) )then
            deallocate(self%pb, stat=alloc_stat)
            if(alloc_stat.ne.0)call allocchk("In: kill_simplex_opt pb ",alloc_stat)
        end if
        self%exists = .false.
    end subroutine kill_opt_simplex

end module simple_opt_simplex
