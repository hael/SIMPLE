! continuous particle swarm optimisation

module simple_opt_particle_swarm
include 'simple_lib.f08'
use simple_optimizer, only: optimizer
use simple_opt_spec,  only: opt_spec
implicit none

public :: opt_particle_swarm
private
#include "simple_local_flags.inc"

type, extends(optimizer) :: opt_particle_swarm
    private
    real, allocatable :: swarm(:,:)       !< particle positions
    real, allocatable :: velocities(:,:)  !< particle velocities
    real              :: yb               !< best cost
    integer           :: best             !< index of best particle position
    real              :: inertia = -0.2   !< velocity memory [0.2,0.3]
    real              :: gamma   = 1.5    !< control parameter [1.5,2.5]
    logical           :: exists = .false. !< to indicate existence
  contains
    procedure :: new          => new_particle_swarm
    procedure :: minimize     => particle_swarm_minimize
    procedure :: kill         => kill_particle_swarm
end type opt_particle_swarm

contains

    !> \brief  is a constructor
    subroutine new_particle_swarm( self, spec )
        class(opt_particle_swarm), intent(inout) :: self !< instance
        class(opt_spec), intent(inout)           :: spec !< specification
        ! destruct if exists
        call self%kill
        ! allocate
        allocate(self%swarm(spec%npop,spec%ndim), self%velocities(spec%npop,spec%ndim), stat=alloc_stat)
        if(alloc_stat.ne.0)call allocchk("In: new_particle_swarm, 1",alloc_stat)
        self%exists = .true. ! indicates existence
    end subroutine new_particle_swarm

    !> \brief  is the particle swarm minimize minimization routine
    subroutine particle_swarm_minimize( self, spec, fun_self, lowest_cost )
        class(opt_particle_swarm), intent(inout) :: self        !< instance
        class(opt_spec),           intent(inout) :: spec        !< specification
        class(*),                  intent(inout) :: fun_self    !< self-pointer for cost function
        real, intent(out)                        :: lowest_cost !< lowest cost
        integer :: t                !< iteration counter
        integer :: npeaks           !< number of local minima
        !real    :: rtol             !< relative tolerance
        real    :: costs(spec%npop) !< particle costs
        integer :: loc(1), nworse
        if( .not. associated(spec%costfun) )then
            THROW_HARD('cost function not associated in opt_spec; particle_swarm_minimize')
        endif
        ! initialize
        call init
        spec%nevals = 0
        npeaks      = 0
        nworse      = 0
        do t=1,spec%maxits ! generations loop
            ! update a randomly selected particle
            call update_particle(irnd_uni(spec%npop))
            if( spec%npeaks > 0 )then
                ! exit when we have idetified spec%npeaks local optima
                if( npeaks == spec%npeaks ) exit
            endif
            if( nworse == 5*spec%npop ) exit ! if no better solutions have been found in the last 5*npop iterations
        end do
        lowest_cost = self%yb
        spec%x = self%swarm(self%best,:)

      contains

        !> \brief  initialize the particle positions & velocities
        subroutine init
            integer :: i, j
            real :: L
            ! obtain particle positions by randomized bounds
            do i=1,spec%npop
                do j=1,spec%ndim
                    L = spec%limits(j,2)-spec%limits(j,1)
                    self%swarm(i,j) = spec%limits(j,1)+ran3()*L
                end do
                if( i == 1 )then
                    if( .not. all(spec%x == 0.) )then ! set one point in the swarm to best point in spec
                        self%swarm(1,:)= spec%x
                    endif
                endif
            end do
            ! init velocities
            do i=1,spec%npop
                do j=1,spec%ndim
                    L = spec%limits(j,2)-spec%limits(j,1)
                    self%velocities(i,j) = ran3()*2*L-L
                end do
            end do
            ! calculate initial costs
            do i=1,spec%npop
                costs(i) = spec%costfun(fun_self, self%swarm(i,:), spec%ndim)
                spec%nevals = spec%nevals+1
            end do
            loc = minloc(costs)
            self%best = loc(1)
            self%yb   = costs(self%best)
        end subroutine init

        subroutine update_particle( i )
            integer, intent(in) :: i
            integer :: j
            real    :: y, L
            ! update velocity & enforce boundaries
            do j=1,spec%ndim
                L = spec%limits(j,2)-spec%limits(j,1)
                self%velocities(i,j) = self%inertia*self%velocities(i,j)+&
                self%gamma*ran3()*(self%swarm(self%best,j)-self%swarm(i,j))
                self%velocities(i,j) = min(L,self%velocities(i,j))
                self%velocities(i,j) = max(-L,self%velocities(i,j))
            end do
            ! move & enforce boundaries
            do j=1,spec%ndim
                self%swarm(i,j) = self%swarm(i,j)+self%velocities(i,j)
                self%swarm(i,j) = min(spec%limits(j,2),self%swarm(i,j))
                self%swarm(i,j) = max(spec%limits(j,1),self%swarm(i,j))
            end do
            ! calculate cost
            y = spec%costfun(fun_self, self%swarm(i,:), spec%ndim)
            costs(i) = y
            spec%nevals = spec%nevals+1
            if( y < self%yb )then
                nworse = 0
                self%best = i
                self%yb   = y
                if( spec%npeaks > 0 )then
                    ! we will output local optima
                    if( y < spec%yb )then
                        npeaks = npeaks+1
                        spec%peaks(npeaks,:spec%ndim)  = self%swarm(self%best,:)
                        spec%peaks(npeaks,spec%ndim+1) = y
                    endif
                endif
            else
                nworse = nworse+1
            endif
        end subroutine update_particle

    end subroutine particle_swarm_minimize

    !> \brief  is a destructor
    subroutine kill_particle_swarm( self )
        class(opt_particle_swarm), intent(inout) :: self !< instance
        if( self%exists )then
            deallocate(self%swarm, self%velocities, stat=alloc_stat)
            if(alloc_stat.ne.0)call allocchk("simple_opt_particle_swarm::kill",alloc_stat)
            self%exists = .false.
        endif
    end subroutine kill_particle_swarm

end module simple_opt_particle_swarm
