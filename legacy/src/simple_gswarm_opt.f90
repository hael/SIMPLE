module simple_gswarm_opt
use simple_optimizer, only: optimizer
use simple_opt_spec,  only: opt_spec
use simple_rnd,       only: ran3, irnd_uni
use simple_defs       ! singleton
implicit none

public :: gswarm_opt
private

type, extends(optimizer) :: gswarm_opt
    private
    real, allocatable    :: swarm(:,:)       !< particle positions
    real, allocatable    :: velocities(:,:)  !< particle velocities
    real, allocatable    :: M(:)             !< masses
    real, allocatable    :: F(:,:,:)         !< forces
    real, allocatable    :: MK(:)            !< masses of the K best
    integer, allocatable :: K(:)             !< the K best solutions
    real                 :: yb               !< best cost
    real                 :: yw               !< worst cost
    integer              :: best             !< index of best particle position
    integer              :: worst            !< index of worst particle position
    real                 :: G=1.             !< gravitational constant
    logical              :: exists = .false. !< to indicate existence
  contains
    procedure :: new      => new_gswarm
    procedure :: minimize => gswarm_minimize
    procedure :: kill     => kill_gswarm
end type gswarm_opt

logical, parameter :: debug=.false.
real, parameter    :: eps=0.001

contains
    
    !> \brief  is a constructor
    subroutine new_gswarm( self, spec )
        use simple_jiffys, only: alloc_err
        class(gswarm_opt), intent(inout) :: self !< instance
        class(opt_spec), intent(inout)   :: spec !< specification
        integer :: alloc_stat
        real    :: x
        ! destruct if exists
        call self%kill
        ! allocate
        allocate(self%swarm(spec%npop,spec%ndim), self%velocities(spec%npop,spec%ndim),&
        self%M(spec%npop), self%F(spec%npop,spec%npop,spec%ndim), self%MK(spec%nnn),&
        self%K(spec%nnn), stat=alloc_stat)
        self%swarm = 0.
        self%M     = 0.
        self%F     = 0.
        self%MK    = 0.
        self%K     = 0
        self%yb    = -huge(x)
        self%yw    = huge(x)
        call alloc_err("In: new_gswarm_opt, 1", alloc_stat)
        self%exists = .true. ! indicates existence
        if( spec%debug ) write(*,*) 'created new particle swarm'
    end subroutine new_gswarm
    
    !> \brief  is the particle swarm minimize minimization routine
    subroutine gswarm_minimize( self, spec, lowest_cost )
        use simple_math, only: euclid, hpsel
        class(gswarm_opt), intent(inout) :: self        !< instance     
        class(opt_spec), intent(inout)           :: spec        !< specification
        real, intent(out)                        :: lowest_cost !< lowest cost
        integer :: t                ! iteration counter
        integer :: npeaks           ! number of local minima
        real    :: rtol             ! relative tolerance
        real    :: costs(spec%npop) ! particle costs
        integer :: loc(1), nworse, i
        if( .not. associated(spec%costfun) )then
            stop 'cost function not associated in opt_spec; gswarm_minimize; simple_gswarm_opt'
        endif
        
        
        
        
        ! initialize
        call init
        spec%nevals = 0
        npeaks      = 0
        nworse      = 0
        do t=1,spec%maxits ! generations loop
            call update_masses
            call update_forces
            call find_kbest
            do i=1,spec%npop
                call update_particle( i )
            end do
            ! update G
            self%G = 0.9*self%G
            
            if( spec%npeaks > 0 )then
                ! exit when we have identified spec%npeaks local optima
                if( npeaks == spec%npeaks ) exit
            endif
            if( nworse == 5*spec%npop ) exit ! if no better solutions have been found in the last 5*npop iterations
        end do
        lowest_cost = self%yb
        spec%x = self%swarm(self%best,:)
        
      contains
        
        !> \brief  initialize the particle positions & velocities
        !!         (1)
        subroutine init
            integer :: i, j
            real    :: L
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
                costs(i) = spec%costfun(self%swarm(i,:), spec%ndim)
                spec%nevals = spec%nevals+1
            end do
            ! set best
            loc = minloc(costs)
            self%best = loc(1)
            self%yb   = costs(self%best)
            ! set worst
            loc = maxloc(costs)
            self%worst = loc(1)
            self%yw    = costs(self%worst)
        end subroutine init
        
        !> \brief  update the masses
        subroutine update_masses
            integer :: i
            real    :: s
            s = 0.
            do i=1,spec%npop
                self%M(i)  = (costs(i)-costs(self%worst))/(costs(self%best)-costs(self%worst))
                s = s+self%M(i)
            end do
            self%M = self%M/s ! array operation
        end subroutine update_masses
        
        !> \brief  update the forces
        subroutine update_forces
            integer :: i, j
            real    :: dist
            do i=1,spec%npop-1
                do j=i+1,spec%npop
                    self%F(i,j,:) = self%G*self%M(i)*self%M(j)*(self%swarm(j,:)-self%swarm(i,:))
                    dist = euclid(self%swarm(i,:),self%swarm(j,:))
                    self%F(i,j,:) = self%F(i,j,:)/(dist+eps)
                end do
            end do
        end subroutine update_forces
        
        !> \brief  find the spec%nnn best solutions
        subroutine find_kbest
            real :: M_copy(spec%npop)
            M_copy = self%M
            call hpsel(self%M, self%MK, self%K)
        end subroutine
        
        !> \brief  update a particle's position    
        subroutine update_particle( i )
            integer, intent(in) :: i
            integer :: j
            real    :: y, L, ran
            real :: actF(spec%ndim) ! acting force
            real :: acc(spec%ndim)  ! acceleration
            ! calculate acting force
            actF = 0.
            do j=1,spec%nnn
                if( self%K(j) /= i )then
                    actF = actF*ran3()*self%F(i,self%K(j),:) 
                endif
            end do
            ! calculate acceleration
            acc = actF/self%M(i)
            ! update velocity & enforce boundaries
            ran = ran3()
            do j=1,spec%ndim
                L = spec%limits(j,2)-spec%limits(j,1)
                self%velocities(i,j) = self%velocities(i,j)*ran+acc(j)
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
            y = spec%costfun(self%swarm(i,:), spec%ndim)
            spec%nevals = spec%nevals+1
            costs(i) = y
            if( i == self%worst )then
                ! set new worst
                loc = maxloc(costs)
                self%worst = loc(1)
                self%yw    = costs(self%worst)
            endif            
            if( y < self%yb )then
                ! update best
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
    
    end subroutine gswarm_minimize
    
    !> \brief  is a destructor
    subroutine kill_gswarm( self )
        class(gswarm_opt), intent(inout) :: self !< instance
        if( self%exists )then
            deallocate(self%swarm, self%velocities,&
            self%M, self%F, self%MK, self%K)
            self%exists = .false.
        endif
    end subroutine kill_gswarm
    
end module simple_gswarm_opt