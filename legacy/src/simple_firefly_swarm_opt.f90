module simple_firefly_swarm_opt
use simple_optimizer, only: optimizer
use simple_opt_spec,  only: opt_spec
use simple_rnd,       only: ran3
use simple_defs       ! singleton
implicit none

public :: firefly_swarm_opt, test_firefly_swarm_opt
private

type :: firefly 
    real :: brightness          !< brightness of firefly (f(x)) 
    real, allocatable :: vec(:) !< solution vector
end type firefly

type, extends(optimizer) :: firefly_swarm_opt
    private
    type(firefly), allocatable :: flies(:)          !< population of fireflies
    real, allocatable          :: dmat_eulers(:,:)  !< distance matrix for Euler angles
    real, allocatable          :: dmat_shifts(:,:)  !< distance matrix for shifts
    real, allocatable          :: alphas(:)         !< randomization parameters
    real, allocatable          :: gammas(:)         !< light absorption coefficients
    type(opt_spec)             :: spec_linmin       !< specification of linmin optimizer
    logical                    :: exists = .false.  !< to indicate existence
 contains
    procedure :: new      => new_firefly_swarm
    procedure :: minimize => firefly_swarm_minimize
    procedure :: kill     => kill_firefly_swarm
end type firefly_swarm_opt

real, parameter    :: beta0 = 1.
logical, parameter :: debug=.false.

contains
    
    !> \brief  is a constructor
    subroutine new_firefly_swarm( self, spec )
        use simple_jiffys,   only: alloc_err
        class(firefly_swarm_opt), intent(inout) :: self !< instance
        class(opt_spec), intent(in)             :: spec !< specification
        integer :: alloc_stat, i
        ! destruct if exists
        call self%kill
        ! allocate
        allocate(self%flies(spec%npop), self%dmat_eulers(spec%npop,spec%npop),&
        self%dmat_shifts(spec%npop,spec%npop), self%alphas(spec%ndim), self%gammas(spec%ndim), stat=alloc_stat)
        call alloc_err("In: new_firefly_swarm_opt, 1", alloc_stat)
        ! make flies
        do i=1,spec%npop
            allocate(self%flies(i)%vec(spec%ndim), stat=alloc_stat)
            call alloc_err("In: new_firefly_swarm_opt, 2", alloc_stat)
        end do
        ! make line minimizer
        call self%spec_linmin%specify('linmin', spec%ndim, maxits=spec%nsample, ftol=spec%ftol)
        call self%spec_linmin%set_costfun(spec%costfun) 
        self%exists = .true. ! indicates existence
        if( spec%debug ) write(*,*) 'created new firefly swarm'
    end subroutine new_firefly_swarm
    
    !> \brief  is the firefly swarm minimize minimization routine
    !!         (1) generate an initial swarm of fireflies
    !!         (2) set brightnesses (I:s) to -cost
    !!         (3) define light absorption coefficients
    !!         while(t <= spec%maxits)
    !!             for i=1,n-1
    !!                 for j=i+1,n
    !!                     if(I(i) > I(j)), move fly j towards i in d-dimension; end if
    !!                     Attractiveness varies with distance d(ij) via exp{-gamma*d(ij)}
    !!                     Evaluate new solutions and update brightnesses (I:s)
    !!                 end for j
    !!             end for i
    !!         end while
    subroutine firefly_swarm_minimize( self, spec, lowest_cost )
        use simple_opt_subs, only: linmin
        class(firefly_swarm_opt), intent(inout) :: self        !< instance     
        class(opt_spec), intent(inout)          :: spec        !< specification
        real, intent(out)                       :: lowest_cost !< lowest cost
        integer :: t      ! iteration counter
        integer :: bright ! brightest (best)  solution
        integer :: dark   ! darkest   (worst) solution
        integer :: npeaks ! number of local minima
        real    :: rtol   ! relative tolerance
        if( .not. associated(spec%costfun) )then
            stop 'cost function not associated in opt_spec; firefly_swarm_minimize; simple_firefly_swarm_opt'
        endif
        ! initialize
        call init
        spec%nevals = 0
        npeaks      = 0
        call calc_brightnesses
        select case(trim(spec%str_mode))
            case('euler')
                call calc_dmat_euler
                if( spec%ndim > 3 ) call calc_dmat_shifts
            case('euclid')
                call calc_dmat
            case DEFAULT
                stop 'unknown mode descriptor; firefly_swarm_minimize; simple_firefly_swarm_opt'
        end select
        if( debug )then
            call print_flies
            call print_dists
        endif
        do t=1,spec%maxits ! generations loop (default maxits=100)
            call move_flies
            call update_alphas
            call find_bright_and_dark
            if( spec%verbose ) write(*,*) 'Iteration: ', t, 'Brightness: ', self%flies(bright)%brightness, 'Solution: ', bright
            rtol=2.*abs(self%flies(bright)%brightness-self%flies(dark)%brightness)/&
            (abs(self%flies(bright)%brightness)+abs(self%flies(dark)%brightness)+TINY)
            if( spec%npeaks > 0 )then
                ! exit when we have identified spec%npeaks local optima
                if( npeaks == spec%npeaks ) exit
            endif
            if(rtol < spec%ftol) exit   
        end do
        call generate_ouput
        
        
      contains
          
          !> \brief  initializes the fly positions, brightnesses, and control parameters
          subroutine init
              integer :: i, j
              real    :: L
              ! obtain flies by randomized bounds
              do i=1,spec%npop
                  do j=1,spec%ndim
                      self%flies(i)%vec(j) = spec%limits(j,1)+ran3()*(spec%limits(j,2)-spec%limits(j,1))
                  end do
                  if( i == 1 )then
                      if( .not. all(spec%x == 0.) )then ! set one point in the pop to best point in spec
                          self%flies(1)%vec = spec%x
                      endif
                  endif
              end do
              ! initialize randomization and attractiveness parameters
              do j=1,spec%ndim
                  L = spec%limits(j,2)-spec%limits(j,1)
                  self%alphas(j) = 0.01*L
                  self%gammas(j) = 0.5/(L*L)
              end do
          end subroutine init
          
          !> \brief  calculates the brightnessses
          subroutine calc_brightnesses
              integer :: i
              do i=1,spec%npop
                  self%flies(i)%brightness = -spec%costfun(self%flies(i)%vec, spec%ndim)
                  spec%nevals = spec%nevals+1
              end do
          end subroutine
          
          !> \brief  calculates the distance matrix for mode=euclid
          subroutine calc_dmat
              use simple_math, only: euclid
              integer :: i, j
              self%dmat_eulers = 0.
              do i=1,spec%npop-1
                  do j=i+1,spec%npop
                      self%dmat_eulers(i,j) = euclid(self%flies(i)%vec,self%flies(j)%vec)
                      self%dmat_eulers(j,i) = self%dmat_eulers(i,j) 
                  end do
              end do
          end subroutine calc_dmat
          
          !> \brief  dynamic update of the distance matrix for mode=euclid
          subroutine update_dmat( fly_moved )
              use simple_math, only: euclid
              integer, intent(in) :: fly_moved
              integer :: i
              do i=1,spec%npop
                  if( i /= fly_moved )then
                      self%dmat_eulers(i,fly_moved) = euclid(self%flies(i)%vec,self%flies(fly_moved)%vec)
                      self%dmat_eulers(fly_moved,i) = self%dmat_eulers(i,fly_moved)
                  endif
              end do
          end subroutine update_dmat
          
          !> \brief  calculates the Euler distance matrix for mode=euler
          subroutine calc_dmat_euler
              use simple_ori,  only: ori
              use simple_math, only: rad2deg
              integer   :: i, j
              type(ori) :: oi, oj
              oi = ori()
              oj = ori()
              self%dmat_eulers = 0.
              do i=1,spec%npop-1
                  do j=i+1,spec%npop
                      select case(spec%ndim)
                          case(2)
                              call oi%e1set(self%flies(i)%vec(1))
                              call oj%e1set(self%flies(j)%vec(1))
                              call oi%e2set(self%flies(i)%vec(2))
                              call oj%e2set(self%flies(j)%vec(2))
                              self%dmat_eulers(i,j) = rad2deg(oi.euldist.oj)
                              self%dmat_eulers(j,i) = self%dmat_eulers(i,j)
                          case(3,5)
                              call oi%set_euler(self%flies(i)%vec(1:3))
                              call oj%set_euler(self%flies(j)%vec(1:3))
                              self%dmat_eulers(i,j) = rad2deg(oi.inpldist.oj)
                              self%dmat_eulers(j,i) = self%dmat_eulers(i,j)
                          case DEFAULT
                              stop 'invalid dimension; calc_dmat_euler; firefly_swarm_minimize'
                      end select
                  end do
              end do
              call oi%kill
              call oj%kill
          end subroutine calc_dmat_euler
          
          !> \brief  dynamic update of the Euler distance matrix for mode=euler
          subroutine update_dmat_euler( fly_moved )
              use simple_ori,  only: ori
              use simple_math, only: rad2deg
              integer, intent(in) :: fly_moved
              integer :: i
              type(ori) :: oi, omoved
              oi = ori()
              omoved = ori()
              do i=1,spec%npop
                  if( i /= fly_moved )then
                      select case(spec%ndim)
                          case(2)
                              call oi%e1set(self%flies(i)%vec(1))
                              call omoved%e1set(self%flies(fly_moved)%vec(1))
                              call oi%e2set(self%flies(i)%vec(2))
                              call omoved%e2set(self%flies(fly_moved)%vec(2))
                              self%dmat_eulers(i,fly_moved) = rad2deg(oi.euldist.omoved)
                              self%dmat_eulers(fly_moved,i) = self%dmat_eulers(i,fly_moved)
                          case(3,5)
                              call oi%set_euler(self%flies(i)%vec(1:3))
                              call omoved%set_euler(self%flies(fly_moved)%vec(1:3))
                              self%dmat_eulers(i,fly_moved) = rad2deg(oi.inpldist.omoved)
                              self%dmat_eulers(fly_moved,i) = self%dmat_eulers(i,fly_moved)
                          case DEFAULT
                              stop 'invalid dimension; calc_dmat_euler; firefly_swarm_minimize'
                      end select
                  endif
              end do
              call oi%kill
              call omoved%kill
          end subroutine update_dmat_euler
          
          !> \brief  calculates the origin shift distance matrix for mode=euler
          subroutine calc_dmat_shifts
              use simple_math, only: euclid
              integer :: i, j
              self%dmat_shifts = 0.
              do i=1,spec%npop-1
                  do j=i+1,spec%npop
                      self%dmat_shifts(i,j) = euclid(self%flies(i)%vec(4:5),self%flies(j)%vec(4:5))
                      self%dmat_shifts(j,i) = self%dmat_shifts(i,j)
                  end do
              end do
          end subroutine calc_dmat_shifts
          
          !> \brief  dynamic update of the shift distance matrix for mode=euler
          subroutine update_dmat_shifts( fly_moved )
              use simple_math, only: euclid
              integer, intent(in) :: fly_moved
              integer :: i
              do i=1,spec%npop
                  if( i /= fly_moved )then
                      self%dmat_shifts(i,fly_moved) = euclid(self%flies(i)%vec(4:5),self%flies(fly_moved)%vec(4:5))
                      self%dmat_shifts(fly_moved,i) = self%dmat_shifts(i,fly_moved)
                  endif
              end do
          end subroutine update_dmat_shifts
          
          !> \brief  moves the population of fireflies according to the attractiveness
          subroutine move_flies
              use simple_ori, only: ori
              integer   :: i, j, k, fly2move, fly2keep
              real      :: beta, ranperm, brightness, y, old_fly(spec%ndim)
              type(ori) :: o
              logical   :: outofrange
              do i=1,spec%npop-1
                  do j=i+1,spec%npop
                      ! select which fly to move
                      if( self%flies(i)%brightness > self%flies(j)%brightness )then
                          fly2move = j
                          fly2keep = i
                      else
                          fly2move = i
                          fly2keep = j
                      endif
                      ! stored the olf fly
                      old_fly = self%flies(fly2move)%vec
                      ! execute the move
                      select case(trim(spec%str_mode))
                          case('euler')
                              o = ori()
                              outofrange = .false.
                              do k=1,spec%ndim
                                  ! calculate random permutation
                                  ranperm = self%alphas(k)*(ran3()-0.5)
                                  if( i <= 3 )then
                                      ! calculate attraction using dmat_eulers
                                      beta = beta0*exp(-self%gammas(k)*self%dmat_eulers(i,j)*self%dmat_eulers(i,j))
                                      self%flies(fly2move)%vec(k) = self%flies(fly2move)%vec(k)+&
                                      beta*(self%flies(fly2keep)%vec(k)-self%flies(fly2move)%vec(k))+ranperm
                                      if( self%flies(fly2move)%vec(k) < spec%limits(k,1) ) outofrange = .true.
                                      if( self%flies(fly2move)%vec(k) > spec%limits(k,2) ) outofrange = .true.
                                  else if( i <= 5 )then
                                      ! calculate attraction using dmat_shifts
                                      ! remeber to enforce the shift constraint using the barrier method
                                      beta = beta0*exp(-self%gammas(k)*self%dmat_shifts(i,j)*self%dmat_shifts(i,j))
                                      self%flies(fly2move)%vec(k) = self%flies(fly2move)%vec(k)+&
                                      beta*(self%flies(fly2keep)%vec(k)-self%flies(fly2move)%vec(k))+ranperm
                                  else
                                      stop 'invalid dimension for mode euler; move_flies; firefly_swarm_minimize'
                                  endif
                              end do
                              if( outofrange )then ! shift the Euler angle
                                  call o%set_euler(self%flies(fly2move)%vec(1:3))
                                  self%flies(fly2move)%vec(1:3) = o%get_euler()
                              endif
                              ! update the distance matrices according to the move
                              call update_dmat_euler(fly2move)
                              if( spec%ndim > 3 ) call update_dmat_shifts(fly2move)
                              call o%kill
                          case('euclid')
                              do k=1,spec%ndim
                                  ! calculate random permutation
                                  ranperm = self%alphas(k)*(ran3()-0.5)
                                  ! calculate attraction using dmat_eulers
                                  beta = beta0*exp(-self%gammas(k)*self%dmat_eulers(i,j)*self%dmat_eulers(i,j))
                                  self%flies(fly2move)%vec(k) = self%flies(fly2move)%vec(k)+&
                                  beta*(self%flies(fly2keep)%vec(k)-self%flies(fly2move)%vec(k))+ranperm
                              end do
                              ! update the distance matrix according to the move
                              call update_dmat(fly2move)
                      end select
                      ! calculate the new brightness
                      brightness = -spec%costfun(self%flies(fly2move)%vec, spec%ndim)
                      spec%nevals = spec%nevals+1
                      if( brightness > self%flies(fly2move)%brightness )then
                          ! set the fly in line minimizer
                          self%spec_linmin%x = self%flies(fly2move)%vec
                          ! determine the direction of improvement
                          self%spec_linmin%xi = old_fly-self%spec_linmin%x
                          ! perform line minimization along this direction
                          self%spec_linmin%nevals = 0
                          call linmin(self%spec_linmin,y)
                          spec%nevals = spec%nevals+self%spec_linmin%nevals
                          self%flies(fly2move)%brightness = -y
                          self%flies(fly2move)%vec = self%spec_linmin%x
                          if( spec%npeaks > 0 )then
                              ! we will output local optima
                              if( y < spec%yb )then
                                  npeaks = npeaks+1
                                  spec%peaks(npeaks,:spec%ndim)  = self%spec_linmin%x
                                  spec%peaks(npeaks,spec%ndim+1) = y
                              endif
                          endif
                      else
                          self%flies(fly2move)%brightness = brightness
                      endif
                      spec%nevals = spec%nevals+1
                  end do
              end do
          end subroutine move_flies
          
          !> \brief  annealing
          subroutine update_alphas
              self%alphas = self%alphas*(0.9**t)
          end subroutine update_alphas
          
          !> \brief  to find brightest and darkest flies
          subroutine find_bright_and_dark
              integer :: i
              real :: brightness, darkness, x
              bright     = 1
              dark       = 1
              brightness = -huge(x)
              darkness   =  huge(x)
              do i=1,spec%npop
                  if( self%flies(i)%brightness > brightness )then
                      brightness = self%flies(i)%brightness
                      bright  = i
                  endif
                  if( self%flies(i)%brightness < darkness )then
                      darkness = self%flies(i)%brightness
                      dark  = i
                  endif
              end do
          end subroutine
          
          !> \brief  to generate output
          subroutine generate_ouput
              use simple_jiffys, only: alloc_err
              integer :: alloc_stat, i
              if( allocated(spec%population) ) deallocate(spec%population)
              allocate(spec%population(spec%npop,spec%ndim+1), stat=alloc_stat)
              call alloc_err("In: generate_output; firefly_swarm_minimize; simple_firefly_swarm_opt", alloc_stat)
              do i=1,spec%npop
                  spec%population(i,1:spec%ndim) = self%flies(i)%vec
                  spec%population(i,spec%ndim+1) = -self%flies(i)%brightness
              end do
              spec%x      = self%flies(bright)%vec
              lowest_cost = -self%flies(bright)%brightness
          end subroutine
          
          subroutine print_flies
              integer :: i
              write(*,*) '*********** PRINTING FLIES ***********'
              do i=1,spec%npop
                  write(*,*) self%flies(i)%vec, self%flies(i)%brightness
              end do
          end subroutine print_flies
          
          subroutine print_dists
              integer :: i, j 
              write(*,*) '*********** PRINTING DISTS ***********'
              do i=1,spec%npop-1
                  do j=i+1,spec%npop
                      write(*,*) i, j, self%dmat_eulers(i,j), self%dmat_eulers(j,i)
                  end do 
                  
              end do
          end subroutine print_dists
        
    end subroutine firefly_swarm_minimize
    
    !> \brief  is a destructor
    subroutine kill_firefly_swarm( self )
        class(firefly_swarm_opt), intent(inout) :: self !< instance
        integer :: i
        if( self%exists )then
            do i=1,size(self%flies)
                deallocate(self%flies(i)%vec)
            end do
            deallocate(self%flies, self%dmat_eulers,&
            self%dmat_shifts, self%alphas, self%gammas)
            call self%spec_linmin%kill
            self%exists = .false.
        endif
    end subroutine kill_firefly_swarm
    
    subroutine test_firefly_swarm_opt
        use simple_testfuns
        type(firefly_swarm_opt) :: fso
        type(opt_spec) :: spec
        real :: gmin, range(2), lims(2,2), lowest_cost
        procedure(testfun), pointer :: costfun_ptr
        costfun_ptr = get_testfun(1, 2, gmin, range)
        lims(:,1) = range(1)
        lims(:,2) = range(2)
        call spec%specify('firefly', 2, 'euclid', limits=lims, verbose=.true.)
        call spec%set_costfun(costfun_ptr)
        call fso%new(spec)
        call fso%minimize(spec, lowest_cost)
        write(*,*) 'spec%x: ', spec%x
        write(*,*) 'lowest_cost: ', lowest_cost
    end subroutine

end module simple_firefly_swarm_opt