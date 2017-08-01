!==Class simple_em_opt
!
! Minimization of an externally defined function by an electromagnetism-like mechanism for
! global optimisation. The code is distributed with the hope that it will be useful, but 
! _WITHOUT_ _ANY_ _WARRANTY_.  Redistribution or modification is regulated by the GNU General 
! Public License. 
! *Author:* Hans Elmlund, 2015-09-29
module simple_em_opt
use simple_optimizer,  only: optimizer
use simple_online_var, only: online_var
use simple_opt_spec,   only: opt_spec
use simple_defs        ! singleton
implicit none

public :: em_opt, test_em_opt
private

type, extends(optimizer) :: em_opt
    private
    real, allocatable :: pop(:,:)       !< solution vector population
    real, allocatable :: costs(:)       !< costs
    real, allocatable :: F(:,:)         !< forces
    real              :: yb=0.          !< best (lowest) cost function val
    integer           :: best           !< index of best
    logical           :: exists=.false. !< to indicate existence
  contains
    procedure :: new      => new_em_opt
    procedure :: minimize => em_minimize
    procedure :: kill     => kill_em_opt
end type

contains

    !> \brief  is a constructor
    subroutine new_em_opt( self, spec )
        use simple_jiffys, only: alloc_err
        class(em_opt), intent(inout)   :: self !< instance
        class(opt_spec), intent(inout) :: spec !< specification
        integer                        :: alloc_stat
        real                           :: x
        call self%kill
        ! allocate & initialize
        allocate(self%pop(spec%npop,spec%ndim), self%costs(spec%npop), self%F(spec%npop,spec%ndim), stat=alloc_stat)
        call alloc_err("In: new_em_opt", alloc_stat)
        self%yb = huge(x)    ! initialize best cost to huge number
        self%exists = .true. ! indicates existence
        if( spec%debug ) write(*,*) 'created em_opt object'
    end subroutine new_em_opt
    
    !> \brief  em minimization
    subroutine em_minimize( self, spec, lowest_cost )
        use simple_rnd, only: gasdev, seed_rnd, ran3
        class(em_opt), intent(inout)   :: self        !< instance
        class(opt_spec), intent(inout) :: spec        !< specification
        real, intent(out)              :: lowest_cost !< lowest cost
        integer                        :: i, jiter, nbetter, conv, niter, npeaks
        if( .not. associated(spec%costfun) )then
            stop 'cost function not associated in opt_spec; em_minimize; simple_em_opt'
        endif
        call init
        call update_best
        ! initialize counters
        conv  = 0
        niter = 0
        if( spec%debug ) write(*,*) 'finalized initialization'
        do jiter=1,spec%maxits  ! iterations loop
            if( spec%verbose )then
                if( jiter == 1 .or. mod(jiter,30) == 0 ) write(*,'(a,1x,i5,1x,a,1x,f12.4)' )&
                'Iteration:', jiter, 'Cost:', self%yb
            endif
            niter = niter+1 ! 4 counting the actual nr of iterations run
            nbetter = 0     ! 4 counting the number of improving solutions
            ! search
            call greedy_local_search
            call calcF
            call moveF
            call update_best
            if( nbetter == 0 )then
                conv = conv+1
            else
                conv = 0
            endif
            if( spec%npeaks > 0 )then
                ! exit when we have identified spec%npeaks local optima
                if( npeaks == spec%npeaks ) exit
            endif
            if( conv == 4*spec%ndim ) exit ! if no better solutions have been found in the last 4*ndim iterations
        end do
        spec%niter = spec%niter+niter
        ! report nr of function evals, nr of iterations, and lowest cost
        spec%nevals = spec%nevals
        spec%niter  = spec%niter
        lowest_cost = self%yb
        
        contains
            
            !> \brief  initialize the population & set best
            subroutine init
                integer :: i, j
                real :: L
                ! initialize counters
                spec%niter = 0
                npeaks = 0
                ! obtain initial solutions by randomized bounds
                do i=1,spec%npop
                    do j=1,spec%ndim
                        L = spec%limits(j,2)-spec%limits(j,1)
                        self%pop(i,j) = spec%limits(j,1)+ran3()*L
                    end do
                    if( i == 1 )then
                        if( .not. all(spec%x == 0.) )then ! set one point in the pop to best point in spec
                            self%pop(1,:)= spec%x
                        endif
                    endif
                end do
            end subroutine init
            
            !> \brief  update the best after move
            subroutine update_best
                integer :: loc(1)
                ! calculate costs
                do i=1,spec%npop
                    self%costs(i) = spec%costfun(self%pop(i,:), spec%ndim)
                    spec%nevals = spec%nevals+1
                end do
                ! update best
                loc        = minloc(self%costs)
                if( loc(1) /= self%best )then
                    nbetter = nbetter+1
                    npeaks  = npeaks+1
                endif
                self%best  = loc(1)
                self%yb    = self%costs(self%best)
                spec%x     = self%pop(self%best,:)
            end subroutine update_best
            
            !> \brief  a greedy local search procedure for improving best
            subroutine greedy_local_search
                integer :: cnt, k
                real    :: len, lambda1, lambda2, y(spec%ndim), cost
                cnt = 1
                len = spec%delta*maxval(spec%limits(:,2)-spec%limits(:,1))
                do k=1,spec%ndim
                    lambda1 = ran3()
                    do while(cnt < spec%nsample)
                        y = self%pop(self%best,:)
                        lambda2 = ran3()
                        if( lambda1 > 0.5 )then
                            y(k) = y(k)+lambda2*len
                        else
                            y(k) = y(k)-lambda2*len
                        endif
                        cost = spec%costfun(y, spec%ndim)
                        if(cost < self%yb)then
                            self%yb = cost
                            self%pop(self%best,:) = y
                            cnt = spec%nsample-1
                            nbetter = nbetter+1
                            npeaks  = npeaks+1
                            if( spec%npeaks > 0 )then
                                if( npeaks == spec%npeaks ) return
                            endif
                        endif
                        cnt = cnt+1
                    end do
                end do
            end subroutine greedy_local_search
            
            !> \brief  calculate the force field
            subroutine calcF
                real    :: sdiff, dist
                integer :: i, j
                real    :: q(spec%npop), diffvec(spec%ndim)
                ! accumulate sum of differences
                sdiff = sum(self%costs-self%costs(self%best))
                ! calculate charges
                q = exp(-spec%ndim*(self%costs-self%costs(self%best))/sdiff)
                ! initialize forces
                self%F = 0.
                ! calculate forces
                do i=1,spec%npop
                    do j=1,spec%npop
                        diffvec = self%pop(j,:)-self%pop(i,:)
                        dist    = sqrt(sum(diffvec**2.))
                        if( self%costs(j) < self%costs(i) )then
                            ! attraction
                            self%F(i,:) = self%F(i,:)+diffvec*((q(i)*q(j))/dist)
                        else
                            ! repulsion
                            self%F(i,:) = self%F(i,:)-diffvec*((q(i)*q(j))/dist)
                        endif
                    end do
                end do
            end subroutine calcF
            
            !> \brief  for moving the field
            subroutine moveF
                use simple_math, only: normvec
                integer :: i, k
                real    :: lambda
                do i=1,spec%npop
                    if( i /= self%best )then
                        lambda = ran3()
                        call normvec(self%F(i,:))
                        do k=1,spec%ndim
                            if( self%F(i,k) > 0. )then
                                self%pop(i,k) = self%pop(i,k)+lambda*self%F(i,k)*(spec%limits(k,2)-self%pop(i,k))
                            else
                                self%pop(i,k) = self%pop(i,k)+lambda*self%F(i,k)*(self%pop(i,k)-spec%limits(k,1))
                            endif
                        end do    
                    endif
                end do
            end subroutine moveF
            
    end subroutine
    
    !> \brief  is a destructor
    subroutine kill_em_opt( self )
        class(em_opt), intent(inout) :: self
        if(self%exists)then
            deallocate(self%pop,self%costs,self%F)
            self%exists = .false.
        endif
    end subroutine
    
    !> \brief  unit test
    subroutine test_em_opt
        use simple_rnd,      only: ran3
        use simple_math,     only: euclid
        use simple_testfuns
        type(opt_spec)              :: ospec
        type(em_opt)                :: illusion
        procedure(testfun), pointer :: costfun_ptr  !< pointer 2 test function
        integer, parameter          :: TSTFUN=18,DIMS=2
        real                        :: limits(DIMS,2),lowest_cost, gmin, range(2)
        integer                     :: i
        ! generate the cost function (here: 1)
        costfun_ptr = get_testfun(TSTFUN, DIMS, gmin, range) ! get testfun, gmin is the global min, range is limits
        limits(:,1) = range(1)
        limits(:,2) = range(2)
        ! specify the optimizer
        call ospec%specify('em',DIMS,limits=limits,nrestarts=1) ! make optimizer spec
        call ospec%set_costfun(costfun_ptr)                     ! set pointer to costfun
        ! initialize with randomized bounds
        do i=1,DIMS
            ospec%x(i) = limits(i,1)+ran3()*(limits(i,2)-limits(i,1))
        end do
        ! generate the optimizer
        call illusion%new(ospec)
        call illusion%minimize(ospec,lowest_cost)
        ! minimize
        print *, 'minimum found:', lowest_cost
        print *, 'correct global minimum:', gmin
        print *, 'solution:', ospec%x
    end subroutine
   
end module simple_em_opt
