!==Class simple_oasis_opt
!
! Minimization of an externally defined function by the OASIS method of Hans Elmlund.
! OASIS stands for Optimization by Adaptive Sampling around Important Solutions.
! The code is distributed with the hope that it will be useful, but _WITHOUT_ _ANY_ _WARRANTY_.
! Redistribution or modification is regulated by the GNU General Public License.
! *Author:* Hans Elmlund, 2013-10-15
module simple_oasis_opt
use simple_defs 
use simple_optimizer,   only: optimizer
use simple_online_var,  only: online_var
use simple_opt_spec,    only: opt_spec
implicit none

public :: oasis_opt!, test_oasis_opt
private

type, extends(optimizer) :: oasis_opt
    private
    type(online_var), allocatable :: vars(:)        !< objects for online variance estimation (per variable)
    type(opt_spec)                :: spec_linmin    !< specification of linmin optimizer
    real, allocatable             :: sdevs(:)       !< standard devations of the variables
    real                          :: yb=0.          !< best (lowest) cost function val
    logical                       :: exists=.false. !< to indicate existence
  contains
    procedure :: new      => new_oasis_opt
    procedure :: minimize => oasis_minimize
    procedure :: kill     => kill_oasis_opt
end type

contains

    !> \brief  is a constructor
    subroutine new_oasis_opt( self, spec )
        use simple_jiffys, only: alloc_err
        class(oasis_opt), intent(inout) :: self !< instance
        class(opt_spec), intent(inout)  :: spec !< specification
        integer                         :: alloc_stat
        real                            :: x
        call self%kill
        ! allocate & initialize
        allocate(self%vars(spec%ndim),self%sdevs(spec%ndim), stat=alloc_stat)
        call alloc_err("In: new_oasis_opt", alloc_stat)
        self%sdevs  = 0.
        self%yb     = huge(x) ! initialize best cost to huge number
        ! make line minimizer
        call self%spec_linmin%specify('linmin', spec%ndim, maxits=spec%nsample, ftol=spec%ftol)
        call self%spec_linmin%set_costfun(spec%costfun)
        self%exists = .true.  ! indicates existence
        if( spec%debug ) write(*,*) 'created oasis object'
    end subroutine

    !> \brief  restarted oasis minimization
    subroutine oasis_minimize( self, spec, lowest_cost )
        use simple_rnd,      only: gasdev, seed_rnd, ran3
        use simple_opt_subs, only: linmin
        class(oasis_opt), intent(inout) :: self        !< instance
        class(opt_spec), intent(inout)  :: spec        !< specification
        real, intent(out)               :: lowest_cost !< lowest cost
        logical                         :: arezero(spec%ndim)
        integer                         :: i, j, k, nbetter, conv, niter
        real                            :: y
        if( .not. associated(spec%costfun) )then
            stop 'cost function not associated in opt_spec; oasis_minimize; simple_oasis_opt'
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
        ! test if sdevs in spec are set
        arezero = .false.
        do i=1,spec%ndim
            if( spec%sdevs(i) == 0. ) arezero(i) = .true.
        end do
        ! generate initial sdevs
        if( all(arezero) )then
            ! by setting all the sdevs to the half the interval
            do i=1,spec%ndim
                spec%sdevs(i) = (spec%limits(i,2)-spec%limits(i,1))/2.
            end do
        endif
        ! set best cost
        spec%nevals = 0
        self%yb     = spec%costfun(spec%x, spec%ndim)
        spec%nevals = spec%nevals+1
        if( spec%npeaks > 0 )then
            spec%peakcnt = 0
            spec%peakcnt = spec%peakcnt+1
            spec%peaks(spec%peakcnt,:spec%ndim)  = spec%x
            spec%peaks(spec%peakcnt,spec%ndim+1) = self%yb
        endif
        ! initialize counters
        spec%niter = 0
        ! run nrestarts
        do i=1,spec%nrestarts
            if( spec%verbose ) write(*,'(a,1x,i5)') '>>> RESTART ROUND:', i
            ! re-seed the random number generator
            call seed_rnd
            ! set sdevs to sdevs in spec
            self%sdevs = spec%sdevs
            ! zero the online_var objects and
            ! communicate the current best point
            do j=1,spec%ndim
                self%vars(j) = online_var()
                call self%vars(j)%add(spec%x(j))
            end do
            conv  = 0
            niter = 0
            if( spec%debug ) write(*,*) 'finalized initialization'
            do j=1,spec%maxits  ! iterations loop
                if( spec%verbose )then
                    if( j == 1 .or. mod(j,30) == 0 ) write(*,'(a,1x,i5,1x,a,1x,f12.4,1x,a,1x,f8.2)' )&
                    'Iteration:', j, 'Cost:', self%yb, 'Avg Sdev:', sum(self%sdevs)/real(spec%ndim)
                endif
                niter = niter+1   ! 4 counting the actual nr of iterations run
                nbetter = 0       ! 4 counting the number of improving solutions
                ! set best point in line minimizer to point in spec
                self%spec_linmin%x = spec%x
                ! OASIS search
                call oasis_srch   ! the sampling/scoring/online updating procedure
                if( nbetter > 0 )then
                    ! determine direction of improvement
                    self%spec_linmin%xi = spec%x-self%spec_linmin%x
                    ! perform line minimization along this direction
                    self%spec_linmin%nevals = 0
                    call linmin(self%spec_linmin,y)
                    spec%nevals = spec%nevals+self%spec_linmin%nevals
                    if( y < self%yb )then
                        self%yb = y
                        spec%x  = self%spec_linmin%x
                        ! add to probabilistic model
                        do k=1,spec%ndim
                            call self%vars(k)%add(spec%x(k))
                        end do
                    endif
                    if( spec%npeaks > 0 )then
                        ! we will output local optima
                        spec%peakcnt = spec%peakcnt+1
                        spec%peaks(spec%peakcnt,:spec%ndim)  = self%spec_linmin%x
                        spec%peaks(spec%peakcnt,spec%ndim+1) = self%yb
                    endif
                endif
                if( spec%debug ) write(*,*) 'finalized line minimization '
                if( nbetter == 0 )then
                    conv = conv+1
                else
                    conv = 0
                endif
                if( spec%npeaks > 0 )then
                    ! exit when we have identified spec%npeaks local optima
                    if( spec%peakcnt == spec%npeaks ) exit
                endif
                if( conv == 4*spec%ndim ) exit ! if no better solutions have been found in the last 4*ndim iterations
            end do
            spec%niter = spec%niter+niter
        end do
        ! report nr of function evals, nr of iterations, and lowest cost
        spec%nevals = spec%nevals/spec%nrestarts
        spec%niter  = spec%niter/spec%nrestarts
        lowest_cost = self%yb

        contains

            subroutine oasis_srch
                use simple_rnd, only: gasdev
                integer :: i, ii
                real    :: pi_old, y, var
                ! generate Gaussian deviates around the best point
                i = 0
                do ii=1,spec%nsample
                    i = i+1
                    if( i > spec%ndim ) i = 1
                    ! store the previous solution element
                    pi_old = spec%x(i)
                    ! change the component using Gaussian random sampling
                    spec%x(i) = gasdev(self%vars(i)%get_mean(), self%sdevs(i)) ! NO LIMITS, SINCE LIMITS GIVE INFINITE LOOP AT TIMES
                    ! score the new solution vector
                    y = spec%costfun(spec%x, spec%ndim)
                    spec%nevals = spec%nevals+1 ! increment the number of costfun evals counter
                    ! update the model if a better solution is found
                    if( y <= self%yb )then
                        nbetter    = nbetter+1           ! increasing the number-of-better-solutions-found-counter
                        self%yb    = y                   ! updating the best cost
                        call self%vars(i)%add(spec%x(i)) ! updating the online variance calculator
                        ! update the standard deviations
                        var = self%vars(i)%get_var()
                        if( var /= 0. )then
                            self%sdevs(i) = max(spec%ftol,sqrt(var))
                        endif
                    else
                        spec%x(i) = pi_old               ! put back the old point
                    endif
                end do
                if( spec%debug ) write(*,*) 'finalized oasis search'
            end subroutine

    end subroutine

    !> \brief  is a destructor
    subroutine kill_oasis_opt( self )
        class(oasis_opt), intent(inout) :: self
        if(self%exists)then
            deallocate(self%vars,self%sdevs)
            call self%spec_linmin%kill
            self%exists = .false.
        endif
    end subroutine

    !> \brief  unit test
    ! subroutine test_oasis_opt
    !     use simple_rnd,      only: ran3
    !     use simple_math,     only: euclid
    !     use simple_testfuns
    !     type(opt_spec)              :: ospec
    !     type(oasis_opt)             :: illusion
    !     procedure(testfun), pointer :: costfun_ptr  !< pointer 2 test function
    !     integer, parameter          :: TSTFUN=3,DIMS=8
    !     real                        :: limits(DIMS,2),lowest_cost, gmin, range(2)
    !     integer                     :: i
    !     ! generate the cost function (here: 1)
    !     costfun_ptr = get_testfun(TSTFUN, DIMS, gmin, range) ! get testfun, gmin is the global min, range is limits
    !     limits(:,1) = range(1)
    !     limits(:,2) = range(2)
    !     ! specify the optimizer
    !     call ospec%specify('oasis',DIMS,limits=limits,nrestarts=1) ! make optimizer spec
    !     call ospec%set_costfun(costfun_ptr)                        ! set pointer to costfun
    !     ! initialize with randomized bounds
    !     do i=1,DIMS
    !         ospec%x(i) = limits(i,1)+ran3()*(limits(i,2)-limits(i,1))
    !     end do
    !     ! generate the optimizer
    !     call illusion%new(ospec)
    !     call illusion%minimize(ospec,lowest_cost)
    !     ! minimize
    !     print *, 'minimum found:', lowest_cost
    !     print *, 'correct global minimum:', gmin
    !     print *, 'solution:', ospec%x
    ! end subroutine

end module simple_oasis_opt
