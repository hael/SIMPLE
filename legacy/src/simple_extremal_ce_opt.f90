module simple_extremal_ce_opt
use simple_optimizer,  only: optimizer
use simple_opt_spec,   only: opt_spec
use simple_rnd,        only: ran3, gasdev
use simple_online_var, only: online_var
use simple_defs        ! singleton
implicit none

public :: extremal_ce_opt
private

type, extends(optimizer) :: extremal_ce_opt
    private
    real, allocatable :: pop(:,:)         !< solution population
!     real, allocatable :: means(:)         !< means of probabilistic model
    real, allocatable :: sdevs(:)         !< standard deviations of probabilistic model
    type(online_var), allocatable :: vars(:)   !< objects for online variance estimation (per variable)
    real              :: yw               !< worst cost
    real              :: yb               !< best cost 
    integer           :: worst            !< index of best particle position
    integer           :: best             !< index of worst particle position
    logical           :: exists = .false. !< to indicate existence
  contains
    procedure :: new      => new_extremal_ce
    procedure :: minimize => extremal_ce_minimize
    procedure :: kill     => kill_extremal_ce
end type extremal_ce_opt

logical, parameter :: debug=.false.

contains
    
    !> \brief  is a constructor
    subroutine new_extremal_ce( self, spec )
        use simple_jiffys,   only: alloc_err
        class(extremal_ce_opt), intent(inout) :: self !< instance
        class(opt_spec), intent(inout)        :: spec !< specification
        integer :: alloc_stat
        ! destruct if exists
        call self%kill
        ! allocate
        allocate(self%pop(spec%npop,spec%ndim),&
        self%sdevs(spec%ndim), self%vars(spec%ndim), stat=alloc_stat) !self%means(spec%ndim)
        call alloc_err("In: new_extremal_ce_opt, 1", alloc_stat)
        self%exists = .true. ! indicates existence
        if( spec%debug ) write(*,*) 'created new extremal_ce obj'
    end subroutine new_extremal_ce
    
    !> \brief  is the extremal cross-entropy routine
    subroutine extremal_ce_minimize( self, spec, lowest_cost )
        use simple_opt_subs, only: linmin
        class(extremal_ce_opt), intent(inout) :: self        !< instance
        class(opt_spec), intent(inout)        :: spec        !< specification
        real, intent(out)                     :: lowest_cost !< lowest cost
        integer :: t                ! iteration counter
        integer :: npeaks           ! number of local minima
        integer :: nworse           ! number of iterations that did not improve the solution pop
        real    :: rtol             ! relative tolerance
        real    :: costs(spec%npop) ! particle costs
        integer :: loc(1), nreplaced
        logical :: err = .false. ! if true variance is zero and the process has converged 
        if( .not. associated(spec%costfun) )then
            stop 'cost function not associated in opt_spec; extremal_ce_minimize; simple_extremal_ce_opt'
        endif
        ! initialize
        call init
        nreplaced   = 0
        spec%nevals = 0
        npeaks      = 0
        do t=1,spec%maxits ! generations loop
            ! generate and evaluate sample from the probabilistic model
            call gen_and_eval_sample
!             if( nreplaced == spec%npop )then
!                 ! update the model
!                 call update_model
!             endif
            if( spec%npeaks > 0 )then
                ! exit when we have identified spec%npeaks local optima
                if( npeaks == spec%npeaks ) exit
            endif
            if( nworse == 4*spec%npop) exit ! if no better solutions have been found in the last 4*npop iterations
            if( err ) exit ! because the variance is zero
        end do
        ! output best
        loc         = minloc(costs)
        self%best   = loc(1)
        lowest_cost = costs(self%best)
        spec%x      = self%pop(self%best,:)
        
      contains
    
            !> \brief  initialize the solution population and the probabilistic model
            subroutine init
                integer :: i, j
                real :: L
                ! obtain particle positions by randomized bounds
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
                ! calculate initial costs
                do i=1,spec%npop
                    costs(i) = spec%costfun(self%pop(i,:), spec%ndim)
                    spec%nevals = spec%nevals+1
                end do
                ! set best/worst
                loc = minloc(costs)
                self%best = loc(1)
                self%yb   = costs(self%best)
                loc = maxloc(costs)
                self%worst = loc(1)
                self%yw    = costs(self%worst)
                ! init sdevs & online model
                do i=1,spec%ndim
                    self%sdevs(i) = (spec%limits(i,2)-spec%limits(i,1))/2.
                    do j=1,spec%npop
                        call self%vars(i)%add(self%pop(j,i))
                    end do
                end do
                ! init means to best solution in pop
!                 do i=1,spec%ndim
!                     self%means(i) = self%pop(self%best,i)
!                 end do      
            end subroutine init
    
            !> \brief  generate and evaluate sample from the probabilistic model
            subroutine gen_and_eval_sample
                integer :: i
                real :: trial(spec%ndim), cost_trial, var
                ! sample the model
                do i=1,spec%ndim
!                     trial(i) = gasdev(self%means(i), self%sdevs(i))
                    trial(i) = gasdev(self%vars(i)%get_mean(), self%sdevs(i))
                end do
                ! calculate the cost of the trial solution
                cost_trial = spec%costfun(trial, spec%ndim)
                spec%nevals = spec%nevals+1
                if( cost_trial < self%yw )then
                    nworse = 0
                    ! add to probabilistic model
                    do i=1,spec%ndim
                        call self%vars(i)%add(trial(i))
                    end do
                    ! update the standard deviations
                    do i=1,spec%ndim
                        var = self%vars(i)%get_var()
                        if( var /= 0. )then
                            self%sdevs(i) = max(spec%ftol,sqrt(var))
                        endif
                    end do
                    ! replace the worst
                    costs(self%worst) = cost_trial
                    self%pop(self%worst,:) = trial
                    ! update worst
                    loc = maxloc(costs)
                    self%worst = loc(1)
                    self%yw    = costs(self%worst)
                    ! update replacement counter
                    nreplaced = nreplaced+1
                else
                    nworse = nworse+1
                endif
            end subroutine gen_and_eval_sample
    
            !> \brief  update the probabilistic model
            !!         may have to implement a learning rate for the parameters of the model
            !!         However, it might not be needed since the process is extremal
!             subroutine update_model
!                 use simple_stat, only: moment
!                 integer :: i
!                 real    :: var, prev_best(spec%ndim), y
!                 ! set best
!                 loc = minloc(costs)
!                 self%best = loc(1)
!                 self%yb   = costs(self%best)
!                 if( spec%npeaks > 0 )then
!                     ! we will output local optima
!                     if( self%yb < spec%yb )then
!                         npeaks = npeaks+1
!                         spec%peaks(npeaks,:spec%ndim)  = self%pop(self%best,:)
!                         spec%peaks(npeaks,spec%ndim+1) = self%yb
!                     endif
!                 endif
!                 ! straight moment calc
!                 do i=1,spec%ndim
!                     call moment(self%pop(:,i), self%means(i), self%sdevs(i), var, err)
!                 end do
!                 ! zero nreplaced
!                 nreplaced = 0
!             end subroutine
            
        end subroutine extremal_ce_minimize
        
        !> \brief  is a destructor
        subroutine kill_extremal_ce( self )
            class(extremal_ce_opt), intent(inout) :: self !< instance
            if( self%exists )then
                deallocate(self%pop, self%sdevs, self%vars) ! self%means
                self%exists = .false.
            endif
        end subroutine kill_extremal_ce

end module simple_extremal_ce_opt   