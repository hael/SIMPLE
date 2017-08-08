! continous optimisation of polar Fourier cross-correlation using differential evolution
module simple_de_pftcc_opt
use simple_pftcc_opt, only: pftcc_opt
use simple_opt_spec,  only: opt_spec
use simple_rnd,       only: ran3, irnd_uni
implicit none

public :: de_pftcc_opt
private

integer, parameter :: N_GENERAL = 136,    N_VALLEY = 126,    N_MULTIMODAL = 103,    N_FLAT = 106
real,    parameter :: F_GENERAL = 0.2790, F_VALLEY = 0.4027, F_MULTIMODAL = 0.3976, F_FLAT = 0.5860
real,    parameter :: X_GENERAL = 0.9813, X_VALLEY = 0.9211, X_MULTIMODAL = 0.9794, X_FLAT = 0.3345
logical, parameter :: DEBUG = .false.

type :: de_pftcc_opt
    private
    real, allocatable :: pop(:,:)          !< solution vector population
    real, allocatable :: costs(:)          !< costs
    integer           :: best              !< index of best population member
    real              :: F  = F_MULTIMODAL !< amplification factor
    real              :: CR = X_MULTIMODAL !< cross-over rate
    logical           :: exists = .false.  !< to indicate existence
  contains
    procedure          :: new          => new_de
    procedure, private :: init_pop
    procedure          :: minimize     => de_minimize
    procedure          :: kill         => kill_de
end type de_pftcc_opt

contains
    
    !> \brief  is a constructor
    subroutine new_de( self, spec )
        use simple_jiffys, only: alloc_err
        class(de_pftcc_opt), intent(inout) :: self !< instance
        class(opt_spec),     intent(inout) :: spec !< specification
        integer :: alloc_stat
        ! destruct if exists
        call self%kill
        ! adjust control parameters according to mode
        select case(trim(spec%str_mode))
            case('general')
                spec%npop = N_GENERAL
                self%F    = F_GENERAL
                self%CR   = X_GENERAL
            case('valley')
                spec%npop = N_VALLEY
                self%F    = F_VALLEY
                self%CR   = X_VALLEY
            case('multimodal')
                spec%npop = N_MULTIMODAL
                self%F    = F_MULTIMODAL
                self%CR   = X_MULTIMODAL
            case('flat')
                spec%npop = N_FLAT
                self%F    = F_FLAT
                self%CR   = X_FLAT
            case DEFAULT
                spec%npop = N_MULTIMODAL
                self%F    = F_MULTIMODAL
                self%CR   = X_MULTIMODAL
        end select
        ! allocate
        allocate(self%pop(spec%npop,spec%ndim), self%costs(spec%npop), stat=alloc_stat)
        call alloc_err("In: new_de; simple_de_opt", alloc_stat)
        self%exists = .true. ! indicates existence
        if( spec%DEBUG ) write(*,*) 'created new differential evolution population'
    end subroutine new_de

    !> \brief  is the particle swarm population initialization routine
    subroutine init_pop( self, ospec )
        class(de_pftcc_opt), intent(inout) :: self        !< instance     
        class(opt_spec),     intent(inout) :: ospec        !< specification
        integer :: i, n_ini
        ! obtain initial solutions by randomized bounds
        do i = 1, ospec%npop
            self%pop(i,:) = gen_individual()
        end do
        if( allocated(ospec%inipopulation) )then
            ! provided init population
            n_ini = size(ospec%inipopulation, dim=1)
            self%pop(1:n_ini,:) = ospec%inipopulation
        endif
        ! set one point in the swarm to best point in spec
        if( .not. all(ospec%x == 0.) )self%pop(1,:) = ospec%x
        contains

            function gen_individual( )result( individual )
                real    :: L, individual(ospec%ndim)
                integer :: j
                individual = 0.
                do j = 1, ospec%ndim
                    L = ospec%limits(j,2) - ospec%limits(j,1)
                    individual(j) = ospec%limits(j,1) + L*ran3()
                end do
            end function gen_individual
    end subroutine init_pop

    !> \brief  is the particle swarm minimize minimization routine
    subroutine de_minimize( self, spec, funcontainer, lowest_cost )
        use simple_math, only: hpsort
        class(de_pftcc_opt), intent(inout) :: self        !< instance     
        class(opt_spec),     intent(inout) :: spec        !< specification
        class(pftcc_opt),    intent(inout) :: funcontainer !< container for the cost function
        real,                intent(out)   :: lowest_cost !< lowest cost
        integer, allocatable :: inds(:)
        integer :: t      ! iteration counter
        integer :: npeaks ! number of local minima
        integer :: loc(1), nworse, X, i
        ! initialization
        call self%init_pop(spec)
        ! calculate initial costs
        do i=1,spec%npop
            self%costs(i) = funcontainer%costfun(self%pop(i,:), spec%ndim)
            spec%nevals   = spec%nevals+1
        end do
        loc       = minloc(self%costs)
        self%best = loc(1)
        ! initialize the remains
        spec%nevals = 0
        npeaks      = 0
        nworse      = 0
        do t=1,spec%maxits ! generations loop
            ! select solution to modify
            X  = irnd_uni(spec%npop)
            call update_agent( X )
            if( spec%npeaks > 0 )then
                ! exit when we have identified spec%npeaks local optima
                if( npeaks == spec%npeaks ) exit
            endif
            if( nworse == spec%npop ) exit ! if no better solutions have been found in the last 4*ndim iterations
        end do
        ! report
        lowest_cost = self%costs(self%best)
        spec%x = self%pop(self%best,:)
        if( spec%npeaks == 0 )then
            ! peaks <- population
            if(allocated(spec%peaks))deallocate(spec%peaks)
            allocate(spec%peaks(spec%npop, 4), source=1.)
            spec%peaks(:,  4) = self%costs 
            spec%peaks(:,1:3) = self%pop
            where(spec%peaks(:,4) >= 0.) spec%peaks(:,4) = 1.
        else
            if( npeaks == spec%npeaks )then
                ! alles klar
            else
                ! reports peaks as top-ranking individuals
                if(allocated(spec%peaks))deallocate(spec%peaks)
                allocate(spec%peaks(spec%npeaks, 4), source=1.)
                allocate(inds(spec%npop))
                inds  = (/ (i, i=1,spec%npop) /)
                call hpsort(spec%npeaks, self%costs, inds)
                spec%peaks(:,  4) = self%costs(inds(spec%npop-spec%npeaks+1:)) 
                spec%peaks(:,1:3) = self%pop(inds(spec%npop-spec%npeaks+1:),:)
                where(spec%peaks(:,4) >= 0.) spec%peaks(:,4) = 1.
                deallocate(inds)
            endif
        endif

        contains

            subroutine update_agent( X )
                integer, intent(in) :: X
                integer :: a, rb, b, i
                real :: trial(spec%ndim), cost_trial
                ! select random disjoint pair
                a  = irnd_uni(spec%npop)
                rb = irnd_uni(spec%npop-1)
                b  = a+rb 
                if( b <= spec%npop )then
                else
                    b = a+rb-spec%npop
                endif
                ! create a trial solution
                do i=1,spec%ndim
                    if( i == X .or. ran3() < self%CR )then
                        trial(i) = self%pop(self%best,i)+self%F*(self%pop(a,i)-self%pop(b,i))
                        ! enforce limits 
                        trial(i) = min(spec%limits(i,2),trial(i))
                        trial(i) = max(spec%limits(i,1),trial(i))
                    else
                        trial(i) = self%pop(X,i) 
                    endif
                end do
                ! calculate cost 
                cost_trial  = funcontainer%costfun(trial, spec%ndim)
                spec%nevals = spec%nevals + 1
                ! update pop if better solution is found
                if( cost_trial < self%costs(X) )then
                    nworse = 0
                    self%pop(X,:) = trial
                    self%costs(X) = cost_trial
                    ! update global best if needed
                    if( cost_trial < self%costs(self%best) ) self%best = X
                    if( spec%npeaks > 0 )then
                        ! we will output local optima
                        if( cost_trial < spec%yb )then
                            npeaks = npeaks+1
                            spec%peaks(npeaks,:spec%ndim)  = trial
                            spec%peaks(npeaks,spec%ndim+1) = cost_trial
                        endif
                    endif
                else
                    nworse = nworse + 1
                endif
            end subroutine update_agent
            
    end subroutine de_minimize

    ! GETTERS

    !> \brief  is a destructor
    subroutine kill_de( self )
        class(de_pftcc_opt), intent(inout) :: self !< instance
        if( self%exists )then
            deallocate(self%pop,self%costs)
            self%exists = .false.
        endif
    end subroutine kill_de

end module simple_de_pftcc_opt
