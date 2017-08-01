!==Class simple_de_opt
!
! constrained minimization of an externally defined function by advanced differential evolution optimization
! A ring-topological neighborhood structure balances between intensification and diversification of the search.
! The code is distributed with the hope that it will be useful, but _WITHOUT_ _ANY_ _WARRANTY_. 
! Redistribution or modification is regulated by the GNU General Public License. 
! *Author:* Hans Elmlund, 2009-05-25.
!
module simple_de_opt
use simple_optimizer, only: optimizer
use simple_ran_tabu,  only: ran_tabu
use simple_rnd,       only: ran3
use simple_defs       ! singleton
implicit none

public :: de_opt
private

!> \brief  is the individual containing, vector length, vector cost, self-adaptive
!!         amplification factor, vector values, best vec in neigh and its cost 
type :: individual
    real              :: cost=0.          !< vector cost
    real              :: F=0.             !< amplification factor
    real              :: w=0.             !< weight factor
    real, allocatable :: vec(:)           !< vector values
    integer           :: best_loc=0       !< best vec in neigh
    real              :: best_loc_cost=0. !< cost of best vec in neigh
end type

type, extends(optimizer) :: de_opt
    private
    type(individual), allocatable :: pop(:)                       !< solution vector population
    type(individual)              :: trial, trial_loc, trial_glob !< trial individuals
    type(ran_tabu)                :: rt                           !< for random number generation
    integer                       :: best_glob=0                  !< globally best vec
    integer                       :: worst_glob=0                 !< globally worst vec
    real                          :: best_glob_cost=0.            !< cost of globally best vec
    real                          :: worst_glob_cost=0.           !< cost of globally worst vec
    logical                       :: exists=.false.               !< to indicate exists 
  contains
    procedure, pass(self) :: new      => new_de_opt
    procedure, pass(self) :: minimize => de_minimize
    procedure, pass(self) :: kill     => kill_de_opt
end type

contains

    !> \brief  is a constructor
    subroutine new_de_opt( self, spec )
        use simple_opt_spec, only: opt_spec
        use simple_jiffys,   only: alloc_err
        class(de_opt), intent(inout) :: self !< instance
        class(opt_spec), intent(in)  :: spec !< specification
        integer                      :: alloc_stat, i
        call self%kill
        allocate(self%pop(spec%npop), self%trial%vec(spec%ndim),&
        self%trial_loc%vec(spec%ndim), self%trial_glob%vec(spec%ndim), stat=alloc_stat )
        call alloc_err( 'new; simple_de_opt, 1', alloc_stat )
        self%rt = ran_tabu(spec%npop)
        do i=1,spec%npop
            self%trial%F  = 0.1+0.9*ran3()
            self%pop(i)%w = ran3()
            allocate(self%pop(i)%vec(spec%ndim), stat=alloc_stat)
            call alloc_err( 'new; simple_de_opt, 2', alloc_stat)
        end do
        self%exists = .true.
    end subroutine
    
    !> \brief  is the bounded DE minimization
    subroutine de_minimize( self, spec, lowest_cost )
        use simple_opt_spec, only: opt_spec
        use simple_math,     only: cyci_1d
        class(de_opt), intent(inout)   :: self        !< instance     
        class(opt_spec), intent(inout) :: spec        !< specification
        real, intent(out)              :: lowest_cost !< lowest cost
        integer                        :: j, i, ir, order(spec%npop)
        real                           :: rtol
        ! initialize
        call init
        call find_best_glob
        call find_best_loc
        call find_worst
        do j=1,spec%maxits ! iterations loop
            ! randomize order
            call self%rt%reset
            call self%rt%ne_ran_iarr(order)
            do ir=1,spec%npop ! individuals loop
                i = order(ir) ! randomized ring search      
                call m_and_x_loc(i)
                call m_and_x_glob(i)
                call comb_and_check(i)
                self%trial%cost = spec%costfun(self%trial%vec, spec%ndim)
                call selec_trial(i)
            end do
            ! calculate relative tolerance
            call find_worst
            rtol=2.0*abs(self%worst_glob_cost-self%best_glob_cost)/&
            (abs(self%worst_glob_cost)+abs(self%best_glob_cost)+TINY)
            if(rtol < spec%ftol) exit
        end do
        spec%x      = self%pop(self%best_glob)%vec
        lowest_cost = self%best_glob_cost
        
        contains
        
            !> \brief  initializes the population using randomized bounds
            subroutine init
                integer :: i, j
                ! first member is the best-so-far solution
                self%pop(1)%vec(:) = spec%x
                ! the others are obtained by randomized bounds
                do i=2,spec%npop
                    self%pop(i)%w = ran3()
                    do j=1,spec%ndim
                        self%pop(i)%vec(j) = spec%limits(j,1)+ran3()*(spec%limits(j,2)-spec%limits(j,1))
                    end do
                end do
                ! calculate costs
                do i=1,spec%npop
                    self%pop(i)%cost = spec%costfun(self%pop(i)%vec, spec%ndim)
                end do
            end subroutine
        
            !> \brief  is for finding the fittest solution globally
            subroutine find_best_glob
                integer :: minpos(1)
                minpos = minloc(self%pop(:)%cost)
                self%best_glob = minpos(1)
                self%best_glob_cost = self%pop(self%best_glob)%cost
            end subroutine
            
            !> \brief  is for finding the worst solution
            subroutine find_worst
                integer :: maxpos(1)
                maxpos = maxloc(self%pop(:)%cost)
                self%worst_glob = maxpos(1)
                self%worst_glob_cost = self%pop(self%worst_glob)%cost
            end subroutine
            
            !> \brief  is for finding the fittest solutions in the local neighborhood model 
            subroutine find_best_loc
                integer :: minpos(1), i, j, counter
                integer :: neigh(2*spec%neigh+1)
                real    :: costs(2*spec%neigh+1)
                do i=1,spec%npop
                    counter = 0
                    do j=i-spec%neigh,i+spec%neigh
                        counter        = counter+1
                        neigh(counter) = cyci_1d([1,spec%npop], j)
                        costs(counter) = self%pop(neigh(counter))%cost
                    end do
                    minpos                    = minloc(costs)
                    self%pop(i)%best_loc      = neigh(minpos(1))
                    self%pop(i)%best_loc_cost = self%pop(neigh(minpos(1)))%cost
                end do
            end subroutine
            
            !> \brief  is for selecting trial solutions for replacement
            subroutine selec_trial( i )
                integer, intent(in) :: i
                integer :: j, ix
                if( self%trial%cost <= self%pop(i)%cost ) then
                    self%pop(i)%cost = self%trial%cost
                    self%pop(i)%F    = self%trial%F
                    self%pop(i)%vec  = self%trial%vec
                    self%pop(i)%w    = self%trial%w
                    if( self%trial%cost <= self%best_glob_cost ) then
                        self%best_glob = i
                        self%best_glob_cost = self%trial%cost
                    endif
                    do j=i-spec%neigh,i+spec%neigh
                        ix = cyci_1d([1,spec%npop], j)
                        if( self%trial%cost <= self%pop(ix)%best_loc_cost ) then
                            self%pop(ix)%best_loc = i
                            self%pop(ix)%best_loc_cost = self%trial%cost
                        endif
                    end do
                end if        
            end subroutine
            
            !> \brief  generates a trial solution vector using a local neighborhood model
            subroutine m_and_x_loc( i )
                use simple_rnd, only: irnd_uni
                integer, intent(in) :: i
                real                :: k, f2
                integer             :: j, a(2), neigh(2*spec%neigh+1), counter, excl
                excl = 0
                counter = 0
                do j=i-spec%neigh,i+spec%neigh
                    counter = counter+1
                    neigh(counter) = cyci_1d([1,spec%npop], j)
                    if( neigh(counter) == self%best_glob )then
                        excl = counter
                    else
                        excl = 0
                    endif
                end do
                ! generate the non-equal two random numbers, not equal to best_glob  
                a(1) = excl
                do while( a(1) == excl ) 
                    a(1) = irnd_uni(2*spec%neigh+1)
                end do
                a(2) = excl
                do while( a(2) == a(1) .or. a(2) == excl )
                    a(2) = irnd_uni(2*spec%neigh+1)
                end do
                ! self-adaption of the amplification factor
                if( ran3() < 0.1 ) then
                    self%trial%F = 0.1+0.9*ran3()
                else
                    self%trial%F = self%pop(i)%F
                endif
                ! this linear combination replaces mutation and crossover
                k = ran3()
                f2 = self%trial%F*k
                self%trial_loc%vec = self%pop(i)%vec+k*(self%pop(self%pop(i)%best_loc)%vec-&
                self%pop(i)%vec)+f2*(self%pop(neigh(a(1)))%vec-self%pop(neigh(a(2)))%vec)
            end subroutine
            
            !> \brief  generates a trial solution vector using a global neighborhood model
            subroutine m_and_x_glob( i )
                integer, intent(in) :: i
                integer             :: a(4)
                real                :: k, f2
                ! generate the non-equal four random numbers, not equal to the global best
                call self%rt%reset
                call self%rt%insert(self%best_glob) 
                call self%rt%ne_ran_iarr(a)         
                ! this 'randomized' linear combination replaces mutation and crossover
                k = ran3()
                f2 = self%trial%F*k
                self%trial_glob%vec = self%pop(i)%vec + k*(self%pop(self%best_glob)%vec-self%pop(i)%vec) +&
                f2*(self%pop(a(1))%vec+self%pop(a(2))%vec-self%pop(a(3))%vec-self%pop(a(4))%vec)
                ! now m&x the weight, note the change of F (no crossover performed)
                self%trial%w = self%pop(i)%w + self%trial%F*(self%pop(self%best_glob)%w-self%pop(i)%w) +&
                self%trial%F*(self%pop(a(1))%w+self%pop(a(2))%w-self%pop(a(3))%w-self%pop(a(4))%w)
                if( self%trial%w > 0.95) self%trial%w = 0.95
                if( self%trial%w < 0.05) self%trial%w = 0.05
            end subroutine
            
            !> \brief  combines the solutions from the two mutation models and checks the result
            subroutine comb_and_check( i )
                use simple_opt_subs,  only: check_and_correct_vec
                integer, intent(in)          :: i
                self%trial%w   = self%pop(i)%w
                self%trial%vec = self%trial%w*self%trial_glob%vec+(1.-self%trial%w)*self%trial_loc%vec
                call check_and_correct_vec(spec, self%trial%vec)
            end subroutine
            
    end subroutine
    
    !> \brief  is a destructor
    subroutine kill_de_opt( self )
        class(de_opt), intent(inout) :: self
        integer :: i
        if( self%exists )then
            call self%rt%kill
            deallocate( self%trial%vec )
            deallocate( self%trial_loc%vec )
            deallocate( self%trial_glob%vec )
            do i=1,size(self%pop)
                deallocate( self%pop(i)%vec )    
            end do    
            deallocate( self%pop )
            self%exists = .false.
       endif
    end subroutine   
    
end module simple_de_opt
