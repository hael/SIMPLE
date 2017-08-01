!==Class simple_sgd_opt
!
! Minimization of an externally defined function by stochastic gradient descent
! The code is distributed with the hope that it will be useful, but _WITHOUT_ _ANY_ _WARRANTY_. 
! Redistribution or modification is regulated by the GNU General Public License. 
! *Author:* Hans Elmlund, 2013-10-15
!
module simple_sgd_opt
use simple_optimizer, only: optimizer
use simple_jiffys     ! singleton
implicit none

public :: sgd_opt
private

type, extends(optimizer) :: sgd_opt
    private
    real, allocatable :: gradient(:)
    logical :: exists=.false.
  contains
    procedure :: new      => new_sgd_opt
    procedure :: minimize => sgd_minimize
    procedure :: kill     => kill_sgd_opt
end type

contains
    
    !> \brief  is a constructor
    subroutine new_sgd_opt( self, spec )
        use simple_opt_spec, only: opt_spec
        class(sgd_opt), intent(inout) :: self !< instance
        class(opt_spec), intent(in)   :: spec !< specification
        integer :: alloc_stat
        call self%kill
        allocate(self%gradient(spec%ndim),stat=alloc_stat)
        call alloc_err('In: new_sgd_opt; simple_sgd_opt', alloc_stat)
        self%gradient = 0.
        self%exists   = .true.
    end subroutine
    
    !>  \brief  stochastic gradient descent minimizer
    subroutine sgd_minimize( self, spec, lowest_cost )
        use simple_opt_spec, only: opt_spec
        class(sgd_opt), intent(inout)  :: self        !< instance
        class(opt_spec), intent(inout) :: spec        !< specification
        real, intent(out)              :: lowest_cost !< minimum function value
        real :: prev_lowest_cost, x
        integer :: i
        if( .not. associated(spec%costfun) )then 
            stop 'cost function not associated in opt_spec; sgd_minimize; simple_sgd_opt'
        endif
        if( .not. associated(spec%gcostfun) )then
            stop 'gradient of cost function not associated in opt_spec; sgd_minimize; simple_sgd_opt'
        endif
        ! calculate initial cost
        spec%nevals = 0 
        lowest_cost = spec%costfun(spec%x, spec%ndim)
        prev_lowest_cost = huge(x)
        spec%nevals = spec%nevals+1
        ! iterate
        do i=1,spec%maxits
            ! calculate gradient
            self%gradient = spec%gcostfun(spec%x, spec%ndim)
            ! update
            spec%x = spec%x-spec%eps*self%gradient
            spec%nevals = spec%nevals+1
            if( mod(i,10)==0 )then
                prev_lowest_cost = lowest_cost
                lowest_cost = spec%costfun(spec%x, spec%ndim)
                print *, 'cost=', lowest_cost
            endif
        end do
        lowest_cost = spec%costfun(spec%x, spec%ndim)
    end subroutine
    
    !> \brief  is a destructor
    subroutine kill_sgd_opt( self )
        class(sgd_opt), intent(inout) :: self !< instance
        if( self%exists )then
            deallocate(self%gradient) 
            self%exists = .false.
        endif
    end subroutine

end module simple_sgd_opt
