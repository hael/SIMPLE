module simple_prime_srch
use simple_defs
use simple_math
use simple_params, only: params
use simple_jiffys, only: alloc_err
use simple_math
implicit none

type prime_srch
    private
    logical                :: exists = .false. !< 2 indicate existence
  contains
    ! CONSTRUCTOR 
    procedure :: new
    ! GETTERS

    ! CALCULATORS

    ! DESTRUCTOR
    procedure :: kill
end type prime_srch

interface prime_srch
    module procedure constructor
end interface prime_srch

contains
    
    !>  \brief  is a constructor
    function constructor( p ) result( self )
        class(params), target, intent(in) :: p !< parameters
        type(prime_srch) :: self
        call self%new(p)
    end function constructor
    
    !>  \brief  is a constructor
    subroutine new( self, p )
        class(prime_srch),     intent(inout) :: self  !< instance
        class(params), target, intent(in)    :: p     !< parameters
        self%exists = .true.
    end subroutine new

    !>  \brief  is a destructor
    subroutine kill( self )
        class(prime_srch), intent(inout) :: self !< instance
        if( self%exists )then
            self%exists = .false.
        endif
    end subroutine kill

end module simple_prime_srch
