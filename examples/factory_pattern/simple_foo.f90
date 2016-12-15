module simple_foo
implicit none

public :: foo
private

type, abstract :: foo
  contains
    procedure(generic_constr),       deferred    :: construct
    procedure(generic_print),        deferred    :: print  
    procedure(generic_destr),        deferred    :: destruct
    procedure(generic_assign), private, deferred :: assign
    generic :: assignment(=) => assign
    procedure(generic_add),    private, deferred :: add
    generic :: operator(+)   => add
end type

abstract interface

    subroutine generic_constr( num, n ) 
        import :: foo
        class(foo), intent(inout) :: num
        integer, intent(in) :: n
    end subroutine
    
    subroutine generic_print( num )
        import :: foo
        class(foo), intent(in) :: num
    end subroutine
    
    subroutine generic_destr( num )
        import :: foo
        class(foo), intent(inout) :: num
    end subroutine
    
    subroutine generic_assign( num1, num2 )
        import :: foo
        class(foo), intent(inout) :: num1
        class(foo), intent(in)  :: num2
    end subroutine
    
    function generic_add( num1, num2 ) result( num )
        import :: foo
        class(foo), intent(in) :: num1, num2
        class(foo), allocatable :: num
    end function
    
end interface

end module simple_foo