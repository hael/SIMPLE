module simple_foor
use simple_foo, only: foo
implicit none

type, extends(foo) :: foor
    integer :: n
    real, allocatable :: arr(:)
  contains
    procedure :: construct       => construct_foor
    procedure :: print           => print_foor
    procedure :: destruct        => destruct_foor
    procedure, private :: assign => assign_foor
    procedure, private :: add    => add_foor
end type

contains

    subroutine construct_foor( num, n )
        class(foor), intent(inout) :: num
        integer, intent(in) :: n
        num%n = n
        allocate( num%arr(n) )
        num%arr = 6.
    end subroutine
    
    subroutine print_foor( num )
        class(foor), intent(in) :: num
        write(*,*) num%arr
    end subroutine
    
    subroutine destruct_foor( num )
        class(foor), intent(inout) :: num
        if( allocated(num%arr) ) deallocate(num%arr)
    end subroutine
    
    subroutine assign_foor( num1, num2 )
        class(foor), intent(inout) :: num1
        class(foo), intent(in)     :: num2
        select type(num2)
            type is (foor)
                num1%arr = num2%arr
        end select
    end subroutine
    
    function add_foor( num1, num2 ) result( num )
        class(foor), intent(in) :: num1
        class(foo),  intent(in) :: num2
        class(foo), allocatable :: num
        select type(num2)
            type is (foor)
                allocate( foor :: num )
                select type (num)
                    type is (foor)
                        num%arr = num1%arr+num2%arr
                end select
        end select
    end function

end module simple_foor