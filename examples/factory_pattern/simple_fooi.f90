module simple_fooi
use simple_foo, only: foo
implicit none

type, extends(foo) :: fooi
    integer :: n
    integer, allocatable :: arr(:)
  contains
    procedure :: construct       => construct_fooi
    procedure :: print           => print_fooi
    procedure :: destruct        => destruct_fooi
    procedure, private :: assign => assign_fooi
    procedure, private :: add    => add_fooi
end type

contains

    subroutine construct_fooi( num, n )
        class(fooi), intent(inout) :: num
        integer, intent(in) :: n
        num%n = n
        allocate( num%arr(n) )
        num%arr = 3
    end subroutine
    
    subroutine print_fooi( num )
        class(fooi), intent(in) :: num
        write(*,*) num%arr
    end subroutine
    
    subroutine destruct_fooi( num )
        class(fooi), intent(inout) :: num
        if( allocated(num%arr) ) deallocate(num%arr)
    end subroutine
    
    subroutine assign_fooi( num1, num2 )
        class(fooi), intent(inout) :: num1
        class(foo), intent(in)     :: num2
        select type(num2)
            type is (fooi)
                num1%arr = num2%arr
        end select
    end subroutine
    
    function add_fooi( num1, num2 ) result( num )
        class(fooi), intent(in) :: num1
        class(foo),  intent(in) :: num2
        class(foo), allocatable :: num
        select type(num2)
            type is (fooi)
                allocate( fooi :: num )
                select type (num)
                    type is (fooi)
                        num%arr = num1%arr+num2%arr
                end select
        end select
    end function

end module simple_fooi