module simple_foofactory
use simple_foo,  only: foo
use simple_foor, only: foor
use simple_fooi, only: fooi
implicit none

public :: foofactory
private

type :: foofactory
    private
    class(foo), allocatable :: foo_type ! Which type of foo? Note 'class' not 'type'
  contains
    procedure :: construct ! Construct                                      
    procedure :: destruct  ! Destruct
end type

contains

    function construct( num, which, n ) result( ptr )
        class(foofactory), intent(inout), target :: num
        character(len=*), intent(in) :: which
        integer, intent(in) :: n
        class(foo), pointer :: ptr
        select case( which )
            case('foor')
                call num%destruct
                allocate(foor :: num%foo_type)
                call num%foo_type%construct(n)
                ptr => num%foo_type
            case('fooi')
                call num%destruct
                allocate(fooi :: num%foo_type)
                call num%foo_type%construct(n)
                ptr => num%foo_type
            case DEFAULT
                ptr => null()
                stop 'Unsupported class in foofactory constructor'
        end select
    end function

    subroutine destruct( num )
        class(foofactory), intent(inout) :: num
        if( allocated(num%foo_type) )then
            call num%foo_type%destruct
            deallocate(num%foo_type)
        endif
    end subroutine

end module simple_foofactory