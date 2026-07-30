!@descr: shared presentation-only value types for the command descriptor layer
module simple_ui_descriptor_types
use simple_string, only: string
implicit none
private

type, public :: ui_choice
    type(string) :: value
    type(string) :: label
    type(string) :: help
end type ui_choice

public :: ui_choices, ui_choices_are_valid

contains

    function ui_choices( values ) result( choices )
        character(len=*), intent(in) :: values(:)
        type(ui_choice), allocatable :: choices(:)
        integer :: i
        allocate(choices(size(values)))
        do i = 1, size(values)
            choices(i)%value = trim(values(i))
            choices(i)%label = trim(values(i))
            choices(i)%help  = ''
        enddo
    end function ui_choices

    pure logical function ui_choices_are_valid( keytype, choices )
        character(len=*), intent(in) :: keytype
        type(ui_choice),  intent(in) :: choices(:)
        integer :: i, j
        ui_choices_are_valid = .false.
        select case( trim(keytype) )
            case('binary')
                if( size(choices) /= 2 ) return
            case('multi')
                if( size(choices) < 2 ) return
            case default
                if( size(choices) == 0 ) ui_choices_are_valid = .true.
                return
        end select
        do i = 1, size(choices)
            if( len_trim(choices(i)%value%to_char()) == 0 ) return
            do j = 1, i - 1
                if( choices(i)%value%to_char() == choices(j)%value%to_char() ) return
            enddo
        enddo
        ui_choices_are_valid = .true.
    end function ui_choices_are_valid

end module simple_ui_descriptor_types
