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

public :: ui_choices_from_legacy_placeholder

contains

    subroutine ui_choices_from_legacy_placeholder( keytype, placeholder, choices, valid )
        character(len=*),              intent(in)  :: keytype, placeholder
        type(ui_choice), allocatable,  intent(out) :: choices(:)
        logical,                       intent(out) :: valid
        integer :: first_delimiter, last_delimiter, choice_count, choice_index, start, finish
        character(len=:), allocatable :: values

        valid = .true.
        if( trim(keytype) /= 'binary' .and. trim(keytype) /= 'multi' )then
            allocate(choices(0))
            return
        endif

        first_delimiter = index(placeholder, '(')
        if( first_delimiter == 0 )then
            valid = .false.
            allocate(choices(0))
            return
        endif
        last_delimiter = index(placeholder(first_delimiter+1:), ')')
        if( last_delimiter == 0 )then
            valid = .false.
            allocate(choices(0))
            return
        endif
        last_delimiter = last_delimiter + first_delimiter
        values = placeholder(first_delimiter+1:last_delimiter-1)
        if( len_trim(values) == 0 )then
            valid = .false.
            allocate(choices(0))
            return
        endif

        choice_count = 1
        do choice_index = 1, len_trim(values)
            if( values(choice_index:choice_index) == '|' ) choice_count = choice_count + 1
        enddo
        if( (trim(keytype) == 'binary' .and. choice_count /= 2) .or. &
            &(trim(keytype) == 'multi' .and. choice_count < 2) )then
            valid = .false.
            allocate(choices(0))
            return
        endif

        allocate(choices(choice_count))
        start = 1
        do choice_index = 1, choice_count
            if( choice_index == choice_count )then
                finish = len_trim(values)
            else
                finish = start + index(values(start:), '|') - 2
            endif
            choices(choice_index)%value = trim(values(start:finish))
            choices(choice_index)%label = trim(values(start:finish))
            choices(choice_index)%help  = ''
            start = finish + 2
        enddo
    end subroutine ui_choices_from_legacy_placeholder

end module simple_ui_descriptor_types
