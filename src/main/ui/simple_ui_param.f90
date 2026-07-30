!@descr: module defining the ui_param type, which is used to define input parameters for the simple_ui_program interface
module simple_ui_param
use simple_string, only: string
use simple_error,  only: simple_exception
use simple_ui_descriptor_types, only: ui_choice, ui_choices_are_valid
implicit none
#include "simple_local_flags.inc"

integer, parameter, public :: UI_PLACEHOLDER_MAX_LEN = 40
public :: ui_placeholder_is_standard

! common input parameter type
type ui_param
    ! deliberately made public (close entaglement with simple_ui_program)
    type(string) :: key
    type(string) :: keytype
    type(string) :: label
    type(string) :: help
    type(string) :: placeholder
    type(string) :: units
    type(string) :: cval_default
    real         :: rval_default = 0.
    logical      :: required = .true.
    logical      :: has_default = .false.
    type(ui_choice), allocatable :: choices(:)
contains
    procedure, private :: set_param_1
    procedure, private :: set_param_2
    generic            :: set_param => set_param_1, set_param_2
    procedure          :: set_choices
    procedure          :: set_generated_default
    final              :: finalize
end type ui_param

contains

    subroutine set_param_1( self, key, keytype, label, help, placeholder, required, default_value, choices )
        class(ui_param), intent(inout) :: self
        character(len=*),       intent(in)    :: key, keytype, label, help, placeholder
        logical,                intent(in)    :: required
        real,                   intent(in)    :: default_value
        type(ui_choice), optional, intent(in)  :: choices(:)
        self%key               = trim(key)
        self%keytype           = trim(keytype)
        self%label       = trim(label)
        self%help        = trim(help)
        self%placeholder = trim(placeholder)
        self%required = required
        self%has_default = .not. self%required
        if( self%has_default ) self%rval_default = default_value
        call self%set_choices(choices)
    end subroutine set_param_1

    subroutine set_param_2( self, key, keytype, label, help, placeholder, required, default_value, choices )
        class(ui_param), intent(inout) :: self
        character(len=*),       intent(in)    :: key, keytype, label, help, placeholder
        logical,                intent(in)    :: required
        character(len=*),       intent(in)    :: default_value
        type(ui_choice), optional, intent(in)  :: choices(:)
        self%key               = trim(key)
        self%keytype           = trim(keytype)
        self%label       = trim(label)
        self%help        = trim(help)
        self%placeholder = trim(placeholder)
        self%required = required
        self%has_default = .not. self%required
        if( self%has_default ) self%cval_default = trim(default_value)
        call self%set_choices(choices)
        if( self%has_default .and. (self%keytype%to_char() == 'binary' .or. &
            &self%keytype%to_char() == 'multi') )then
            if( .not. choice_is_available(self%choices, self%cval_default%to_char()) )then
                THROW_HARD('UI parameter default is not a declared choice: '//self%key%to_char())
            endif
        endif
    end subroutine set_param_2

    subroutine set_choices( self, choices )
        class(ui_param), intent(inout) :: self
        type(ui_choice), optional, intent(in) :: choices(:)
        if( allocated(self%choices) ) deallocate(self%choices)
        if( present(choices) )then
            allocate(self%choices(size(choices)))
            self%choices = choices
        else
            allocate(self%choices(0))
        endif
        if( .not. ui_choices_are_valid(self%keytype%to_char(), self%choices) )then
            THROW_HARD('Invalid choices for UI parameter: '//self%key%to_char())
        endif
        self%placeholder = standardize_placeholder(self%keytype%to_char())
    end subroutine set_choices

    subroutine set_generated_default( self, value, found )
        class(ui_param), intent(inout) :: self
        character(len=*), intent(in)   :: value
        logical,          intent(in)   :: found
        integer :: ios
        logical :: had_default
        if( .not. found ) return
        had_default = self%has_default
        select case( trim(self%keytype%to_char()) )
            case('num', 'int', 'float')
                read(value, *, iostat=ios) self%rval_default
                if( ios /= 0 )then
                    THROW_HARD('Generated UI default is not numeric: '//self%key%to_char())
                endif
            case default
                if( self%keytype%to_char() == 'binary' .or. self%keytype%to_char() == 'multi' )then
                    if( .not. choice_is_available(self%choices, trim(value)) )then
                        ! A baseline parameter value may be valid globally but outside
                        ! this program's declared subset. Retain its validated local
                        ! default; an input without one remains a descriptor error.
                        if( had_default ) return
                        THROW_HARD('Generated UI default is not a declared choice: '//self%key%to_char())
                    endif
                endif
                self%cval_default = trim(value)
        end select
        self%has_default = .true.
    end subroutine set_generated_default

    pure logical function choice_is_available( choices, value )
        type(ui_choice),  intent(in) :: choices(:)
        character(len=*), intent(in) :: value
        integer :: i
        choice_is_available = .false.
        do i = 1, size(choices)
            if( choices(i)%value%to_char() == trim(value) )then
                choice_is_available = .true.
                return
            endif
        enddo
    end function choice_is_available

    pure function standardize_placeholder( keytype ) result( placeholder )
        character(len=*), intent(in) :: keytype
        character(len=UI_PLACEHOLDER_MAX_LEN) :: placeholder
        select case( trim(keytype) )
            case('binary', 'multi', 'hidden_dir', 'hidden_float', 'hidden_int')
                placeholder = ''
            case('dir')
                placeholder = 'e.g. /path/to/folder'
            case('file')
                placeholder = 'e.g. input.mrc'
            case('num', 'int', 'float')
                placeholder = 'e.g. 10'
            case('str', 'string')
                placeholder = 'e.g. value'
            case default
                placeholder = 'e.g. value'
        end select
    end function standardize_placeholder

    pure logical function ui_placeholder_is_standard( keytype, placeholder )
        character(len=*), intent(in) :: keytype, placeholder
        ui_placeholder_is_standard = trim(placeholder) == trim(standardize_placeholder(keytype))
    end function ui_placeholder_is_standard

    subroutine finalize(self)
        type(ui_param), intent(inout) :: self
        call self%key%kill()
        call self%keytype%kill()
        call self%label%kill()
        call self%help%kill()
        call self%placeholder%kill()
        call self%units%kill()
        call self%cval_default%kill()
    end subroutine finalize

end module simple_ui_param
