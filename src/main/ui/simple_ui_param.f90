!@descr: module defining the ui_param type, which is used to define input parameters for the simple_ui_program interface
module simple_ui_param
use simple_string, only: string
use simple_error,  only: simple_exception
use simple_ui_descriptor_types, only: ui_choice, ui_choices_from_legacy_placeholder
use simple_ui_visibility, only: UI_VIS_DEVELOPER, ui_visibility_is_valid
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
    type(string) :: gui_submenu
    type(string) :: active_flags
    type(string) :: exclusive_group
    type(string) :: cval_default
    real         :: rval_default = 0.
    logical      :: required = .true.
    logical      :: has_default = .false.
    integer      :: visibility = UI_VIS_DEVELOPER
    logical      :: online   = .false.
    type(ui_choice), allocatable :: choices(:)
contains
    procedure, private :: set_param_1
    procedure, private :: set_param_2
    generic            :: set_param => set_param_1, set_param_2
    procedure          :: apply_gui_overrides
    procedure          :: refresh_legacy_choices
    final              :: finalize
end type ui_param

contains

    subroutine set_param_1( self, key, keytype, label, help, placeholder, required, default_value )
        class(ui_param), intent(inout) :: self
        character(len=*),       intent(in)    :: key, keytype, label, help, placeholder
        logical,                intent(in)    :: required
        real,                   intent(in)    :: default_value
        self%key               = trim(key)
        self%keytype           = trim(keytype)
        self%label       = trim(label)
        self%help        = trim(help)
        self%placeholder = trim(placeholder)
        self%required = required
        self%has_default = .not. self%required
        self%visibility = UI_VIS_DEVELOPER
        if( self%has_default ) self%rval_default = default_value
        call self%refresh_legacy_choices()
    end subroutine set_param_1

    subroutine set_param_2( self, key, keytype, label, help, placeholder, required, default_value )
        class(ui_param), intent(inout) :: self
        character(len=*),       intent(in)    :: key, keytype, label, help, placeholder
        logical,                intent(in)    :: required
        character(len=*),       intent(in)    :: default_value
        self%key               = trim(key)
        self%keytype           = trim(keytype)
        self%label       = trim(label)
        self%help        = trim(help)
        self%placeholder = trim(placeholder)
        self%required = required
        self%has_default = .not. self%required
        self%visibility = UI_VIS_DEVELOPER
        if( self%has_default ) self%cval_default = trim(default_value)
        call self%refresh_legacy_choices()
    end subroutine set_param_2

    subroutine apply_gui_overrides(p, gui_submenu, gui_exclusive_group, gui_active_flags, gui_online, gui_visibility)
        class(ui_param),            intent(inout) :: p
        character(len=*), optional, intent(in)    :: gui_submenu, gui_exclusive_group, gui_active_flags
        logical,          optional, intent(in)    :: gui_online
        integer,          optional, intent(in)    :: gui_visibility
        if( present(gui_submenu))         p%gui_submenu     = trim(gui_submenu)
        if( present(gui_exclusive_group)) p%exclusive_group = trim(gui_exclusive_group)
        if( present(gui_active_flags))    p%active_flags    = trim(gui_active_flags)
        if( present(gui_online))          p%online          = gui_online
        if( present(gui_visibility) )then
            if( .not. ui_visibility_is_valid(gui_visibility) )then
                THROW_HARD('ui_param%apply_gui_overrides received an invalid visibility level')
            endif
            p%visibility = gui_visibility
        endif
    end subroutine apply_gui_overrides

    subroutine refresh_legacy_choices( self )
        class(ui_param), intent(inout) :: self
        character(len=:), allocatable :: legacy_placeholder
        logical :: valid
        legacy_placeholder = self%placeholder%to_char()
        call ui_choices_from_legacy_placeholder(self%keytype%to_char(), legacy_placeholder, &
            &self%choices, valid)
        self%placeholder = standardize_placeholder(self%keytype%to_char())
    end subroutine refresh_legacy_choices

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
        call self%gui_submenu%kill()
        call self%active_flags%kill()
        call self%exclusive_group%kill()
        call self%cval_default%kill()
    end subroutine finalize

end module simple_ui_param
