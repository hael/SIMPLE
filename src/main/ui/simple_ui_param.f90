!@descr: module defining the ui_param type, which is used to define input parameters for the simple_ui_program interface
module simple_ui_param
use simple_string, only: string
use simple_error,  only: simple_exception
use simple_ui_descriptor_types, only: ui_choice, ui_choices_from_legacy_placeholder
use simple_ui_visibility, only: UI_VIS_DEVELOPER, ui_visibility_from_advanced, ui_visibility_is_valid
implicit none
#include "simple_local_flags.inc"

! common input parameter type
type ui_param
    ! deliberately made public (close entaglement with simple_ui_program)
    type(string) :: key
    type(string) :: keytype
    type(string) :: descr_short
    type(string) :: descr_long
    type(string) :: descr_placeholder
    type(string) :: units
    type(string) :: gui_submenu
    type(string) :: active_flags
    type(string) :: exclusive_group
    type(string) :: cval_default
    real         :: rval_default = 0.
    logical      :: required = .true.
    logical      :: has_default = .false.
    logical      :: advanced = .true.
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

    subroutine set_param_1( self, key, keytype, descr_short, descr_long, descr_placeholder, required, default_value )
        class(ui_param), intent(inout) :: self
        character(len=*),       intent(in)    :: key, keytype, descr_short, descr_long, descr_placeholder
        logical,                intent(in)    :: required
        real,                   intent(in)    :: default_value
        self%key               = trim(key)
        self%keytype           = trim(keytype)
        self%descr_short       = trim(descr_short)
        self%descr_long        = trim(descr_long)
        self%descr_placeholder = trim(descr_placeholder)
        self%required = required
        self%has_default = .not. self%required
        if( self%has_default ) self%rval_default = default_value
        call self%refresh_legacy_choices()
    end subroutine set_param_1

    subroutine set_param_2( self, key, keytype, descr_short, descr_long, descr_placeholder, required, default_value )
        class(ui_param), intent(inout) :: self
        character(len=*),       intent(in)    :: key, keytype, descr_short, descr_long, descr_placeholder
        logical,                intent(in)    :: required
        character(len=*),       intent(in)    :: default_value
        self%key               = trim(key)
        self%keytype           = trim(keytype)
        self%descr_short       = trim(descr_short)
        self%descr_long        = trim(descr_long)
        self%descr_placeholder = trim(descr_placeholder)
        self%required = required
        self%has_default = .not. self%required
        if( self%has_default ) self%cval_default = trim(default_value)
        call self%refresh_legacy_choices()
    end subroutine set_param_2

    subroutine apply_gui_overrides(p, gui_submenu, gui_exclusive_group, gui_active_flags, gui_advanced, gui_online, gui_visibility)
        class(ui_param),            intent(inout) :: p
        character(len=*), optional, intent(in)    :: gui_submenu, gui_exclusive_group, gui_active_flags
        logical,          optional, intent(in)    :: gui_advanced, gui_online
        integer,          optional, intent(in)    :: gui_visibility
        if( present(gui_submenu))         p%gui_submenu     = trim(gui_submenu)
        if( present(gui_exclusive_group)) p%exclusive_group = trim(gui_exclusive_group)
        if( present(gui_active_flags))    p%active_flags    = trim(gui_active_flags)
        if( present(gui_online))          p%online          = gui_online
        if( present(gui_advanced))        p%advanced        = gui_advanced
        if( present(gui_visibility) )then
            if( .not. ui_visibility_is_valid(gui_visibility) )then
                THROW_HARD('ui_param%apply_gui_overrides received an invalid visibility level')
            endif
            p%visibility = gui_visibility
        else if( present(gui_advanced) )then
            p%visibility = ui_visibility_from_advanced(gui_advanced)
        endif
    end subroutine apply_gui_overrides

    subroutine refresh_legacy_choices( self )
        class(ui_param), intent(inout) :: self
        logical :: valid
        call ui_choices_from_legacy_placeholder(self%keytype%to_char(), self%descr_placeholder%to_char(), &
            &self%choices, valid)
        if( .not. valid ) return
    end subroutine refresh_legacy_choices

    subroutine finalize(self)
        type(ui_param), intent(inout) :: self
        call self%key%kill()
        call self%keytype%kill()
        call self%descr_short%kill()
        call self%descr_long%kill()
        call self%descr_placeholder%kill()
        call self%units%kill()
        call self%gui_submenu%kill()
        call self%active_flags%kill()
        call self%exclusive_group%kill()
        call self%cval_default%kill()
    end subroutine finalize

end module simple_ui_param
