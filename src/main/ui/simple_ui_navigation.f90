!@descr: generic suite and program-group navigation descriptors
module simple_ui_navigation
use simple_string, only: string
implicit none
private

integer, parameter, public :: UI_LAYOUT_GROUPED = 1
integer, parameter, public :: UI_LAYOUT_FLAT    = 2

type, public :: ui_program_group
    type(string)              :: id
    type(string)              :: title
    integer                   :: display_order = 0
    type(string), allocatable :: program_names(:)
contains
    procedure :: new => new_ui_program_group
    procedure :: nprograms
end type ui_program_group

type, public :: ui_suite
    type(string)                         :: id
    type(string)                         :: executable
    type(string)                         :: display_name
    integer                              :: layout = UI_LAYOUT_GROUPED
    type(ui_program_group), allocatable  :: groups(:)
contains
    procedure :: new => new_ui_suite
    procedure :: ngroups
end type ui_suite

public :: ui_layout_is_valid
public :: ui_layout_name

contains

    pure logical function ui_layout_is_valid( layout )
        integer, intent(in) :: layout
        ui_layout_is_valid = layout == UI_LAYOUT_GROUPED .or. layout == UI_LAYOUT_FLAT
    end function ui_layout_is_valid

    pure function ui_layout_name( layout ) result( name )
        integer, intent(in) :: layout
        character(len=7)    :: name

        select case( layout )
            case( UI_LAYOUT_GROUPED )
                name = 'grouped'
            case( UI_LAYOUT_FLAT )
                name = 'flat'
            case default
                name = 'invalid'
        end select
    end function ui_layout_name

    subroutine new_ui_program_group( self, id, title, display_order, program_names )
        class(ui_program_group),      intent(inout) :: self
        character(len=*),             intent(in)    :: id, title
        integer,                      intent(in)    :: display_order
        type(string),       optional, intent(in)    :: program_names(:)

        self%id            = trim(id)
        self%title         = trim(title)
        self%display_order = display_order
        if( allocated(self%program_names) ) deallocate(self%program_names)
        if( present(program_names) )then
            allocate(self%program_names(size(program_names)), source=program_names)
        else
            allocate(self%program_names(0))
        endif
    end subroutine new_ui_program_group

    pure integer function nprograms( self )
        class(ui_program_group), intent(in) :: self
        if( allocated(self%program_names) )then
            nprograms = size(self%program_names)
        else
            nprograms = 0
        endif
    end function nprograms

    subroutine new_ui_suite( self, id, executable, display_name, layout, groups )
        class(ui_suite),                 intent(inout) :: self
        character(len=*),                intent(in)    :: id, executable, display_name
        integer,                         intent(in)    :: layout
        type(ui_program_group), optional, intent(in)    :: groups(:)

        self%id           = trim(id)
        self%executable   = trim(executable)
        self%display_name = trim(display_name)
        self%layout       = layout
        if( allocated(self%groups) ) deallocate(self%groups)
        if( present(groups) )then
            allocate(self%groups(size(groups)), source=groups)
        else
            allocate(self%groups(0))
        endif
    end subroutine new_ui_suite

    pure integer function ngroups( self )
        class(ui_suite), intent(in) :: self
        if( allocated(self%groups) )then
            ngroups = size(self%groups)
        else
            ngroups = 0
        endif
    end function ngroups

end module simple_ui_navigation
