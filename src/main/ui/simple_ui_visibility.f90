!@descr: typed visibility levels for command-descriptor presentation
module simple_ui_visibility
implicit none
private

integer, parameter, public :: UI_VIS_STANDARD  = 1
integer, parameter, public :: UI_VIS_ADVANCED  = 2
integer, parameter, public :: UI_VIS_DEVELOPER = 3

public :: ui_visibility_is_valid
public :: ui_visibility_name

contains

    pure logical function ui_visibility_is_valid( visibility )
        integer, intent(in) :: visibility
        ui_visibility_is_valid = visibility == UI_VIS_STANDARD .or. &
                                 visibility == UI_VIS_ADVANCED .or. &
                                 visibility == UI_VIS_DEVELOPER
    end function ui_visibility_is_valid

    pure function ui_visibility_name( visibility ) result( name )
        integer, intent(in) :: visibility
        character(len=9)    :: name
        select case( visibility )
            case( UI_VIS_STANDARD )
                name = 'standard'
            case( UI_VIS_ADVANCED )
                name = 'advanced'
            case( UI_VIS_DEVELOPER )
                name = 'developer'
            case default
                name = 'invalid'
        end select
    end function ui_visibility_name

end module simple_ui_visibility
