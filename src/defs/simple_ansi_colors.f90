module simple_ansi_colors
implicit none
character(len=1), parameter :: C_ESC     = achar(27)
character(len=2), parameter :: C_START   = C_ESC // '['
character(len=1), parameter :: C_END     = 'm'
character(len=*), parameter :: C_BLACK   = '30'
character(len=*), parameter :: C_RED     = '31'
character(len=*), parameter :: C_GREEN   = '32'
character(len=*), parameter :: C_YELLOW  = '33'
character(len=*), parameter :: C_BLUE    = '34'
character(len=*), parameter :: C_MAGENTA = '35'
character(len=*), parameter :: C_CYAN    = '36'
character(len=*), parameter :: C_WHITE   = '37'
character(len=*), parameter :: C_CLEAR   = C_START // '0' // C_END

contains

    function color( str, code ) result( out )
        character(len=*), intent(in) :: str
        character(len=*), intent(in) :: code
        character(len=:), allocatable :: out
        out = C_START // CODE // C_END // trim(str) // C_CLEAR
    end function color

end module simple_ansi_colors
