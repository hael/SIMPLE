module simple_ansi_ctrls
implicit none
character(len=*), parameter :: C_ESC            = achar(27)
character(len=*), parameter :: C_START          = C_ESC // '['
character(len=*), parameter :: C_END            = 'm'
character(len=*), parameter :: C_BOLD           = '1'
character(len=*), parameter :: C_GRAY           = '2'
character(len=*), parameter :: C_ITALIC         = '3'
character(len=*), parameter :: C_UNDERLINED     = '4'
character(len=*), parameter :: C_INVERTED       = '7'
character(len=*), parameter :: C_STRIKED        = '9'
character(len=*), parameter :: C_BLACK          = '30'
character(len=*), parameter :: C_RED            = '31'
character(len=*), parameter :: C_GREEN          = '32'
character(len=*), parameter :: C_YELLOW         = '33'
character(len=*), parameter :: C_BLUE           = '34'
character(len=*), parameter :: C_MAGENTA        = '35'
character(len=*), parameter :: C_CYAN           = '36'
character(len=*), parameter :: C_WHITE          = '37'
character(len=*), parameter :: C_MARKED_BLACK   = '40'
character(len=*), parameter :: C_MARKED_RED     = '41'
character(len=*), parameter :: C_MARKED_GREEN   = '42'
character(len=*), parameter :: C_MARKED_YELLOW  = '43'
character(len=*), parameter :: C_MARKED_BLUE    = '44'
character(len=*), parameter :: C_MARKED_MAGENTA = '45'
character(len=*), parameter :: C_MARKED_CYAN    = '46'
character(len=*), parameter :: C_MARKED_WHITE   = '46'
character(len=*), parameter :: C_CLEAR          = C_START // '0' // C_END

contains

    function format_str( str, code ) result( out )
        character(len=*), intent(in) :: str
        character(len=*), intent(in) :: code
        character(len=:), allocatable :: out
        out = C_START // CODE // C_END // str // C_CLEAR
    end function format_str

end module simple_ansi_ctrls
