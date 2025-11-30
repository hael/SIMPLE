module simple_defs_string
use, intrinsic :: iso_c_binding
implicit none
character(len=*),             parameter :: LOWER_CASE_LETTERS = 'abcdefghijklmnopqrstuvwxyz'
character(len=*),             parameter :: UPPER_CASE_LETTERS = 'ABCDEFGHIJKLMNOPQRSTUVWXYZ'
character(len=*),             parameter :: INTEGERS = '1234567890'
character(kind=c_char,len=*), parameter :: NEW_LINES_C = C_FORM_FEED // C_NEW_LINE // C_CARRIAGE_RETURN // C_VERTICAL_TAB
character(kind=c_char,len=*), parameter :: BLANK_C_CHARACTERS = C_NULL_CHAR // C_HORIZONTAL_TAB
character(len=*),             parameter :: BLANK_CHARACTERS =  ' '//BLANK_C_CHARACTERS
character(len=*),             parameter :: BLANKS_AND_NEW_LINES = BLANK_CHARACTERS // NEW_LINES_C
end module simple_defs_string