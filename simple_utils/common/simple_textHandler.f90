! 
!*******************************************************************************
!     Author: Frederic D.R. Bonnet date: 4th of MarchOctober 2015.
!
! Name:
! simple_textHandler - Various utilities and handling routines and functions
!                       to handle text input-output on CPU for the front end
!                       and data treatement.
!
! Description:
! simple_textHandler provides routine for handling text from input-output 
! of data sets on CPU.
!*******************************************************************************
!
module simple_textHandler

  use simple_defs

  implicit none

  interface
     !TODO: insert the interface routine here.
  end interface

contains
!*******************************************************************************
! DESCRIPTION
! Convert the input string to lower case.
!*******************************************************************************
! SOURCE
  function lower_case(string) result (output)
    implicit none

    !global variables
    character(len = *), intent(in) :: string
    character(len = len(string))   :: output
    ! Local variables.
    integer                        :: n, length
    !counters
    integer                        ::  i

    !start of the execution commands

    length = len(string)

    do i = 1, length
       n = ichar(string(i:i))            !character to integer value
       if ((n >= 65).and.(n <= 91)) n = n+32
       output(i:i) = char(n)
    end do
  end function lower_case
!*******************************************************************************
! DESCRIPTION
! Convert the input string to upper case.
!*******************************************************************************
! SOURCE
  function upper_case(string) result (output)
    implicit none

    !global variables
    character(len = *), intent(in) :: string
    character(len = len(string))   :: output
    ! Local variables.
    integer                        :: n, length
    !counters
    integer                        :: i

    !start of the execution commands

    length = len(string)

    do i = 1, length
      n = ichar(string(i:i))    !character to integer value
      if ((n >= 97).and.(n <= 123)) n = n-32
      output(i:i) = char(n)
    end do
  end function upper_case

end module simple_textHandler
