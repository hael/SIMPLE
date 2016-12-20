!
!     cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!     Returns the significant length of a string.
!     Every character is significant with the 
!     exception of:
!     blank     (32)
!     null      (0)
!     reqd. routines - NONE
!
      function strlen(string)
      implicit none

      character(len=*) string
      integer strlen
      integer i, blank

      blank = ichar(' ')

      strlen = len(string)
      i = ichar(string(strlen:strlen))
      do while ((i.eq.blank .or. i.eq.0) .and. strlen.gt.0)
         strlen = strlen - 1
         i = ichar(string(strlen:strlen))
      end do

      return
      end function strlen
