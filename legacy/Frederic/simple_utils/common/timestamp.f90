!*******************************************************************************
!     Author: Frederic D.R. Bonnet date: 08th of August 2013.
!
! Name:
! timestamp - Various utilities and matrix getter for other modules.
!
! DESCRIPTION
! subroutine timestamp prints the current date and time.
!*******************************************************************************
! SOURCE
subroutine timestamp()
  implicit none

  character(len=8)      :: date
  character(len=10)     :: time
  character(len=5)      :: zone
  integer,dimension(8)  :: values
  integer               :: year, month, day
  integer               :: UTC, hour, min, sec, msec
  character(len=9),parameter,dimension(12) :: all_month = (/ &
       'January  ', 'February ', 'March    ', 'April    ', &
       'May      ', 'June     ', 'July     ', 'August   ', &
       'September', 'October  ', 'November ', 'December ' /)

  call date_and_time ( date, time, zone, values )

  year  = values(1)
  month = values(2)
  day   = values(3)
  UTC   = values(4)
  hour  = values(5)
  min   = values(6)
  sec   = values(7)
  msec  = values(8)

  write (*,'(1x,a,1x,a,1x,i2,1x,i4,2x,i2,a1,i2.2,a1,i2.2,a1,i3.3)') &
       'current date and time:',trim (all_month(month)), day, year, &
       hour,':', min, ':', sec, '.', msec

  return
end subroutine timestamp
