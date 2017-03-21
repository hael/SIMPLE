!==Module simple_timer
!
! simple_timer is a little module for calculating the relative and actual CPU-time.
! The code is distributed with the hope that it will be useful, but _WITHOUT_ _ANY_ _WARRANTY_.
! Redistribution or modification is regulated by the GNU General Public License.
! *Author:* Hans Elmlund, 2009-10-01.
! 
!==Changes are documented below
!
!* incorporated in the _SIMPLE_ library, HE 2009-10-01
!
module simple_timer
use simple_jiffys ! singleton
use simple_defs   ! singleton
implicit none

private :: raise_sys_error

contains

  double precision  :: start = 0
  double precision  :: finish = 0


  double precision function tic( time )
    double precision                  :: time(2)
    double precision,save :: last_time = 0
    double precision      :: this_time
        intrinsic SYSTEM_CLOCK_TIME
        call SYSTEM_CLOCK_TIME(this_time)
        time(1) = real(this_time-last_time)
        time(2) = 0.
        dtime = time(1)
        last_time = this_time
    end function tic


    !
    real function toc( time )
        real :: time(2)
        call SYSTEM_CLOCK_TIME(etime)
        time(1) = etime
        time(2) = 0
    end function etime

    !> is for getting the actual cputime
    function getabscpu( lprint ) result( actual )
        logical, intent(in) :: lprint
        real                :: tarray(2)
        real                :: actual
        actual = etime( tarray )
        if( lprint )then
            write(*,'(A,2X,F9.2)') 'Actual cpu-time:', actual
        endif
    end function getabscpu



    PROGRAM simple_timer
      INTEGER , INTENT(OUT):: count, count_rate, count_max
      CALL SYSTEM_CLOCK(count, count_rate, count_max)
      WRITE(*,*) count, count_rate, count_max
    END PROGRAM

    
end module simple_timer

