program speed
!    program to time 3 different methods of calculating
!    the dot product of long vectors
      intrinsic system_clock
      integer nrep,i , iii
      parameter (nrep=50000000)
      real   b,c,cfac
      common b,c,cfac
      integer,parameter :: tmax=8
      real, dimension(tmax) :: timings
      integer :: iti
#ifdef DOUBLE_CLOCK
      integer(8) :: icount1,irate,icmax,icount2
#endif
      character(len=160) :: t_char
      interface saxpy
         function saxpy (a,x,y)
            real x,y,a,saxpy
!    Declare x,y,a as input only to saxpy
            intent (in) x,y,a
         end function saxpy
      end interface saxpy
!    The following is a Statement function
      slaxpy(x,y,z)= x*y+z
      iti=8
      cfac=.25
      b=1.
      c=.1
!      allocate(timings(iti))
      iti=1
!    Use the Fortran 90 Intrinsic Subroutine to determine
!    the current time in clock ticks (icount2), the clock
!    rate in clicks per second, and the largest possible
!    count before the clock resets.

!    CAUTION:  On this and probably most Unix work stations
!    this clock is measuring real time, not time your program
!    spends running.  If you are sharing the machine with others
!    you will count the time they have the CPU also.  USE the
!    unix "users" command to check for a situation when you are
!    the only user on the machine before running the program.
!    To filter out system activity run the program many times
!    and select results with the lowest total times.

      call system_clock(icount1,irate,icmax)

      !      print *, 'clock rate = ',irate, ' ticks per second'
      timings(iti) = REAL(irate)/REAL(1.0E+06); iti=iti+1;
      call system_clock(icount1,irate,icmax)

!     The "1000" loop just makes sure that I do lots of work
!     to get good statistics.  Note that coding is set so
!     results of each pass through the loop are a little
!     different.  Without the "+.00001*dotpro" optimization
!     features on many compilers will only execute the loop
!     once.

!     Do a simple multiply and add many times
!     I am intentially suppressing any vectorization and forcing
!     a result that won't let the optimizer throw away unused
!     calculations

      do 1999 iii=1,4
         iti=2
      call system_clock(icount1,irate,icmax)
      do 1000 i = 1,nrep
         c = cfac*c +b
 1000    continue
      call system_clock(icount2,irate,icmax)
      time = real(icount2-icount1)/real(irate)
    !  print *, 'Time for local evaluation = ', time, ' seconds'
      timings(iti) = time; iti=iti+1;

!    Repeat the same work using a standard function

      call system_clock(icount1,irate,icmax)
      do 1001 i = 1,nrep
         c = saxpy (cfac, c,b)
 1001    continue
      call system_clock(icount2,irate,icmax)
      time = real(icount2-icount1)/real(irate)
      timings(iti) = time; iti=iti+1; !print *, 'Time for function evaluation = ', time, ' seconds'

!    Repeat the same work using a Subroutine

      call system_clock(icount1,irate,icmax)
      do 1002 i = 1,nrep
         call ssaxpy(cfac,c,b,c)
 1002    continue
      call system_clock(icount2,irate,icmax)
      time = real(icount2-icount1)/real(irate)
      timings(iti) = time; iti=iti+1; !print *, 'Time for Subroutine = ', time, ' seconds'

!    Repeat the same work using a standard function

      call system_clock(icount1,irate,icmax)
      do 1003 i = 1,nrep
         c = slaxpy(cfac,c,b)
 1003    continue
      call system_clock(icount2,irate,icmax)
      time = real(icount2-icount1)/real(irate)
      timings(iti) = time; iti=iti+1; !print *, 'Time for Statement Function = ', time, ' seconds'

!    Repeat the same work using an internal function

      call system_clock(icount1,irate,icmax)
      do 1004 i = 1,nrep
         c = siaxpy(cfac,c,b)
 1004    continue
      call system_clock(icount2,irate,icmax)
      time = real(icount2-icount1)/real(irate)
      timings(iti) = time; iti=iti+1; !print *, 'Time for Internal Function = ', time, ' seconds'

!    Repeat the same work using a Subroutine with
!    passing through COMMON

      call system_clock(icount1,irate,icmax)
      do 1005 i = 1,nrep
         call scsaxpy
 1005    continue
      call system_clock(icount2,irate,icmax)
      time = real(icount2-icount1)/real(irate)
      timings(iti) = time; iti=iti+1; !print *, 'Time for Subroutine with COMMON = ', time, ' seconds'

!    Do a baseline with just the DO loop and increment

      call system_clock(icount1,irate,icmax)
      do 1006 i = 1,nrep
 1006    continue
      call system_clock(icount2,irate,icmax)
      time = real(icount2-icount1)/real(irate)
      timings(iti) = time; iti=iti+1; !print *, 'Time for bare Loop = ', time, ' seconds'
1999 continue

      write(*,"(8A20)") "Ticks/s",  "Local fn (s)", "Standard fn (s)", &
           "Subroutine (s)",  "Statement (s)", "Internal fn (s)", &
           "COMMON (s)", "Empty Loop (s)"
      

      write(t_char,'(8es20.8)') timings
      print '(a160)', adjustr(t_char)
      !i=2
      !write(*,"(g12.0)",advance="no") timings(1)
      !write(*,"(6g10.10)",advance="no") (timings(i),i=2,7)
!      timings(2), timings(3), timings(4), &
!           timings(5), timings(6), timings(7)



      write(*,*)'-----------------------------------------------'

      !      print *, c,i

      stop

    contains
      function siaxpy(a1,x1,y1)
         real a1,x1,y1,siaxpy
         siaxpy =  a1*x1 + y1
         return
         end function
      end 

      function saxpy(a,x,y)

!   Multply all contents of  "x" by the scalar "a"
!   then add the result to  "y"

!   John Mahaffy    4/3/96

      implicit none
      real a,x,y,saxpy
      intent (in) a,x,y
      saxpy =  a*x + y
      return
      end
      subroutine ssaxpy(a,x,y,z)

!   Multply all contents of  "x" by the scalar "a"
!   then add the result to  "y"

!   John Mahaffy    4/3/96

      implicit none
      real a,x,y,z
      z     =  a*x + y
      return
    end 
      subroutine scsaxpy

!   Multply all contents of  "c" by the scalar "cfac"
!   then add the result to  "b"

!   John Mahaffy    4/3/96

      implicit none
      real   b,c,cfac
      common b,c,cfac
      c     =  cfac*c + b
      return
    end 

      
