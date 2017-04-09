
#if defined  _WIN32
#define DEV_NULL "nul"
#else
#define DEV_NULL "/dev/null"
#endif

  
#define TBLOCK(XX) \
XX = tic();
!  write(*,'(A)') "Begin TBLOCK"

#define TSTOP(XX) \
  write(*,'(A,A,1i4,A,1d20.10)') __FILE__,":",__LINE__,":Elapsed time (s) ", toc()

program simple_timer_basic_test
     use precision_m
     use simple_timer
     implicit none
     integer(dp),parameter :: nrep = 500000000
     real(dp) :: xx, c, cfac, b, etime
     integer(dp) ::  t1, t2
     integer(dp) :: i

     c = .1
     cfac = .25
     b = 1.
     xx = 12.0_dp



    
!tstart = call tic(timer)
!do  i = 1,nrep
!  c = cfac*c +b
!end do
!tend = toc(timer,tstart)
!etime =  elapsedTime(timer,tend,tstart)
!print *, 'Time for local evaluation = ', etime, ' seconds'
! tstart = functic(timer)

     write (*, '(A)') "Welcome to the HighRes Fortran Timer\n"
     t1= get_clock_rate()
     print *, 'Clock rate: ', t1
     write (*, '(A)') "1. Simple timestamp and diff "
     t1 = tic()
     write (*, '(A,1d20.10)') "   t1 = ", real(t1, dp)
     do i = 1, nrep
        c = cfac*c + b
     end do
     t2 = tic()
     write (*, '(A,1d20.10)') "   t2 = ", real(t2, dp)
     etime = tdiff( t2, t1)
     write (*, '(A,1d20.10)') 'Time for simple evaluation (s) = ', etime

          call reset_timer()
     write(*,"(A)") "  "
     write (*, '(A)') "2. Simple tic toc usage "
     write (*, '(A)') "2. t1=TIC; etime=TOC(T1)"
     write (*, '(A)') "2. Calling tic(timer)"
     t1 = tic()
     write (*, '(A,1d20.10)') "2. t1 = ", real(t1, dp)
     c=.1
     do i = 1, nrep
        c = cfac*c + b
     end do
     etime = toc(t1)
     write (*, '(A,1d20.10)') "2.   toc(t1) = ", etime


          call reset_timer()
     write (*, '(A)') "3. Simple tic toc "
     t1 = tic()
     write (*, '(A,1d20.10)') "3.   t1 = ", real(t1, dp)
     c=.1
     do i = 1, nrep
        c = cfac*c + b
     end do
     etime = toc()
     write (*, '(A,1d20.10)') "3.   toc() = ", etime


     call reset_timer()
     write (*, '(A)') "4.  Testing return of toc in write cmd "
     t1 = tic()
     write (*, '(A,1d20.10)') "4.  t1 = ", real(t1, dp)
     c=.1
     do i = 1, nrep
        c = cfac*c + b
     end do
     write (*, '(A,1d20.10)') "4.  toc in write ", toc(t1)

         call reset_timer()
     write (*, '(A)') "5.  Testing empty tic() lhs "
     open(10,file=DEV_NULL)
     write (10,'(A,1i20)') "5.   result tic ", tic()
    ! write (*, '(A,1d20.10)') "2. t1 = ", real(t1, dp)
     c=.1
     do i = 1, nrep
        c = cfac*c + b
       
     end do
     write (*, '(A,1d20.10)') "5.   tic/toc in write ", toc()

     call reset_timer()
     write (*, '(A)') "6.  Testing tic toc in preprocessor macro "

     TBLOCK(t1)
     c=.1
     do i = 1, nrep
        c = cfac*c + b
     end do
     TSTOP(t1)

  end program simple_timer_basic_test
