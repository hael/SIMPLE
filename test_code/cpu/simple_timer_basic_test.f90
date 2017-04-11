
#if defined  _WIN32
#define DEV_NULL "nul"
#else
#define DEV_NULL "/dev/null"
#endif

#define NREP_MAX 10000000


#define CAT(prefix, suffix)            prefix ## suffix
#define _UNIQUE_LABEL(prefix, suffix)  CAT(prefix, suffix)
#define UNIQUE_LABEL(prefix)           _UNIQUE_LABEL(prefix, __LINE__)


#define TBLOCK() \
  t1 = tic();

#define TSTOP() \
  write(*,'(A,A,1i4,A,1d20.10)') __FILE__,":",__LINE__,":Elapsed time (s) ", toc(); 

#define TIMER_LOOP_START()    \
  integer(dp) :: timestamp; \
  real, dimension(3):: elapsed; \
  integer :: ii \
  do ii=1,3 \
     timestamp = tic();


#define TIMER_LOOP_END()           \
     elapsed(ii) = toc(timestamp); \
end do; \
write(*,'(A,A,1i4,A)') __FILE__,":",__LINE__, ' *** Timed loop *** '; \
write(*,'(A,1d20.10)') "  Average over 3 (sec):",  SUM(elapsed,DIM=1) / REAL(3.,dp); \
write(*,'(A,1d20.10)') "  Min time (sec) ", MINVAL(elapsed);


program simple_timer_basic_test
     use precision_m
     use simple_timer
     implicit none
     integer(dp),parameter :: nrep = NREP_MAX
     real(dp)    :: xx, c, cfac, b
     real(dp)    :: etime
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

     write (*, '(A)') "HighRes Fortran Timer"
     t1= get_clock_rate()
     print *, 'Clock rate: ', t1
     print *,"   "
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
     write (*, '(A,1d20.10)') '2.   toc(t1) = ', etime
     write(*,"(A)") "  "


          call reset_timer()
     write (*, '(A)') "3. Simple tic toc "
     t1 = tic()
     write (*, '(A,1d20.10)') "3.   t1 = ", real(t1, dp)
     c=.1
     do i = 1, nrep
        c = cfac*c + b
     end do
     etime = toc()
     write (*, '(A,1d20.10)') '3.   toc() = ', etime
     write(*,"(A)") ' '


     call reset_timer()
     write (*, '(A)') "4.  Testing return of toc in write cmd "
     t1 = tic()
     write (*, '(A,1d20.10)') "4.  t1 = ", real(t1, dp)
     c=.1
     do i = 1, nrep
        c = cfac*c + b
     end do
     write (*, '(A,1d20.10)') '4.  toc in write ', toc(t1)
     write(*,"(A)") "  "

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
     write(*,"(A)") ' '
     write (*, '(A)') '6.  Testing tic toc in preprocessor macro  '
     TBLOCK()
     c=.1
     do i = 1, nrep
        c = cfac*c + b
     end do
     TSTOP()
     write(*,"(A)") ' '

     call reset_timer()
     write(*,"(A)") ' '
     write (*, '(A)') '7.  Testing tic toc in preprocessor macro with subroutine  '
     TBLOCK()
     c = saxy(c)
     TSTOP()
     write(*,"(A)") ' '

     call reset_timer()
     write(*,"(A)") ' '
     write (*, '(A)') '8.  Testing timed_loop macro - must be in own subroutine '
     c = test_loop()




   contains



     function saxy(c_in) result(c)
       real(dp), intent(in) :: c_in
       integer(dp),parameter :: nrep = NREP_MAX
       real(dp)              :: xx, c, cfac, b
       integer(dp)           :: i

       c = c_in
       c=.1
       do i = 1, nrep
          c = cfac*c + b
       end do
     end function saxy

     ! Loop macro has to declare an array and an int
     ! -- this means it has to be at the top or in
     !    its own function
     function test_loop() result(cloop)
       real(dp) :: cloop
       integer(dp) :: timestamp
       real(dp), dimension(3):: elapsed
       integer :: ii
       do ii=1,3 
          timestamp = tic();

          cloop= REAL(.1,dp)
          cloop = saxy(cloop)


          elapsed(ii) = toc(timestamp)
       end do
       write(*,'(A,A,1i4,A)') __FILE__,":",__LINE__, ' *** Timed loop *** '
       write(*,'(A,1d20.10)') "  Average over 3 (sec):",  SUM(elapsed,DIM=1) / REAL(3.,dp)
       write(*,'(A,1d20.10)') "  Min time (sec) ", MINVAL(elapsed)
       write(*,"(A)") ' '

     end function test_loop


  end program simple_timer_basic_test
