!------------------------------------------------------------------------------!
! SIMPLE v3.0         Elmlund & Elmlund Lab          simplecryoem.com          !
!------------------------------------------------------------------------------!
!> Simple timer module: High resolution (nanoseconds) timer in Fortran
!! See simple_timer.h for useful macros
!!
!!     use simple_timer
!!     integer(8)  :: tn;
!!     tn = tic()
!!     ...
!!     write(*,'(A,F20.10)') ">>> Elapsed time (sec) ", toc()
!! 
!! \author Michael Eager 2017
!! REVISIONS:
!! Version 0.1:  64 bit INT implementation of system_clock
!!  present in gfortran and pgfortran
!! Version 0.2: Special loop and profiling timers added May 2017
!!
! The SIMPLE code is distributed with the hope that it will be
! useful, but WITHOUT ANY WARRANTY. Redistribution and modification is regulated
! by the GNU General Public License.
! -----------------------------------------------------------------------------!

module simple_timer
!     use simple_jiffys ! singleton
use simple_defs   ! singleton, fp_kind declared
!   use precision_m
implicit none
!  private :: raise_sys_error
#ifndef TIMER_I4
integer, parameter,public :: timer_int_kind = fp_kind
#else
integer, parameter,public :: timer_int_kind = sp
#endif
#include "simple_local_flags.inc"
private
   integer(timer_int_kind), public   :: clock_ticks_per_second = INT(0, timer_int_kind) !< Number of counts per second
   integer(timer_int_kind), public   :: last_time_point = INT(0, timer_int_kind)        !< Current timesamp
   integer, public       :: idx_elapsed = 0, num_elapsed = 3
   integer, public       :: num_profile_loops, num_profile_vars
   logical, public       :: in_loop = .false.
   integer, public       :: profile_counter
   integer, parameter, public :: MAX_TOKENS = 10                          !< number of entries in token dictionary
   integer, parameter, public :: MAX_TOKEN_CHARSIZE = 30                  !< max char length of tokens
   real(timer_int_kind), allocatable, public :: elapsed_times(:)                 !< timer loop storage of elapsed times
   real(timer_int_kind), allocatable, public :: profile_matrix(:, :)             !< profiler storage of elapsed times for each token
   integer(timer_int_kind), dimension(MAX_TOKENS), public :: profile_last_timerstamp !< storage of time stamps for profiler
   character(len=MAX_TOKEN_CHARSIZE), dimension(MAX_TOKENS), public :: profile_labels = ""

   public :: gettime,cast_time_char
   public :: tic, tickrate
   public :: toc, tdiff, tocprint
   public :: now, reset_timer
   public :: timer_loop_start, in_timer_loop, timer_loop_end
   public :: timer_profile_setup, timer_profile_start, timer_profile_break, timer_profile_report
   public ::  tic_i4, toc_i4, tdiff_i4, tickrate_i4  !! Only for testing 32-bit system_clock (1 ms resolution)
contains

integer function gettime ()
#ifdef PGI
    include 'lib3f.h'          ! time
    gettime=time()
#elif defined(INTEL)
    use ifport
    !integer :: gettime
    call time(gettime)
#else
    gettime= time()
#endif
end function gettime

function cast_time_char (arg)
    character(len=24) :: cast_time_char
    integer, intent(in) :: arg
#ifdef PGI
    include 'lib3f.h'          ! time
    cast_time_char=ctime(arg)
#elif defined(INTEL)
    !! TODO fix intel ctime
    use ifport
    integer :: now
    call time(now)
#else
    cast_time_char=ctime(arg)
#endif
    

end function cast_time_char


!< Force timestamps and clock rate to zero
   subroutine reset_timer
      last_time_point = INT(0, timer_int_kind)
      if (allocated(elapsed_times)) deallocate (elapsed_times)
   end subroutine reset_timer

!< Get system_clock timestamp
   integer(timer_int_kind) function tic()
      call system_clock(count=tic)
      last_time_point = tic
   end function tic


!< Get the clock tick count per second
   integer(timer_int_kind) function tickrate()
      tickrate = INT(0, timer_int_kind)
      if (clock_ticks_per_second .eq. 0) call system_clock(count_rate=tickrate)
      clock_ticks_per_second = tickrate
#ifdef _DEBUG
      write (*, '(A,1d20.10)') " CLOCK_RATE(ticks/sec) ", REAL(clock_ticks_per_second, timer_int_kind)
#endif
   end function tickrate

!< Calculate the time from two timestamps
   real(timer_int_kind) function tdiff(tfinal, tstart)
      integer(timer_int_kind), intent(in) :: tfinal
      integer(timer_int_kind), intent(in) :: tstart
! integer(timer_int_kind)                       ::  end_point
! if(present(tstart)) last_time_point = tstart
! if(.not. present(tfinal)) then
!    call system_clock(count=end_point)
!    tfinal = end_point
! end if
      if (clock_ticks_per_second .eq. INT(0, timer_int_kind)) call system_clock(count_rate=clock_ticks_per_second)
! Calulate the time difference
      tdiff = REAL(tfinal - tstart, timer_int_kind)/REAL(clock_ticks_per_second, timer_int_kind)
!      last_time_point = tfinal
   end function tdiff

!< Complete the timing regime using a reference timestamp or the one
!  in last_time_point
   real(timer_int_kind) function toc(tstart)
      integer(timer_int_kind), intent(in), optional ::  tstart
      integer(timer_int_kind)                       ::  end_point
      if (present(tstart)) last_time_point = tstart
      call system_clock(count=end_point)
      toc = tdiff(end_point, last_time_point)
      last_time_point = end_point
   end function toc

!< Complete the timing regime using a reference timestamp or the one
!  in last_time_point
   subroutine tocprint(tstart, comment)
      character(len=*), intent(inout), optional :: comment
      integer(timer_int_kind), intent(in), optional ::  tstart
      integer(timer_int_kind)                     ::  end_point
      real(timer_int_kind)                        :: elapsed
      if (.not. present(comment)) comment = " Simple timer "

      if (present(tstart)) last_time_point = tstart
      call system_clock(count=end_point)
#ifdef _DEBUG
      write (*, '(A,1ES20.10)') " TOC Time stamp ", REAL(end_point, timer_int_kind)
#endif
      elapsed = tdiff(end_point, last_time_point)
      last_time_point = end_point
      write (*, '(A,A,1ES20.10)') trim(comment), " Elapsed time ", elapsed
   end subroutine tocprint

!> print current time and date
   subroutine now
      character(len=8)  :: date
      character(len=10) :: time
      print *, "System_clock: ", tic()
      call date_and_time(date, time)
      write (*, '(A,A,A,A,A,A,A)') 'Date: ', date(7:8), '-', date(5:6), '-', date(1:4), '\n'
      write (*, '(A,A,A,A,A,A,A)') 'Time: ', time(1:2), ':', time(3:4), ':', time(5:10), '\n'
   end subroutine now

!< in_timer_loop checks the time within a timer loop
! It does not start the timer or set the in_loop variable
! \return  Returns false on the final loop- so that timer_loop_end can finish
   logical function in_timer_loop()
      in_timer_loop = .false.
      if (.not. in_loop) then
         print *, "Failed timer_loop: Timer loop did not start"
      else
         if (idx_elapsed .lt. num_elapsed) then
            elapsed_times(idx_elapsed) = toc()
            idx_elapsed = idx_elapsed + 1
            in_timer_loop = .true.
         end if
      end if
   end function in_timer_loop
!<Begin timer loop
   subroutine timer_loop_start(num)
      integer, intent(in), optional :: num
      integer(timer_int_kind)  :: dummytimestamp
      num_elapsed = 3
      if (present(num)) num_elapsed = num
      call reset_timer()
      if (num_elapsed .gt. 1) then
         allocate (elapsed_times(num_elapsed))
         idx_elapsed = 1
         dummytimestamp = tic()
         in_loop = .true.
      end if
DebugPrint 'Size of elapsed array ', size(elapsed_times)

   end subroutine timer_loop_start
!< end timer loop
   subroutine timer_loop_end(COMMENT)
      character(len=*), intent(in), optional :: COMMENT
      character(len=128) :: strcomment
      if (.not. present(COMMENT)) then
         strcomment = ' no comment '
      else
         if (len_trim(COMMENT) .le. 128) then
            strcomment = trim(adjustl(COMMENT))
         else
            stop "Timer loop error - comment string must be less than 128 characters"
         end if
      end if
      if (.not. in_loop) then
         print *, "Failed timer_loop_end: Timer loop did not start"
      else
         write (*, '(A,A)') "******* TIMER LOOP ", trim(strcomment)
         write (*, '(A,1i8)') '*** Iterations:  ', num_elapsed
         if (idx_elapsed .eq. num_elapsed) then
            elapsed_times(idx_elapsed) = toc()
            write (*, '(A,1ES20.10)') "*** Average (sec):", &
               SUM(elapsed_times, DIM=1)/REAL(num_elapsed, timer_int_kind)
            write (*, '(A,1ES20.10,A,1i3)') "*** Longest run(sec) ", &
               MAXVAL(elapsed_times, DIM=1), &
               '    at ', MAXLOC(elapsed_times, DIM=1)
            write (*, '(A,1ES20.10,A,1i3)') "*** Shortest run(sec) ", &
               MINVAL(elapsed_times, DIM=1), &
               '   at ', MINLOC(elapsed_times, DIM=1)
         else
            write (*, '(A,1i8)') '*** Failed at iteration ', idx_elapsed
         end if
         write (*, '(A)') "******* TIMER LOOP **************"
         in_loop = .false.
         if (allocated(elapsed_times)) then
            deallocate (elapsed_times)
         end if
      end if
   end subroutine timer_loop_end

   !< Setup profiling
   subroutine timer_profile_setup(nLoops, nVars, vin)
   use simple_strings
      integer,          intent(in)    :: nLoops
      integer,          intent(in)    :: nVars
      character(len=*), intent(inout) :: vin
      character(len=20), allocatable  :: v(:)
      integer :: nargs_parse, ind, tmp_nargs=0
      DebugPrint " timer_profile_setup ", char(nLoops), "  ", char(nVars)
      DebugPrint vin
      ! define number of tokens
      if(nVars .le. 0)then
          call removepunct(vin)
          nargs_parse = cntRecsPerLine(vin)
      else
          nargs_parse = nVars
      end if

      if (nLoops .lt. 1) then

          DebugPrint "timer_profile_setup error -- must have more than 1 loop"
         stop
      elseif (nargs_parse .gt. MAX_TOKENS .or. nargs_parse .le. 0) then
          print*, "timer_profile_setup  nargs_parse error -- outside range"
          stop
      elseif (len_trim(vin)==0 )then
          DebugPrint "timer_profile_setup error -- token string is empty"
          stop
      else
          call removepunct(vin)

          DebugPrint " timer_profile_setup remove punct"
          DebugPrint vin


          allocate(v(nargs_parse))
          ind = index(vin, ',')
          if( ind == 0 ) then
              call parse(vin,' ', v,nargs_parse)
              DebugPrint " timer_profile_setup no-comma token input"
              DebugPrint vin
          else 
              call parse(vin,',',v,tmp_nargs)
          end if
          DebugPrint " timer_profile_setup parsed tokens"
          DebugPrint v
          if (nargs_parse .ne. tmp_nargs .or. size(v,1) .ne. nargs_parse) then
              DebugPrint "timer_profile_setup error -- parsing token string error ", nargs_parse, size(v,1)
              stop
          end if
          DebugPrint  " timer_profile_setup parsed tokens OK"
          DebugPrint  v

         ! profile_labels are a fixed size
         if (nVars .ge. 1 .and. v(1) .ne. "") then
            if (len_trim(v(1)) .le. MAX_TOKEN_CHARSIZE) then
               profile_labels(1) = trim(adjustl(v(1)))
               num_profile_vars = 1
            else
               stop 'Error: Timer profile token 1 too long'
            end if
         end if
         if (nVars .ge. 2 .and. v(2) .ne. "") then
            if (len_trim(v(2)) .le. MAX_TOKEN_CHARSIZE) then
               profile_labels(2) = trim(adjustl(v(2)))
               num_profile_vars = 2
            else
               stop 'Error: Timer profile token 2 too long'
            end if
         end if
         if (nVars .ge. 3 .and. v(3) .ne. "") then
            if (len_trim(v(3)) .le. MAX_TOKEN_CHARSIZE) then
               profile_labels(3) = trim(adjustl(v(3)))
               num_profile_vars = 3
            else
               stop 'Error: Timer profile token 3 too long'
            end if
         end if
         if (nVars .ge. 4 .and. v(4) .ne. "") then
            if (len_trim(v(4)) .le. MAX_TOKEN_CHARSIZE) then
               profile_labels(4) = trim(adjustl(v(4)))
               num_profile_vars = 4
            else
               stop 'Error: Timer profile token 4 too long'
            end if
         end if
         if (nVars .ge. 5 .and. v(5) .ne. "") then
            if (len_trim(v(5)) .le. MAX_TOKEN_CHARSIZE) then
               profile_labels(5) = trim(adjustl(v(5)))
               num_profile_vars = 5
            else
               stop 'Error: Timer profile token 5 too long'
            end if
        end if
        if (nVars .ge. 6 .and. v(6) .ne. "") then
            if (len_trim(v(6)) .le. MAX_TOKEN_CHARSIZE) then
               profile_labels(6) = trim(adjustl(v(6)))
               num_profile_vars = 6
            else
               stop 'Error: Timer profile token 6 too long'
            end if
         end if
         if (nVars .ge. 7 .and. v(7) .ne. "") then
            if (len_trim(v(7)) .le. MAX_TOKEN_CHARSIZE) then
               profile_labels(7) = trim(adjustl(v(7)))
               num_profile_vars = 7
            else
               stop 'Error: Timer profile token 7 too long'
            end if
         end if
         if (nVars .ge. 8 .and. v(8) .ne. "") then
            if (len_trim(v(8)) .le. MAX_TOKEN_CHARSIZE) then
               profile_labels(8) = trim(adjustl(v(8)))
               num_profile_vars = 8
            else
               stop 'Error: Timer profile token 8 too long'
            end if
         end if
         if (nVars .ge. 9 .and. v(9) .ne. "") then
            if (len_trim(v(9)) .le. MAX_TOKEN_CHARSIZE) then
               profile_labels(9) = trim(adjustl(v(9)))
               num_profile_vars = 9
            else
               stop 'Error: Timer profile token 9 too long'
            end if
         end if
         if (nVars .ge. 10 .and. v(10) .ne. "") then
            if (len_trim(v(10)) .le. MAX_TOKEN_CHARSIZE) then
               profile_labels(10) = trim(adjustl(v(10)))
               num_profile_vars = 10
            else
               stop 'Error: Timer profile token 10 too long'
            end if
         end if
      end if
      if (nVars .ne. num_profile_vars) then
         stop 'timer profile setup error: vars input > internal num_profile_vars'
      end if
      num_profile_loops = nLoops
      if (allocated(profile_matrix)) deallocate (profile_matrix)
!if (allocated(profile_last_timerstamp)) deallocate (profile_last_timerstamp)
! allocate (profile_last_timerstamp(num_profile_loops))
      allocate (profile_matrix(num_profile_loops, num_profile_vars))
      profile_matrix = REAL(0.0, timer_int_kind)
      profile_last_timerstamp = INT(0, timer_int_kind)
      DebugPrint  " Profile matrix size ", size(profile_matrix, 1), size(profile_matrix, 2)
      DebugPrint profile_matrix(1:10, 1:2)
   end subroutine timer_profile_setup

   !< Within profile loop - start timer with token 
   subroutine timer_profile_start(token)
      character(len=*), intent(in) :: token
      integer ::  ival
      do ival = 1, num_profile_vars
         if (.not. (INDEX(profile_labels(ival), trim(adjustl(token))) == 0)) then
            profile_last_timerstamp(ival) = tic()
            exit
         end if
      end do

      if (ival .gt. num_profile_vars) then
         write (*, '(A,A,A,1i10)') "Error Timer_Profile_start:", &
            trim(adjustl(token)), " label index outside range ", ival
#ifdef _DEBUG
      else
         DebugPrint "Label: ", profile_labels(ival), " time stamp "
#endif
      end if
   end subroutine timer_profile_start

   !< Within profile loop - get elapsed time for token 'token' and reset
   subroutine timer_profile_break(token)
      character(len=*), intent(in) :: token
      integer(timer_int_kind) :: tmp_tstamp = INT(0, timer_int_kind)
      integer ::  ival = 0, iloop = 0
!
! Need bounds checking of timestamps and matrix
      do ival = 1, num_profile_vars
         iloop = 0
         if (.not. (INDEX(profile_labels(ival), trim(adjustl(token))) == 0)) then
             DebugPrint  'Timer profile break: Found label ', profile_labels(ival)
            do iloop = 1, num_profile_loops
               if (profile_matrix(iloop, ival) .eq. 0) then
                  tmp_tstamp = tic()
                  profile_matrix(iloop, ival) = tdiff(tmp_tstamp, profile_last_timerstamp(ival))
                  profile_last_timerstamp(ival) = tmp_tstamp
                  exit
               end if
            end do
            exit
         end if
      end do
      if (tmp_tstamp .eq. 0) then
         write (*, '(A,2i8,A)') "Error profile: No time stamp created. loop,val,label:", &
            iloop, ival, trim(adjustl(token))
      end if
      if (ival .gt. num_profile_vars + 1) then
         write (*, '(A,1i8)') "Error profile: label index outside range ", ival
      end if
      if (iloop .gt. num_profile_loops + 1) then
         write (*, '(A,1i8)') "Error profile: loop index outside range ", iloop
      end if
#ifdef _DEBUG

      if ((ival .gt. num_profile_vars) .or. (iloop .gt. num_profile_loops)) then
         write (*, '(A,2i10)') "Timer_Profile_break: label/loop index outside range ", ival, iloop
      else
         DebugPrint "Label: ", profile_labels(ival), " time ", profile_matrix(iloop, ival)
      end if

#endif
   end subroutine timer_profile_break

!< Profile report
   subroutine timer_profile_report(COMMENT, totaltime)
      character(len=*), intent(in) :: COMMENT
      real(timer_int_kind), intent(in):: totaltime
      real(timer_int_kind)            :: total_avg
      real(timer_int_kind), allocatable:: avgtime_token(:)
      integer :: ival, iloop
!    if (.not.present(COMMENT)) COMMENT="PROFILE"

      if (allocated(profile_matrix)) then
         total_avg = totaltime/REAL(num_profile_loops, timer_int_kind)
         allocate (avgtime_token(num_profile_vars))
         do ival = 1, num_profile_vars
            avgtime_token(ival) = SUM(profile_matrix(:, ival), &
                                      MASK=(profile_matrix(:, ival) /= 0.), DIM=1)/REAL(num_profile_loops, timer_int_kind)
         end do

         write (*, '(A,A,A)') "** PROFILE REPORT : ", trim(adjustl(COMMENT))
         write (*, '(A,A,A,1i4)') '** FILE:LINE: ', __FILE__, ":", __LINE__
         write (*, '(A,1i8,A)') '** Iterations: ', num_profile_loops, ' timed loops'

         do ival = 1, num_profile_vars
            write (*, '(A,A)') '**** Label name: ', trim(profile_labels(ival))
            write (*, '(A,1ES20.6,2x,1ES20.5,A)') "**** Average (sec):", &
               avgtime_token(ival), avgtime_token(ival)/total_avg, '%'
            if (num_profile_loops .gt. 1) then
               write (*, '(A,1ES20.6,A,1i10)') "**** Longest run(sec) ", &
                  MAXVAL(profile_matrix(:, ival), MASK=(profile_matrix(:, ival) /= 0.), DIM=1), &
                  '    at ', MAXLOC(profile_matrix(:, ival), MASK=(profile_matrix(:, ival) /= 0.), DIM=1)
               write (*, '(A,1ES20.6,A,1i10)') "**** Shortest run(sec) ", &
                  &REAL(MINVAL(profile_matrix(:, ival), MASK=(profile_matrix(:, ival) /= 0.0_dp), DIM=1),dp), &
                  &'   at ', INT(MINLOC(profile_matrix(:, ival), MASK=(profile_matrix(:, ival) /= 0.), DIM=1),sp)
            end if
         end do
         write (*, '(A,1d20.6)') "** Total time (sec):", totaltime
         write (*, '(A,1d20.6)') "** Average iteration (sec):", total_avg

         write (*, '(A,A,A)') "******* END ", trim(COMMENT), " REPORT **************"
         deallocate (profile_matrix,avgtime_token)
      end if
! unset labels
      profile_labels = ""
      profile_last_timerstamp = INT(0, timer_int_kind)
   end subroutine timer_profile_report


!!! Testing I4 system_clock


   !< Get system_clock timestamp
   integer(kind=4) function tic_i4()
       call system_clock(count=tic_i4)
        last_time_point = INT(tic_i4, timer_int_kind)
   end function tic_i4


   !< Get the clock tick count per second
   integer(kind=4) function tickrate_i4()
      tickrate_i4 = INT(0, kind=4)
      call system_clock(count_rate=tickrate_i4)
      clock_ticks_per_second = INT(tickrate_i4, timer_int_kind)
   end function tickrate_i4

   !< Calculate the time from two timestamps
   real function tdiff_i4(tfinal, tstart)
      integer(kind=4), intent(in) :: tfinal
      integer(kind=4), intent(in) :: tstart
      integer(kind=4) :: ticks_per_second
      call system_clock(count_rate=ticks_per_second)
      ! Calulate the time difference
      tdiff_i4 = REAL(tfinal - tstart, kind=4)/REAL(ticks_per_second, kind=4)
   end function tdiff_i4

   !< Complete the timing regime using a reference timestamp or the one
   !  in last_time_point
   real function toc_i4(tstart)
      integer(kind=4), intent(inout), optional ::  tstart
      integer(kind=4)                          ::  end_point
      if (.not.present(tstart))  tstart = INT(last_time_point, kind=4)
      call system_clock(count=end_point)
      toc_i4 = tdiff_i4(end_point, tstart)
      last_time_point = INT(end_point, timer_int_kind)
   end function toc_i4


end module simple_timer

