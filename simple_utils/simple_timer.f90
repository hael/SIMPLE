!==Module simple_timer
!
!>\brief High resolution timer in fortran
!
!  64 bit INT implementation of system_clock
!  present in gfortran and pgfortran
!<------------------------------------
module simple_timer
!     use simple_jiffys ! singleton
!     use simple_defs   ! singleton
  use precision_m
  implicit none

!  private :: raise_sys_error
  private
  integer(dp),public   :: clock_ticks_per_second=INT(0,dp) !< Number of counts per second
  integer(dp),public   :: last_time_point=INT(0,dp) !< Current timesamp
  integer,public       :: idx_elapsed=0,num_elapsed=3
  integer,public       :: num_profile_loops,num_profile_vars
  logical,public       :: inloop=.false.
  integer,parameter,public :: ic_profile=5 ! number of entries in language dictionary
  integer,parameter,public :: lc_profile=30 ! max length of tokens
  real(dp),allocatable,public :: elapsed_times(:)
  real(dp),allocatable,public :: profile_matrix(:,:)
  integer(dp),dimension(ic_profile),public :: profile_last_timerstamp
  character(len=lc_profile),dimension(ic_profile),public :: profile_labels=" "

  public :: tic,tickrate
  public :: toc,tdiff,tocprint
  public :: now,reset_timer
  public :: timer_loop_start,in_timer_loop,timer_loop_end
  public :: timer_profile_setup,timer_profile_start,timer_profile_break

contains

!< Force timestamps and clock rate to zero
  subroutine reset_timer
    last_time_point=INT(0,dp)
    if (allocated(elapsed_times)) deallocate (elapsed_times)
  end subroutine reset_timer

!< Get system_clock timestamp
  integer(dp) function tic()
    call system_clock(count=tic)
    last_time_point=tic
  end function tic

!< Get the clock tick count per second
  integer(dp) function tickrate()
    tickrate=INT(0,dp)
    if (clock_ticks_per_second.eq.0) call system_clock(count_rate=tickrate)
    clock_ticks_per_second=tickrate
#ifdef _DEBUG
    write (*,'(A,1d20.10)') " CLOCK_RATE(ticks/sec) ",REAL(clock_ticks_per_second,dp)
#endif
  end function tickrate

!< Calculate the time from two timestamps
  real(dp) function tdiff(tfinal,tstart)
    integer(dp),intent(in) :: tfinal
    integer(dp),intent(in) :: tstart
! integer(dp)                       ::  end_point
! if(present(tstart)) last_time_point = tstart
! if(.not. present(tfinal)) then
!    call system_clock(count=end_point)
!    tfinal = end_point
! end if
    if (clock_ticks_per_second.eq.INT(0,dp)) call system_clock(count_rate=clock_ticks_per_second)
! Calulate the time difference
    tdiff=REAL(tfinal-tstart,dp)/REAL(clock_ticks_per_second,dp)
!      last_time_point = tfinal
  end function tdiff

!< Complete the timing regime using a reference timestamp or the one
!  in last_time_point
  real(dp) function toc(start_optional)
    integer(dp),intent(in),optional ::  start_optional
    integer(dp)                       ::  end_point
    if (present(start_optional)) last_time_point=start_optional
    call system_clock(count=end_point)
    toc=tdiff(end_point,last_time_point)
    last_time_point=end_point
  end function toc

!< Complete the timing regime using a reference timestamp or the one
!  in last_time_point
  real(dp) function tocprint(start_optional)
    integer(dp),intent(in),optional ::  start_optional
    integer(dp)                       ::  end_point
    if (present(start_optional)) last_time_point=start_optional
    call system_clock(count=end_point)
#ifdef _DEBUG
    write (*,'(A,1d20.10)') " TOC Time stamp ",REAL(end_point,dp)
#endif
    tocprint=tdiff(end_point,last_time_point)
    last_time_point=end_point
    write (*,'(A,1d20.10)') " Elapsed time ",tocprint
  end function tocprint

!> print current time and date
  subroutine now
    character(len=8)  :: date
    character(len=10) :: time

    print*,"System_clock: ",tic()
!***********************************************************************************
    call date_and_time(date,time)
    write (*,'(A,A,A,A,A,A,A)') 'Date: ',date(7:8),'-',date(5:6),'-',date(1:4),'\n'
    write (*,'(A,A,A,A,A,A,A)') 'Time: ',time(1:2),':',time(3:4),':',time(5:10),'\n'
  end subroutine now

!< in_timer_loop checks the time within a timer loop
! It does not start the timer or set the inloop variable
! It returns false on the final loop- so that timer_loop_end can finish
  logical function in_timer_loop()
    in_timer_loop=.false.
    if (.not.inloop) then
      print*,"Failed timer_loop: Timer loop did not start"
    else
      if (idx_elapsed.lt.num_elapsed) then
        elapsed_times(idx_elapsed)=toc()
        idx_elapsed=idx_elapsed+1
        in_timer_loop=.true.
      end if
    end if
  end function in_timer_loop

  subroutine timer_loop_start(num)
    integer,intent(in),optional :: num
    integer(dp)  :: dummytimestamp
    num_elapsed=3
    if (present(num)) num_elapsed=num
    call reset_timer()
    if (num_elapsed.gt.1) then
      allocate (elapsed_times(num_elapsed))
      idx_elapsed=1
      dummytimestamp=tic()
      inloop=.true.
    end if
#ifdef _DEBUG
    print*,'Size of elapsed array ',size(elapsed_times)
#endif
  end subroutine timer_loop_start

  subroutine timer_loop_end
  if (.not.inloop) then
    print*,"Failed timer_loop_end: Timer loop did not start"
  else
    write (*,'(A)') "******* TIMER LOOP **************"
    write (*,'(A,A,A,1i4)') '*** FILE:LINE: ',__FILE__,":",__LINE__
    write (*,'(A,1i8,A)') '*** Iterations: ',num_elapsed,' timed loops'

    if (idx_elapsed.eq.num_elapsed) then
      elapsed_times(idx_elapsed)=toc()
      write (*,'(A,1d20.10)') "*** Average (sec):",SUM(elapsed_times,DIM=1)/REAL(num_elapsed,dp)
      write (*,'(A,1d20.10,A,1i3)') "*** Longest run(sec) ",MAXVAL(elapsed_times,DIM=1), &
        '    at ',MAXLOC(elapsed_times,DIM=1)
      write (*,'(A,1d20.10,A,1i3)') "*** Shortest run(sec) ",MINVAL(elapsed_times,DIM=1), &
        '   at ',MINLOC(elapsed_times,DIM=1)
    else
      write (*,'(A,1i8)') '*** Failed at iteration ',idx_elapsed
    end if
    write (*,'(A)') "******* TIMER LOOP **************"
    inloop=.false.
    if (allocated(elapsed_times)) then
      deallocate (elapsed_times)
    end if
  end if
  end subroutine timer_loop_end

  subroutine timer_profile_setup(nLoops,v1,v2,v3,v4,v5)
    integer,intent(in) :: nLoops
    character(len=20),intent(in),optional :: v1
    character(len=20),intent(in),optional :: v2
    character(len=20),intent(in),optional :: v3
    character(len=20),intent(in),optional :: v4
    character(len=20),intent(in),optional :: v5
    !   integer :: nVars=0
    if (nLoops.lt.1) then
      print*,"timer_profile_setup error -- must have more than 1 loop"
      return
    end if

    if (.not.present(v1)) then
      print*,"timer_profile_setup error -- maximum 5 labels, minimum 1"
      return
    else
      !      allocate (character(len=20) :: profile_labels(5))
      if (present(v1)) then
        if (len_trim(v1).le.lc_profile) then
          profile_labels(1)=adjustl(v1)
          num_profile_vars=1
        end if
      end if
      if (present(v2)) then
        if (len_trim(v1).le.lc_profile) then
          profile_labels(2)=adjustl(v2)
          num_profile_vars=1
        end if
      end if
      if (present(v3)) then
        if (len_trim(v1).le.lc_profile) then
          profile_labels(3)=adjustl(v3)
          num_profile_vars=3
        end if
      end if
      if (present(v4)) then
        if (len_trim(v1).le.lc_profile) then
          profile_labels(4)=adjustl(v4)
          num_profile_vars=4
        end if
      end if
      if (present(v5)) then
        if (len_trim(v1).le.lc_profile) then
          profile_labels(5)=adjustl(v5)
          num_profile_vars=5
        end if
      end if
    end if
    num_profile_loops=nLoops
    if (allocated(profile_matrix)) deallocate (profile_matrix)
    !if (allocated(profile_last_timerstamp)) deallocate (profile_last_timerstamp)
    ! allocate (profile_last_timerstamp(num_profile_loops))
    allocate (profile_matrix(num_profile_loops,num_profile_vars))
    profile_matrix=REAL(0.0,dp)
    profile_last_timerstamp=INT(0,dp)

  end subroutine timer_profile_setup

!< Within loop - start timer
  subroutine timer_profile_start(LABEL)
    character(len=lc_profile),intent(inout) :: LABEL
    integer ::  ival
    do ival=1,num_profile_vars
      if (.not. (INDEX(trim(LABEL),profile_labels(ival))==0)) then
        profile_last_timerstamp(ival)=tic()
        exit
      end if
    end do
#ifdef _DEBUG
    if (ival.gt.num_profile_vars) then
      write (*,'(A,1i4,A,1i10)') "Timer_Profile_start:", &
        __LINE__," label index outside range ",ival
    else
      print*,"Label: ",profile_labels(ival)," time stamp "
    end if
#endif
  end subroutine timer_profile_start

!< Within loop - get elapsed time for label
  subroutine timer_profile_break(LABEL)
    character(len=lc_profile),intent(inout) :: LABEL
    integer(dp) :: tmp_tstamp=INT(0,dp)
    integer ::  ival=0,iloop=0
!
! Need bounds checking of timestamps and matrix
    do ival=1,num_profile_vars
      iloop=0
      if (.not. (INDEX(LABEL,profile_labels(ival))==0)) then
        do iloop=1,num_profile_loops
          if (profile_matrix(iloop,ival).ne.0.0) then
            tmp_tstamp=tic()
            profile_matrix(iloop,ival)= &
              tdiff(profile_last_timerstamp(ival),tmp_tstamp)
            profile_last_timerstamp(ival)=tmp_tstamp
            exit
          end if
        end do
        exit
      end if
    end do
    if (tmp_tstamp.eq.0) then
      write (*,'(A,A,A,A)') "Timer_Profile_break:",__LINE__," No time stamp created for ",label
    end if
    if (ival.gt.num_profile_vars) then
      write (*,'(A,1i4,A,1i8)') "Timer_Profile_break:",__LINE__," label index outside range ",ival
    end if
    if (iloop.gt.num_profile_loops) then
      write (*,'(A,1i4,A,1i8)') "Timer_Profile_break:",__LINE__," loop index outside range ",iloop
    end if
#ifdef _DEBUG
    if ((ival.gt.num_profile_vars).or.(iloop.gt.num_profile_loops)) then
      write (*,'(A,1i4,A,2i10)') "Timer_Profile_break:", &
        __LINE__," label/loop index outside range ",ival,iloop
    else
      print*,"Label: ",profile_labels(ival)," time ",profile_matrix(iloop,ival)
    end if
#endif

  end subroutine timer_profile_break

!< Profile report
  subroutine timer_profile_report(COMMENT)
    character(len=*),intent(inout),optional :: COMMENT
    integer :: ival,iloop
    if (.not.present(COMMENT)) COMMENT="PROFILE"
    write (*,'(A,A,A)') "*** LOOP REPORT : ",trim(COMMENT)
    write (*,'(A,A,A,1i4)') '*** FILE:LINE: ',__FILE__,":",__LINE__
    write (*,'(A,1i8,A)') '*** Iterations: ',num_profile_loops,' timed loops'

    do ival=1,num_profile_vars
      write (*,'(A,A)') '**** Label name: ',profile_labels(ival)
#ifdef _DEBUG
      do iloop=1,num_profile_loops
        if (profile_matrix(iloop,ival).eq.0) write (*,'(A,2i10)') "    Zero elapsed time at ",ival,iloop
      end do
#endif
      write (*,'(A,1d20.10)') "**** Average (sec):", &
        SUM(profile_matrix(:,ival),DIM=1)/REAL(num_profile_loops,dp)
      write (*,'(A,1d20.10,A,1i3)') "**** Longest run(sec) ", &
        MAXVAL(profile_matrix(:,ival),DIM=1),'    at ',MAXLOC(profile_matrix(:,ival),DIM=1)
      write (*,'(A,1d20.10,A,1i3)') "**** Shortest run(sec) ", &
        MINVAL(profile_matrix(:,ival),DIM=1),'   at ',MINLOC(profile_matrix(:,ival),DIM=1)
    end do
    write (*,'(A,A,A)') "******* END ",trim(COMMENT)," REPORT **************"
    deallocate (profile_matrix)
!    deallocate (profile_last_timerstamp)
!    deallocate (profile_labels)
    ! unset label
    profile_labels=" "
    profile_last_timerstamp=INT(0,dp)
  end subroutine timer_profile_report

end module simple_timer
