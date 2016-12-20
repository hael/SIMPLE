! 
!
!*******************************************************************************
!     Author: Frederic D.R. Bonnet date: 9th of October 2014.
!
! Name:
! simple_timing - Various utilities and timing routine for both CPU and GPU
!                 for other modules.
!
! Description:
! The routines finishing with an _c denotes a routine written as extern C code.
! these have been interfaced here via the interface module.
! The interfaced C code for the timming routine written in C are highly
! portable and very light which can be called from anywhere. These may be 
! called outside from any programs simultaneously and independently.
!
! timing_simple provides timing routine for the CPU wall clocks and GPU clocks
! to be used in other modules for benchmarking and peformance testing.
!
! Usage:
! To use the module ones just needs to call once the routine 
! start_Alltimers_cpu() at the strat of the main program then wrapp around the
! area of code that needs to be times with call start_timer_cpu("routine_name")
! and call stop_timer_cpu("routine_name"). It is important that the name are
! the same as the name array is sorted using shell sorted algorithm. Then at
! the end of the main program call stop_Alltimers_cpu(). The call start/stop
! routine can be spread out through the modules and routine everywhere in the
! code up to 50 timers. running simutaneously with a total cpu time in seconds
! and wall time in seconds. A cpu and Wall (% and sec) is reported for each
! timer
!
! Example:
!
! program myprog
! use simple_timing
! implicit none
!
! call start_Alltimers_cpu()
!
! call start_timer_cpu("timer #1")
!       ....do work....
! call stop_timer_cpu("timer #1")
!
!      
!      ....untimed work....
!
!
! call start_timer_cpu("timer #2")
!       ....do work....
! call stop_timer_cpu("timer #2")
!
!
! call stop_Alltimers_cpu()
! end program myprog
!*******************************************************************************
!
module simple_timing

  use simple_defs
  use simple_sorting
  use simple_error_handling

  implicit none

  !global variables
  logical                     :: enable_timing = .false.!global timer
  !local variables
  integer, private            :: ntimers = -1     !Number of named timers
  !The parameters
  integer, parameter, private :: LEN_TIMERS  = 32 !Length of the name of timer.
  integer, parameter, private :: NMAX_TIMERS = 100 !Maximum # of named timers.
  integer, parameter, private :: STACKSIZE_TIMERS = 16!Size of stack of timers.
  !error handlers
  integer, save, public :: TIMING_INVALID
  
#ifdef CUDA
  !TODO: make these global variables
  !cudaEvent_t start_event, stop_event; 
#endif

  !abstract type timer_type
  type timer_type
     double precision :: cpu_start, wall_start
     double precision :: cpu = 0.0d0, wall = 0.0d0
     double precision :: cpu_tot = 0.0d0, wall_tot = 0.0d0 
  end type timer_type

  !timers definition set a global and private variables 
  integer, private :: istacktimers = 0 !Pointer to the top of stack of timers.
  integer,dimension(STACKSIZE_TIMERS), private :: stacktimers !Stack of timers.
  integer,dimension(0:NMAX_TIMERS),private :: sorted_timers!Index timers_name
  character(len=LEN_TIMERS),dimension(0:NMAX_TIMERS), private :: timers_name
  type (timer_type),dimension(0:NMAX_TIMERS), private, save :: timers!max 51

  interface

     subroutine gettimeofday_c(Time)
       implicit none
       double precision,dimension(2) :: Time
     end subroutine gettimeofday_c

     subroutine elapsed_time_c(start_Time, end_Time, elapsed_time)
       implicit none
       double precision              :: elapsed_time
       double precision,dimension(2) :: start_Time
       double precision,dimension(2) :: end_Time
     end subroutine elapsed_time_c

     subroutine getcputime_c(cpu_time, elapsed_cpu_time)
       implicit none
       double precision,dimension(2) :: cpu_time
       double precision              :: elapsed_cpu_time
     end subroutine getcputime_c

  end interface
  !public methods
  public :: simple_timing_errors
  
contains
!*******************************************************************************
!Name: timers for routine profing initialisers
!
! DESCRIPTION
! subroutine to start/stop all the timers on CPU.
!
!*******************************************************************************
! SOURCE
subroutine simple_timing_errors()
  use simple_eglossary
  implicit none
  call file_err_define(err_name='TIMING_INVALID',&
       err_msg='Error in timing routines',&
       err_id=TIMING_INVALID,&
       err_action='Control the running conditions of f_timing routines called')
  return
end subroutine simple_timing_errors
!*******************************************************************************
!Name: simple_file_timing_initialise
!
! DESCRIPTION
! timers for routine profing initialisers on CPU only.
!
!*******************************************************************************
! SOURCE
subroutine simple_file_timing_initialise()
  implicit none
  !TODO: here needs to be implemented
  return
end subroutine simple_file_timing_initialise
!*******************************************************************************
!Name: reset_timers
!
! DESCRIPTION
! subroutine to start/stop all the timers on CPU.
!
!*******************************************************************************
! SOURCE
subroutine reset_timers_cpu(timer)
  implicit none

  !global variables
  type (timer_type), intent(inout) :: timer

  !start of the execution commands.
  timer%cpu = 0.0d0
  timer%wall = 0.0d0
  timer%cpu_tot = 0.0d0
  timer%wall_tot = 0.0d0

  return
end subroutine reset_timers_cpu
!*******************************************************************************
!     Author: Frederic D.R. Bonnet date: 08th of August 2013.
!
! Name:
! timer_timeStamp_cpu - Timestamper for code starting and finishing code for
! other modules.
!
! DESCRIPTION
! subroutine timestamp prints the current date and time.
!*******************************************************************************
! SOURCE
subroutine timer_timeStamp_cpu(values)
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

  call date_and_time(date, time, zone, values)

  year  = values(1)
  month = values(2)
  day   = values(3)
  UTC   = values(4)
  hour  = values(5)
  min   = values(6)
  sec   = values(7)
  msec  = values(8)

  !  write (*,'(1x,a,1x,a,1x,i2,1x,i4,2x,i2,a1,i2.2,a1,i2.2,a1,i3.3)') &
  !       'current date and time:',trim (all_month(month)), day, year, &
  !       hour,':', min, ':', sec, '.', msec

  return
end subroutine timer_timeStamp_cpu
!*******************************************************************************
! DESCRIPTION
! subroutine to start/stop all the timers on CPU.
!
!*******************************************************************************
! SOURCE
  subroutine start_Alltimers_cpu()
    implicit none

    !start of the execution commandds
    call hello_timers_cpu()

    enable_timing = .true.

    call start_timer_cpu("Other timers")

    return
  end subroutine start_Alltimers_cpu
  !stop all the timers.
  subroutine stop_Alltimers_cpu()
    implicit none

    !start of the execution commandds

    !stopping the timers
    call stop_timer_cpu("Other timers")

    !printing the timing results to screen
    call print_timer_cpu()

    !greeting off the timers users
    call bye_timers_cpu()

    return
  end subroutine stop_Alltimers_cpu
!*******************************************************************************
! DESCRIPTION
! subroutine to start the timer on CPU.
!
!*******************************************************************************
! SOURCE
  subroutine start_timer_cpu(name)
    implicit none

    !global variables
    character(len=*),intent(in) :: name
    !local variables
    integer                     :: ctimer
    double precision            :: cpu,wall
    !counters
    integer                     :: itimer

    !start of the execution commands
    !write(*,*)"start of timing on CPU the routine: ",name

    if (.not. enable_timing) return
    
    !get the wall and cpu time
    call get_wall_time_cpu(cpu,wall)

    !suspend the current timer on top of stack
    ctimer = -1
    if (istacktimers > 0) then
       ctimer = stacktimers(istacktimers)
       timers(ctimer)%cpu  = timers(ctimer)%cpu + cpu - timers(ctimer)%cpu_start
       timers(ctimer)%wall = timers(ctimer)%wall + wall - timers(ctimer)%wall_start
    end if

    !create the input timer on the top of the stack
    itimer = create_timer(name)

    !push input on the stack

    istacktimers = istacktimers+1
    if (istacktimers > STACKSIZE_TIMERS) then
       write (*, '(a)') "start_timer : Timers stack overflow."
       write (*, '(a)') "start_timer : Stack ="
       do itimer = STACKSIZE_TIMERS, 1, -1
          write (*, '(i3," : ",a)' ) itimer, timers_name(stacktimers(itimer))
       end do
       stop "start_timer : Error, timers stack overflow."
    end if
    stacktimers(istacktimers) = itimer
    
    !Start timers.

    timers(itimer)%cpu = 0.d0
    timers(itimer)%wall = 0.d0
    timers(itimer)%cpu_start  = cpu
    timers(itimer)%wall_start = wall

    return
  end subroutine start_timer_cpu

  !create the timer with the argument name

  integer function create_timer(name)
    implicit none
    !global variables
    character(len = *), intent(in) :: name
    !local variables
    integer                        :: increment
    !counters
    integer                        :: i,j,k

    !start of the execution commands

    !check if the input timer already exist
    create_timer = get_timer_index(name)
    if (create_timer > 0 ) return
    !if not create it
    if (ntimers >= NMAX_TIMERS ) then
       create_timer = 0
       return
    end if

    !increment the ntimers index
    ntimers = ntimers + 1
    timers_name(ntimers) = name(1:min(LEN_TIMERS,len(name)))
    
    call reset_timers_cpu(timers(ntimers))

    !Sort the timers by the shellsort method.
    !Worst case O(n^2), best case O(n*log(n)) to speed up access.

    call shellSort_cpu(sorted_timers, timers_name, ntimers)

    !return the index of the new timer
    create_timer = ntimers

    return
  end function create_timer

  integer function get_timer_index(index_name)
    implicit none
    !global variables
    character(len=*),intent(in)   :: index_name       
    !local variables 
    integer :: itimer, jtimer, jtimer1, jtimer2
    character(len = LEN_TIMERS) :: name

    !start of the execution commands
    get_timer_index = -1 

    if (ntimers < 0) return
    name = index_name(1:min(LEN_TIMERS, len(index_name)))

    jtimer1 = 0
    jtimer2 = ntimers

    if (timers_name(sorted_timers(jtimer1)) == name) then
       get_timer_index = sorted_timers(jtimer1)
       return
    end if

    if (timers_name(sorted_timers(jtimer2)) == name) then
       get_timer_index = sorted_timers(jtimer2)
       return
    end if

    do while (jtimer2 - jtimer1 > 1)

       jtimer = (jtimer1+jtimer2)/2
       itimer = sorted_timers(jtimer)

       if (timers_name(itimer) == name) then
          get_timer_index = itimer
          return
       end if

       if (timers_name(itimer)  > name) then
          jtimer2 = jtimer
       else
          jtimer1 = jtimer
       end if

    end do

    return
  end function get_timer_index
!*******************************************************************************
! DESCRIPTION
! subroutine to stop the timer on CPU.
!
!*******************************************************************************
! SOURCE
  subroutine stop_timer_cpu(name)
    implicit none
    !global variables
    character(len=*), intent(in) :: name
    !local variables
    integer                     :: ctimer
    double precision            :: cpu,wall
    !counters
    integer                     :: itimer

    !start of the execution commands
    !write(*,*)"stoping of timing on CPU the routine: ",name

    if (.not.enable_timing) return

    ! Check if the input timer exists.

    itimer = max(0, get_timer_index(name)) !Defaults timer #0 ("Others timers").

    !Check that the input timer matches the current timer(the one on the top 
    !stack).

    if (istacktimers <= 0) &
         &         stop "stop_timer : Error, the stack of timers is empty."

    if (stacktimers(istacktimers) /= itimer) then
       write (*, '("stop_timer : Error, the input timer (",a,") does not match the current timer (",a,").")') &
            &       trim(name), trim(timers_name(stacktimers(istacktimers)))
       write (*, '(a)') "stop_timer : Stack ="
       do itimer = istacktimers, 1, -1
          write (*, '(i3," : ",a)' ) itimer, timers_name(stacktimers(itimer))
       end do
       stop "stop_timer: Error, the input timer do not match"
    end if

    !Get the CPU and wall time.

    call get_wall_time_cpu(cpu, wall)

    !Stop the current timer

    timers(itimer)%cpu  = timers(itimer)%cpu  + cpu - timers(itimer)%cpu_start
    timers(itimer)%wall = timers(itimer)%wall + wall - timers(itimer)%wall_start
    timers(itimer)%cpu_tot  = timers(itimer)%cpu_tot  + timers(itimer)%cpu
    timers(itimer)%wall_tot = timers(itimer)%wall_tot + timers(itimer)%wall
    
    !pop it from the stack.

    istacktimers = istacktimers - 1

    ! Restart the new "current" timer (the one now on the top of the stack).

    if (istacktimers > 0) then
       ctimer = stacktimers(istacktimers)
       timers(ctimer)%cpu_start  = cpu
       timers(ctimer)%wall_start = wall
    end if

    return
  end subroutine stop_timer_cpu
!*******************************************************************************
! DESCRIPTION
! subroutine to calculate the cpu and wall time on CPU.
!
!*******************************************************************************
! SOURCE
  subroutine get_wall_time_cpu(cpu,wall)
    implicit none
    !global variables
    double precision              :: cpu,wall
    !local variables
    integer                       :: init
    integer                       :: year_init, month_init
    integer                       :: month, month_now
    integer,dimension(8)          :: values
    double precision,dimension(2) :: cpu_time
    real(kind(1.e0))              :: tt(2)
    ! Parameters.
    integer,parameter             :: &
         &      ndays(24)=(/31, 28, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31, &
         &                  31, 28, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31/)

    interface external_c_functions
       function get_cpu_time_c(tt)
         use simple_defs
         implicit none
         real(sp), intent(out) :: tt(2)
         real(dp) :: get_cpu_time_c
       end function get_cpu_time_c
    end interface external_c_functions

    !start of the execution commands

    !store the month and the year at first call
    if (init == 0 ) then
       init = 1
       call timer_timeStamp_cpu(values)
       year_init  = values(1)
       month_init = values(2)
    end if

    !now get the CPU time from the interfaced C code getcputime_c(cpu_time)
#if defined (MACOSX)
    cpu = get_cpu_time_c(tt)
#else
    call getcputime_c(cpu_time, cpu)
#endif

    !now get the wall time for the 
    ! First compute the number of seconds from the beginning of the month.

    wall = (values(3) * 24.d0 + values(5)) * 3600.d0 + &
&           values(6) * 60.d0 + values(7) + values(8) * 0.001d0

    !Compute the number of seconds to be added if the month has changed since
    !first call. This fails if the program ran one whole year!

    month_now = values(2)
    if (month_now /= month_init) then
      if (values(1) == year_init+1) then
        month_now = month_now+12
      end if
      do month = month_init, month_now-1
        wall = wall+86400.d0*ndays(month)
      end do
    end if

    ! Now take into account bissextile years.

    if ((mod(year_init, 4) == 0).and.(month_init <=  2).and.(month_now >  2)) &
&     wall = wall+3600.d0
    if ((mod(values(1), 4) == 0).and.(month_init <= 14).and.(month_now > 14)) &
&     wall = wall+3600.d0

    return
  end subroutine get_wall_time_cpu
!*******************************************************************************
! DESCRIPTION
! subroutine to print the timer on CPU.
!
!*******************************************************************************
! SOURCE
  subroutine print_timer_cpu()
    implicit none

    !local variables
    double precision               :: cputotal,walltotal
    character(len = 3)             :: fmtname
    character(len = 9)             :: hline
    !Counters
    integer                        :: itimer,jtimer
    !start of the execution commands

    !if ( ntimers < 0 ) return

    !compute the total CPU and wall time
    cputotal = 0.0d0
    walltotal = 0.0d0
    cputotal  = sum(timers(:)%cpu_tot)
    walltotal = sum(timers(:)%wall_tot)

    !printing the statistics to screen

    write (fmtname, '("a",i2)') LEN_TIMERS

    write (*, '(/'//fmtname//',4x,a7,1x,a6,3x,a9,1x,a6)') &
&     "", "cpu (s)", "%cpu", "wall (s)", "%wall"

    do jtimer = 0, ntimers
       itimer = sorted_timers(jtimer)
       if (itimer == 0) cycle
       write (*, '('//fmtname//',1x,f9.2,3x,f6.2,1x,f9.2,2x,f6.2)') &
            &  timers_name(itimer), &
            &  timers(itimer)%cpu_tot, 100.d0*timers(itimer)%cpu_tot/cputotal, &
            &  timers(itimer)%wall_tot, 100.d0*timers(itimer)%wall_tot/walltotal
    end do

    write (*, '('//fmtname//',1x,f9.2,3x,f6.2,1x,f9.2,2x,f6.2)') &
         &     timers_name(0), &
         &     timers(0)%cpu_tot,  100.d0*timers(0)%cpu_tot/cputotal, &
         &     timers(0)%wall_tot, 100.d0*timers(0)%wall_tot/walltotal

    write (hline, '("(",i2,"(""-""))")') LEN_TIMERS+42
    write (*, hline)

    write (*, '('//fmtname//',1x,f12.2,1x,6x,3x,f12.2)')"Total", &
         & cputotal, walltotal

    return
  end subroutine print_timer_cpu
!*******************************************************************************
! DESCRIPTION
!
!                    ********** GPU timers ***********
!
! subroutines to create start and stop events calculate the elapsed time
! and synchronize the events for coordinated timing. Routine to also destroy
! the event to clear the memory properly for the timer on GPU.
!
!*******************************************************************************
! SOURCE
  subroutine cuda_EventCreate_start_gpu()!TODO: insert the event creator
    implicit none
    
    !start of the execution commands

#ifdef CUDA
    !TODO:
    !    cudaEventCreate(&start_event);
    !    cudaEventCreate(&stop_event);
#else
    write(*,*)"!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!"
    write(*,*)"!!!You need to compile with the -CUDA  !!!"
    write(*,*)"                                          "
    write(*,*)"!!! normal timer will be called instead   "
    write(*,*)"!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!"
#endif

    return
  end subroutine cuda_EventCreate_start_gpu
!*******************************************************************************
  subroutine cuda_EventDestroy_gpu()!TODO: insert the event destroyer
    implicit none

    !start of the execution commands

#ifdef CUDA
    !TODO:
    !  cuEventDestroy(start_event);
    !  cuEventDestroy(stop_event);
#else
    write(*,*)"!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!"
    write(*,*)"!!!You need to compile with the -CUDA  !!!"
    write(*,*)"                                          "
    write(*,*)"!!! normal timer will be called instead   "
    write(*,*)"!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!"
#endif

    return
  end subroutine cuda_EventDestroy_gpu
!*******************************************************************************
  subroutine cuda_EventRecord_start_timer_gpu() !TODO: (start) event recorder
    implicit none

#ifdef CUDA
    !TODO: insert the
    !
    !    cudaEventRecord(start_event,0); //start of the cudaEventRecord
    !
    !    cudaEventRecord(stop_event,0);      //stop of the cudaEventRecord
    !    cudaEventSynchronize(stop_event);   //wait until the event is finished.
    !    cudaEventElapsedTime(&elapse, start_event, stop_event);
    !
    !    //making sure that the GPU is synchronized with CPU before proceeding
    !    cuCtxSynchronize();
#else
    write(*,*)"!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!"
    write(*,*)"!!!You need to compile with the -CUDA  !!!"
    write(*,*)"                                          "
    write(*,*)"!!! normal timer will be called instead   "
    write(*,*)"!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!"
#endif

    return
  end subroutine cuda_EventRecord_start_timer_gpu
!*******************************************************************************
  subroutine cuda_EventRecord_stop_timer_gpu()!TODO: (stop) event recorder
    implicit none

#ifdef CUDA
    !TODO: insert the 
    !
    !    cudaEventRecord(start_event,0); //start of the cudaEventRecord
    !
    !    cudaEventRecord(stop_event,0);      //stop of the cudaEventRecord
    !    cudaEventSynchronize(stop_event);   //wait until the event is finished.
    !    cudaEventElapsedTime(&elapse, start_event, stop_event);
    !
    !    //making sure that the GPU is synchronized with CPU before proceeding
    !    cuCtxSynchronize();
#else
    write(*,*)"!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!"
    write(*,*)"!!!You need to compile with the -CUDA  !!!"
    write(*,*)"                                          "
    write(*,*)"!!! normal timer will be called instead   "
    write(*,*)"!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!"
#endif

    return
  end subroutine cuda_EventRecord_stop_timer_gpu
!*******************************************************************************
  subroutine cuda_EventSynchronize_timer_gpu()!TODO: event Synchronizer
    implicit none

#ifdef CUDA
    !TODO: insert the 
    !
    !    cudaEventRecord(start_event,0); //start of the cudaEventRecord
    !
    !    cudaEventRecord(stop_event,0);      //stop of the cudaEventRecord
    !    cudaEventSynchronize(stop_event);   //wait until the event is finished.
    !    cudaEventElapsedTime(&elapse, start_event, stop_event);
    !
    !    //making sure that the GPU is synchronized with CPU before proceeding
    !    cuCtxSynchronize();
#else
    write(*,*)"!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!"
    write(*,*)"!!!You need to compile with the -CUDA  !!!"
    write(*,*)"                                          "
    write(*,*)"!!! normal timer will be called instead   "
    write(*,*)"!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!"
#endif

    return
  end subroutine cuda_EventSynchronize_timer_gpu
!*******************************************************************************
  subroutine cuda_EventElapsedTime_timer_gpu()!TODO: event ElapsedTimer
    implicit none

#ifdef CUDA
    !TODO: insert the 
    !
    !    cudaEventRecord(start_event,0); //start of the cudaEventRecord
    !
    !    cudaEventRecord(stop_event,0);      //stop of the cudaEventRecord
    !    cudaEventSynchronize(stop_event);   //wait until the event is finished.
    !    cudaEventElapsedTime(&elapse, start_event, stop_event);
    !
    !    //making sure that the GPU is synchronized with CPU before proceeding
    !    cuCtxSynchronize();
#else
    write(*,*)"!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!"
    write(*,*)"!!!You need to compile with the -CUDA  !!!"
    write(*,*)"                                          "
    write(*,*)"!!! normal timer will be called instead   "
    write(*,*)"!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!"
#endif

    return
  end subroutine cuda_EventElapsedTime_timer_gpu
!*******************************************************************************
! DESCRIPTION
! subroutine to greet and good bye the timers on CPU.
!
!*******************************************************************************
! SOURCE

  !Helloer of the timers
  subroutine hello_timers_cpu()
    implicit none
    write(*,*)"Hello timers on CPU"
    return
  end subroutine hello_timers_cpu
  !goodbyer of the timers
  subroutine bye_timers_cpu()
    implicit none
    write(*,*)"Good bye timers on CPU"
    return
  end subroutine bye_timers_cpu

end module simple_timing
