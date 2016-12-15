! 
!*******************************************************************************
!     Author: Frederic D.R. Bonnet date: 28th of June 2016.
!
! Name:
! simple_dynamic_memory - Module that manages the dynamic memory allocation
!
!*******************************************************************************
!
module simple_dynamic_memory
  use simple_memory_profiling
  use simple_file_utils
  use simple_eglossary
  use simple_yaml_strings
  use simple_error_handling
  use simple_err_defs
  use simple_module_file_malloc
  implicit none

  private 
  logical, parameter :: track_origins=.true.!< keeps track of all allocation
  integer, parameter :: max_ctrl = 5 !Maximum number of nested levels
  integer, parameter :: namelen=file_malloc_namelen !Character length (simple_file_malloc module var)
  integer :: ictrl=0                 !Id of active control structure (<=max_ct

  character(len=*), parameter :: processid='Process Id'
  character(len=*), parameter :: main='Main_program'
  character(len=*), parameter :: subprograms='Subroutines'
  character(len=*), parameter :: no_of_calls='No. of calls'
  character(len=*), parameter :: prof_enabled='Profiling Enabled'
  character(len=*), parameter :: t0_time='Time of last opening'
  character(len=*), parameter :: tot_time='Total time (s)'

  !to be initialized in the dynamic_memory module
  integer, save :: ERR_INVALID_MALLOC
  integer, save :: ERR_MALLOC_INTERNAL

  public :: file_malloc_initialise
  public :: file_malloc_set_status
  !error handlers
  public :: dynamic_memory_errors

  !control structure of file_lib library
  !contains all the global variables needed
  type :: mem_ctrl
     logical :: profile_initialized  !global variables for initialization
     logical :: routine_opened       !global variable
     logical :: profile_routine      !profiling decider
     character(len=namelen) :: present_routine !name of the active routine 
     character(len=256) :: logfile !output log file
     integer :: logfile_unit !logfile unit for the stream
     integer :: output_level !output level for login
     !eglossaries needed for profiling storage
     type(eglossary), pointer :: egloss_global!Mem status at higher level
     type(eglossary), pointer :: egloss_routine!Mem status inside the routine
     type(eglossary), pointer :: egloss_calling_sequence !routine profiling
     type(eglossary), pointer :: egloss_codepoint!points to location in the previous eglossary
  end type mem_ctrl

  !if the library routines are called without initialization
  type(mem_ctrl), dimension(0:max_ctrl), save :: mems
  !memeory profiler 
  type(memory_state), save :: memstate 

contains

  !Error messages associated to the dynamic memory
  subroutine dynamic_memory_errors()
    implicit none

    call file_err_define(err_name='ERR_ALLOCATE',&
         err_msg='Allocation error',err_id=ERR_ALLOCATE,&
         err_action='Control the order of the allocation or if the memory limit has been reached')
    call file_err_define(err_name='ERR_DEALLOCATE',&
         err_msg='Deallocation error',err_id=ERR_DEALLOCATE,&
         err_action='Control the order of the allocation or if the memory limit has been reached')
    call file_err_define(err_name='ERR_MEMLIMIT',&
         err_msg='Memory limit reached',err_id=ERR_MEMLIMIT,&
         err_action='Control the size of the arrays needed for this run with bigdft-tool program')
    call file_err_define(err_name='ERR_INVALID_COPY',&
         err_msg='Copy not allowed',&
         err_id=ERR_INVALID_COPY,&
         err_action=&
         'A f_memcpy command failed, probably invalid sizes: check sizes of arrays at runtime')
    call file_err_define(err_name='ERR_INVALID_MALLOC',&
         err_msg='Invalid specification of f_malloc',&
         err_id=ERR_INVALID_MALLOC,&
         err_action='Put coherent data for the memory space allocation')
    call file_err_define(err_name='ERR_MALLOC_INTERNAL',&
         err_msg='Internal error of memory profiler',&
         err_id=ERR_MALLOC_INTERNAL,&
         err_action='An invalid operation occurs, submit bug report to developers')
    return
  end subroutine dynamic_memory_errors

  !Nullifies the pointers and initialises the ctrl mem
  pure subroutine nullify_mem_ctrl(mem)
    implicit none
    type(mem_ctrl),intent(out) :: mem
    mem%profile_initialized=.false.
    mem%routine_opened=.false.
    mem%profile_routine=.true.
    mem%present_routine=repeat(' ',namelen)
    mem%logfile=repeat(' ',len(mem%logfile))
    mem%logfile_unit=-1 !not initiliased
    mem%output_level=0
    nullify(mem%egloss_global)
    nullify(mem%egloss_routine)
    nullify(mem%egloss_calling_sequence)
    nullify(mem%egloss_codepoint)
    return
  end subroutine nullify_mem_ctrl
  !initialises the mem ctrl
  subroutine initialise_mem_ctrl(mem)
    implicit none
    type(mem_ctrl) :: mem
    mem%profile_initialized=.true.
    !initalize the eglossary with the allocation information
    nullify(mem%egloss_routine)
    call egloss_init(mem%egloss_global)
    call egloss_init(mem%egloss_calling_sequence)
    !in principle the calling sequence starts from the main
    mem%egloss_codepoint => mem%egloss_calling_sequence
    call set_routine_info(mem%present_routine,mem%profile_routine)
    return
  end subroutine initialise_mem_ctrl

  !> Transfer to the f_malloc_module the information of the routine
  subroutine set_routine_info(name,profile)
    !use simple_module_file_malloc
    implicit none
    logical, intent(in) :: profile
    character(len=*), intent(in) :: name
    !function calls
    integer :: strlen
    !TODO: need to fix the undefined reference to file_malloc_routine_name
!    file_malloc_routine_name(1:len(file_malloc_routine_name))=name
!    file_malloc_default_profiling=profile
    return
  end subroutine set_routine_info

  !initilaises the ctrl memory
  function mem_ctrl_init() result(mem)
    type(mem_ctrl) :: mem
    call nullify_mem_ctrl(mem)
    call initialise_mem_ctrl(mem)
  end function mem_ctrl_init
    
  !Set to 0 counters
  pure subroutine memstate_init(memstate)
    implicit none
    type(memory_state), intent(out) :: memstate

    memstate%memtot%memory=int(0,kind=8)
    memstate%memtot%peak=int(0,kind=8)
    memstate%memtot%routine=''
    memstate%memtot%array=''

    memstate%memalloc=0
    memstate%memdealloc=0

    memstate%memloc%routine='routine'
    memstate%memloc%array='array'
    memstate%memloc%memory=int(0,kind=8) !fake initialisation to print the first routine
    memstate%memloc%peak=int(0,kind=8)
  end subroutine memstate_init

  !subroutine to initioalise the file memory
  subroutine file_malloc_initialise()
    implicit none
    !local variables
    integer :: rc
    character(len=2) :: char_out
    !Functions calls
    integer :: convert_int2char_pos_c
    ictrl = ictrl + 1
    rc = convert_int2char_pos_c(char_out,abs(max_ctrl))
    if (file_err(ictrl > max_ctrl,&
         'The number of active instances cannot exceed'//trim(char_out),&
         ERR_MALLOC_INTERNAL)) return
    mems(ictrl)=mem_ctrl_init()

    if (ictrl==1) call memstate_init(memstate)

    write(*,*) "memstate: ",memstate
    !initialize the memprofiling counters
    !TODO: need to pu the time stamp of the started profile using the timestamp
    call set(mems(ictrl)%egloss_global//'Timestamp of Profile initialization',&
         trim(yaml_date_and_time_toa()))
    
    !Process Id (used to dump)
    call set(mems(ictrl)%egloss_global//processid,0)

    !TODO: need a profiling starting method for the main program
    call file_routine(id=main)
    
    call file_malloc_set_status(0.e0)

    return
  end subroutine file_malloc_initialise

  !This routine adds the corresponding subroutine name to the eglossary
  subroutine file_routine(id,profile)
    use simple_yaml_output
    implicit none
    logical, intent(in), optional :: profile
    character(len=*), intent(in), optional :: id
    !local variables
    integer(kind=8) :: itime
    integer :: lgt,ncalls
    !fiunction calls
    integer(kind=8) :: getcpunanotime_c
    
    if (file_err(ictrl == 0,&
      'ERROR (f_routine): the routine f_malloc_initialize has not been called',&
      ERR_MALLOC_INTERNAL)) return
    
    if (.not. present(id)) return !no effect

    itime = getcpunanotime_c(itime)

    write(*,*) "in file_routine, itime: ",itime

    !TODO: need a proper p[rofiling methods to integrate into here from
    !      simple_timing.f90 module. This module will with file_timer_interrupt

    !desactivate profile_routine if the mother routine has desactivated it
    if (present(profile)) mems(ictrl)%profile_routine= &
                          mems(ictrl)%profile_routine .and. profile

    if (track_origins) then
       if(associated(mems(ictrl)%egloss_routine)) then
          !appending the routine to the global eglossary
          call egloss_appender(mems(ictrl)%egloss_global, &
                               mems(ictrl)%egloss_routine)
          nullify(mems(ictrl)%egloss_routine)
       end if
       !previous routine has not been closed yet
       if (mems(ictrl)%routine_opened) then
          mems(ictrl)%egloss_codepoint=>mems(ictrl)%egloss_codepoint//subprograms
       end if
       mems(ictrl)%routine_opened=.true.
       !see if the key existed in the codepoint
       if (has_key(mems(ictrl)%egloss_codepoint,trim(id))) then
          !retrieve number of calls and increase it
          ncalls=mems(ictrl)%egloss_codepoint//trim(id)//no_of_calls
          call set(mems(ictrl)%egloss_codepoint//trim(id)//no_of_calls,ncalls+1)
          !write the starting point for the time
          call set(mems(ictrl)%egloss_codepoint//trim(id)//t0_time,itime)
          call set(mems(ictrl)%egloss_codepoint//trim(id)//prof_enabled,mems(ictrl)%profile_routine)
       else
          !create a new eglossary
          call set(mems(ictrl)%egloss_codepoint//trim(id),&
               egloss_new((/no_of_calls .is. yaml_toa(1), &
                                t0_time .is. yaml_toa(itime),&
                               tot_time .is. yaml_toa(0.d0,fmt='(f4.1)'), &
                           prof_enabled .is. yaml_toa(mems(ictrl)%profile_routine)/)))
       end if
       !then fix the new codepoint from this one
       mems(ictrl)%egloss_codepoint=>mems(ictrl)%egloss_codepoint//trim(id)

       lgt=min(len(id),namelen)
       mems(ictrl)%present_routine=repeat(' ',namelen)
       mems(ictrl)%present_routine(1:lgt)=id(1:lgt)
    end if !end of track_origins if block
    
    !TODO: need to finish the method implementation

    !need to set the routine info
    !need to resume the timer if it has been interrupted
    
    return

  end subroutine file_routine
  
  !subroutine to set the status
  subroutine file_malloc_set_status(memory_limit,output_level,&
       logfile_name,iproc)
    !use simple_file_utils
    use simple_error_handling
    !use simple_yaml_strings
    !use simple_eglossary
    use simple_yaml_output
    implicit none
    real(kind=4), intent(in), optional :: memory_limit     !Memory limit
    integer, intent(in), optional :: output_level !Level of output for memocc
                                                  !0 no file, 1 light, 2 full
    character(len=*), intent(in), optional :: logfile_name !Name of the logfile
    integer, intent(in), optional :: iproc !Proc Id (used to dump, by default 0)
    !local variables
    integer :: unit_in
    integer :: jproc
    !indexers
    integer :: jctrl
    !TODO(file_malloc_set_status): need to finish implementation of the method  
    write(*,*) "in file_malloc_set_status: "
    write(*,*) "memory_limit: ",memory_limit
    write(*,*) "output_level: ",output_level
    write(*,*) "logfile_name: ",logfile_name
    write(*,*) "iproc: ",iproc

    if (file_err(ictrl == 0,&
         'ERROR (file_malloc_set_status): the routine file_malloc_initialize has not been called',&
         ERR_MALLOC_INTERNAL)) return

    if (present(output_level)) then
       if (output_level>0) then
          jproc=0
          !TODO: need to check if we already know which proc we are
          !jproc=mems(ictrl)%egloss_global .get. processid
          if (present(iproc)) jproc=iproc
          !if the logfile_name is not there throw an error
          if ( .not. present(logfile_name)) then
             call file_err_throw('Error, file_malloc_set_status needs logfile_name for nontrivial output level',&
                  err_id=ERR_INVALID_MALLOC)
             return
          end if
          !TODO: need to close the previous stream
          write(*,*)"mems(ictrl)%logfile_unit: ",mems(ictrl)%logfile_unit
          if (mems(ictrl)%logfile_unit > 0 .and. jproc==0) then
             !TODO: need to implement the close_stream method
             !call yaml_close_stream(unit=mems(ictrl)%logfile_unit)
          end if
          do jctrl=ictrl-1,1,-1
             if (trim(logfile_name)==mems(jctrl)%logfile) &
                  call file_err_throw('Logfile name "'//trim(logfile_name)//&
                  '" in f_malloc_set_status invalid, already in use for instance No.'//jctrl&
                  ,err_id=ERR_INVALID_MALLOC)
             exit
          end do
          unit_in=-1 !SM: unt otherwise not defined for jproc/=0
          if (jproc==0) then
             !check if the file is opened
             call file_unit(trim(logfile_name),unit_in)
             !after this check an opened filename may now be closed
             call file_close(unit_in)
             !we can now open the file
             !get a free unit and start from there
             unit_in=file_get_free_unit(98)
             write(*,*) unit_in      
             call yaml_set_stream(unit=unit_in,filename=trim(logfile_name), &
                  position='rewind',setdefault=.false., &
                  record_length=131)
             
             if (output_level==2) then
                call yaml_comment('Present Array','yes',unit_in)
                call yaml_sequence_open('List of allocations',unit=unit_in)
             end if
          end if
          write(*,*) " file_malloc_set_status line 269 unit_in: ",unit_in
          !store the found unit in the structure
          mems(ictrl)%logfile_unit=unit_in
          call file_strcpy(dest=mems(ictrl)%logfile,src=logfile_name)
       end if
       mems(ictrl)%output_level=output_level
    end if
       
    !TODO: set memory limit
    if (present(memory_limit)) call file_set_memory_limit(memory_limit)
    
    return
  end subroutine file_malloc_set_status

  
end module simple_dynamic_memory
