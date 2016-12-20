!> \brief SIMPLE systemQuery class for GPU
module simple_systemQuery_cpu
use, intrinsic :: iso_c_binding
use simple_defs

implicit none

public :: systemQuery_cpu, hello_systemQuery_cpu, bye_systemQuery_cpu
!public :: get_warning_dataStruct

type :: systemQuery_cpu

   type(systemDetails) :: hstD
   integer             :: n_CPU_cores
#if defined (MACOSX)
   integer*8           :: h_memSize
   integer             :: h_memUser
#elif defined (LINUX)   
   integer*8           :: h_memSize
   integer*8           :: avail_Mem
#endif
   logical             :: existence_hstQ_cpu=.false. !< objects exist or not

contains
  !constructor
  procedure :: new_systemQuery_cpu
  !destructor
  procedure :: kill_systemQuery_cpu
  !Terminators
  procedure :: terminate
  !setters
  procedure :: set_ncpu_cores
  procedure :: set_memory_size_host
  !getters
  procedure :: get_ncpu_cores
  procedure :: get_SysTotalMem_size
  procedure :: get_SysAvailMem_size
  !warners
  procedure :: get_warning_dataStruct
end type systemQuery_cpu

interface systemQuery_cpu
   module procedure constructor_systemQuery_cpu
end interface systemQuery_cpu

interface
#if defined (LINUX)
   function get_CPU_cores(nCPU_cores, hstD)
     import :: systemDetails
     type(systemDetails) :: hstD
     integer :: nCPU_cores
     integer :: get_CPU_cores
   end function get_CPU_cores

   function get_memorySize_host(memSize, hstD)
     import :: systemDetails
     type(systemDetails) :: hstD
     integer*8 :: memSize
     integer :: get_memorySize_host
   end function get_memorySize_host

   function get_available_memory_host(availMem, hstD)
     import :: systemDetails
     type(systemDetails) :: hstD
     integer*8 :: availMem
     integer :: get_available_memory_host
   end function get_available_memory_host
#endif
end interface

contains 

  !CONSTRUCTORS

  !> \brief is a systemQuery constructor
  function constructor_systemQuery_cpu(hstD) result(hstQ_cpu)
    type(systemQuery_cpu) :: hstQ_cpu
    type(systemDetails) :: hstD
    call hstQ_cpu%new_systemQuery_cpu(hstD)
  end function constructor_systemQuery_cpu

  subroutine new_systemQuery_cpu(hstQ_cpu,hstD)
    class(systemQuery_cpu) :: hstQ_cpu
    type(systemDetails) :: hstD
    !kill the existing object before allocating a new one
    call hstQ_cpu%kill_systemQuery_cpu

    call gather_systemDetails(hstQ_cpu,hstD)

    !set the existence of the object to true
    hstQ_cpu%existence_hstQ_cpu = .true.

    return
  end subroutine new_systemQuery_cpu

  !WORKERS METHODS

  subroutine gather_systemDetails(hstQ_cpu,hstD)
    class(systemQuery_cpu), intent(inout) :: hstQ_cpu
    type(systemDetails) :: hstD

    call set_ncpu_cores(hstQ_cpu,hstD)
    call set_memory_size_host(hstQ_cpu,hstD)

    return
  end subroutine gather_systemDetails

  !SETTERS

  !settiong the memory size of the system.
  subroutine set_memory_size_host(hstQ_cpu,hstD)
    class(systemQuery_cpu), intent(inout) :: hstQ_cpu
    type(systemDetails) :: hstD
#if defined (MACOSX)
    integer*8 :: memSize
    integer ::  userMem
#elif defined (LINUX)
    integer*8 :: memSize
    integer*8 :: availMem
#endif
    integer :: rc !return code
    interface external_c_function_get_memorySize_host
       function get_memorySize_host_c(memSize, hstD)
         import :: systemDetails
         type(systemDetails) :: hstD
         integer*8 :: memSize
         integer :: get_memorySize_host_c
       end function get_memorySize_host_c
    end interface external_c_function_get_memorySize_host
    interface external_c_function_get_userMem_host
       function get_userMem_host_c(userMem, hstD)
         import :: systemDetails
         type(systemDetails) :: hstD
         integer :: userMem
         integer :: get_userMem_host_c
       end function get_userMem_host_c
    end interface external_c_function_get_userMem_host
    
#if defined (MACOSX)
    rc = get_memorySize_host_c(memSize, hstD)
    hstQ_cpu%h_memSize = hstD%mem_Size
    rc = get_userMem_host_c(userMem, hstD)
    hstQ_cpu%h_memUser = hstD%mem_User
    write(*,*)"memSize= ",memSize,userMem
#elif defined (LINUX)
    rc = get_memorySize_host(memSize, hstD)
    hstQ_cpu%h_memSize = hstD%mem_Size
    rc = get_available_memory_host(availMem, hstD)
    hstQ_cpu%avail_Mem = hstD%avail_Mem
#endif
    return
  end subroutine set_memory_size_host

  !setting the number of cores on the system
  subroutine set_ncpu_cores(hstQ_cpu,hstD)
    class(systemQuery_cpu), intent(inout) :: hstQ_cpu
    type(systemDetails) :: hstD
    integer :: nCPU_cores
    integer :: rc !return code
    interface external_c_function_get_CPU_cores
       function get_CPU_cores_c(nCPU_cores, hstD)
         import :: systemDetails
         type(systemDetails) :: hstD
         integer :: nCPU_cores
         integer :: get_CPU_cores_c
       end function get_CPU_cores_c
    end interface external_c_function_get_CPU_cores
    
#if defined (MACOSX)
    rc = get_CPU_cores_c(nCPU_cores, hstD)
    hstQ_cpu%n_CPU_cores = hstD%nCPUcores
    write(*,*)"nCPU_cores =",nCPU_cores
#elif defined (LINUX)
    rc = get_CPU_cores(nCPU_cores, hstD)
    hstQ_cpu%n_CPU_cores = hstD%nCPUcores
#endif

    return
  end subroutine set_ncpu_cores

  !GETTERS

  function n_CPU_cores(hstQ_cpu) result(nCPU_cores_out)
    class(systemQuery_cpu), intent(in) :: hstQ_cpu
    integer :: nCPU_cores_out
    nCPU_cores_out = hstQ_cpu%n_CPU_cores
  end function n_CPU_cores

  function get_ncpu_cores(hstQ_cpu) result(n_CPU_cores_out)
    class(systemQuery_cpu), intent(in) :: hstQ_cpu
    integer :: n_CPU_cores_out
    n_CPU_cores_out = hstQ_cpu%n_CPU_cores
  end function get_ncpu_cores
  
  !getter for the total memory size on the system
  function get_SysTotalMem_size(hstQ_cpu) result(h_TotalMemSize_out)
    class(systemQuery_cpu), intent(inout) :: hstQ_cpu
    integer*8           :: h_TotalMemSize_out
#if defined (MACOSX)
    h_TotalMemSize_out = hstQ_cpu%h_memSize 
#elif defined (LINUX)
    h_TotalMemSize_out = hstQ_cpu%h_memSize
#endif
  end function get_SysTotalMem_size
  !getter for the available memory size on the system
  function get_SysAvailMem_size(hstQ_cpu) result(h_AvailMemSize_out)
    class(systemQuery_cpu), intent(inout) :: hstQ_cpu
    integer*8           :: h_AvailMemSize_out
#if defined (MACOSX)
    h_AvailMemSize_out = hstQ_cpu%h_memUser
#elif defined (LINUX)
    h_AvailMemSize_out = hstQ_cpu%avail_Mem
#endif
  end function get_SysAvailMem_size
  
  !GREETERS

  !> \brief hello greeting routine for the object
  subroutine hello_systemQuery_cpu(err)
    implicit none
    integer :: err
    !start of the execution commands
    write(*,*) "Hello systemQuery CPU world"
    write(*,*)
    return
  end subroutine hello_systemQuery_cpu

  !> \brief bye greeting routine for the object
  subroutine bye_systemQuery_cpu()
    implicit none
    !start of the execution commands
    write(*,*) "Bye systemQuery CPU world"
    write(*,*)
    return
  end subroutine bye_systemQuery_cpu

  !TERMINATORS because of incompatiblilties errors
  subroutine terminate(hstQ_cpu)
    class(systemQuery_cpu) :: hstQ_cpu
    call hstQ_cpu%kill_systemQuery_cpu

    call bye_systemQuery_cpu()
    stop
    return
  end subroutine terminate

  !WARNERS

  subroutine get_warning_dataStruct(hstQ_cpu)
    !implicit none
    class(systemQuery_cpu) :: hstQ_cpu
    write(*,*)
    write(*,*)"***********data structure and getter do not match***************"
    write(*,*)
    return
  end subroutine get_warning_dataStruct
  
  !DESTRUCTOR

  !> \brief is a systemQuery_cpu destructor
  subroutine kill_systemQuery_cpu(hstQ_cpu)
    class(systemQuery_cpu), intent(inout) :: hstQ_cpu
    if ( hstQ_cpu%existence_hstQ_cpu) then

       hstQ_cpu%existence_hstQ_cpu = .false.

    end if
    return
  end subroutine kill_systemQuery_cpu

end module simple_systemQuery_cpu
