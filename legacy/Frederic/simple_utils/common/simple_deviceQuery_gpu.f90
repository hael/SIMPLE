!> \brief SIMPLE deviceQuery class for GPU
module simple_deviceQuery_gpu
use, intrinsic :: iso_c_binding
use simple_defs
use simple_cuda_defs

implicit none

public :: deviceQuery_gpu, hello_deviceQuery_gpu, bye_deviceQuery_gpu
public :: get_devD

type :: deviceQuery_gpu

   type(deviceDetails)  :: devD
   integer              :: devCnt
   character,pointer    :: devname(:)
   real(sp)             :: devVer
   real(sp)             :: runVer
   real(sp),allocatable :: totGlobalMem_MB(:)
   integer*8,allocatable:: totGlobalMem_bytes(:)
   integer,allocatable  :: n_MultiProc(:)
   integer,allocatable  :: n_cudacor_per_MultiProc(:)
   integer,allocatable  :: n_cudacor(:)
   integer,allocatable  :: SM_suitable(:)
   integer,allocatable  :: n_registers_per_blk(:)
   integer,allocatable  :: n_warp(:)
   integer,allocatable  :: n_threads_per_mp(:)
   integer,allocatable  :: n_threads_per_blk(:)
   integer,allocatable  :: dev_is_ecc(:)
   integer,allocatable  :: dev_is_p2p(:,:)
   logical              :: existence_devQ_gpu=.false. !< objects exist or not

contains
  !constructor
  procedure :: new_deviceQuery_gpu
  !destructor
  procedure :: kill_deviceQuery_gpu
  !Terminators
  procedure :: terminate
  !setters
  procedure :: set_devCnt
  procedure :: set_p2p
  procedure :: set_ecc_support
  procedure :: set_thread_details
  procedure :: set_nregisters
  procedure :: set_cuda_cores
  procedure :: set_tot_global_mem_MB
  procedure :: set_tot_global_mem_bytes
  procedure :: set_cuda_rnt_version
  procedure :: set_cuda_drv_version
  procedure :: set_devName
  !getters
  procedure :: get_devCnt
  procedure :: get_devname
  procedure :: get_cuda_drv_version
  procedure :: get_cuda_rnt_version
  procedure :: get_totGlobalMem_MB
  procedure :: get_totGlobalMem_bytes
  procedure :: get_dev_cuda_mp
  procedure :: get_dev_cuda_core_per_mp
  procedure :: get_dev_cuda_core
  procedure :: get_dev_SM_suitable
  procedure :: get_dev_nregisters_per_blk
  procedure :: get_dev_nwarp
  procedure :: get_dev_n_threads_per_mp
  procedure :: get_dev_n_threads_per_blk
  procedure :: get_eecsupport
  procedure :: get_p2p
  procedure :: get_devD
  !warners
  procedure :: get_warning_dataStruct_gpu
  !initialisers
  procedure :: initialiase_deviceDetails
end type deviceQuery_gpu

interface deviceQuery_gpu
   module procedure constructor_deviceQuery_gpu
end interface deviceQuery_gpu

interface
#if defined (LINUX)
   function get_dev_count(devD)
     import :: deviceDetails
     type(deviceDetails) :: devD
     integer :: get_dev_count
   end function get_dev_count

   function deviceDetails_init(idev,devD)
     import :: deviceDetails
     type(deviceDetails) :: devD
     integer :: idev
     integer :: deviceDetails_init
   end function deviceDetails_init

   function get_dev_Name(idev, devD)
     import :: deviceDetails
     type(deviceDetails) :: devD
     integer :: idev
     integer :: get_dev_Name
   end function get_dev_Name

   function get_cuda_driverVersion(driverVersion, devD)
     import :: deviceDetails
     type(deviceDetails) :: devD
     integer :: idev
     integer :: driverVersion
     integer :: get_cuda_driverVersion
   end function get_cuda_driverVersion
 
   function get_cuda_runtimeVersion(runtimeVersion, devD)
     import :: deviceDetails
     type(deviceDetails) :: devD
     integer :: runtimeVersion
     integer :: get_cuda_runtimeVersion
   end function get_cuda_runtimeVersion

   function get_tot_global_mem_MB(totalGlobalMem, devD)
     import :: deviceDetails
     type(deviceDetails) :: devD
     real    :: totalGlobalMem
     integer :: get_tot_global_mem_MB
   end function get_tot_global_mem_MB

   function get_tot_global_mem_bytes(totalGlobalMem, devD)
     import :: deviceDetails
     type(deviceDetails) :: devD
     integer*8 :: totalGlobalMem
     integer   :: get_tot_global_mem_bytes
   end function get_tot_global_mem_bytes

   function get_CUDA_cores(idev, major, minor, nmp, &
        cuda_cores_per_mp, ncuda_cores, devD)
     import :: deviceDetails
     type(deviceDetails) :: devD
     integer :: idev
     integer :: major,minor
     integer :: nmp, cuda_cores_per_mp
     integer :: ncuda_cores
     integer :: get_CUDA_cores_c
   end function get_CUDA_cores

   function get_thread_details(warpSize, max_threads_per_mp, &
        max_threads_per_blk, devD)
     import :: deviceDetails
     type(deviceDetails) :: devD
     integer :: warpSize, max_threads_per_mp, max_threads_per_blk
     integer :: get_thread_details
   end function get_thread_details
   
   function get_nregisters(nregisters_per_blk, devD)
     import :: deviceDetails
     type(deviceDetails) :: devD
     integer :: nregisters_per_blk
     integer :: get_nregisters
   end function get_nregisters

   function get_eec_support(is_ecc, idev, devD)
     import :: deviceDetails
     type(deviceDetails) :: devD
     integer :: idev
     integer :: is_ecc
     integer :: get_eec_support
   end function get_eec_support

   function get_peer_to_peer_capabilities(devD)
     import :: deviceDetails
     type(deviceDetails) :: devD
     integer :: get_peer_to_peer_capabilities
   end function get_peer_to_peer_capabilities
#endif
end interface

contains 
  !CONSTRUCTORS
  !> \brief is a deviceQuery constructor
  function constructor_deviceQuery_gpu(devD,idev) result(devQ_gpu)
    type(deviceQuery_gpu) :: devQ_gpu
    type(deviceDetails) :: devD
    integer,optional :: idev
    call devQ_gpu%new_deviceQuery_gpu(devD,idev)
  end function constructor_deviceQuery_gpu

  subroutine new_deviceQuery_gpu(devQ_gpu,devD,idev)
    class(deviceQuery_gpu) :: devQ_gpu
    type(deviceDetails) :: devD
    integer,optional :: idev
    !kill the existing object before allocating a new one
    call devQ_gpu%kill_deviceQuery_gpu
    !get the total number of devices avaible on sys
    call set_devCnt(devQ_gpu,devD)
    if (devQ_gpu%devCnt == 0 ) call devQ_gpu%terminate()

    !allocating arrays for details for n devices
    allocate(        devQ_gpu%totGlobalMem_MB(0:(get_devCnt(devQ_gpu)-1)))
    allocate(     devQ_gpu%totGlobalMem_bytes(0:(get_devCnt(devQ_gpu)-1)))
    allocate(            devQ_gpu%n_MultiProc(0:(get_devCnt(devQ_gpu)-1)))
    allocate(devQ_gpu%n_cudacor_per_MultiProc(0:(get_devCnt(devQ_gpu)-1)))
    allocate(              devQ_gpu%n_cudacor(0:(get_devCnt(devQ_gpu)-1)))
    allocate(            devQ_gpu%SM_suitable(0:(get_devCnt(devQ_gpu)-1)))
    allocate(    devQ_gpu%n_registers_per_blk(0:(get_devCnt(devQ_gpu)-1)))
    allocate(                 devQ_gpu%n_warp(0:(get_devCnt(devQ_gpu)-1)))
    allocate(       devQ_gpu%n_threads_per_mp(0:(get_devCnt(devQ_gpu)-1)))
    allocate(      devQ_gpu%n_threads_per_blk(0:(get_devCnt(devQ_gpu)-1)))
    allocate(             devQ_gpu%dev_is_ecc(0:(get_devCnt(devQ_gpu)-1)))
    allocate(             devQ_gpu%dev_is_p2p(0:MAX_N_GPU-1,0:MAX_N_GPU-1))

    !allocate(                devQ_gpu%devname(0:(get_devCnt(devQ_gpu)-1)))

    call gather_deviceDetails(devQ_gpu,devD,idev)

    !set the existence of the object to true
    devQ_gpu%existence_devQ_gpu = .true.

    return
  end subroutine new_deviceQuery_gpu

  !WORKERS METHODS
  subroutine gather_deviceDetails(devQ_gpu,devD,idev)
    class(deviceQuery_gpu), intent(inout) :: devQ_gpu
    type(deviceDetails) :: devD
    integer,optional :: idev
    integer  :: jdev
    !initialising the data structure and allocating in C

    call set_cuda_drv_version(devQ_gpu,devD)
    call set_cuda_rnt_version(devQ_gpu,devD)

    if ( present(idev) ) then
       call initialiase_deviceDetails(devQ_gpu,idev,devD)
       call               set_devName(devQ_gpu,idev,devD)
       call     set_tot_global_mem_MB(devQ_gpu,idev,devD)
       call  set_tot_global_mem_bytes(devQ_gpu,idev,devD)
       call            set_cuda_cores(devQ_gpu,idev,devD)
       call            set_nregisters(devQ_gpu,idev,devD)
       call        set_thread_details(devQ_gpu,idev,devD)
       call           set_ecc_support(devQ_gpu,idev,devD)
    else
       do jdev = 0,(get_devCnt(devQ_gpu)-1)
          call initialiase_deviceDetails(devQ_gpu,jdev,devD)
          call               set_devName(devQ_gpu,jdev,devD)
          call     set_tot_global_mem_MB(devQ_gpu,jdev,devD)
          call  set_tot_global_mem_bytes(devQ_gpu,jdev,devD)
          call            set_cuda_cores(devQ_gpu,jdev,devD)
          call            set_nregisters(devQ_gpu,jdev,devD)
          call        set_thread_details(devQ_gpu,jdev,devD)
          call           set_ecc_support(devQ_gpu,jdev,devD)
       end do
    end if
       
    !TODO: need to add the logic for multi GPU selection process

    call set_p2p(devQ_gpu,devD)

    return
  end subroutine gather_deviceDetails

  !SETTERS
  !setting the peer to peer
  subroutine set_p2p(devQ_gpu,devD)
    class(deviceQuery_gpu) :: devQ_gpu
    type(deviceDetails) :: devD
    integer :: idev
    integer :: is_ecc
    integer :: rc !return code
    interface  external_c_function_get_peer_to_peer_capabilities
       function get_peer_to_peer_capabilities_c(devD)
         import :: deviceDetails
         type(deviceDetails) :: devD
         integer :: get_peer_to_peer_capabilities_c
       end function get_peer_to_peer_capabilities_c
    end interface external_c_function_get_peer_to_peer_capabilities
#if defined (MACOSX) && defined (CUDA)
    rc = get_peer_to_peer_capabilities_c(devD)
    devQ_gpu%dev_is_p2p = devD%is_p2p
#elif defined (LINUX) && defined (CUDA)
    rc = get_peer_to_peer_capabilities(devD)
    devQ_gpu%dev_is_p2p = devD%is_p2p
#else
    call devQ_gpu%terminate()
#endif
    return
  end subroutine set_p2p
  !setting the ecc support 
  subroutine set_ecc_support(devQ_gpu,idev,devD)
    class(deviceQuery_gpu) :: devQ_gpu
    type(deviceDetails) :: devD
    integer :: idev
    integer :: is_ecc
    integer :: rc !return code
    interface  external_c_function_get_eec_support
       function get_eec_support_c(is_ecc, idev, devD)
         import :: deviceDetails
         type(deviceDetails) :: devD
         integer :: idev
         integer :: is_ecc
         integer :: get_eec_support_c
       end function get_eec_support_c
    end interface external_c_function_get_eec_support
#if defined (MACOSX) && defined (CUDA)
    rc = get_eec_support_c(is_ecc, idev, devD)
    devQ_gpu%dev_is_ecc(idev) = devD%is_ecc
#elif defined (LINUX) && defined (CUDA)
    rc = get_eec_support(is_ecc, idev, devD)
    devQ_gpu%dev_is_ecc(idev) = devD%is_ecc
#else
    call devQ_gpu%terminate()
#endif
    return
  end subroutine set_ecc_support
  !setting up the thread details 
  subroutine set_thread_details(devQ_gpu,idev,devD)
    class(deviceQuery_gpu) :: devQ_gpu
    type(deviceDetails) :: devD
    integer :: idev
    integer :: warpSize, max_threads_per_mp, max_threads_per_blk
    integer :: rc !return code
    interface external_c_function_get_thread_details
       function get_thread_details_c(warpSize, max_threads_per_mp, max_threads_per_blk, devD)
         import :: deviceDetails
         type(deviceDetails) :: devD
         integer :: warpSize, max_threads_per_mp, max_threads_per_blk
         integer :: get_thread_details_c
       end function get_thread_details_c
    end interface external_c_function_get_thread_details
#if defined (MACOSX) && defined (CUDA)
    rc = get_thread_details_c(warpSize, max_threads_per_mp, max_threads_per_blk, devD)
    devQ_gpu%n_warp (idev) = devD%warpSze
    devQ_gpu%n_threads_per_mp (idev) = devD%maxthreads_per_mp
    devQ_gpu%n_threads_per_blk (idev) = devD%maxthreads_per_blk
#elif defined (LINUX) && defined (CUDA)
    rc = get_thread_details(warpSize, max_threads_per_mp, max_threads_per_blk, devD)
    devQ_gpu%n_warp (idev) = devD%warpSze
    devQ_gpu%n_threads_per_mp (idev) = devD%maxthreads_per_mp
    devQ_gpu%n_threads_per_blk (idev) = devD%maxthreads_per_blk
#else
    call devQ_gpu%terminate()
#endif
    return
  end subroutine set_thread_details
  !setting up the registers
  subroutine set_nregisters(devQ_gpu,idev,devD)
    class(deviceQuery_gpu) :: devQ_gpu
    type(deviceDetails) :: devD
    integer :: idev
    integer :: nregisters_per_blk
    integer :: rc !return code
    interface  external_c_function_get_nregisters
       function get_nregisters_c(nregisters_per_blk, devD)
         import :: deviceDetails
         type(deviceDetails) :: devD
         integer :: nregisters_per_blk
         integer :: get_nregisters_c
       end function get_nregisters_c
    end interface external_c_function_get_nregisters
#if defined (MACOSX) && defined (CUDA)
    rc = get_nregisters_c(nregisters_per_blk, devD)
    devQ_gpu%n_registers_per_blk(idev) = devD%nregisters_per_blk
#elif defined (LINUX) && defined (CUDA)
    rc = get_nregisters(nregisters_per_blk, devD)
    devQ_gpu%n_registers_per_blk(idev) = devD%nregisters_per_blk
#else
    call devQ_gpu%terminate()
#endif
    return
  end subroutine set_nregisters

  !setting the CUDA cores
  subroutine set_cuda_cores(devQ_gpu,idev,devD)
    class(deviceQuery_gpu) :: devQ_gpu
    type(deviceDetails) :: devD
    integer :: idev
    integer :: nmp, cuda_cores_per_mp, ncuda_cores, SMsuitable
    integer :: rc !return code
    interface  external_c_function_CUDA_cores_c
       function get_CUDA_cores_c(idev, major, minor, nmp, cuda_cores_per_mp, ncuda_cores, devD)
         import :: deviceDetails
         type(deviceDetails) :: devD
         integer :: idev
         integer :: major,minor
         integer :: nmp, cuda_cores_per_mp
         integer :: ncuda_cores
         integer :: get_CUDA_cores_c
       end function get_CUDA_cores_c
    end interface external_c_function_CUDA_cores_c

#if defined (MACOSX) && defined (CUDA)
    rc = get_CUDA_cores_c(idev, GPU_CARD_MAJOR, GPU_CARD_MINOR, nmp, cuda_cores_per_mp, ncuda_cores, devD)
    devQ_gpu%n_MultiProc(idev) = devD%nMultiProc
    devQ_gpu%n_cudacor_per_MultiProc(idev) = devD%ncc_per_mp
    devQ_gpu%n_cudacor(idev) = devD%ncc
    devQ_gpu%SM_suitable(idev) = devD%is_SMsuitable
#elif defined (LINUX) && defined (CUDA)
    rc = get_CUDA_cores(idev, GPU_CARD_MAJOR, GPU_CARD_MINOR, nmp, cuda_cores_per_mp, ncuda_cores, devD)
    devQ_gpu%n_MultiProc(idev) = devD%nMultiProc
    devQ_gpu%n_cudacor_per_MultiProc(idev) = devD%ncc_per_mp
    devQ_gpu%n_cudacor(idev) = devD%ncc
    devQ_gpu%SM_suitable(idev) = devD%is_SMsuitable
#else
    call devQ_gpu%terminate()
#endif
    return
  end subroutine set_cuda_cores

  !setting the total memory for that device
  subroutine set_tot_global_mem_MB(devQ_gpu,idev,devD)
    class(deviceQuery_gpu) :: devQ_gpu
    type(deviceDetails) :: devD
    integer  :: idev
    real(sp) :: totalGlobalMem
    integer  :: rc !return code
    interface  external_c_function_tot_global_mem_MB
       function get_tot_global_mem_MB_c(totalGlobalMem, devD)
         import :: deviceDetails
         type(deviceDetails) :: devD
         real    :: totalGlobalMem
         integer :: get_tot_global_mem_MB_c
       end function get_tot_global_mem_MB_c
    end interface external_c_function_tot_global_mem_MB
#if defined (MACOSX) && defined (CUDA)
    rc = get_tot_global_mem_MB_c(totalGlobalMem, devD)
    devQ_gpu%totGlobalMem_MB(idev) = devD%tot_global_mem_MB
#elif defined (LINUX) && defined (CUDA)
    rc = get_tot_global_mem_MB(totalGlobalMem, devD)
    devQ_gpu%totGlobalMem_MB(idev) = devD%tot_global_mem_MB
#else
    call devQ_gpu%terminate()
#endif
    return
  end subroutine set_tot_global_mem_MB

  subroutine set_tot_global_mem_bytes(devQ_gpu,idev,devD)
    class(deviceQuery_gpu) :: devQ_gpu
    type(deviceDetails) :: devD
    integer   :: idev
    integer*8 :: totalGlobalMem
    integer   :: rc !return code
    interface  external_c_function_tot_global_mem_bytes
       function get_tot_global_mem_bytes_c(totalGlobalMem, devD)
         import :: deviceDetails
         type(deviceDetails) :: devD
         integer*8 :: totalGlobalMem
         integer   :: get_tot_global_mem_bytes_c
       end function get_tot_global_mem_bytes_c
    end interface external_c_function_tot_global_mem_bytes
#if defined (MACOSX) && defined (CUDA)
    rc = get_tot_global_mem_bytes_c(totalGlobalMem, devD)
    devQ_gpu%totGlobalMem_bytes(idev) = devD%tot_global_mem_bytes
#elif defined (LINUX) && defined (CUDA)
    rc = get_tot_global_mem_bytes(totalGlobalMem, devD)
    devQ_gpu%totGlobalMem_bytes(idev) = devD%tot_global_mem_bytes
#else
    call devQ_gpu%terminate()
#endif
    return
  end subroutine set_tot_global_mem_bytes

  !setting the runtime version
  subroutine set_cuda_rnt_version(devQ_gpu,devD)
    class(deviceQuery_gpu) :: devQ_gpu
    type(deviceDetails) :: devD
    integer  :: runtimeVersion
    integer  :: rc !return code
    interface external_c_function_runtimeVersion
       function get_cuda_runtimeVersion_c(runtimeVersion, devD)
         import :: deviceDetails
         type(deviceDetails) :: devD
         integer :: runtimeVersion
         integer :: get_cuda_runtimeVersion_c
       end function get_cuda_runtimeVersion_c
    end interface external_c_function_runtimeVersion
#if defined (MACOSX) && defined (CUDA)
    rc = get_cuda_runtimeVersion_c(runtimeVersion, devD)
    devQ_gpu%runVer = devD%d_runver
#elif defined (LINUX) && defined (CUDA)
    rc = get_cuda_runtimeVersion(runtimeVersion, devD)
    devQ_gpu%runVer = devD%d_runver
#else
    call devQ_gpu%terminate()
#endif
    return
  end subroutine set_cuda_rnt_version

  !setter the driver version
  subroutine set_cuda_drv_version(devQ_gpu,devD)
    class(deviceQuery_gpu) :: devQ_gpu
    type(deviceDetails) :: devD
    integer  :: driverVersion
    integer  :: rc !return code
    interface external_c_function_driverVersion
       function get_cuda_driverVersion_c(driverVersion, devD)
         import :: deviceDetails
         type(deviceDetails) :: devD
         integer :: driverVersion
         integer :: get_cuda_driverVersion_c
       end function get_cuda_driverVersion_c
    end interface external_c_function_driverVersion
#if defined (MACOSX) && defined (CUDA)
    rc = get_cuda_driverVersion_c(driverVersion, devD)
    devQ_gpu%devVer = devD%d_ver
#elif defined (LINUX) && defined (CUDA)
    rc = get_cuda_driverVersion(driverVersion, devD)
    devQ_gpu%devVer = devD%d_ver
#else
    call devQ_gpu%terminate()
#endif
    return
  end subroutine set_cuda_drv_version

  !setter for the devices name
  subroutine set_devName(devQ_gpu,idev,devD)
    class(deviceQuery_gpu) :: devQ_gpu
    type(deviceDetails) :: devD
    integer :: idev
    integer :: rc !return code
    interface external_c_function_devname
       function get_dev_Name_c(idev, devD)
         import :: deviceDetails
         type(deviceDetails) :: devD
         integer :: idev
         integer :: get_dev_Name_c
       end function get_dev_Name_c
    end interface external_c_function_devname
#if defined (MACOSX) && defined (CUDA)
    rc = get_dev_Name_c(idev,devD)
    call c_f_pointer(devD%dev_name, devQ_gpu%devname, [DEVNAME_STRING_LENGTH])
#elif defined (LINUX) && defined (CUDA)
    rc = get_dev_Name(idev,devD)
    call c_f_pointer(devD%dev_name, devQ_gpu%devname, [DEVNAME_STRING_LENGTH])
#else
    call devQ_gpu%terminate()
#endif
    return
  end subroutine set_devName

  !setter for the device counts these are the available GPU on sys
  subroutine set_devCnt(devQ_gpu,devD)
    class(deviceQuery_gpu), intent(inout) :: devQ_gpu
    type(deviceDetails) :: devD
    integer :: rc !return code
    interface external_c_functions
       function get_dev_count_c(devD)
         import :: deviceDetails
         type(deviceDetails) :: devD
         integer :: get_dev_count_c
       end function get_dev_count_c
    end interface external_c_functions
#if defined (MACOSX) && defined (CUDA)
    rc = get_dev_count_c(devD)
    devQ_gpu%devCnt = devD%ndev
#elif defined (LINUX) && defined (CUDA)
    rc = get_dev_count(devD)
    devQ_gpu%devCnt = devD%ndev
#else
    call devQ_gpu%terminate()
#endif
    return
  end subroutine set_devCnt

  !GETTERS
  function get_devCnt(devQ_gpu) result(devCnt_out)
    class(deviceQuery_gpu), intent(in) :: devQ_gpu
    integer                     :: devCnt_out
    devCnt_out = devQ_gpu%devCnt
  end function get_devCnt

  function get_devname(devQ_gpu,idev) result(devName_out)
    class(deviceQuery_gpu), intent(in) :: devQ_gpu
    character(len=1),allocatable :: devName_out(:)
    integer :: idev
    allocate(devName_out(0:DEVNAME_STRING_LENGTH))
    devName_out = devQ_gpu%devname
  end function get_devname

  function get_cuda_drv_version(devQ_gpu) result(devVer_out)
    class(deviceQuery_gpu), intent(in) :: devQ_gpu
    real(sp)             :: devVer_out
    devVer_out = devQ_gpu%devVer
  end function get_cuda_drv_version 

  function get_cuda_rnt_version(devQ_gpu) result(runVer_out)
    class(deviceQuery_gpu), intent(in) :: devQ_gpu
    real(sp)             :: runVer_out
    runVer_out = devQ_gpu%runVer
  end function get_cuda_rnt_version

  function get_totGlobalMem_MB(devQ_gpu,idev) result(mem_MB_out)
    class(deviceQuery_gpu), intent(in) :: devQ_gpu
    real(sp)             :: mem_MB_out
    integer :: idev
    mem_MB_out = devQ_gpu%totGlobalMem_MB(idev)
  end function get_totGlobalMem_MB

  function get_totGlobalMem_bytes(devQ_gpu,idev) result(mem_bytes_out)
    class(deviceQuery_gpu), intent(in) :: devQ_gpu
    integer*8             :: mem_bytes_out
    integer :: idev
    mem_bytes_out = devQ_gpu%totGlobalMem_bytes(idev)
  end function get_totGlobalMem_bytes

  function get_dev_cuda_mp(devQ_gpu,idev) result(mp_out)
    class(deviceQuery_gpu), intent(in) :: devQ_gpu
    integer             :: mp_out
    integer :: idev
    mp_out = devQ_gpu%n_MultiProc(idev)
  end function get_dev_cuda_mp

  function get_dev_cuda_core_per_mp(devQ_gpu,idev) result(core_per_mp_out)
    class(deviceQuery_gpu), intent(in) :: devQ_gpu
    integer             :: core_per_mp_out
    integer :: idev
    core_per_mp_out = devQ_gpu%n_cudacor_per_MultiProc(idev)
  end function get_dev_cuda_core_per_mp

  function get_dev_cuda_core(devQ_gpu,idev) result(core_out)
    class(deviceQuery_gpu), intent(in) :: devQ_gpu
    integer             :: core_out
    integer :: idev
    core_out = devQ_gpu%n_cudacor(idev)
  end function get_dev_cuda_core

  function get_dev_SM_suitable(devQ_gpu,idev) result(SMsuit_out)
    class(deviceQuery_gpu), intent(in) :: devQ_gpu
    integer             :: SMsuit_out
    integer :: idev
    SMsuit_out = devQ_gpu%SM_suitable(idev)
  end function get_dev_SM_suitable

  function get_dev_nregisters_per_blk(devQ_gpu,idev) result(nregisters_out)
    class(deviceQuery_gpu), intent(in) :: devQ_gpu
    integer             :: nregisters_out
    integer :: idev
    nregisters_out = devQ_gpu%n_registers_per_blk(idev)
  end function get_dev_nregisters_per_blk

  function get_dev_nwarp(devQ_gpu,idev) result(nwarp_out)
    class(deviceQuery_gpu), intent(in) :: devQ_gpu
    integer             :: nwarp_out
    integer :: idev
    nwarp_out = devQ_gpu%n_warp(idev)
  end function get_dev_nwarp

  function get_dev_n_threads_per_mp(devQ_gpu,idev) result(nthreads_per_mp_out)
    class(deviceQuery_gpu), intent(in) :: devQ_gpu
    integer             :: nthreads_per_mp_out
    integer :: idev
    nthreads_per_mp_out = devQ_gpu%n_threads_per_mp(idev)
  end function get_dev_n_threads_per_mp

  function get_dev_n_threads_per_blk(devQ_gpu,idev) result(nthreads_per_blk_out)
    class(deviceQuery_gpu), intent(in) :: devQ_gpu
    integer             :: nthreads_per_blk_out
    integer :: idev
    nthreads_per_blk_out = devQ_gpu%n_threads_per_blk(idev)
  end function get_dev_n_threads_per_blk

  function get_eecsupport(devQ_gpu,idev) result(ecc_out)
    class(deviceQuery_gpu), intent(in) :: devQ_gpu
    integer             :: ecc_out
    integer :: idev
    ecc_out = devQ_gpu%dev_is_ecc(idev)
  end function get_eecsupport

  function get_p2p(devQ_gpu) result(p2p_out)
    class(deviceQuery_gpu), intent(in) :: devQ_gpu
    integer,dimension(0:MAX_N_GPU-1,0:MAX_N_GPU-1) :: p2p_out
    p2p_out = devQ_gpu%dev_is_p2p
  end function get_p2p

  function get_devD(devQ_gpu) result(devD_out)
    class(deviceQuery_gpu), intent(in) :: devQ_gpu
    type(deviceDetails) :: devD_out
    devD_out  = devQ_gpu%devD
  end function get_devD
  
  !INITIALISERS
  subroutine initialiase_deviceDetails(devQ_gpu,idev,devD)
    class(deviceQuery_gpu), intent(inout) :: devQ_gpu
    type(deviceDetails) :: devD
    integer :: idev
    integer :: rc !return code
    interface external_c_function_init
       function deviceDetails_init_c(idev,devD)
         import :: deviceDetails
         type(deviceDetails) :: devD
         integer :: idev
         integer :: deviceDetails_init_c
       end function deviceDetails_init_c
    end interface external_c_function_init
#if defined (MACOSX) && defined (CUDA)
    rc = deviceDetails_init_c(idev,devD)
#elif defined (LINUX) && defined (CUDA)
    rc = deviceDetails_init(idev,devD)
#else
    call devQ_gpu%terminate()
#endif
    return
  end subroutine initialiase_deviceDetails

  !GREETERS
  !> \brief hello greeting routine for the object
  subroutine hello_deviceQuery_gpu(err)
    implicit none
    integer :: err
    !start of the execution commands
    write(*,*) "Hello deviceQuery GPU world"
    write(*,*)
    return
  end subroutine hello_deviceQuery_gpu

  !> \brief bye greeting routine for the object
  subroutine bye_deviceQuery_gpu()
    implicit none
    !start of the execution commands
    write(*,*) "Bye deviceQuery GPU world"
    write(*,*)
    return
  end subroutine bye_deviceQuery_gpu

  !TERMINATORS because of incompatiblilties errors
  subroutine terminate(devQ_gpu)
    class(deviceQuery_gpu) :: devQ_gpu
    call devQ_gpu%kill_deviceQuery_gpu
    write(*,*)"**************************WARNING*******************************"
    write(*,*)"*There are no GPU devices available for computation            *"
    write(*,*)"*You need to check that your hardware is suitable and has GPU  *"
    write(*,*)"*computational capacities or that you have compiled with       *"
    write(*,*)"*-DCUDA to acces the CUDA environment computation.             *"
    write(*,*)"****************************************************************"
    call bye_deviceQuery_gpu()
    stop
    return
  end subroutine terminate

  !WARNERS

  subroutine get_warning_dataStruct_gpu(devQ_gpu)
    !implicit none
    class(deviceQuery_gpu) :: devQ_gpu
    write(*,*)
    write(*,*)"***********data structure and getter do not match***************"
    write(*,*)
    return
  end subroutine get_warning_dataStruct_gpu

  !DESTRUCTOR
  !> \brief is a deviceQuery_gpu destructor
  subroutine kill_deviceQuery_gpu(devQ_gpu)
    class(deviceQuery_gpu), intent(inout) :: devQ_gpu
    if ( devQ_gpu%existence_devQ_gpu) then

if(allocated(devQ_gpu%totGlobalMem_MB    )) deallocate(devQ_gpu%totGlobalMem_MB)
if(allocated(devQ_gpu%totGlobalMem_bytes )) deallocate(devQ_gpu%totGlobalMem_bytes)
if(allocated(devQ_gpu%n_MultiProc        )) deallocate(devQ_gpu%n_MultiProc)
if(allocated(devQ_gpu%n_cudacor_per_MultiProc)) deallocate(devQ_gpu%n_cudacor_per_MultiProc)
if(allocated(devQ_gpu%n_cudacor          )) deallocate(devQ_gpu%n_cudacor)
if(allocated(devQ_gpu%SM_suitable        )) deallocate(devQ_gpu%SM_suitable)
if(allocated(devQ_gpu%n_registers_per_blk)) deallocate(devQ_gpu%n_registers_per_blk)
if(allocated(devQ_gpu%n_warp             )) deallocate(devQ_gpu%n_warp)
if(allocated(devQ_gpu%n_threads_per_mp   )) deallocate(devQ_gpu%n_threads_per_mp)
if(allocated(devQ_gpu%n_threads_per_blk  )) deallocate(devQ_gpu%n_threads_per_blk)
if(allocated(devQ_gpu%n_registers_per_blk)) deallocate(devQ_gpu%n_registers_per_blk)
if(allocated(devQ_gpu%dev_is_ecc         )) deallocate(devQ_gpu%dev_is_ecc)
if(allocated(devQ_gpu%dev_is_p2p         )) deallocate(devQ_gpu%dev_is_p2p)

       devQ_gpu%existence_devQ_gpu = .false.

    end if
    return
  end subroutine kill_deviceQuery_gpu

end module simple_deviceQuery_gpu
