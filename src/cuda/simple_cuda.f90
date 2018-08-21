module simple_cuda
include 'simple_lib.f08'
use, intrinsic :: ISO_C_BINDING
use CUDA
#include "simple_cuda_handle.inc"
implicit none

 type( cudaDeviceProp ) :: cuda_properties
 integer                :: cuda_deviceCount
 integer                :: cuda_currentDevice

contains

    subroutine check_cuda_device
        use CUDA
        implicit none
        integer :: device
        integer :: deviceCount
        integer :: gpuDeviceCount
        type( cudaDeviceProp ) :: properties
        integer(KIND(cudaSuccess)) :: cudaResultCode, err
        logical :: error_found
        deviceCount=0
        gpuDeviceCount=0
        err = cudaThreadSynchronize()

        ! Obtain the device count from the system and check for errors
        cudaResultCode = cudaGetDeviceCount(deviceCount)
        if (cudaResultCode /= cudaSuccess) then
            deviceCount = 0
        end if
        call my_cudaErrorCheck(cudaResultCode, error_found)
        if(error_found )&
        THROW_HARD("CUDA device count failed")
        ! Check for the device to be an emulator. If not, increment the counter
        do device = 0, deviceCount - 1, 1
            cudaResultCode = cudaGetDeviceProperties(properties, device)
            if (properties%major /= 9999) gpuDeviceCount = gpuDeviceCount + 1
            call my_cudaErrorCheck(cudaResultCode, error_found)
            if(error_found) THROW_HARD("CUDA device property failed "//int2str(device))
        end do
        cuda_deviceCount= gpuDeviceCount
        print*, gpuDeviceCount, " GPU CUDA device(s) found"

        do device = 0, deviceCount - 1, 1
            print*, "    Cuda device ", device
            cudaResultCode = cudaGetDeviceProperties(properties,device)
             call my_cudaErrorCheck(cudaResultCode, error_found)
             if( error_found) THROW_HARD("CUDA device property failed "//int2str(device))
             if(device==0) cuda_properties = properties
             call print_cuda_properties(properties)
         end do
         err = cudaThreadSynchronize()
         print *, " Check CUDA completed "
    end subroutine check_cuda_device

    subroutine set_cuda_device(d)
        integer , intent(in) :: d
        integer :: device
        integer :: deviceCount
        type( cudaDeviceProp ) :: properties
        integer(KIND(cudaSuccess)) :: cudaResultCode
        logical :: error_found
        deviceCount=0
        ! Obtain the device count from the system and check for errors
        cudaResultCode = cudaGetDeviceCount(deviceCount)
        call my_cudaErrorCheck(cudaResultCode, error_found)
        if( error_found) THROW_HARD("cudaGetDeviceCount failed")

        if (d >= deviceCount) then
            print *, "simple_cuda::set_cuda_device index too high "
            return
        endif

        cudaResultCode = cudaGetDeviceProperties(properties, d)
        call my_cudaErrorCheck(cudaResultCode, error_found)
        if( error_found) THROW_HARD("cudaGetDeviceProperties failed")
        cuda_properties = properties

        cudaResultCode = cudaSetDevice(d)
        call my_cudaErrorCheck(cudaResultCode, error_found)
        if( error_found) THROW_HARD("cudaSetDevice failed")
        cuda_currentDevice = d
        print *, " New CUDA device set successfully ", cuda_currentDevice
        call print_cuda_properties(cuda_properties)
    end subroutine set_cuda_device

    function get_cuda_property(d,properties) result(error_found)
         integer,                intent(in)  :: d
         type( cudaDeviceProp ), intent(out) :: properties
         integer :: deviceCount
         integer(KIND(cudaSuccess)) :: cudaResultCode
         logical :: error_found
         deviceCount=0

        ! Obtain the device count from the system and check for errors
        cudaResultCode = cudaGetDeviceCount(deviceCount)
        call my_cudaErrorCheck(cudaResultCode, error_found)
        if( error_found) THROW_HARD("cudaGetDeviceCount failed")

        if (d >= deviceCount) then
            print *, "simple_cuda::get_cuda_device index too high "
            error_found=.false.
            return
        end if

        cudaResultCode = cudaGetDeviceProperties(properties, d)
        call my_cudaErrorCheck(cudaResultCode, error_found)
        if( error_found) THROW_HARD("cudaGetDeviceProperties failed")

    end function get_cuda_property

    subroutine print_cuda_properties(properties, flags)
        type( cudaDeviceProp ), intent(in) :: properties
        integer, intent(in), optional :: flags
        integer :: flags_here
        flags_here=0
        print*," GPU Device Name             : ", properties%name
        print"('  CUDA Major Version          : ', i0)", properties%major
        print"('  CUDA Minor Version          : ', i0)", properties%minor
        print"('  Total Global Memory     (MB): ',f9.3)", real(properties%totalGlobalMem,dp)/real(1024._dp**2._dp,dp)
        print"('  Shared Memory per Block (kB): ',f9.3)", real(properties%sharedMemPerBlock)/real(1024)
        print"('  Regs Memory per Block   (kB): ',f9.3)", real(properties%regsPerBlock)/real(1024)
        print"('  Warp Size                   : ', i0)", properties%warpSize
        print"('  Memory Pitch            (GB): ',f9.3)",real(properties%memPitch)/real(1024._dp**3._dp,dp)
        print"('  Max Threads per block       : ', i0)", properties%maxThreadsPerBlock
        print"('  Max Threads per dim         : ', i0,1x,i0,1x,i0)", properties%maxThreadsDim
        print"('  Max Grid Size               : ', i0,1x,i0,1x,i0)", properties%maxGridSize
        print"('  GPU Clock Rate              : ', i0)", properties%clockRate
        print"('  Total Const Memory      (kB): ',f9.3)", real(properties%totalConstMem)/real(1024)
        print"('  Texture Alignment           : ', i0)", properties%textureAlignment
        print"('  Device Overlap              : ', i0)", properties%deviceOverlap
        print"('  MultiProcessorCount         : ', i0)", properties%multiProcessorCount
        print"('  Kernel Exec TimeoutEnabled  : ', i0)", properties%kernelExecTimeoutEnabled
        print"('  Integrated                  : ', i0)", properties%integrated
        print"('  CanMapHostMemory            : ', i0)", properties%canMapHostMemory
        print"('  Compute Mode                : ', i0)", properties%computeMode
        print"('  MaxTexture1D                : ', i0)", properties%maxTexture1D
        print"('  MaxTexture2D                : ', i0,1x,i0)", properties%maxTexture2D
        print"('  MaxTexture3D                : ', i0,1x,i0,1x,i0)", properties%maxTexture3D
        !print*, "    maxTexture2DArray(3)     :", properties%maxTexture2DArray
        print"('  Surface Alignment           : ', i0)", properties%surfaceAlignment
        print"('  Concurrent Kernels          : ', i0)", properties%concurrentKernels
        print"('  ECCEnabled                  : ', i0)", properties%ECCEnabled
        print"('  pciBusID                    : ', i0)", properties%pciBusID
        print"('  pciDeviceID                 : ', i0)", properties%pciDeviceID ! PCI device (sometimes called slot) identifier of the device
        if(flags_here > 0) then
            print*, "    pci Domain ID                               :", properties%pciDomainID
            print*, "    TCC Driver                                  :", properties%tccDriver
            print*, "    Async Engine Count                          :", properties%asyncEngineCount
            print*, "    Unified Addressing                          :", properties%unifiedAddressing
            print*, "    memory Clock Rate                           :", properties%memoryClockRate
            print*, "    Memory Bus Width                            :", properties%memoryBusWidth
            print*, "    L2 Cache Size                    (kB)       :", real(properties%l2CacheSize)/real(1024)
            print*, "    Max Threads Per Multi Processor             :", properties%maxThreadsPerMultiProcessor
            print*, "    Stream Priorities Supported                 :", properties%streamPrioritiesSupported
            print*, "    Global L1 Cache Supported                   :", properties%globalL1CacheSupported
            print*, "    Local L1 Cache Supported                    :", properties%localL1CacheSupported
            print*, "    Shared Mem Per Multiprocessor    (kB)       :", real(properties%sharedMemPerMultiprocessor)/real(1024)
            print*, "    Regs Per Multiprocessor          (kB)       :", real(properties%regsPerMultiprocessor)/real(1024)
            print*, "    Managed Memory                              :", properties%managedMemory !! in CUDA 8.0 changed from managedMemorySupported
            print*, "    is Multi Gpu Board                          :", properties%isMultiGpuBoard
            print*, "    Multi Gpu Board Group ID                    :", properties%multiGpuBoardGroupID
            print*, "    Single To Double Precision Perf Ratio       :", properties%singleToDoublePrecisionPerfRatio
            print*, "    Pageable Memory Access                      :", properties%pageableMemoryAccess
            print*, "    Concurrent Managed Access                   :", properties%concurrentManagedAccess
            print*, "    ComputeP reemption Supported                :", properties%computePreemptionSupported
            print*, "    Can Use Host Pointer For RegisteredMem      :", properties%canUseHostPointerForRegisteredMem
            print*, "    Cooperative Launch                          :", properties%cooperativeLaunch
            print*, "    Cooperative Multi Device Launch             :", properties%cooperativeMultiDeviceLaunch
            print*, "    Pageable Memory Access Uses Host PageTables :", properties%pageableMemoryAccessUsesHostPageTables
            print*, "    Direct Managed Mem Access From Host         :",properties%directManagedMemAccessFromHost
        end if
    end subroutine print_cuda_properties


    subroutine cuda_query_version
        integer (KIND(cudaSuccess)) ::  ierr
        integer(c_int) ::    cuVer
        logical :: error_found
        error_found=.false.
        ierr = cudaRuntimeGetVersion(cuVer)
        call my_cudaErrorCheck(ierr, error_found)
        if( error_found) THROW_HARD("cudaGetDeviceProperties failed ")
        write (*,'(A,I0)') 'CUDA Runtime Version: ', cuVer
    end subroutine cuda_query_version
    subroutine cuda_query_driver_version
        integer (KIND(cudaSuccess)) ::  ierr
        integer(c_int):: driverVersion
        logical :: error_found
        error_found=.false.
        ierr = cudaDriverGetVersion(driverVersion)
        call my_cudaErrorCheck(ierr,error_found)
        if (error_found  ) THROW_HARD ("cudaDriverGetVersion failed")
        write (*,'(A,I0)') 'CUDA Driver Version: ', driverVersion
     end subroutine cuda_query_driver_version
     subroutine cuda_query_device_count
         integer (KIND(cudaSuccess)) ::  ierr
         integer(c_int)::  deviceCount
        logical :: error_found
        error_found=.false.
        ierr = cudaGetDeviceCount(deviceCount)
        call my_cudaErrorCheck(ierr,error_found)
        if (error_found  ) THROW_HARD("cudaGetDeviceCount failed")
        write (*,'(A,I0)') 'CUDA Devices Count: ', deviceCount
    end subroutine cuda_query_device_count
     subroutine cuda_query_thread_limit
         integer (KIND(cudaLimitStackSize))::limitStackSize= cudaLimitStackSize
          integer(c_int):: pValue
         integer (KIND(cudaSuccess)) ::  ierr
         logical :: error_found
         error_found=.false.
         ierr = cudaThreadGetLimit(pValue,limitStackSize)
         call my_cudaErrorCheck(ierr,error_found)
         if (error_found  ) THROW_HARD ("cudaThreadGetLimit failed")
         write (*,'(A,I0)') 'CUDA Thread Limit : ', pValue
     end subroutine cuda_query_thread_limit
     subroutine cuda_thread_synchronize
         integer (KIND(cudaSuccess)) ::  ierr
         logical :: error_found
         error_found=.false.
         ierr = cudaThreadSynchronize()
         call my_cudaErrorCheck(ierr,error_found)
         if (error_found  ) THROW_HARD ( "cudaThreadSynchronize failed")
         print *," cudaThreadSynchronized "
     end subroutine cuda_thread_synchronize
     subroutine cuda_print_mem_info
         integer(c_int):: free, total
         integer (KIND(cudaSuccess)) ::  ierr
         logical :: error_found
         error_found=.false.
         ierr = cudaMemGetInfo(free,total)
         call my_cudaErrorCheck(ierr,error_found)
         if (error_found  ) THROW_HARD ("cudaMemGetInfo failed")
         write (*,'(A,I0)') 'CUDA Memory Info (free/total) : ', free,total
     end subroutine cuda_print_mem_info
 end module simple_cuda
