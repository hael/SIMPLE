!@descr: contains utilities for gpu offloading with OpenMP and
! bindings to CUDAToolkit
module simple_gpu_utils
use simple_core_module_api
use, intrinsic :: iso_c_binding, only: c_int, c_ptr
implicit none
#include "simple_local_flags.inc"

#ifdef USE_OPENMP_OFFLOAD
interface

    integer(c_int) function cudaSetDevice(device) bind(C, name='cudaSetDevice')
       import :: c_int
       integer(c_int), value :: device
    end function cudaSetDevice

    integer(c_int) function cudaDeviceSynchronize() bind(C, name='cudaDeviceSynchronize')
       import :: c_int
    end function cudaDeviceSynchronize

    integer(c_int) function cufftPlan1d(plan, nx, ffttype, batch) bind(C, name='cufftPlan1d')
       import :: c_int
       integer(c_int) :: plan
       integer(c_int), value :: nx, ffttype, batch
    end function cufftPlan1d

    integer(c_int) function cufftPlan2d(plan, nx, ny, ffttype) bind(C, name='cufftPlan2d')
       import :: c_int
       integer(c_int) :: plan
       integer(c_int), value :: nx, ny, ffttype
    end function cufftPlan2d

    integer(c_int) function cufftExecR2C(plan, idata, odata) bind(C, name='cufftExecR2C')
       import :: c_int, c_ptr
       integer(c_int), value :: plan
       type(c_ptr), value :: idata
       type(c_ptr), value :: odata
    end function cufftExecR2C

    integer(c_int) function cufftExecC2R(plan, idata, odata) bind(C, name='cufftExecC2R')
       import :: c_int, c_ptr
       integer(c_int), value :: plan
       type(c_ptr), value :: idata
       type(c_ptr), value :: odata
    end function cufftExecC2R

    integer(c_int) function cufftDestroy(plan) bind(C, name='cufftDestroy')
       import :: c_int
       integer(c_int), value :: plan
    end function cufftDestroy
end interface
#endif

interface set_offload_device
    module procedure set_offload_device_1, set_offload_device_2
end interface

contains

    subroutine set_offload_device_1( cline, id )
        use simple_cmdline, only: cmdline
        class(cmdline), intent(in)  :: cline
        integer,        intent(out) :: id
        integer :: ndevices
        id = -1
        if( cline%defined('device') )then
            ! device id provided
            id = cline%get_iarg('device')
            call set_offload_device_2( id )
        else
            ! use default device
#ifdef USE_OPENMP_OFFLOAD
            ndevices = omp_get_num_devices()
            if( ndevices < 1 ) then
                THROW_WARN('A device could not be set: No device identified')
                return
            endif
            id = omp_get_default_device()
            call set_offload_device_2( id )
#endif            
        endif
    end subroutine set_offload_device_1

    subroutine set_offload_device_2( id )
        integer, intent(in) :: id
        integer :: ndevices, host_device
#ifdef USE_OPENMP_OFFLOAD
        ndevices = omp_get_num_devices()
        if( ndevices < 1 ) then
            THROW_WARN('Device could not be set: No device identified')
            return
        endif
        if( (id > ndevices) .or. (id < 0) ) then
            THROW_WARN('Device could not be set: Invalid device ID provided')
            return
        endif
        host_device = omp_get_initial_device()
        if( host_device == id ) then
            THROW_WARN('Device could not be set: device ID provided is the host')
        endif
        call omp_set_default_device( id )
#endif
    end subroutine set_offload_device_2

end module simple_gpu_utils
