!@descr: contains utilities for gpu offloading with OpenMP and
! bindings to CUDAToolkit
module simple_gpu_utils
use simple_core_module_api
use, intrinsic :: iso_c_binding, only: c_char, c_float, c_int, c_null_char, c_ptr, c_size_t
implicit none
#include "simple_local_flags.inc"

integer(c_int), parameter :: CUFFT_R2C             = int(z'2a', c_int)
integer(c_int), parameter :: CUFFT_R2C_MANY        = int(z'2a', c_int)
integer(c_int), parameter :: CUFFT_C2R             = int(z'2c', c_int)
integer(c_int), parameter :: CUDA_SUCCESS          = 0_c_int
integer(c_int), parameter :: CUFFT_SUCCESS         = 0_c_int
integer(c_int), parameter :: CUBLAS_STATUS_SUCCESS = 0_c_int
integer(c_int), parameter :: CUBLAS_OP_N           = 0_c_int
integer(c_int), parameter :: CUBLAS_OP_T           = 1_c_int
integer(c_int), parameter :: CUBLAS_OP_C           = 2_c_int

#ifdef USE_OPENMP_OFFLOAD
    ! Mirrors the C cudaDeviceProp struct layout (CUDA 11/12, 64-bit Linux).
    ! Explicit pad1 aligns totalGlobalMem to 8 bytes after luidDeviceNodeMask.
    ! pad2 extends the type to 2048 bytes
    type, bind(C) :: cuda_device_prop
        character(c_char)  :: name(256)
        character(c_char)  :: uuid(16)
        character(c_char)  :: luid(8)
        integer(c_int)     :: luidDeviceNodeMask
        character(c_char)  :: pad1(4)                ! alignment padding
        integer(c_size_t)  :: totalGlobalMem
        integer(c_size_t)  :: sharedMemPerBlock
        integer(c_int)     :: regsPerBlock
        integer(c_int)     :: warpSize
        integer(c_size_t)  :: memPitch
        integer(c_int)     :: maxThreadsPerBlock
        integer(c_int)     :: maxThreadsDim(3)
        integer(c_int)     :: maxGridSize(3)
        integer(c_int)     :: clockRate
        integer(c_size_t)  :: totalConstMem
        integer(c_int)     :: major
        integer(c_int)     :: minor
        integer(c_size_t)  :: textureAlignment
        integer(c_size_t)  :: texturePitchAlignment
        integer(c_int)     :: deviceOverlap
        integer(c_int)     :: multiProcessorCount
        character(c_char)  :: pad2(1656)             ! pad to 2048 bytes total
    end type cuda_device_prop

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

    integer(c_int) function cufftPlan3d(plan, nx, ny, nz, ffttype) bind(C, name='cufftPlan3d')
       import :: c_int
       integer(c_int) :: plan
       integer(c_int), value :: nx, ny, nz, ffttype
    end function cufftPlan3d

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

    integer(c_int) function cufftPlanMany(plan, rank, n, inembed, istride, idist, onembed,&
                            &ostride, odist, ffttype, batch) bind(C, name='cufftPlanMany')
       import :: c_int
       integer(c_int) :: plan
       integer(c_int), value :: rank
       integer(c_int), dimension(*), intent(in) :: n
       integer(c_int), dimension(*), intent(in) :: inembed
       integer(c_int), value :: istride, idist
       integer(c_int), dimension(*), intent(in) :: onembed
       integer(c_int), value :: ostride, odist
       integer(c_int), value :: ffttype, batch
    end function cufftPlanMany

    integer(c_int) function cufftDestroy(plan) bind(C, name='cufftDestroy')
       import :: c_int
       integer(c_int), value :: plan
    end function cufftDestroy

    integer(c_int) function cublasCreate(handle) bind(C, name='cublasCreate_v2')
       import :: c_int, c_ptr
       type(c_ptr) :: handle
    end function cublasCreate

    integer(c_int) function cublasDestroy(handle) bind(C, name='cublasDestroy_v2')
       import :: c_int, c_ptr
       type(c_ptr), value :: handle
    end function cublasDestroy

    integer(c_int) function cublasSgemm(handle, transa, transb, m, n, k, alpha, a, lda,&
                            &b, ldb, beta, c, ldc) bind(C, name='cublasSgemm_v2')
       import :: c_float, c_int, c_ptr
       type(c_ptr),    value :: handle
       integer(c_int), value :: transa, transb, m, n, k, lda, ldb, ldc
       real(c_float),  intent(in) :: alpha, beta
       type(c_ptr),    value :: a, b, c
    end function cublasSgemm

    integer(c_int) function cudaGetDeviceProperties(prop, device) bind(C, name='cudaGetDeviceProperties')
       import :: c_int, cuda_device_prop
       type(cuda_device_prop), intent(out) :: prop
       integer(c_int), value :: device
    end function cudaGetDeviceProperties

    integer(c_int) function cudaMemGetInfo(free_mem, total_mem) bind(C, name='cudaMemGetInfo')
       import :: c_int, c_size_t
       integer(c_size_t), intent(out) :: free_mem, total_mem
    end function cudaMemGetInfo
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
#ifdef USE_OPENMP_OFFLOAD
        ! Project built with GPU offloading
        if( cline%defined('device') )then
            ! device id provided
            id = cline%get_iarg('device')
            call set_offload_device_2( id )
        else
            ! use default device
            ndevices = omp_get_num_devices()
            if( ndevices < 1 ) then
                THROW_HARD('A device could not be set: No device identified')
                return
            endif
            id = omp_get_default_device()
            call set_offload_device_2( id )
        endif
        if( id < 0 ) then
            THROW_HARD('A device could not be set')
        endif
#else
        ! CPU only branch
        if( cline%defined('device') )then
            THROW_WARN('Ignoring DEVICE input. Build with USE_OPENMP_OFFLOAD=ON for device offloading')
        endif
#endif
    end subroutine set_offload_device_1

    subroutine set_offload_device_2( id )
        integer, intent(in) :: id
        integer :: ndevices, host_device
        integer(c_int) :: ierr
#ifdef USE_OPENMP_OFFLOAD
        ndevices = omp_get_num_devices()
        if( ndevices < 1 ) then
            THROW_WARN('Device could not be set: No device identified')
            return
        endif
        if( (id >= ndevices) .or. (id < 0) ) then
            THROW_WARN('Device could not be set: Invalid device ID provided')
            return
        endif
        host_device = omp_get_initial_device()
        if( host_device == id ) then
            THROW_WARN('Device could not be set: device ID provided is the host')
        endif
        call omp_set_default_device( id )
        ierr = cudaSetDevice(int(id, c_int))
        if( ierr /= 0_c_int ) then
            write(logfhandle, '(A,I0,A,I0)') &
                '>>> WARNING: cudaSetDevice failed for OpenMP device ID ', id, ', CUDA status = ', ierr
        endif
        write(logfhandle, '(A,I0)') '>>> OpenMP offload device ID: ', id
#endif
    end subroutine set_offload_device_2

    subroutine print_gpu_specs( id )
        integer, intent(in) :: id
#ifdef USE_OPENMP_OFFLOAD
        type(cuda_device_prop) :: prop
        integer(c_int)         :: ierr
        integer(c_size_t)      :: free_mem, total_mem
        character(len=256)     :: devname
        integer                :: i
        ierr = cudaGetDeviceProperties(prop, int(id, c_int))
        if( ierr /= CUDA_SUCCESS ) then
            write(logfhandle, '(A,I0,A,I0)') &
                '>>> WARNING: cudaGetDeviceProperties failed for device ', id, ', status = ', ierr
            return
        endif
        ierr = cudaMemGetInfo(free_mem, total_mem)
        if( ierr /= CUDA_SUCCESS ) then
            write(logfhandle, '(A,I0,A,I0)') &
                '>>> WARNING: cudaMemGetInfo failed for device ', id, ', status = ', ierr
            free_mem  = 0_c_size_t
            total_mem = prop%totalGlobalMem
        endif
        ! Parse null-terminated device name into a Fortran string
        devname = ' '
        do i = 1, 256
            if( prop%name(i) == c_null_char ) exit
            devname(i:i) = prop%name(i)
        end do
        write(logfhandle, '(/,A)')          '>>> GPU DEVICE PROPERTIES'
        write(logfhandle, '(A,I0)')         '    Device ID              : ', id
        write(logfhandle, '(A,A)')          '    Name                   : ', trim(devname)
        write(logfhandle, '(A,I0,A,I0)')   '    Compute capability     : ', prop%major, '.', prop%minor
        write(logfhandle, '(A,I0)')         '    Multiprocessors        : ', prop%multiProcessorCount
        write(logfhandle, '(A,F8.3,A)')    '    Total global memory    : ', &
            dble(prop%totalGlobalMem) / (1024.d0**3), ' GB'
        write(logfhandle, '(A,F8.3,A,F8.3,A)') '    Memory free / in use   : ', &
            dble(free_mem)            / (1024.d0**3), ' GB / ', &
            dble(total_mem - free_mem)/ (1024.d0**3), ' GB'
        write(logfhandle, '(A,I0,A)')      '    Shared mem / block     : ', &
            prop%sharedMemPerBlock / 1024, ' KB'
        write(logfhandle, '(A,I0,A)')      '    Constant memory        : ', &
            prop%totalConstMem     / 1024, ' KB'
        write(logfhandle, '(A,I0)')        '    Registers / block      : ', prop%regsPerBlock
        write(logfhandle, '(A,I0)')        '    Warp size              : ', prop%warpSize
        write(logfhandle, '(A,I0)')        '    Max threads / block    : ', prop%maxThreadsPerBlock
        write(logfhandle, '(A,3(I0,1X))') '    Max block dimensions   : ', prop%maxThreadsDim
        write(logfhandle, '(A,3(I0,1X))') '    Max grid dimensions    : ', prop%maxGridSize
        write(logfhandle, '(A,F8.3,A)')   '    Clock rate             : ', &
            real(prop%clockRate) / 1.0e6, ' GHz'
#else
        write(logfhandle, '(A)') '>>> GPU specs unavailable: build with USE_OPENMP_OFFLOAD=ON'
#endif
    end subroutine print_gpu_specs

end module simple_gpu_utils
