!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
  module cuda_h
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
    use, intrinsic :: ISO_C_BINDING
    use cuda_unknowns

    implicit none

    enum, bind(C) !:: CUctx_flags_enum
      enumerator :: CU_CTX_SCHED_AUTO=0
      enumerator :: CU_CTX_SCHED_SPIN=1
      enumerator :: CU_CTX_SCHED_YIELD=2
      enumerator :: CU_CTX_SCHED_MASK=3
      enumerator :: CU_CTX_BLOCKING_SYNC=4
      enumerator :: CU_CTX_MAP_HOST=8
      enumerator :: CU_CTX_LMEM_RESIZE_TO_MAX=16
      enumerator :: CU_CTX_FLAGS_MASK=31
    end enum ! CUctx_flags_enum

    enum, bind(C) !:: CUevent_flags_enum
      enumerator :: CU_EVENT_DEFAULT=0
      enumerator :: CU_EVENT_BLOCKING_SYNC=1
    end enum ! CUevent_flags_enum

    enum, bind(C) !:: CUarray_format_enum
      enumerator :: CU_AD_FORMAT_UNSIGNED_INT8=1
      enumerator :: CU_AD_FORMAT_UNSIGNED_INT16=2
      enumerator :: CU_AD_FORMAT_UNSIGNED_INT32=3
      enumerator :: CU_AD_FORMAT_SIGNED_INT8=8
      enumerator :: CU_AD_FORMAT_SIGNED_INT16=9
      enumerator :: CU_AD_FORMAT_SIGNED_INT32=10
      enumerator :: CU_AD_FORMAT_HALF=16
      enumerator :: CU_AD_FORMAT_FLOAT=32
    end enum ! CUarray_format_enum

    enum, bind(C) !:: CUaddress_mode_enum
      enumerator :: CU_TR_ADDRESS_MODE_WRAP=0
      enumerator :: CU_TR_ADDRESS_MODE_CLAMP=1
      enumerator :: CU_TR_ADDRESS_MODE_MIRROR=2
    end enum ! CUaddress_mode_enum

    enum, bind(C) !:: CUfilter_mode_enum
      enumerator :: CU_TR_FILTER_MODE_POINT=0
      enumerator :: CU_TR_FILTER_MODE_LINEAR=1
    end enum ! CUfilter_mode_enum

    enum, bind(C) !:: CUdevice_attribute_enum
      enumerator :: CU_DEVICE_ATTRIBUTE_MAX_THREADS_PER_BLOCK=1
      enumerator :: CU_DEVICE_ATTRIBUTE_MAX_BLOCK_DIM_X=2
      enumerator :: CU_DEVICE_ATTRIBUTE_MAX_BLOCK_DIM_Y=3
      enumerator :: CU_DEVICE_ATTRIBUTE_MAX_BLOCK_DIM_Z=4
      enumerator :: CU_DEVICE_ATTRIBUTE_MAX_GRID_DIM_X=5
      enumerator :: CU_DEVICE_ATTRIBUTE_MAX_GRID_DIM_Y=6
      enumerator :: CU_DEVICE_ATTRIBUTE_MAX_GRID_DIM_Z=7
      enumerator :: CU_DEVICE_ATTRIBUTE_MAX_SHARED_MEMORY_PER_BLOCK=8
      enumerator :: CU_DEVICE_ATTRIBUTE_SHARED_MEMORY_PER_BLOCK=8
      enumerator :: CU_DEVICE_ATTRIBUTE_TOTAL_CONSTANT_MEMORY=9
      enumerator :: CU_DEVICE_ATTRIBUTE_WARP_SIZE=10
      enumerator :: CU_DEVICE_ATTRIBUTE_MAX_PITCH=11
      enumerator :: CU_DEVICE_ATTRIBUTE_MAX_REGISTERS_PER_BLOCK=12
      enumerator :: CU_DEVICE_ATTRIBUTE_REGISTERS_PER_BLOCK=12
      enumerator :: CU_DEVICE_ATTRIBUTE_CLOCK_RATE=13
      enumerator :: CU_DEVICE_ATTRIBUTE_TEXTURE_ALIGNMENT=14
      enumerator :: CU_DEVICE_ATTRIBUTE_GPU_OVERLAP=15
      enumerator :: CU_DEVICE_ATTRIBUTE_MULTIPROCESSOR_COUNT=16
      enumerator :: CU_DEVICE_ATTRIBUTE_KERNEL_EXEC_TIMEOUT=17
      enumerator :: CU_DEVICE_ATTRIBUTE_INTEGRATED=18
      enumerator :: CU_DEVICE_ATTRIBUTE_CAN_MAP_HOST_MEMORY=19
      enumerator :: CU_DEVICE_ATTRIBUTE_COMPUTE_MODE=20
      enumerator :: CU_DEVICE_ATTRIBUTE_MAXIMUM_TEXTURE1D_WIDTH=21
      enumerator :: CU_DEVICE_ATTRIBUTE_MAXIMUM_TEXTURE2D_WIDTH=22
      enumerator :: CU_DEVICE_ATTRIBUTE_MAXIMUM_TEXTURE2D_HEIGHT=23
      enumerator :: CU_DEVICE_ATTRIBUTE_MAXIMUM_TEXTURE3D_WIDTH=24
      enumerator :: CU_DEVICE_ATTRIBUTE_MAXIMUM_TEXTURE3D_HEIGHT=25
      enumerator :: CU_DEVICE_ATTRIBUTE_MAXIMUM_TEXTURE3D_DEPTH=26
      enumerator :: CU_DEVICE_ATTRIBUTE_MAXIMUM_TEXTURE2D_ARRAY_WIDTH=27
      enumerator :: CU_DEVICE_ATTRIBUTE_MAXIMUM_TEXTURE2D_ARRAY_HEIGHT=28
      enumerator :: CU_DEVICE_ATTRIBUTE_MAXIMUM_TEXTURE2D_ARRAY_NUMSLICES=29
      enumerator :: CU_DEVICE_ATTRIBUTE_SURFACE_ALIGNMENT=30
      enumerator :: CU_DEVICE_ATTRIBUTE_CONCURRENT_KERNELS=31
      enumerator :: CU_DEVICE_ATTRIBUTE_ECC_ENABLED=32
      enumerator :: CU_DEVICE_ATTRIBUTE_PCI_BUS_ID=33
      enumerator :: CU_DEVICE_ATTRIBUTE_PCI_DEVICE_ID=34
    end enum ! CUdevice_attribute_enum

    enum, bind(C) !:: CUfunction_attribute_enum
      enumerator :: CU_FUNC_ATTRIBUTE_MAX_THREADS_PER_BLOCK=0
      enumerator :: CU_FUNC_ATTRIBUTE_SHARED_SIZE_BYTES=1
      enumerator :: CU_FUNC_ATTRIBUTE_CONST_SIZE_BYTES=2
      enumerator :: CU_FUNC_ATTRIBUTE_LOCAL_SIZE_BYTES=3
      enumerator :: CU_FUNC_ATTRIBUTE_NUM_REGS=4
      enumerator :: CU_FUNC_ATTRIBUTE_PTX_VERSION=5
      enumerator :: CU_FUNC_ATTRIBUTE_BINARY_VERSION=6
      enumerator :: CU_FUNC_ATTRIBUTE_MAX
    end enum ! CUfunction_attribute_enum

    enum, bind(C) !:: CUfunc_cache_enum
      enumerator :: CU_FUNC_CACHE_PREFER_NONE=0
      enumerator :: CU_FUNC_CACHE_PREFER_SHARED=1
      enumerator :: CU_FUNC_CACHE_PREFER_L1=2
    end enum ! CUfunc_cache_enum

    enum, bind(C) !:: CUmemorytype_enum
      enumerator :: CU_MEMORYTYPE_HOST=1
      enumerator :: CU_MEMORYTYPE_DEVICE=2
      enumerator :: CU_MEMORYTYPE_ARRAY=3
    end enum ! CUmemorytype_enum

    enum, bind(C) !:: CUcomputemode_enum
      enumerator :: CU_COMPUTEMODE_DEFAULT=0
      enumerator :: CU_COMPUTEMODE_EXCLUSIVE=1
      enumerator :: CU_COMPUTEMODE_PROHIBITED=2
    end enum ! CUcomputemode_enum

    enum, bind(C) !:: CUjit_option_enum
      enumerator :: CU_JIT_MAX_REGISTERS=0
      enumerator :: CU_JIT_THREADS_PER_BLOCK
      enumerator :: CU_JIT_WALL_TIME
      enumerator :: CU_JIT_INFO_LOG_BUFFER
      enumerator :: CU_JIT_INFO_LOG_BUFFER_SIZE_BYTES
      enumerator :: CU_JIT_ERROR_LOG_BUFFER
      enumerator :: CU_JIT_ERROR_LOG_BUFFER_SIZE_BYTES
      enumerator :: CU_JIT_OPTIMIZATION_LEVEL
      enumerator :: CU_JIT_TARGET_FROM_CUCONTEXT
      enumerator :: CU_JIT_TARGET
      enumerator :: CU_JIT_FALLBACK_STRATEGY
    end enum ! CUjit_option_enum

    enum, bind(C) !:: CUjit_target_enum
      enumerator :: CU_TARGET_COMPUTE_10=0
      enumerator :: CU_TARGET_COMPUTE_11
      enumerator :: CU_TARGET_COMPUTE_12
      enumerator :: CU_TARGET_COMPUTE_13
      enumerator :: CU_TARGET_COMPUTE_20
    end enum ! CUjit_target_enum

    enum, bind(C) !:: CUjit_fallback_enum
      enumerator :: CU_PREFER_PTX=0
      enumerator :: CU_PREFER_BINARY
    end enum ! CUjit_fallback_enum

    enum, bind(C) !:: CUgraphicsRegisterFlags_enum
      enumerator :: CU_GRAPHICS_REGISTER_FLAGS_NONE=0
    end enum ! CUgraphicsRegisterFlags_enum

    enum, bind(C) !:: CUgraphicsMapResourceFlags_enum
      enumerator :: CU_GRAPHICS_MAP_RESOURCE_FLAGS_NONE=0
      enumerator :: CU_GRAPHICS_MAP_RESOURCE_FLAGS_READ_ONLY=1
      enumerator :: CU_GRAPHICS_MAP_RESOURCE_FLAGS_WRITE_DISCARD=2
    end enum ! CUgraphicsMapResourceFlags_enum

    enum, bind(C) !:: CUarray_cubemap_face_enum
      enumerator :: CU_CUBEMAP_FACE_POSITIVE_X=0
      enumerator :: CU_CUBEMAP_FACE_NEGATIVE_X=1
      enumerator :: CU_CUBEMAP_FACE_POSITIVE_Y=2
      enumerator :: CU_CUBEMAP_FACE_NEGATIVE_Y=3
      enumerator :: CU_CUBEMAP_FACE_POSITIVE_Z=4
      enumerator :: CU_CUBEMAP_FACE_NEGATIVE_Z=5
    end enum ! CUarray_cubemap_face_enum

    enum, bind(C) !:: CUlimit_enum
      enumerator :: CU_LIMIT_STACK_SIZE=0
      enumerator :: CU_LIMIT_PRINTF_FIFO_SIZE=1
    end enum ! CUlimit_enum

    enum, bind(C) !:: cudaError_enum
      enumerator :: CUDA_SUCCESS=0
      enumerator :: CUDA_ERROR_INVALID_VALUE=1
      enumerator :: CUDA_ERROR_OUT_OF_MEMORY=2
      enumerator :: CUDA_ERROR_NOT_INITIALIZED=3
      enumerator :: CUDA_ERROR_DEINITIALIZED=4
      enumerator :: CUDA_ERROR_NO_DEVICE=100
      enumerator :: CUDA_ERROR_INVALID_DEVICE=101
      enumerator :: CUDA_ERROR_INVALID_IMAGE=200
      enumerator :: CUDA_ERROR_INVALID_CONTEXT=201
      enumerator :: CUDA_ERROR_CONTEXT_ALREADY_CURRENT=202
      enumerator :: CUDA_ERROR_MAP_FAILED=205
      enumerator :: CUDA_ERROR_UNMAP_FAILED=206
      enumerator :: CUDA_ERROR_ARRAY_IS_MAPPED=207
      enumerator :: CUDA_ERROR_ALREADY_MAPPED=208
      enumerator :: CUDA_ERROR_NO_BINARY_FOR_GPU=209
      enumerator :: CUDA_ERROR_ALREADY_ACQUIRED=210
      enumerator :: CUDA_ERROR_NOT_MAPPED=211
      enumerator :: CUDA_ERROR_NOT_MAPPED_AS_ARRAY=212
      enumerator :: CUDA_ERROR_NOT_MAPPED_AS_POINTER=213
      enumerator :: CUDA_ERROR_ECC_UNCORRECTABLE=214
      enumerator :: CUDA_ERROR_UNSUPPORTED_LIMIT=215
      enumerator :: CUDA_ERROR_INVALID_SOURCE=300
      enumerator :: CUDA_ERROR_FILE_NOT_FOUND=301
      enumerator :: CUDA_ERROR_SHARED_OBJECT_SYMBOL_NOT_FOUND=302
      enumerator :: CUDA_ERROR_SHARED_OBJECT_INIT_FAILED=303
      enumerator :: CUDA_ERROR_INVALID_HANDLE=400
      enumerator :: CUDA_ERROR_NOT_FOUND=500
      enumerator :: CUDA_ERROR_NOT_READY=600
      enumerator :: CUDA_ERROR_LAUNCH_FAILED=700
      enumerator :: CUDA_ERROR_LAUNCH_OUT_OF_RESOURCES=701
      enumerator :: CUDA_ERROR_LAUNCH_TIMEOUT=702
      enumerator :: CUDA_ERROR_LAUNCH_INCOMPATIBLE_TEXTURING=703
      enumerator :: CUDA_ERROR_POINTER_IS_64BIT=800
      enumerator :: CUDA_ERROR_SIZE_IS_64BIT=801
      enumerator :: CUDA_ERROR_UNKNOWN=999
    end enum ! cudaError_enum

    type, bind(C) :: timespec
      integer(c_long) :: tv_sec
      integer(c_long) :: tv_nsec
    end type timespec

    type, bind(C) :: timeval
      integer(c_long) :: tv_sec
      integer(c_long) :: tv_usec
    end type timeval

    type, bind(C) :: pthread_internal_slist__
      type (c_ptr) :: next__
    end type pthread_internal_slist__

    type, bind(C) :: drand48_data
      integer(c_short) :: x__(3)
      integer(c_short) :: old_x__(3)
      integer(c_short) :: c__
      integer(c_short) :: init__
      integer(c_long_long) :: a__
    end type drand48_data

    type, bind(C) :: CUuuid_st
      character(c_char) :: bytes(16)
    end type CUuuid_st

    type, bind(C) :: CUdevprop_st
      integer(c_int) :: maxThreadsPerBlock
      integer(c_int) :: maxThreadsDim(3)
      integer(c_int) :: maxGridSize(3)
      integer(c_int) :: sharedMemPerBlock
      integer(c_int) :: totalConstantMemory
      integer(c_int) :: SIMDWidth
      integer(c_int) :: memPitch
      integer(c_int) :: regsPerBlock
      integer(c_int) :: clockRate
      integer(c_int) :: textureAlign
    end type CUdevprop_st

    type, bind(C) :: CUDA_MEMCPY2D_st
      integer(c_int) :: srcXInBytes
      integer(c_int) :: srcY
      integer (KIND(CU_MEMORYTYPE_HOST)) :: srcMemoryType
      type(c_ptr) :: srcHost
      integer(c_int) :: srcDevice
      type (cudaArray_t) :: srcArray
      integer(c_int) :: srcPitch
      integer(c_int) :: dstXInBytes
      integer(c_int) :: dstY
      integer (KIND(CU_MEMORYTYPE_HOST)) :: dstMemoryType
      type(c_ptr) :: dstHost
      integer(c_int) :: dstDevice
      type (cudaArray_t) :: dstArray
      integer(c_int) :: dstPitch
      integer(c_int) :: WidthInBytes
      integer(c_int) :: Height
    end type CUDA_MEMCPY2D_st

    type, bind(C) :: CUDA_MEMCPY3D_st
      integer(c_int) :: srcXInBytes
      integer(c_int) :: srcY
      integer(c_int) :: srcZ
      integer(c_int) :: srcLOD
      integer (KIND(CU_MEMORYTYPE_HOST)) :: srcMemoryType
      type(c_ptr) :: srcHost
      integer(c_int) :: srcDevice
      type (cudaArray_t) :: srcArray
      type(c_ptr) :: reserved0
      integer(c_int) :: srcPitch
      integer(c_int) :: srcHeight
      integer(c_int) :: dstXInBytes
      integer(c_int) :: dstY
      integer(c_int) :: dstZ
      integer(c_int) :: dstLOD
      integer (KIND(CU_MEMORYTYPE_HOST)) :: dstMemoryType
      type(c_ptr) :: dstHost
      integer(c_int) :: dstDevice
      type (cudaArray_t) :: dstArray
      type(c_ptr) :: reserved1
      integer(c_int) :: dstPitch
      integer(c_int) :: dstHeight
      integer(c_int) :: WidthInBytes
      integer(c_int) :: Height
      integer(c_int) :: Depth
    end type CUDA_MEMCPY3D_st

    interface ! [['int', None], 'posix_openpt', [['int', None, '__oflag']]]
      function posix_openpt(oflag__) result( res ) bind(C, name="posix_openpt")
        use, intrinsic :: ISO_C_BINDING
        implicit none
        integer(c_int), value :: oflag__
        integer(c_int) :: res
      end function posix_openpt
    end interface

    interface ! [['CUresult', None], 'cuInit', [['unsigned', 'int', None, 'Flags']]]
      function cuInit(Flags) result( res ) bind(C, name="cuInit")
        use, intrinsic :: ISO_C_BINDING
        import CUDA_SUCCESS
        implicit none
        integer(c_int), value :: Flags
        integer (KIND(CUDA_SUCCESS)) :: res
      end function cuInit
    end interface

    interface ! [['CUresult', None], 'cuDriverGetVersion', [['int', '*', 'driverVersion']]]
      function cuDriverGetVersion(driverVersion) result( res ) bind(C, name="cuDriverGetVersion")
        use, intrinsic :: ISO_C_BINDING
        import CUDA_SUCCESS
        implicit none
        integer(c_int) :: driverVersion
        integer (KIND(CUDA_SUCCESS)) :: res
      end function cuDriverGetVersion
    end interface

    interface ! [['CUresult', None], 'cuDeviceGet', [['CUdevice', '*', 'device'], ['int', None, 'ordinal']]]
      function cuDeviceGet(device,ordinal) result( res ) bind(C, name="cuDeviceGet")
        use, intrinsic :: ISO_C_BINDING
        import CUDA_SUCCESS
        implicit none
        integer(c_int) :: device
        integer(c_int), value :: ordinal
        integer (KIND(CUDA_SUCCESS)) :: res
      end function cuDeviceGet
    end interface

    interface ! [['CUresult', None], 'cuCtxCreate', [['CUcontext', '*', 'pctx'], ['unsigned', 'int', None, 'flags'], ['CUdevice', None, 'dev']]]
      function cuCtxCreate(pctx,flags,dev) result( res ) bind(C, name="cuCtxCreate")
        use, intrinsic :: ISO_C_BINDING
        import CUDA_SUCCESS
        import cudaContext_t
        implicit none
        type (cudaContext_t) :: pctx
        integer(c_int), value :: flags
        integer(c_int), value :: dev
        integer (KIND(CUDA_SUCCESS)) :: res
      end function cuCtxCreate
    end interface

    interface ! [['CUresult', None], 'cuModuleLoad', [['CUmodule', '*', 'module'], ['const', 'char', '*', 'fname']]]
      function cuModuleLoad(module,fname) result( res ) bind(C, name="cuModuleLoad")
        use, intrinsic :: ISO_C_BINDING
        import CUDA_SUCCESS
        import cudaModule_t
        implicit none
        type (cudaModule_t) :: module
        character(c_char) :: fname
        integer (KIND(CUDA_SUCCESS)) :: res
      end function cuModuleLoad
    end interface

    interface ! [['CUresult', None], 'cuMemGetInfo', [['unsigned', 'int', '*', 'free'], ['unsigned', 'int', '*', 'total']]]
      function cuMemGetInfo(free,total) result( res ) bind(C, name="cuMemGetInfo")
        use, intrinsic :: ISO_C_BINDING
        import CUDA_SUCCESS
        implicit none
        integer(c_int) :: free
        integer(c_int) :: total
        integer (KIND(CUDA_SUCCESS)) :: res
      end function cuMemGetInfo
    end interface

    interface ! [['CUresult', None], 'cuMemAlloc', [['CUdeviceptr', '*', 'dptr'], ['unsigned', 'int', None, 'bytesize']]]
      function cuMemAlloc(dptr,bytesize) result( res ) bind(C, name="cuMemAlloc")
        use, intrinsic :: ISO_C_BINDING
        import CUDA_SUCCESS
        implicit none
        integer(c_int) :: dptr
        integer(c_int), value :: bytesize
        integer (KIND(CUDA_SUCCESS)) :: res
      end function cuMemAlloc
    end interface

    interface ! [['CUresult', None], 'cuMemAllocHost', [['void', '**', 'pp'], ['unsigned', 'int', None, 'bytesize']]]
      function cuMemAllocHost(pp,bytesize) result( res ) bind(C, name="cuMemAllocHost")
        use, intrinsic :: ISO_C_BINDING
        import CUDA_SUCCESS
        implicit none
        type(c_ptr) :: pp
        integer(c_int), value :: bytesize
        integer (KIND(CUDA_SUCCESS)) :: res
      end function cuMemAllocHost
    end interface

    interface ! [['CUresult', None], 'cuMemHostAlloc', [['void', '**', 'pp'], ['size_t', None, 'bytesize'], ['unsigned', 'int', None, 'Flags']]]
      function cuMemHostAlloc(pp,bytesize,Flags) result( res ) bind(C, name="cuMemHostAlloc")
        use, intrinsic :: ISO_C_BINDING
        import CUDA_SUCCESS
        implicit none
        type(c_ptr) :: pp
        integer(c_int), value :: bytesize
        integer(c_int), value :: Flags
        integer (KIND(CUDA_SUCCESS)) :: res
      end function cuMemHostAlloc
    end interface

    interface ! [['CUresult', None], 'cuMemHostGetDevicePointer', [['CUdeviceptr', '*', 'pdptr'], ['void', '*', 'p'], ['unsigned', 'int', None, 'Flags']]]
      function cuMemHostGetDevicePointer(pdptr,p,Flags) result( res ) bind(C, name="cuMemHostGetDevicePointer")
        use, intrinsic :: ISO_C_BINDING
        import CUDA_SUCCESS
        implicit none
        integer(c_int) :: pdptr
        type(c_ptr), value :: p
        integer(c_int), value :: Flags
        integer (KIND(CUDA_SUCCESS)) :: res
      end function cuMemHostGetDevicePointer
    end interface

    interface ! [['CUresult', None], 'cuMemcpyHtoD', [['CUdeviceptr', None, 'dstDevice'], ['const', 'void', '*', 'srcHost'], ['unsigned', 'int', None, 'ByteCount']]]
      function cuMemcpyHtoD(dstDevice,srcHost,ByteCount) result( res ) bind(C, name="cuMemcpyHtoD")
        use, intrinsic :: ISO_C_BINDING
        import CUDA_SUCCESS
        implicit none
        integer(c_int), value :: dstDevice
        type(c_ptr), value :: srcHost
        integer(c_int), value :: ByteCount
        integer (KIND(CUDA_SUCCESS)) :: res
      end function cuMemcpyHtoD
    end interface

    interface ! [['CUresult', None], 'cuMemcpyDtoD', [['CUdeviceptr', None, 'dstDevice'], ['CUdeviceptr', None, 'srcDevice'], ['unsigned', 'int', None, 'ByteCount']]]
      function cuMemcpyDtoD(dstDevice,srcDevice,ByteCount) result( res ) bind(C, name="cuMemcpyDtoD")
        use, intrinsic :: ISO_C_BINDING
        import CUDA_SUCCESS
        implicit none
        integer(c_int), value :: dstDevice
        integer(c_int), value :: srcDevice
        integer(c_int), value :: ByteCount
        integer (KIND(CUDA_SUCCESS)) :: res
      end function cuMemcpyDtoD
    end interface

    interface ! [['CUresult', None], 'cuMemcpyDtoA', [['CUarray', None, 'dstArray'], ['unsigned', 'int', None, 'dstOffset'], ['CUdeviceptr', None, 'srcDevice'], ['unsigned', 'int', None, 'ByteCount']]]
      function cuMemcpyDtoA(dstArray,dstOffset,srcDevice,ByteCount) result( res ) bind(C, name="cuMemcpyDtoA")
        use, intrinsic :: ISO_C_BINDING
        import CUDA_SUCCESS
        import cudaArray_t
        implicit none
        type (cudaArray_t), value :: dstArray
        integer(c_int), value :: dstOffset
        integer(c_int), value :: srcDevice
        integer(c_int), value :: ByteCount
        integer (KIND(CUDA_SUCCESS)) :: res
      end function cuMemcpyDtoA
    end interface

    interface ! [['CUresult', None], 'cuMemcpyHtoA', [['CUarray', None, 'dstArray'], ['unsigned', 'int', None, 'dstOffset'], ['const', 'void', '*', 'srcHost'], ['unsigned', 'int', None, 'ByteCount']]]
      function cuMemcpyHtoA(dstArray,dstOffset,srcHost,ByteCount) result( res ) bind(C, name="cuMemcpyHtoA")
        use, intrinsic :: ISO_C_BINDING
        import CUDA_SUCCESS
        import cudaArray_t
        implicit none
        type (cudaArray_t), value :: dstArray
        integer(c_int), value :: dstOffset
        type(c_ptr), value :: srcHost
        integer(c_int), value :: ByteCount
        integer (KIND(CUDA_SUCCESS)) :: res
      end function cuMemcpyHtoA
    end interface

    interface ! [['CUresult', None], 'cuMemcpyAtoA', [['CUarray', None, 'dstArray'], ['unsigned', 'int', None, 'dstOffset'], ['CUarray', None, 'srcArray'], ['unsigned', 'int', None, 'srcOffset'], ['unsigned', 'int', None, 'ByteCount']]]
      function cuMemcpyAtoA(dstArray,dstOffset,srcArray,srcOffset,ByteCount) result( res ) bind(C, name="cuMemcpyAtoA")
        use, intrinsic :: ISO_C_BINDING
        import CUDA_SUCCESS
        import cudaArray_t
        import cudaArray_t
        implicit none
        type (cudaArray_t), value :: dstArray
        integer(c_int), value :: dstOffset
        type (cudaArray_t), value :: srcArray
        integer(c_int), value :: srcOffset
        integer(c_int), value :: ByteCount
        integer (KIND(CUDA_SUCCESS)) :: res
      end function cuMemcpyAtoA
    end interface

    interface ! [['CUresult', None], 'cuMemcpyHtoDAsync', [['CUdeviceptr', None, 'dstDevice'], ['const', 'void', '*', 'srcHost'], ['unsigned', 'int', None, 'ByteCount'], ['CUstream', None, 'hStream']]]
      function cuMemcpyHtoDAsync(dstDevice,srcHost,ByteCount,hStream) result( res ) bind(C, name="cuMemcpyHtoDAsync")
        use, intrinsic :: ISO_C_BINDING
        import CUDA_SUCCESS
        import cudaStream_t
        implicit none
        integer(c_int), value :: dstDevice
        type(c_ptr), value :: srcHost
        integer(c_int), value :: ByteCount
        type (cudaStream_t), value :: hStream
        integer (KIND(CUDA_SUCCESS)) :: res
      end function cuMemcpyHtoDAsync
    end interface

    interface ! [['CUresult', None], 'cuMemcpyDtoDAsync', [['CUdeviceptr', None, 'dstDevice'], ['CUdeviceptr', None, 'srcDevice'], ['unsigned', 'int', None, 'ByteCount'], ['CUstream', None, 'hStream']]]
      function cuMemcpyDtoDAsync(dstDevice,srcDevice,ByteCount,hStream) result( res ) bind(C, name="cuMemcpyDtoDAsync")
        use, intrinsic :: ISO_C_BINDING
        import CUDA_SUCCESS
        import cudaStream_t
        implicit none
        integer(c_int), value :: dstDevice
        integer(c_int), value :: srcDevice
        integer(c_int), value :: ByteCount
        type (cudaStream_t), value :: hStream
        integer (KIND(CUDA_SUCCESS)) :: res
      end function cuMemcpyDtoDAsync
    end interface

    interface ! [['CUresult', None], 'cuMemcpyHtoAAsync', [['CUarray', None, 'dstArray'], ['unsigned', 'int', None, 'dstOffset'], ['const', 'void', '*', 'srcHost'], ['unsigned', 'int', None, 'ByteCount'], ['CUstream', None, 'hStream']]]
      function cuMemcpyHtoAAsync(dstArray,dstOffset,srcHost,ByteCount,hStream) result( res ) bind(C, name="cuMemcpyHtoAAsync")
        use, intrinsic :: ISO_C_BINDING
        import CUDA_SUCCESS
        import cudaArray_t
        import cudaStream_t
        implicit none
        type (cudaArray_t), value :: dstArray
        integer(c_int), value :: dstOffset
        type(c_ptr), value :: srcHost
        integer(c_int), value :: ByteCount
        type (cudaStream_t), value :: hStream
        integer (KIND(CUDA_SUCCESS)) :: res
      end function cuMemcpyHtoAAsync
    end interface

    interface ! [['CUresult', None], 'cuMemsetD8', [['CUdeviceptr', None, 'dstDevice'], ['unsigned', 'char', None, 'uc'], ['unsigned', 'int', None, 'N']]]
      function cuMemsetD8(dstDevice,uc,N) result( res ) bind(C, name="cuMemsetD8")
        use, intrinsic :: ISO_C_BINDING
        import CUDA_SUCCESS
        implicit none
        integer(c_int), value :: dstDevice
        integer(c_signed_char), value :: uc
        integer(c_int), value :: N
        integer (KIND(CUDA_SUCCESS)) :: res
      end function cuMemsetD8
    end interface

    interface ! [['CUresult', None], 'cuMemsetD2D8', [['CUdeviceptr', None, 'dstDevice'], ['unsigned', 'int', None, 'dstPitch'], ['unsigned', 'char', None, 'uc'], ['unsigned', 'int', None, 'Width'], ['unsigned', 'int', None, 'Height']]]
      function cuMemsetD2D8(dstDevice,dstPitch,uc,Width,Height) result( res ) bind(C, name="cuMemsetD2D8")
        use, intrinsic :: ISO_C_BINDING
        import CUDA_SUCCESS
        implicit none
        integer(c_int), value :: dstDevice
        integer(c_int), value :: dstPitch
        integer(c_signed_char), value :: uc
        integer(c_int), value :: Width
        integer(c_int), value :: Height
        integer (KIND(CUDA_SUCCESS)) :: res
      end function cuMemsetD2D8
    end interface

    interface ! [['CUresult', None], 'cuFuncSetBlockShape', [['CUfunction', None, 'hfunc'], ['int', None, 'x'], ['int', None, 'y'], ['int', None, 'z']]]
      function cuFuncSetBlockShape(hfunc,x,y,z) result( res ) bind(C, name="cuFuncSetBlockShape")
        use, intrinsic :: ISO_C_BINDING
        import CUDA_SUCCESS
        import cudaFunction_t
        implicit none
        type (cudaFunction_t), value :: hfunc
        integer(c_int), value :: x
        integer(c_int), value :: y
        integer(c_int), value :: z
        integer (KIND(CUDA_SUCCESS)) :: res
      end function cuFuncSetBlockShape
    end interface

    interface ! [['CUresult', None], 'cuTexRefCreate', [['CUtexref', '*', 'pTexRef']]]
      function cuTexRefCreate(pTexRef) result( res ) bind(C, name="cuTexRefCreate")
        use, intrinsic :: ISO_C_BINDING
        import CUDA_SUCCESS
        import cudaTexRef_t
        implicit none
        type (cudaTexRef_t) :: pTexRef
        integer (KIND(CUDA_SUCCESS)) :: res
      end function cuTexRefCreate
    end interface

    interface ! [['CUresult', None], 'cuTexRefSetArray', [['CUtexref', None, 'hTexRef'], ['CUarray', None, 'hArray'], ['unsigned', 'int', None, 'Flags']]]
      function cuTexRefSetArray(hTexRef,hArray,Flags) result( res ) bind(C, name="cuTexRefSetArray")
        use, intrinsic :: ISO_C_BINDING
        import CUDA_SUCCESS
        import cudaTexRef_t
        import cudaArray_t
        implicit none
        type (cudaTexRef_t), value :: hTexRef
        type (cudaArray_t), value :: hArray
        integer(c_int), value :: Flags
        integer (KIND(CUDA_SUCCESS)) :: res
      end function cuTexRefSetArray
    end interface

    interface ! [['CUresult', None], 'cuTexRefGetAddress', [['CUdeviceptr', '*', 'pdptr'], ['CUtexref', None, 'hTexRef']]]
      function cuTexRefGetAddress(pdptr,hTexRef) result( res ) bind(C, name="cuTexRefGetAddress")
        use, intrinsic :: ISO_C_BINDING
        import CUDA_SUCCESS
        import cudaTexRef_t
        implicit none
        integer(c_int) :: pdptr
        type (cudaTexRef_t), value :: hTexRef
        integer (KIND(CUDA_SUCCESS)) :: res
      end function cuTexRefGetAddress
    end interface

    interface ! [['CUresult', None], 'cuSurfRefSetArray', [['CUsurfref', None, 'hSurfRef'], ['CUarray', None, 'hArray'], ['unsigned', 'int', None, 'Flags']]]
      function cuSurfRefSetArray(hSurfRef,hArray,Flags) result( res ) bind(C, name="cuSurfRefSetArray")
        use, intrinsic :: ISO_C_BINDING
        import CUDA_SUCCESS
        import cudaSurfRef_t
        import cudaArray_t
        implicit none
        type (cudaSurfRef_t), value :: hSurfRef
        type (cudaArray_t), value :: hArray
        integer(c_int), value :: Flags
        integer (KIND(CUDA_SUCCESS)) :: res
      end function cuSurfRefSetArray
    end interface

    interface ! [['CUresult', None], 'cuParamSetSize', [['CUfunction', None, 'hfunc'], ['unsigned', 'int', None, 'numbytes']]]
      function cuParamSetSize(hfunc,numbytes) result( res ) bind(C, name="cuParamSetSize")
        use, intrinsic :: ISO_C_BINDING
        import CUDA_SUCCESS
        import cudaFunction_t
        implicit none
        type (cudaFunction_t), value :: hfunc
        integer(c_int), value :: numbytes
        integer (KIND(CUDA_SUCCESS)) :: res
      end function cuParamSetSize
    end interface

    interface ! [['CUresult', None], 'cuLaunch', [['CUfunction', None, 'f']]]
      function cuLaunch(f) result( res ) bind(C, name="cuLaunch")
        use, intrinsic :: ISO_C_BINDING
        import CUDA_SUCCESS
        import cudaFunction_t
        implicit none
        type (cudaFunction_t), value :: f
        integer (KIND(CUDA_SUCCESS)) :: res
      end function cuLaunch
    end interface

    interface ! [['CUresult', None], 'cuEventCreate', [['CUevent', '*', 'phEvent'], ['unsigned', 'int', None, 'Flags']]]
      function cuEventCreate(phEvent,Flags) result( res ) bind(C, name="cuEventCreate")
        use, intrinsic :: ISO_C_BINDING
        import CUDA_SUCCESS
        import cudaEvent_t
        implicit none
        type (cudaEvent_t) :: phEvent
        integer(c_int), value :: Flags
        integer (KIND(CUDA_SUCCESS)) :: res
      end function cuEventCreate
    end interface

    interface ! [['CUresult', None], 'cuStreamCreate', [['CUstream', '*', 'phStream'], ['unsigned', 'int', None, 'Flags']]]
      function cuStreamCreate(phStream,Flags) result( res ) bind(C, name="cuStreamCreate")
        use, intrinsic :: ISO_C_BINDING
        import CUDA_SUCCESS
        import cudaStream_t
        implicit none
        type (cudaStream_t) :: phStream
        integer(c_int), value :: Flags
        integer (KIND(CUDA_SUCCESS)) :: res
      end function cuStreamCreate
    end interface

    interface ! [['CUresult', None], 'cuGraphicsUnregisterResource', [['CUgraphicsResource', None, 'resource']]]
      function cuGraphicsUnregisterResource(resource) result( res ) bind(C, name="cuGraphicsUnregisterResource")
        use, intrinsic :: ISO_C_BINDING
        import CUDA_SUCCESS
        import cudaGraphicsResource_t
        implicit none
        type (cudaGraphicsResource_t), value :: resource
        integer (KIND(CUDA_SUCCESS)) :: res
      end function cuGraphicsUnregisterResource
    end interface

    interface ! [['CUresult', None], 'cuCtxSetLimit', [['CUlimit', None, 'limit'], ['size_t', None, 'value']]]
      function cuCtxSetLimit(limit,value) result( res ) bind(C, name="cuCtxSetLimit")
        use, intrinsic :: ISO_C_BINDING
        import CUDA_SUCCESS
        import CU_LIMIT_STACK_SIZE
        implicit none
        integer (KIND(CU_LIMIT_STACK_SIZE)), value :: limit
        integer(c_int), value :: value
        integer (KIND(CUDA_SUCCESS)) :: res
      end function cuCtxSetLimit
    end interface

  end module cuda_h
