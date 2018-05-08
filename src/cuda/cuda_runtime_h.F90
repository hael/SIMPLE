!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
  module cuda_runtime_h
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
    use, intrinsic :: ISO_C_BINDING
    use cuda_unknowns

    implicit none

    enum, bind(C) !:: cudaRoundMode
      enumerator :: cudaRoundNearest
      enumerator :: cudaRoundZero
      enumerator :: cudaRoundPosInf
      enumerator :: cudaRoundMinInf
    end enum ! cudaRoundMode

    enum, bind(C) !:: cudaError
      enumerator :: cudaSuccess=0
      enumerator :: cudaErrorMissingConfiguration=1
      enumerator :: cudaErrorMemoryAllocation=2
      enumerator :: cudaErrorInitializationError=3
      enumerator :: cudaErrorLaunchFailure=4
      enumerator :: cudaErrorPriorLaunchFailure=5
      enumerator :: cudaErrorLaunchTimeout=6
      enumerator :: cudaErrorLaunchOutOfResources=7
      enumerator :: cudaErrorInvalidDeviceFunction=8
      enumerator :: cudaErrorInvalidConfiguration=9
      enumerator :: cudaErrorInvalidDevice=10
      enumerator :: cudaErrorInvalidValue=11
      enumerator :: cudaErrorInvalidPitchValue=12
      enumerator :: cudaErrorInvalidSymbol=13
      enumerator :: cudaErrorMapBufferObjectFailed=14
      enumerator :: cudaErrorUnmapBufferObjectFailed=15
      enumerator :: cudaErrorInvalidHostPointer=16
      enumerator :: cudaErrorInvalidDevicePointer=17
      enumerator :: cudaErrorInvalidTexture=18
      enumerator :: cudaErrorInvalidTextureBinding=19
      enumerator :: cudaErrorInvalidChannelDescriptor=20
      enumerator :: cudaErrorInvalidMemcpyDirection=21
      enumerator :: cudaErrorAddressOfConstant=22
      enumerator :: cudaErrorTextureFetchFailed=23
      enumerator :: cudaErrorTextureNotBound=24
      enumerator :: cudaErrorSynchronizationError=25
      enumerator :: cudaErrorInvalidFilterSetting=26
      enumerator :: cudaErrorInvalidNormSetting=27
      enumerator :: cudaErrorMixedDeviceExecution=28
      enumerator :: cudaErrorCudartUnloading=29
      enumerator :: cudaErrorUnknown=30
      enumerator :: cudaErrorNotYetImplemented=31
      enumerator :: cudaErrorMemoryValueTooLarge=32
      enumerator :: cudaErrorInvalidResourceHandle=33
      enumerator :: cudaErrorNotReady=34
      enumerator :: cudaErrorInsufficientDriver=35
      enumerator :: cudaErrorSetOnActiveProcess=36
      enumerator :: cudaErrorInvalidSurface=37
      enumerator :: cudaErrorNoDevice=38
      enumerator :: cudaErrorECCUncorrectable=39
      enumerator :: cudaErrorSharedObjectSymbolNotFound=40
      enumerator :: cudaErrorSharedObjectInitFailed=41
      enumerator :: cudaErrorUnsupportedLimit=42
      enumerator :: cudaErrorDuplicateVariableName=43
      enumerator :: cudaErrorDuplicateTextureName=44
      enumerator :: cudaErrorDuplicateSurfaceName=45
      enumerator :: cudaErrorDevicesUnavailable=46
      enumerator :: cudaErrorStartupFailure=127
      enumerator :: cudaErrorApiFailureBase=10000
    end enum ! cudaError

    enum, bind(C) !:: cudaChannelFormatKind
      enumerator :: cudaChannelFormatKindSigned=0
      enumerator :: cudaChannelFormatKindUnsigned=1
      enumerator :: cudaChannelFormatKindFloat=2
      enumerator :: cudaChannelFormatKindNone=3
    end enum ! cudaChannelFormatKind

    enum, bind(C) !:: cudaMemcpyKind
      enumerator :: cudaMemcpyHostToHost=0
      enumerator :: cudaMemcpyHostToDevice=1
      enumerator :: cudaMemcpyDeviceToHost=2
      enumerator :: cudaMemcpyDeviceToDevice=3
    end enum ! cudaMemcpyKind

    enum, bind(C) !:: cudaGraphicsRegisterFlags
      enumerator :: cudaGraphicsRegisterFlagsNone=0
    end enum ! cudaGraphicsRegisterFlags

    enum, bind(C) !:: cudaGraphicsMapFlags
      enumerator :: cudaGraphicsMapFlagsNone=0
      enumerator :: cudaGraphicsMapFlagsReadOnly=1
      enumerator :: cudaGraphicsMapFlagsWriteDiscard=2
    end enum ! cudaGraphicsMapFlags

    enum, bind(C) !:: cudaGraphicsCubeFace
      enumerator :: cudaGraphicsCubeFacePositiveX=0
      enumerator :: cudaGraphicsCubeFaceNegativeX=1
      enumerator :: cudaGraphicsCubeFacePositiveY=2
      enumerator :: cudaGraphicsCubeFaceNegativeY=3
      enumerator :: cudaGraphicsCubeFacePositiveZ=4
      enumerator :: cudaGraphicsCubeFaceNegativeZ=5
    end enum ! cudaGraphicsCubeFace

    enum, bind(C) !:: cudaFuncCache
      enumerator :: cudaFuncCachePreferNone=0
      enumerator :: cudaFuncCachePreferShared=1
      enumerator :: cudaFuncCachePreferL1=2
    end enum ! cudaFuncCache

    enum, bind(C) !:: cudaComputeMode
      enumerator :: cudaComputeModeDefault=0
      enumerator :: cudaComputeModeExclusive=1
      enumerator :: cudaComputeModeProhibited=2
    end enum ! cudaComputeMode

    enum, bind(C) !:: cudaLimit
      enumerator :: cudaLimitStackSize=0
      enumerator :: cudaLimitPrintfFifoSize=1
    end enum ! cudaLimit

    enum, bind(C) !:: cudaSurfaceBoundaryMode
      enumerator :: cudaBoundaryModeZero=0
      enumerator :: cudaBoundaryModeClamp=1
      enumerator :: cudaBoundaryModeTrap=2
    end enum ! cudaSurfaceBoundaryMode

    enum, bind(C) !:: cudaSurfaceFormatMode
      enumerator :: cudaFormatModeForced
      enumerator :: cudaFormatModeAuto
    end enum ! cudaSurfaceFormatMode

    enum, bind(C) !:: cudaTextureAddressMode
      enumerator :: cudaAddressModeWrap
      enumerator :: cudaAddressModeClamp
      enumerator :: cudaAddressModeMirror
    end enum ! cudaTextureAddressMode

    enum, bind(C) !:: cudaTextureFilterMode
      enumerator :: cudaFilterModePoint
      enumerator :: cudaFilterModeLinear
    end enum ! cudaTextureFilterMode

    enum, bind(C) !:: cudaTextureReadMode
      enumerator :: cudaReadModeElementType
      enumerator :: cudaReadModeNormalizedFloat
    end enum ! cudaTextureReadMode

    type, bind(C) :: cudaChannelFormatDesc
      integer(c_int) :: x
      integer(c_int) :: y
      integer(c_int) :: z
      integer(c_int) :: w
      integer (KIND(cudaChannelFormatKindSigned)) :: f
    end type cudaChannelFormatDesc

    type, bind(C) :: cudaPitchedPtr
      type(c_ptr) :: ptr
      integer(c_int) :: pitch
      integer(c_int) :: xsize
      integer(c_int) :: ysize
    end type cudaPitchedPtr

    type, bind(C) :: cudaExtent
      integer(c_int) :: width
      integer(c_int) :: height
      integer(c_int) :: depth
    end type cudaExtent

    type, bind(C) :: cudaPos
      integer(c_int) :: x
      integer(c_int) :: y
      integer(c_int) :: z
    end type cudaPos

    type, bind(C) :: cudaMemcpy3DParms
      type (c_ptr) :: srcArray
      type (cudaPos) :: srcPos
      type (cudaPitchedPtr) :: srcPtr
      type (c_ptr) :: dstArray
      type (cudaPos) :: dstPos
      type (cudaPitchedPtr) :: dstPtr
      type (cudaExtent) :: extent
      integer (KIND(cudaMemcpyHostToHost)) :: kind
    end type cudaMemcpy3DParms

    type, bind(C) :: cudaFuncAttributes
      integer(c_int) :: sharedSizeBytes
      integer(c_int) :: constSizeBytes
      integer(c_int) :: localSizeBytes
      integer(c_int) :: maxThreadsPerBlock
      integer(c_int) :: numRegs
      integer(c_int) :: ptxVersion
      integer(c_int) :: binaryVersion
      integer(c_int) :: cudaReserved__(6)
    end type cudaFuncAttributes

    type, bind(C) :: cudaDeviceProp
      character(c_char) :: name(256)
      integer(c_int) :: totalGlobalMem
      integer(c_int) :: sharedMemPerBlock
      integer(c_int) :: regsPerBlock
      integer(c_int) :: warpSize
      integer(c_int) :: memPitch
      integer(c_int) :: maxThreadsPerBlock
      integer(c_int) :: maxThreadsDim(3)
      integer(c_int) :: maxGridSize(3)
      integer(c_int) :: clockRate
      integer(c_int) :: totalConstMem
      integer(c_int) :: major
      integer(c_int) :: minor
      integer(c_int) :: textureAlignment
      integer(c_int) :: deviceOverlap
      integer(c_int) :: multiProcessorCount
      integer(c_int) :: kernelExecTimeoutEnabled
      integer(c_int) :: integrated
      integer(c_int) :: canMapHostMemory
      integer(c_int) :: computeMode
      integer(c_int) :: maxTexture1D
      integer(c_int) :: maxTexture2D(2)
      integer(c_int) :: maxTexture3D(3)
      integer(c_int) :: maxTexture2DArray(3)
      integer(c_int) :: surfaceAlignment
      integer(c_int) :: concurrentKernels
      integer(c_int) :: ECCEnabled
      integer(c_int) :: pciBusID
      integer(c_int) :: pciDeviceID
      integer(c_int) :: cudaReserved__(22)
    end type cudaDeviceProp

    type, bind(C) :: surfaceReference
      type (cudaChannelFormatDesc) :: channelDesc
    end type surfaceReference

    type, bind(C) :: textureReference
      integer(c_int) :: normalized
      integer (KIND(cudaFilterModePoint)) :: filterMode
      integer (KIND(cudaAddressModeWrap)) :: addressMode(3)
      type (cudaChannelFormatDesc) :: channelDesc
      integer(c_int) :: cudaReserved__(16)
    end type textureReference

    type, bind(C) :: uchar1
      integer(c_signed_char) :: x
    end type uchar1

    type, bind(C) :: uchar3
      integer(c_signed_char) :: x
      integer(c_signed_char) :: y
      integer(c_signed_char) :: z
    end type uchar3

    type, bind(C) :: short1
      integer(c_short) :: x
    end type short1

    type, bind(C) :: short3
      integer(c_short) :: x
      integer(c_short) :: y
      integer(c_short) :: z
    end type short3

    type, bind(C) :: int1
      integer(c_int) :: x
    end type int1

    type, bind(C) :: uint1
      integer(c_int) :: x
    end type uint1

    type, bind(C) :: int3
      integer(c_int) :: x
      integer(c_int) :: y
      integer(c_int) :: z
    end type int3

    type, bind(C) :: uint3
      integer(c_int) :: x
      integer(c_int) :: y
      integer(c_int) :: z
    end type uint3

    type, bind(C) :: long1
      integer(c_long) :: x
    end type long1

    type, bind(C) :: long3
      integer(c_long) :: x
      integer(c_long) :: y
      integer(c_long) :: z
    end type long3

    type, bind(C) :: ulong3
      integer(c_long) :: x
      integer(c_long) :: y
      integer(c_long) :: z
    end type ulong3

    type, bind(C) :: float1
      real(c_float) :: x
    end type float1

    type, bind(C) :: float3
      real(c_float) :: x
      real(c_float) :: y
      real(c_float) :: z
    end type float3

    type, bind(C) :: longlong1
      integer(c_long_long) :: x
    end type longlong1

    type, bind(C) :: ulonglong1
      integer(c_long_long) :: x
    end type ulonglong1

    type, bind(C) :: longlong3
      integer(c_long_long) :: x
      integer(c_long_long) :: y
      integer(c_long_long) :: z
    end type longlong3

    type, bind(C) :: ulonglong3
      integer(c_long_long) :: x
      integer(c_long_long) :: y
      integer(c_long_long) :: z
    end type ulonglong3

    type, bind(C) :: double1
      real(c_double) :: x
    end type double1

    type, bind(C) :: double3
      real(c_double) :: x
      real(c_double) :: y
      real(c_double) :: z
    end type double3

    type, bind(C) :: dim3
      integer(c_int) :: x
      integer(c_int) :: y
      integer(c_int) :: z
    end type dim3

    interface ! [['cudaError_t', None], 'cudaMalloc3D', [['struct', 'cudaPitchedPtr', '*', 'pitchedDevPtr'], ['struct', 'cudaExtent', None, 'extent']]]
      function cudaMalloc3D(pitchedDevPtr,extent) result( res ) bind(C, name="cudaMalloc3D")
        use, intrinsic :: ISO_C_BINDING
        import cudaSuccess
        import cudaPitchedPtr
        import cudaExtent
        implicit none
        type (cudaPitchedPtr) :: pitchedDevPtr
        type (cudaExtent), value :: extent
        integer (KIND(cudaSuccess)) :: res
      end function cudaMalloc3D
    end interface

    interface ! [['cudaError_t', None], 'cudaMemset3D', [['struct', 'cudaPitchedPtr', None, 'pitchedDevPtr'], ['int', None, 'value'], ['struct', 'cudaExtent', None, 'extent']]]
      function cudaMemset3D(pitchedDevPtr,value,extent) result( res ) bind(C, name="cudaMemset3D")
        use, intrinsic :: ISO_C_BINDING
        import cudaSuccess
        import cudaPitchedPtr
        import cudaExtent
        implicit none
        type (cudaPitchedPtr), value :: pitchedDevPtr
        integer(c_int), value :: value
        type (cudaExtent), value :: extent
        integer (KIND(cudaSuccess)) :: res
      end function cudaMemset3D
    end interface

    interface ! [['cudaError_t', None], 'cudaMalloc', [['void', '**', 'devPtr'], ['size_t', None, 'size']]]
      function cudaMalloc(devPtr,size) result( res ) bind(C, name="cudaMalloc")
        use, intrinsic :: ISO_C_BINDING
        import cudaSuccess
        implicit none
        type(c_ptr) :: devPtr
        integer(c_int), value :: size
        integer (KIND(cudaSuccess)) :: res
      end function cudaMalloc
    end interface

    interface ! [['cudaError_t', None], 'cudaMallocHost', [['void', '**', 'ptr'], ['size_t', None, 'size']]]
      function cudaMallocHost(ptr,size) result( res ) bind(C, name="cudaMallocHost")
        use, intrinsic :: ISO_C_BINDING
        import cudaSuccess
        implicit none
        type(c_ptr) :: ptr
        integer(c_int), value :: size
        integer (KIND(cudaSuccess)) :: res
      end function cudaMallocHost
    end interface

    interface ! [['cudaError_t', None], 'cudaMallocPitch', [['void', '**', 'devPtr'], ['size_t', '*', 'pitch'], ['size_t', None, 'width'], ['size_t', None, 'height']]]
      function cudaMallocPitch(devPtr,pitch,width,height) result( res ) bind(C, name="cudaMallocPitch")
        use, intrinsic :: ISO_C_BINDING
        import cudaSuccess
        implicit none
        type(c_ptr) :: devPtr
        integer(c_int) :: pitch
        integer(c_int), value :: width
        integer(c_int), value :: height
        integer (KIND(cudaSuccess)) :: res
      end function cudaMallocPitch
    end interface

    interface ! [['cudaError_t', None], 'cudaFree', [['void', '*', 'devPtr']]]
      function cudaFree(devPtr) result( res ) bind(C, name="cudaFree")
        use, intrinsic :: ISO_C_BINDING
        import cudaSuccess
        implicit none
        type(c_ptr), value :: devPtr
        integer (KIND(cudaSuccess)) :: res
      end function cudaFree
    end interface

    interface ! [['cudaError_t', None], 'cudaFreeHost', [['void', '*', 'ptr']]]
      function cudaFreeHost(ptr) result( res ) bind(C, name="cudaFreeHost")
        use, intrinsic :: ISO_C_BINDING
        import cudaSuccess
        implicit none
        type(c_ptr), value :: ptr
        integer (KIND(cudaSuccess)) :: res
      end function cudaFreeHost
    end interface

    interface ! [['cudaError_t', None], 'cudaFreeArray', [['struct', 'cudaArray', '*', 'array']]]
      function cudaFreeArray(array) result( res ) bind(C, name="cudaFreeArray")
        use, intrinsic :: ISO_C_BINDING
        import cudaSuccess
        implicit none
        type (c_ptr), value :: array
        integer (KIND(cudaSuccess)) :: res
      end function cudaFreeArray
    end interface

    interface ! [['cudaError_t', None], 'cudaHostAlloc', [['void', '**', 'pHost'], ['size_t', None, 'bytes'], ['unsigned', 'int', None, 'flags']]]
      function cudaHostAlloc(pHost,bytes,flags) result( res ) bind(C, name="cudaHostAlloc")
        use, intrinsic :: ISO_C_BINDING
        import cudaSuccess
        implicit none
        type(c_ptr) :: pHost
        integer(c_int), value :: bytes
        integer(c_int), value :: flags
        integer (KIND(cudaSuccess)) :: res
      end function cudaHostAlloc
    end interface

    interface ! [['cudaError_t', None], 'cudaHostGetDevicePointer', [['void', '**', 'pDevice'], ['void', '*', 'pHost'], ['unsigned', 'int', None, 'flags']]]
      function cudaHostGetDevicePointer(pDevice,pHost,flags) result( res ) bind(C, name="cudaHostGetDevicePointer")
        use, intrinsic :: ISO_C_BINDING
        import cudaSuccess
        implicit none
        type(c_ptr) :: pDevice
        type(c_ptr), value :: pHost
        integer(c_int), value :: flags
        integer (KIND(cudaSuccess)) :: res
      end function cudaHostGetDevicePointer
    end interface

    interface ! [['cudaError_t', None], 'cudaHostGetFlags', [['unsigned', 'int', '*', 'pFlags'], ['void', '*', 'pHost']]]
      function cudaHostGetFlags(pFlags,pHost) result( res ) bind(C, name="cudaHostGetFlags")
        use, intrinsic :: ISO_C_BINDING
        import cudaSuccess
        implicit none
        integer(c_int) :: pFlags
        type(c_ptr), value :: pHost
        integer (KIND(cudaSuccess)) :: res
      end function cudaHostGetFlags
    end interface

    interface ! [['cudaError_t', None], 'cudaMemGetInfo', [['size_t', '*', 'free'], ['size_t', '*', 'total']]]
      function cudaMemGetInfo(free,total) result( res ) bind(C, name="cudaMemGetInfo")
        use, intrinsic :: ISO_C_BINDING
        import cudaSuccess
        implicit none
        integer(c_int) :: free
        integer(c_int) :: total
        integer (KIND(cudaSuccess)) :: res
      end function cudaMemGetInfo
    end interface

    interface ! [['cudaError_t', None], 'cudaMemcpy', [['void', '*', 'dst'], ['const', 'void', '*', 'src'], ['size_t', None, 'count'], ['enum', 'cudaMemcpyKind', None, 'kind']]]
      function cudaMemcpy(dst,src,count,kind_arg) result( res ) bind(C, name="cudaMemcpy")
        use, intrinsic :: ISO_C_BINDING
        import cudaSuccess
        import cudaMemcpyHostToHost
        implicit none
        type(c_ptr), value :: dst
        type(c_ptr), value :: src
        integer(c_int), value :: count
        integer (KIND(cudaMemcpyHostToHost)), value :: kind_arg
        integer (KIND(cudaSuccess)) :: res
      end function cudaMemcpy
    end interface

    interface ! [['cudaError_t', None], 'cudaMemcpyToArray', [['struct', 'cudaArray', '*', 'dst'], ['size_t', None, 'wOffset'], ['size_t', None, 'hOffset'], ['const', 'void', '*', 'src'], ['size_t', None, 'count'], ['enum', 'cudaMemcpyKind', None, 'kind']]]
      function cudaMemcpyToArray(dst,wOffset,hOffset,src,count,kind_arg) result( res ) bind(C, name="cudaMemcpyToArray")
        use, intrinsic :: ISO_C_BINDING
        import cudaSuccess
        import cudaMemcpyHostToHost
        implicit none
        type (c_ptr), value :: dst
        integer(c_int), value :: wOffset
        integer(c_int), value :: hOffset
        type(c_ptr), value :: src
        integer(c_int), value :: count
        integer (KIND(cudaMemcpyHostToHost)), value :: kind_arg
        integer (KIND(cudaSuccess)) :: res
      end function cudaMemcpyToArray
    end interface

    interface ! [['cudaError_t', None], 'cudaMemcpy2D', [['void', '*', 'dst'], ['size_t', None, 'dpitch'], ['const', 'void', '*', 'src'], ['size_t', None, 'spitch'], ['size_t', None, 'width'], ['size_t', None, 'height'], ['enum', 'cudaMemcpyKind', None, 'kind']]]
      function cudaMemcpy2D(dst,dpitch,src,spitch,width,height,kind_arg) result( res ) bind(C, name="cudaMemcpy2D")
        use, intrinsic :: ISO_C_BINDING
        import cudaSuccess
        import cudaMemcpyHostToHost
        implicit none
        type(c_ptr), value :: dst
        integer(c_int), value :: dpitch
        type(c_ptr), value :: src
        integer(c_int), value :: spitch
        integer(c_int), value :: width
        integer(c_int), value :: height
        integer (KIND(cudaMemcpyHostToHost)), value :: kind_arg
        integer (KIND(cudaSuccess)) :: res
      end function cudaMemcpy2D
    end interface

    interface ! [['cudaError_t', None], 'cudaMemcpy2DToArray', [['struct', 'cudaArray', '*', 'dst'], ['size_t', None, 'wOffset'], ['size_t', None, 'hOffset'], ['const', 'void', '*', 'src'], ['size_t', None, 'spitch'], ['size_t', None, 'width'], ['size_t', None, 'height'], ['enum', 'cudaMemcpyKind', None, 'kind']]]
      function cudaMemcpy2DToArray(dst,wOffset,hOffset,src,spitch,width,height,kind_arg) result( res ) bind(C, name="cudaMemcpy2DToArray")
        use, intrinsic :: ISO_C_BINDING
        import cudaSuccess
        import cudaMemcpyHostToHost
        implicit none
        type (c_ptr), value :: dst
        integer(c_int), value :: wOffset
        integer(c_int), value :: hOffset
        type(c_ptr), value :: src
        integer(c_int), value :: spitch
        integer(c_int), value :: width
        integer(c_int), value :: height
        integer (KIND(cudaMemcpyHostToHost)), value :: kind_arg
        integer (KIND(cudaSuccess)) :: res
      end function cudaMemcpy2DToArray
    end interface

    interface ! [['cudaError_t', None], 'cudaMemcpyToSymbol', [['const', 'char', '*', 'symbol'], ['const', 'void', '*', 'src'], ['size_t', None, 'count'], ['size_t', None, 'offset'], ['enum', 'cudaMemcpyKind', None, 'kind']]]
      function cudaMemcpyToSymbol(symbol,src,count,offset,kind_arg) result( res ) bind(C, name="cudaMemcpyToSymbol")
        use, intrinsic :: ISO_C_BINDING
        import cudaSuccess
        import cudaMemcpyHostToHost
        implicit none
        character(c_char) :: symbol
        type(c_ptr), value :: src
        integer(c_int), value :: count
        integer(c_int), value :: offset
        integer (KIND(cudaMemcpyHostToHost)), value :: kind_arg
        integer (KIND(cudaSuccess)) :: res
      end function cudaMemcpyToSymbol
    end interface

    interface ! [['cudaError_t', None], 'cudaMemcpyFromSymbol', [['void', '*', 'dst'], ['const', 'char', '*', 'symbol'], ['size_t', None, 'count'], ['size_t', None, 'offset'], ['enum', 'cudaMemcpyKind', None, 'kind']]]
      function cudaMemcpyFromSymbol(dst,symbol,count,offset,kind_arg) result( res ) bind(C, name="cudaMemcpyFromSymbol")
        use, intrinsic :: ISO_C_BINDING
        import cudaSuccess
        import cudaMemcpyHostToHost
        implicit none
        type(c_ptr), value :: dst
        character(c_char) :: symbol
        integer(c_int), value :: count
        integer(c_int), value :: offset
        integer (KIND(cudaMemcpyHostToHost)), value :: kind_arg
        integer (KIND(cudaSuccess)) :: res
      end function cudaMemcpyFromSymbol
    end interface

    interface ! [['cudaError_t', None], 'cudaMemcpyAsync', [['void', '*', 'dst'], ['const', 'void', '*', 'src'], ['size_t', None, 'count'], ['enum', 'cudaMemcpyKind', None, 'kind'], ['cudaStream_t', None, 'stream']]]
      function cudaMemcpyAsync(dst,src,count,kind_arg,stream) result( res ) bind(C, name="cudaMemcpyAsync")
        use, intrinsic :: ISO_C_BINDING
        import cudaSuccess
        import cudaMemcpyHostToHost
        import cudaStream_t
        implicit none
        type(c_ptr), value :: dst
        type(c_ptr), value :: src
        integer(c_int), value :: count
        integer (KIND(cudaMemcpyHostToHost)), value :: kind_arg
        type (cudaStream_t), value :: stream
        integer (KIND(cudaSuccess)) :: res
      end function cudaMemcpyAsync
    end interface

    interface ! [['cudaError_t', None], 'cudaMemcpyToArrayAsync', [['struct', 'cudaArray', '*', 'dst'], ['size_t', None, 'wOffset'], ['size_t', None, 'hOffset'], ['const', 'void', '*', 'src'], ['size_t', None, 'count'], ['enum', 'cudaMemcpyKind', None, 'kind'], ['cudaStream_t', None, 'stream']]]
      function cudaMemcpyToArrayAsync(dst,wOffset,hOffset,src,count,kind_arg,stream) result( res ) bind(C, name="cudaMemcpyToArrayAsync")
        use, intrinsic :: ISO_C_BINDING
        import cudaSuccess
        import cudaMemcpyHostToHost
        import cudaStream_t
        implicit none
        type (c_ptr), value :: dst
        integer(c_int), value :: wOffset
        integer(c_int), value :: hOffset
        type(c_ptr), value :: src
        integer(c_int), value :: count
        integer (KIND(cudaMemcpyHostToHost)), value :: kind_arg
        type (cudaStream_t), value :: stream
        integer (KIND(cudaSuccess)) :: res
      end function cudaMemcpyToArrayAsync
    end interface

    interface ! [['cudaError_t', None], 'cudaMemcpy2DAsync', [['void', '*', 'dst'], ['size_t', None, 'dpitch'], ['const', 'void', '*', 'src'], ['size_t', None, 'spitch'], ['size_t', None, 'width'], ['size_t', None, 'height'], ['enum', 'cudaMemcpyKind', None, 'kind'], ['cudaStream_t', None, 'stream']]]
      function cudaMemcpy2DAsync(dst,dpitch,src,spitch,width,height,kind_arg,stream) result( res ) bind(C, name="cudaMemcpy2DAsync")
        use, intrinsic :: ISO_C_BINDING
        import cudaSuccess
        import cudaMemcpyHostToHost
        import cudaStream_t
        implicit none
        type(c_ptr), value :: dst
        integer(c_int), value :: dpitch
        type(c_ptr), value :: src
        integer(c_int), value :: spitch
        integer(c_int), value :: width
        integer(c_int), value :: height
        integer (KIND(cudaMemcpyHostToHost)), value :: kind_arg
        type (cudaStream_t), value :: stream
        integer (KIND(cudaSuccess)) :: res
      end function cudaMemcpy2DAsync
    end interface

    interface ! [['cudaError_t', None], 'cudaMemcpy2DToArrayAsync', [['struct', 'cudaArray', '*', 'dst'], ['size_t', None, 'wOffset'], ['size_t', None, 'hOffset'], ['const', 'void', '*', 'src'], ['size_t', None, 'spitch'], ['size_t', None, 'width'], ['size_t', None, 'height'], ['enum', 'cudaMemcpyKind', None, 'kind'], ['cudaStream_t', None, 'stream']]]
      function cudaMemcpy2DToArrayAsync(dst,wOffset,hOffset,src,spitch,width,height,kind_arg,stream) result( res ) bind(C, name="cudaMemcpy2DToArrayAsync")
        use, intrinsic :: ISO_C_BINDING
        import cudaSuccess
        import cudaMemcpyHostToHost
        import cudaStream_t
        implicit none
        type (c_ptr), value :: dst
        integer(c_int), value :: wOffset
        integer(c_int), value :: hOffset
        type(c_ptr), value :: src
        integer(c_int), value :: spitch
        integer(c_int), value :: width
        integer(c_int), value :: height
        integer (KIND(cudaMemcpyHostToHost)), value :: kind_arg
        type (cudaStream_t), value :: stream
        integer (KIND(cudaSuccess)) :: res
      end function cudaMemcpy2DToArrayAsync
    end interface

    interface ! [['cudaError_t', None], 'cudaMemcpyToSymbolAsync', [['const', 'char', '*', 'symbol'], ['const', 'void', '*', 'src'], ['size_t', None, 'count'], ['size_t', None, 'offset'], ['enum', 'cudaMemcpyKind', None, 'kind'], ['cudaStream_t', None, 'stream']]]
      function cudaMemcpyToSymbolAsync(symbol,src,count,offset,kind_arg,stream) result( res ) bind(C, name="cudaMemcpyToSymbolAsync")
        use, intrinsic :: ISO_C_BINDING
        import cudaSuccess
        import cudaMemcpyHostToHost
        import cudaStream_t
        implicit none
        character(c_char) :: symbol
        type(c_ptr), value :: src
        integer(c_int), value :: count
        integer(c_int), value :: offset
        integer (KIND(cudaMemcpyHostToHost)), value :: kind_arg
        type (cudaStream_t), value :: stream
        integer (KIND(cudaSuccess)) :: res
      end function cudaMemcpyToSymbolAsync
    end interface

    interface ! [['cudaError_t', None], 'cudaMemcpyFromSymbolAsync', [['void', '*', 'dst'], ['const', 'char', '*', 'symbol'], ['size_t', None, 'count'], ['size_t', None, 'offset'], ['enum', 'cudaMemcpyKind', None, 'kind'], ['cudaStream_t', None, 'stream']]]
      function cudaMemcpyFromSymbolAsync(dst,symbol,count,offset,kind_arg,stream) result( res ) bind(C, name="cudaMemcpyFromSymbolAsync")
        use, intrinsic :: ISO_C_BINDING
        import cudaSuccess
        import cudaMemcpyHostToHost
        import cudaStream_t
        implicit none
        type(c_ptr), value :: dst
        character(c_char) :: symbol
        integer(c_int), value :: count
        integer(c_int), value :: offset
        integer (KIND(cudaMemcpyHostToHost)), value :: kind_arg
        type (cudaStream_t), value :: stream
        integer (KIND(cudaSuccess)) :: res
      end function cudaMemcpyFromSymbolAsync
    end interface

    interface ! [['cudaError_t', None], 'cudaMemset', [['void', '*', 'devPtr'], ['int', None, 'value'], ['size_t', None, 'count']]]
      function cudaMemset(devPtr,value,count) result( res ) bind(C, name="cudaMemset")
        use, intrinsic :: ISO_C_BINDING
        import cudaSuccess
        implicit none
        type(c_ptr), value :: devPtr
        integer(c_int), value :: value
        integer(c_int), value :: count
        integer (KIND(cudaSuccess)) :: res
      end function cudaMemset
    end interface

    interface ! [['cudaError_t', None], 'cudaMemset2D', [['void', '*', 'devPtr'], ['size_t', None, 'pitch'], ['int', None, 'value'], ['size_t', None, 'width'], ['size_t', None, 'height']]]
      function cudaMemset2D(devPtr,pitch,value,width,height) result( res ) bind(C, name="cudaMemset2D")
        use, intrinsic :: ISO_C_BINDING
        import cudaSuccess
        implicit none
        type(c_ptr), value :: devPtr
        integer(c_int), value :: pitch
        integer(c_int), value :: value
        integer(c_int), value :: width
        integer(c_int), value :: height
        integer (KIND(cudaSuccess)) :: res
      end function cudaMemset2D
    end interface

    interface ! [['cudaError_t', None], 'cudaGetSymbolAddress', [['void', '**', 'devPtr'], ['const', 'char', '*', 'symbol']]]
      function cudaGetSymbolAddress(devPtr,symbol) result( res ) bind(C, name="cudaGetSymbolAddress")
        use, intrinsic :: ISO_C_BINDING
        import cudaSuccess
        implicit none
        type(c_ptr) :: devPtr
        character(c_char) :: symbol
        integer (KIND(cudaSuccess)) :: res
      end function cudaGetSymbolAddress
    end interface

    interface ! [['cudaError_t', None], 'cudaGetSymbolSize', [['size_t', '*', 'size'], ['const', 'char', '*', 'symbol']]]
      function cudaGetSymbolSize(size,symbol) result( res ) bind(C, name="cudaGetSymbolSize")
        use, intrinsic :: ISO_C_BINDING
        import cudaSuccess
        implicit none
        integer(c_int) :: size
        character(c_char) :: symbol
        integer (KIND(cudaSuccess)) :: res
      end function cudaGetSymbolSize
    end interface

    interface ! [['cudaError_t', None], 'cudaGetDeviceCount', [['int', '*', 'count']]]
      function cudaGetDeviceCount(count) result( res ) bind(C, name="cudaGetDeviceCount")
        use, intrinsic :: ISO_C_BINDING
        import cudaSuccess
        implicit none
        integer(c_int) :: count
        integer (KIND(cudaSuccess)) :: res
      end function cudaGetDeviceCount
    end interface

    interface ! [['cudaError_t', None], 'cudaGetDeviceProperties', [['struct', 'cudaDeviceProp', '*', 'prop'], ['int', None, 'device']]]
      function cudaGetDeviceProperties(prop,device) result( res ) bind(C, name="cudaGetDeviceProperties")
        use, intrinsic :: ISO_C_BINDING
        import cudaSuccess
        import cudaDeviceProp
        implicit none
        type (cudaDeviceProp) :: prop
        integer(c_int), value :: device
        integer (KIND(cudaSuccess)) :: res
      end function cudaGetDeviceProperties
    end interface

    interface ! [['cudaError_t', None], 'cudaSetDevice', [['int', None, 'device']]]
      function cudaSetDevice(device) result( res ) bind(C, name="cudaSetDevice")
        use, intrinsic :: ISO_C_BINDING
        import cudaSuccess
        implicit none
        integer(c_int), value :: device
        integer (KIND(cudaSuccess)) :: res
      end function cudaSetDevice
    end interface

    interface ! [['cudaError_t', None], 'cudaGetDevice', [['int', '*', 'device']]]
      function cudaGetDevice(device) result( res ) bind(C, name="cudaGetDevice")
        use, intrinsic :: ISO_C_BINDING
        import cudaSuccess
        implicit none
        integer(c_int) :: device
        integer (KIND(cudaSuccess)) :: res
      end function cudaGetDevice
    end interface

    interface ! [['cudaError_t', None], 'cudaSetValidDevices', [['int', '*', 'device_arr'], ['int', None, 'len']]]
      function cudaSetValidDevices(device_arr,len) result( res ) bind(C, name="cudaSetValidDevices")
        use, intrinsic :: ISO_C_BINDING
        import cudaSuccess
        implicit none
        integer(c_int) :: device_arr
        integer(c_int), value :: len
        integer (KIND(cudaSuccess)) :: res
      end function cudaSetValidDevices
    end interface

    interface ! [['cudaError_t', None], 'cudaSetDeviceFlags', [['int', None, 'flags']]]
      function cudaSetDeviceFlags(flags) result( res ) bind(C, name="cudaSetDeviceFlags")
        use, intrinsic :: ISO_C_BINDING
        import cudaSuccess
        implicit none
        integer(c_int), value :: flags
        integer (KIND(cudaSuccess)) :: res
      end function cudaSetDeviceFlags
    end interface

    interface ! [['struct', 'cudaChannelFormatDesc', None], 'cudaCreateChannelDesc', [['int', None, 'x'], ['int', None, 'y'], ['int', None, 'z'], ['int', None, 'w'], ['enum', 'cudaChannelFormatKind', None, 'f']]]
      function cudaCreateChannelDesc(x,y,z,w,f) result( res ) bind(C, name="cudaCreateChannelDesc")
        use, intrinsic :: ISO_C_BINDING
        import cudaChannelFormatDesc
        import cudaChannelFormatKindSigned
        implicit none
        integer(c_int), value :: x
        integer(c_int), value :: y
        integer(c_int), value :: z
        integer(c_int), value :: w
        integer (KIND(cudaChannelFormatKindSigned)), value :: f
        type (cudaChannelFormatDesc) :: res
      end function cudaCreateChannelDesc
    end interface

    interface ! [['cudaError_t', None], 'cudaGetLastError', ['void']]
      function cudaGetLastError() result( res ) bind(C, name="cudaGetLastError")
        use, intrinsic :: ISO_C_BINDING
        import cudaSuccess
        implicit none
        integer (KIND(cudaSuccess)) :: res
      end function cudaGetLastError
    end interface

    interface ! [['cudaError_t', None], 'cudaPeekAtLastError', ['void']]
      function cudaPeekAtLastError() result( res ) bind(C, name="cudaPeekAtLastError")
        use, intrinsic :: ISO_C_BINDING
        import cudaSuccess
        implicit none
        integer (KIND(cudaSuccess)) :: res
      end function cudaPeekAtLastError
    end interface

    interface ! [['const', 'char', '*'], 'cudaGetErrorString', [['cudaError_t', None, 'error']]]
      function cudaGetErrorString(error) result( res ) bind(C, name="cudaGetErrorString")
        use, intrinsic :: ISO_C_BINDING
        import cudaSuccess
        implicit none
        integer (KIND(cudaSuccess)), value :: error
        type(c_ptr) :: res
      end function cudaGetErrorString
    end interface

    interface ! [['cudaError_t', None], 'cudaConfigureCall', [['dim3', None, 'gridDim'], ['dim3', None, 'blockDim'], ['size_t', None, 'sharedMem'], ['cudaStream_t', None, 'stream']]]
      function cudaConfigureCall(gridDim,blockDim,sharedMem,stream) result( res ) bind(C, name="cudaConfigureCall")
        use, intrinsic :: ISO_C_BINDING
        import cudaSuccess
        import dim3
        import dim3
        import cudaStream_t
        implicit none
        type (dim3), value :: gridDim
        type (dim3), value :: blockDim
        integer(c_int), value :: sharedMem
        type (cudaStream_t), value :: stream
        integer (KIND(cudaSuccess)) :: res
      end function cudaConfigureCall
    end interface

    interface ! [['cudaError_t', None], 'cudaSetupArgument', [['const', 'void', '*', 'arg'], ['size_t', None, 'size'], ['size_t', None, 'offset']]]
      function cudaSetupArgument(arg,size,offset) result( res ) bind(C, name="cudaSetupArgument")
        use, intrinsic :: ISO_C_BINDING
        import cudaSuccess
        implicit none
        type(c_ptr), value :: arg
        integer(c_int), value :: size
        integer(c_int), value :: offset
        integer (KIND(cudaSuccess)) :: res
      end function cudaSetupArgument
    end interface

    interface ! [['cudaError_t', None], 'cudaFuncSetCacheConfig', [['const', 'char', '*', 'func'], ['enum', 'cudaFuncCache', None, 'cacheConfig']]]
      function cudaFuncSetCacheConfig(func,cacheConfig) result( res ) bind(C, name="cudaFuncSetCacheConfig")
        use, intrinsic :: ISO_C_BINDING
        import cudaSuccess
        import cudaFuncCachePreferNone
        implicit none
        character(c_char) :: func
        integer (KIND(cudaFuncCachePreferNone)), value :: cacheConfig
        integer (KIND(cudaSuccess)) :: res
      end function cudaFuncSetCacheConfig
    end interface

    interface ! [['cudaError_t', None], 'cudaLaunch', [['const', 'char', '*', 'entry']]]
      function cudaLaunch(entry) result( res ) bind(C, name="cudaLaunch")
        use, intrinsic :: ISO_C_BINDING
        import cudaSuccess
        implicit none
        character(c_char) :: entry
        integer (KIND(cudaSuccess)) :: res
      end function cudaLaunch
    end interface

    interface ! [['cudaError_t', None], 'cudaFuncGetAttributes', [['struct', 'cudaFuncAttributes', '*', 'attr'], ['const', 'char', '*', 'func']]]
      function cudaFuncGetAttributes(attr,func) result( res ) bind(C, name="cudaFuncGetAttributes")
        use, intrinsic :: ISO_C_BINDING
        import cudaSuccess
        import cudaFuncAttributes
        implicit none
        type (cudaFuncAttributes) :: attr
        character(c_char) :: func
        integer (KIND(cudaSuccess)) :: res
      end function cudaFuncGetAttributes
    end interface

    interface ! [['cudaError_t', None], 'cudaStreamCreate', [['cudaStream_t', '*', 'pStream']]]
      function cudaStreamCreate(pStream) result( res ) bind(C, name="cudaStreamCreate")
        use, intrinsic :: ISO_C_BINDING
        import cudaSuccess
        import cudaStream_t
        implicit none
        type (cudaStream_t) :: pStream
        integer (KIND(cudaSuccess)) :: res
      end function cudaStreamCreate
    end interface

    interface ! [['cudaError_t', None], 'cudaStreamDestroy', [['cudaStream_t', None, 'stream']]]
      function cudaStreamDestroy(stream) result( res ) bind(C, name="cudaStreamDestroy")
        use, intrinsic :: ISO_C_BINDING
        import cudaSuccess
        import cudaStream_t
        implicit none
        type (cudaStream_t), value :: stream
        integer (KIND(cudaSuccess)) :: res
      end function cudaStreamDestroy
    end interface

    interface ! [['cudaError_t', None], 'cudaStreamSynchronize', [['cudaStream_t', None, 'stream']]]
      function cudaStreamSynchronize(stream) result( res ) bind(C, name="cudaStreamSynchronize")
        use, intrinsic :: ISO_C_BINDING
        import cudaSuccess
        import cudaStream_t
        implicit none
        type (cudaStream_t), value :: stream
        integer (KIND(cudaSuccess)) :: res
      end function cudaStreamSynchronize
    end interface

    interface ! [['cudaError_t', None], 'cudaStreamQuery', [['cudaStream_t', None, 'stream']]]
      function cudaStreamQuery(stream) result( res ) bind(C, name="cudaStreamQuery")
        use, intrinsic :: ISO_C_BINDING
        import cudaSuccess
        import cudaStream_t
        implicit none
        type (cudaStream_t), value :: stream
        integer (KIND(cudaSuccess)) :: res
      end function cudaStreamQuery
    end interface

    interface ! [['cudaError_t', None], 'cudaEventCreate', [['cudaEvent_t', '*', 'event']]]
      function cudaEventCreate(event) result( res ) bind(C, name="cudaEventCreate")
        use, intrinsic :: ISO_C_BINDING
        import cudaSuccess
        import cudaEvent_t
        implicit none
        type (cudaEvent_t) :: event
        integer (KIND(cudaSuccess)) :: res
      end function cudaEventCreate
    end interface

    interface ! [['cudaError_t', None], 'cudaEventCreateWithFlags', [['cudaEvent_t', '*', 'event'], ['int', None, 'flags']]]
      function cudaEventCreateWithFlags(event,flags) result( res ) bind(C, name="cudaEventCreateWithFlags")
        use, intrinsic :: ISO_C_BINDING
        import cudaSuccess
        import cudaEvent_t
        implicit none
        type (cudaEvent_t) :: event
        integer(c_int), value :: flags
        integer (KIND(cudaSuccess)) :: res
      end function cudaEventCreateWithFlags
    end interface

    interface ! [['cudaError_t', None], 'cudaEventRecord', [['cudaEvent_t', None, 'event'], ['cudaStream_t', None, 'stream']]]
      function cudaEventRecord(event,stream) result( res ) bind(C, name="cudaEventRecord")
        use, intrinsic :: ISO_C_BINDING
        import cudaSuccess
        import cudaEvent_t
        import cudaStream_t
        implicit none
        type (cudaEvent_t), value :: event
        type (cudaStream_t), value :: stream
        integer (KIND(cudaSuccess)) :: res
      end function cudaEventRecord
    end interface

    interface ! [['cudaError_t', None], 'cudaEventQuery', [['cudaEvent_t', None, 'event']]]
      function cudaEventQuery(event) result( res ) bind(C, name="cudaEventQuery")
        use, intrinsic :: ISO_C_BINDING
        import cudaSuccess
        import cudaEvent_t
        implicit none
        type (cudaEvent_t), value :: event
        integer (KIND(cudaSuccess)) :: res
      end function cudaEventQuery
    end interface

    interface ! [['cudaError_t', None], 'cudaEventSynchronize', [['cudaEvent_t', None, 'event']]]
      function cudaEventSynchronize(event) result( res ) bind(C, name="cudaEventSynchronize")
        use, intrinsic :: ISO_C_BINDING
        import cudaSuccess
        import cudaEvent_t
        implicit none
        type (cudaEvent_t), value :: event
        integer (KIND(cudaSuccess)) :: res
      end function cudaEventSynchronize
    end interface

    interface ! [['cudaError_t', None], 'cudaEventDestroy', [['cudaEvent_t', None, 'event']]]
      function cudaEventDestroy(event) result( res ) bind(C, name="cudaEventDestroy")
        use, intrinsic :: ISO_C_BINDING
        import cudaSuccess
        import cudaEvent_t
        implicit none
        type (cudaEvent_t), value :: event
        integer (KIND(cudaSuccess)) :: res
      end function cudaEventDestroy
    end interface

    interface ! [['cudaError_t', None], 'cudaEventElapsedTime', [['float', '*', 'ms'], ['cudaEvent_t', None, 'start'], ['cudaEvent_t', None, 'end']]]
      function cudaEventElapsedTime(ms,start,end) result( res ) bind(C, name="cudaEventElapsedTime")
        use, intrinsic :: ISO_C_BINDING
        import cudaSuccess
        import cudaEvent_t
        import cudaEvent_t
        implicit none
        real(c_float) :: ms
        type (cudaEvent_t), value :: start
        type (cudaEvent_t), value :: end
        integer (KIND(cudaSuccess)) :: res
      end function cudaEventElapsedTime
    end interface

    interface ! [['cudaError_t', None], 'cudaSetDoubleForDevice', [['double', '*', 'd']]]
      function cudaSetDoubleForDevice(d) result( res ) bind(C, name="cudaSetDoubleForDevice")
        use, intrinsic :: ISO_C_BINDING
        import cudaSuccess
        implicit none
        real(c_double) :: d
        integer (KIND(cudaSuccess)) :: res
      end function cudaSetDoubleForDevice
    end interface

    interface ! [['cudaError_t', None], 'cudaSetDoubleForHost', [['double', '*', 'd']]]
      function cudaSetDoubleForHost(d) result( res ) bind(C, name="cudaSetDoubleForHost")
        use, intrinsic :: ISO_C_BINDING
        import cudaSuccess
        implicit none
        real(c_double) :: d
        integer (KIND(cudaSuccess)) :: res
      end function cudaSetDoubleForHost
    end interface

    interface ! [['cudaError_t', None], 'cudaThreadExit', ['void']]
      function cudaThreadExit() result( res ) bind(C, name="cudaThreadExit")
        use, intrinsic :: ISO_C_BINDING
        import cudaSuccess
        implicit none
        integer (KIND(cudaSuccess)) :: res
      end function cudaThreadExit
    end interface

    interface ! [['cudaError_t', None], 'cudaThreadSynchronize', ['void']]
      function cudaThreadSynchronize() result( res ) bind(C, name="cudaThreadSynchronize")
        use, intrinsic :: ISO_C_BINDING
        import cudaSuccess
        implicit none
        integer (KIND(cudaSuccess)) :: res
      end function cudaThreadSynchronize
    end interface

    interface ! [['cudaError_t', None], 'cudaThreadSetLimit', [['enum', 'cudaLimit', None, 'limit'], ['size_t', None, 'value']]]
      function cudaThreadSetLimit(limit,value) result( res ) bind(C, name="cudaThreadSetLimit")
        use, intrinsic :: ISO_C_BINDING
        import cudaSuccess
        import cudaLimitStackSize
        implicit none
        integer (KIND(cudaLimitStackSize)), value :: limit
        integer(c_int), value :: value
        integer (KIND(cudaSuccess)) :: res
      end function cudaThreadSetLimit
    end interface

    interface ! [['cudaError_t', None], 'cudaThreadGetLimit', [['size_t', '*', 'pValue'], ['enum', 'cudaLimit', None, 'limit']]]
      function cudaThreadGetLimit(pValue,limit) result( res ) bind(C, name="cudaThreadGetLimit")
        use, intrinsic :: ISO_C_BINDING
        import cudaSuccess
        import cudaLimitStackSize
        implicit none
        integer(c_int) :: pValue
        integer (KIND(cudaLimitStackSize)), value :: limit
        integer (KIND(cudaSuccess)) :: res
      end function cudaThreadGetLimit
    end interface

    interface ! [['cudaError_t', None], 'cudaDriverGetVersion', [['int', '*', 'driverVersion']]]
      function cudaDriverGetVersion(driverVersion) result( res ) bind(C, name="cudaDriverGetVersion")
        use, intrinsic :: ISO_C_BINDING
        import cudaSuccess
        implicit none
        integer(c_int) :: driverVersion
        integer (KIND(cudaSuccess)) :: res
      end function cudaDriverGetVersion
    end interface

    interface ! [['cudaError_t', None], 'cudaRuntimeGetVersion', [['int', '*', 'runtimeVersion']]]
      function cudaRuntimeGetVersion(runtimeVersion) result( res ) bind(C, name="cudaRuntimeGetVersion")
        use, intrinsic :: ISO_C_BINDING
        import cudaSuccess
        implicit none
        integer(c_int) :: runtimeVersion
        integer (KIND(cudaSuccess)) :: res
      end function cudaRuntimeGetVersion
    end interface

    interface ! [['cudaError_t', None], 'cudaGraphicsUnregisterResource', [['struct', 'cudaGraphicsResource', '*', 'resource']]]
      function cudaGraphicsUnregisterResource(resource) result( res ) bind(C, name="cudaGraphicsUnregisterResource")
        use, intrinsic :: ISO_C_BINDING
        import cudaSuccess
        implicit none
        type (c_ptr), value :: resource
        integer (KIND(cudaSuccess)) :: res
      end function cudaGraphicsUnregisterResource
    end interface

    interface ! [['cudaError_t', None], 'cudaGraphicsResourceSetMapFlags', [['struct', 'cudaGraphicsResource', '*', 'resource'], ['unsigned', 'int', None, 'flags']]]
      function cudaGraphicsResourceSetMapFlags(resource,flags) result( res ) bind(C, name="cudaGraphicsResourceSetMapFlags")
        use, intrinsic :: ISO_C_BINDING
        import cudaSuccess
        implicit none
        type (c_ptr), value :: resource
        integer(c_int), value :: flags
        integer (KIND(cudaSuccess)) :: res
      end function cudaGraphicsResourceSetMapFlags
    end interface

    interface ! [['cudaError_t', None], 'cudaGraphicsMapResources', [['int', None, 'count'], ['struct', 'cudaGraphicsResource', '**', 'resources'], ['cudaStream_t', None, 'stream']]]
      function cudaGraphicsMapResources(count,resources,stream) result( res ) bind(C, name="cudaGraphicsMapResources")
        use, intrinsic :: ISO_C_BINDING
        import cudaSuccess
        import cudaStream_t
        implicit none
        integer(c_int), value :: count
        type (c_ptr) :: resources
        type (cudaStream_t), value :: stream
        integer (KIND(cudaSuccess)) :: res
      end function cudaGraphicsMapResources
    end interface

    interface ! [['cudaError_t', None], 'cudaGraphicsUnmapResources', [['int', None, 'count'], ['struct', 'cudaGraphicsResource', '**', 'resources'], ['cudaStream_t', None, 'stream']]]
      function cudaGraphicsUnmapResources(count,resources,stream) result( res ) bind(C, name="cudaGraphicsUnmapResources")
        use, intrinsic :: ISO_C_BINDING
        import cudaSuccess
        import cudaStream_t
        implicit none
        integer(c_int), value :: count
        type (c_ptr) :: resources
        type (cudaStream_t), value :: stream
        integer (KIND(cudaSuccess)) :: res
      end function cudaGraphicsUnmapResources
    end interface

    interface ! [['cudaError_t', None], 'cudaGraphicsResourceGetMappedPointer', [['void', '**', 'devPtr'], ['size_t', '*', 'size'], ['struct', 'cudaGraphicsResource', '*', 'resource']]]
      function cudaGraphicsResourceGetMappedPointer(devPtr,size,resource) result( res ) bind(C, name="cudaGraphicsResourceGetMappedPointer")
        use, intrinsic :: ISO_C_BINDING
        import cudaSuccess
        implicit none
        type(c_ptr) :: devPtr
        integer(c_int) :: size
        type (c_ptr), value :: resource
        integer (KIND(cudaSuccess)) :: res
      end function cudaGraphicsResourceGetMappedPointer
    end interface

    interface ! [['cudaError_t', None], 'cudaGraphicsSubResourceGetMappedArray', [['struct', 'cudaArray', '**', 'arrayPtr'], ['struct', 'cudaGraphicsResource', '*', 'resource'], ['unsigned', 'int', None, 'arrayIndex'], ['unsigned', 'int', None, 'mipLevel']]]
      function cudaGraphicsSubResourceGetMappedArray(arrayPtr,resource,arrayIndex,mipLevel) result( res ) bind(C, name="cudaGraphicsSubResourceGetMappedArray")
        use, intrinsic :: ISO_C_BINDING
        import cudaSuccess
        implicit none
        type (c_ptr) :: arrayPtr
        type (c_ptr), value :: resource
        integer(c_int), value :: arrayIndex
        integer(c_int), value :: mipLevel
        integer (KIND(cudaSuccess)) :: res
      end function cudaGraphicsSubResourceGetMappedArray
    end interface

  end module cuda_runtime_h
