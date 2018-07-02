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

    enum, bind(C) !:: cudaError  driver_types.h
        enumerator  :: cudaSuccess=0
        enumerator  :: cudaErrorMissingConfiguration=1
        enumerator  :: cudaErrorMemoryAllocation=2
        enumerator  :: cudaErrorInitializationError=3
        enumerator  :: cudaErrorLaunchFailure=4
        enumerator  :: cudaErrorPriorLaunchFailure=5
        enumerator  :: cudaErrorLaunchTimeout=6
        enumerator  :: cudaErrorLaunchOutOfResources=7
        enumerator  :: cudaErrorInvalidDeviceFunction=8
        enumerator  :: cudaErrorInvalidConfiguration=9
        enumerator  :: cudaErrorInvalidDevice=10
        enumerator  :: cudaErrorInvalidValue=11
        enumerator  :: cudaErrorInvalidPitchValue=12
        enumerator  :: cudaErrorInvalidSymbol=13
        enumerator  :: cudaErrorMapBufferObjectFailed=14
        enumerator  :: cudaErrorUnmapBufferObjectFailed=15
        enumerator  :: cudaErrorInvalidHostPointer=16
        enumerator  :: cudaErrorInvalidDevicePointer=17
        enumerator  :: cudaErrorInvalidTexture=18
        enumerator  :: cudaErrorInvalidTextureBinding=19
        enumerator  :: cudaErrorInvalidChannelDescriptor=20
        enumerator  :: cudaErrorInvalidMemcpyDirection=21
        enumerator  :: cudaErrorAddressOfConstant=22
        enumerator  :: cudaErrorTextureFetchFailed=23
        enumerator  :: cudaErrorTextureNotBound=24
        enumerator  :: cudaErrorSynchronizationError=25
        enumerator  :: cudaErrorInvalidFilterSetting=26
        enumerator  :: cudaErrorInvalidNormSetting=27
        enumerator  :: cudaErrorMixedDeviceExecution=28
        enumerator  :: cudaErrorCudartUnloading=29
        enumerator  :: cudaErrorUnknown=30
        enumerator  :: cudaErrorNotYetImplemented=31
        enumerator  :: cudaErrorMemoryValueTooLarge=32
        enumerator  :: cudaErrorInvalidResourceHandle=33
        enumerator  :: cudaErrorNotReady=34
        enumerator  :: cudaErrorInsufficientDriver=35
        enumerator  :: cudaErrorSetOnActiveProcess=36
        enumerator  :: cudaErrorInvalidSurface=37
        enumerator  :: cudaErrorNoDevice=38
        enumerator  :: cudaErrorECCUncorrectable=39
        enumerator  :: cudaErrorSharedObjectSymbolNotFound=40
        enumerator  :: cudaErrorSharedObjectInitFailed=41
        enumerator  :: cudaErrorUnsupportedLimit=42
        enumerator  :: cudaErrorDuplicateVariableName=43
        enumerator  :: cudaErrorDuplicateTextureName=44
        enumerator  :: cudaErrorDuplicateSurfaceName=45
        enumerator  :: cudaErrorDevicesUnavailable=46
        enumerator  :: cudaErrorInvalidKernelImage           =     47
        enumerator  :: cudaErrorNoKernelImageForDevice       =     48
        enumerator  :: cudaErrorIncompatibleDriverContext    =     49
        enumerator  :: cudaErrorPeerAccessAlreadyEnabled     =     50
        enumerator  :: cudaErrorPeerAccessNotEnabled         =     51
        enumerator  :: cudaErrorDeviceAlreadyInUse           =     54
        enumerator  :: cudaErrorProfilerDisabled             =     55
        enumerator  :: cudaErrorProfilerNotInitialized       =     56
        enumerator  :: cudaErrorProfilerAlreadyStarted       =     57
        enumerator  :: cudaErrorProfilerAlreadyStopped       =    58
        enumerator  :: cudaErrorAssert                        =    59
        enumerator  :: cudaErrorTooManyPeers                 =     60
        enumerator  :: cudaErrorHostMemoryAlreadyRegistered  =     61
        enumerator  :: cudaErrorHostMemoryNotRegistered      =     62
        enumerator  :: cudaErrorOperatingSystem              =     63
        enumerator  :: cudaErrorPeerAccessUnsupported        =     64
        enumerator  :: cudaErrorLaunchMaxDepthExceeded       =     65
        enumerator  :: cudaErrorLaunchFileScopedTex          =     66
        enumerator  :: cudaErrorLaunchFileScopedSurf         =     67
        enumerator  :: cudaErrorSyncDepthExceeded            =     68
        enumerator  :: cudaErrorLaunchPendingCountExceeded   =     69
        enumerator  :: cudaErrorNotPermitted                 =     70
        enumerator  :: cudaErrorNotSupported                 =     71
        enumerator  :: cudaErrorHardwareStackError           =     72
        enumerator  :: cudaErrorIllegalInstruction           =     73
        enumerator  :: cudaErrorMisalignedAddress            =     74
        enumerator  :: cudaErrorInvalidAddressSpace          =     75
        enumerator  :: cudaErrorInvalidPc                    =     76
        enumerator  :: cudaErrorIllegalAddress               =     77
        enumerator  :: cudaErrorInvalidPtx                   =     78
        enumerator  :: cudaErrorInvalidGraphicsContext       =     79
        enumerator  :: cudaErrorNvlinkUncorrectable          =     80
        enumerator  :: cudaErrorJitCompilerNotFound          =     81
        enumerator  :: cudaErrorCooperativeLaunchTooLarge    =     82
        enumerator  :: cudaErrorStartupFailure=127   ! z'7f'  0x7f
        enumerator  :: cudaErrorApiFailureBase=10000
    end enum ! cudaError



    enum, bind(C) !:: cudaChannelFormatKind
        enumerator :: cudaChannelFormatKindSigned=0
        enumerator :: cudaChannelFormatKindUnsigned=1
        enumerator :: cudaChannelFormatKindFloat=2
        enumerator :: cudaChannelFormatKindNone=3
    end enum ! cudaChannelFormatKind

    enum, bind(C) !:: cudaMemoryType
        enumerator :: cudaMemoryTypeHost=1
        enumerator :: cudaMemoryTypeDevice=2
    end enum ! cudaMemoryType

    enum, bind(C) !:: cudaMemcpyKind
        enumerator :: cudaMemcpyHostToHost=0
        enumerator :: cudaMemcpyHostToDevice=1
        enumerator :: cudaMemcpyDeviceToHost=2
        enumerator :: cudaMemcpyDeviceToDevice=3
        enumerator :: cudaMemcpyDefault=4  ! Direction of the transfer is inferred from the pointer values. Requires unified virtual addressing
    end enum ! cudaMemcpyKind

    enum, bind(C) !:: cudaGraphicsRegisterFlags
        enumerator :: cudaGraphicsRegisterFlagsNone=0
        enumerator :: cudaGraphicsRegisterFlagsReadOnly       = 1 !!  CUDA will not write to this resource
        enumerator :: cudaGraphicsRegisterFlagsWriteDiscard     = 2  !!  CUDA will only write to and will not read from this resource
        enumerator :: cudaGraphicsRegisterFlagsSurfaceLoadStore = 4  !! CUDA will bind this resource to a surface reference
        enumerator :: cudaGraphicsRegisterFlagsTextureGather    = 8  !! CUDA will perform texture gather operations on this resource
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

    enum, bind(C) !:: cudaFuncAttribute
        enumerator ::  cudaFuncAttributeMaxDynamicSharedMemorySize = 8!!  Maximum dynamic shared memory size
        enumerator ::  cudaFuncAttributePreferredSharedMemoryCarveout = 9 !! Preferred shared memory-L1 cache split ratio
        enumerator ::  cudaFuncAttributeMax
    end enum
    enum, bind(C) !:: cudaFuncCache
        enumerator :: cudaFuncCachePreferNone=0
        enumerator :: cudaFuncCachePreferShared=1
        enumerator :: cudaFuncCachePreferL1=2
    end enum ! cudaFuncCache
    enum, bind(C) !::cudaSharedMemConfig
        enumerator ::  cudaSharedMemBankSizeDefault   = 0
        enumerator ::  cudaSharedMemBankSizeFourByte  = 1
        enumerator ::  cudaSharedMemBankSizeEightByte = 2
    end enum
    enum, bind(C) !:: cudaComputeMode
        enumerator :: cudaComputeModeDefault=0
        enumerator :: cudaComputeModeExclusive=1
        enumerator :: cudaComputeModeProhibited=2
        enumerator ::cudaComputeModeExclusiveProcess=3
    end enum ! cudaComputeMode

    enum, bind(C) !:: cudaLimit
        enumerator :: cudaLimitStackSize=0
        enumerator :: cudaLimitPrintfFifoSize=1
        enumerator :: cudaLimitMallocHeapSize               = 2 !! GPU malloc heap size
        enumerator :: cudaLimitDevRuntimeSyncDepth          = 3 !! GPU device runtime synchronize depth
        enumerator :: cudaLimitDevRuntimePendingLaunchCount = 4 !!  GPU device runtime pending launch count
    end enum ! cudaLimit
    enum, bind(C) !:: cudaMemoryAdvise
        enumerator ::  cudaMemAdviseSetReadMostly          = 1 !!  Data will mostly be read and only occassionally be written to
        enumerator :: cudaMemAdviseUnsetReadMostly        = 2 !! Undo the effect of ::cudaMemAdviseSetReadMostly
        enumerator :: cudaMemAdviseSetPreferredLocation   = 3 !! Set the preferred location for the data as the specified device
        enumerator ::  cudaMemAdviseUnsetPreferredLocation = 4 !! Clear the preferred location for the data
        enumerator ::  cudaMemAdviseSetAccessedBy          = 5 !! Data will be accessed by the specified device, so prevent page faults as much as possible
        enumerator :: cudaMemAdviseUnsetAccessedBy        = 6 !! Let the Unified Memory subsystem decide on the page faulting policy for the specified device
    end enum ! cudaMemoryAdvise
    enum, bind(C) !:: cudaMemRangeAttribute
        enumerator ::     cudaMemRangeAttributeReadMostly           = 1  !! Whether the range will mostly be read and only occassionally be written to
        enumerator ::     cudaMemRangeAttributePreferredLocation    = 2  !! The preferred location of the range
        enumerator ::     cudaMemRangeAttributeAccessedBy           = 3  !! Memory range has ::cudaMemAdviseSetAccessedBy set for specified device
        enumerator ::     cudaMemRangeAttributeLastPrefetchLocation = 4  !! The last location to which the range was prefetched
    end enum ! !:: cudaMemRangeAttribute

    enum, bind(C) !::cudaOutputMode
        enumerator :: cudaKeyValuePair =0
        enumerator :: cudaCSV
    end enum !

    !! CUDA device attributes
    enum, bind(C)  ! cudaDeviceAttr
        enumerator :: cudaDeviceAttr=0
        enumerator :: cudaDevAttrMaxThreadsPerBlock             = 1   !!   /**< Maximum number of threads per block */
        enumerator :: cudaDevAttrMaxBlockDimX                   = 2   !!   /**< Maximum block dimension X */
        enumerator :: cudaDevAttrMaxBlockDimY                   = 3   !!   /**< Maximum block dimension Y */
        enumerator :: cudaDevAttrMaxBlockDimZ                   = 4   !!   /**< Maximum block dimension Z */
        enumerator :: cudaDevAttrMaxGridDimX                    = 5   !!   /**< Maximum grid dimension X */
        enumerator :: cudaDevAttrMaxGridDimY                    = 6   !!   /**< Maximum grid dimension Y */
        enumerator :: cudaDevAttrMaxGridDimZ                    = 7   !!   /**< Maximum grid dimension Z */
        enumerator :: cudaDevAttrMaxSharedMemoryPerBlock        = 8   !!   /**< Maximum shared memory available per block in bytes */
        enumerator :: cudaDevAttrTotalConstantMemory            = 9   !!   /**< Memory available on device for __constant__ variables in a CUDA C kernel in bytes */
        enumerator :: cudaDevAttrWarpSize                       = 10   !!  /**< Warp size in threads */
        enumerator :: cudaDevAttrMaxPitch                       = 11   !!  /**< Maximum pitch in bytes allowed by memory copies */
        enumerator :: cudaDevAttrMaxRegistersPerBlock           = 12   !!  /**< Maximum number of 32-bit registers available per block */
        enumerator :: cudaDevAttrClockRate                      = 13   !!  /**< Peak clock frequency in kilohertz */
        enumerator :: cudaDevAttrTextureAlignment               = 14   !!  /**< Alignment requirement for textures */
        enumerator :: cudaDevAttrGpuOverlap                     = 15   !!  /**< Device can possibly copy memory and execute a kernel concurrently */
        enumerator :: cudaDevAttrMultiProcessorCount            = 16   !!  /**< Number of multiprocessors on device */
        enumerator :: cudaDevAttrKernelExecTimeout              = 17   !!  /**< Specifies whether there is a run time limit on kernels */
        enumerator :: cudaDevAttrIntegrated                     = 18   !!  /**< Device is integrated with host memory */
        enumerator :: cudaDevAttrCanMapHostMemory               = 19   !!  /**< Device can map host memory into CUDA address space */
        enumerator :: cudaDevAttrComputeMode                    = 20   !!  /**< Compute mode (See ::cudaComputeMode for details) */
        enumerator :: cudaDevAttrMaxTexture1DWidth              = 21   !!  /**< Maximum 1D texture width */
        enumerator :: cudaDevAttrMaxTexture2DWidth              = 22   !!  /**< Maximum 2D texture width */
        enumerator :: cudaDevAttrMaxTexture2DHeight             = 23   !!  /**< Maximum 2D texture height */
        enumerator :: cudaDevAttrMaxTexture3DWidth              = 24   !!  /**< Maximum 3D texture width */
        enumerator :: cudaDevAttrMaxTexture3DHeight             = 25   !!  /**< Maximum 3D texture height */
        enumerator :: cudaDevAttrMaxTexture3DDepth              = 26   !!  /**< Maximum 3D texture depth */
        enumerator :: cudaDevAttrMaxTexture2DLayeredWidth       = 27   !!  /**< Maximum 2D layered texture width */
        enumerator :: cudaDevAttrMaxTexture2DLayeredHeight      = 28   !!  /**< Maximum 2D layered texture height */
        enumerator :: cudaDevAttrMaxTexture2DLayeredLayers      = 29   !!  /**< Maximum layers in a 2D layered texture */
        enumerator :: cudaDevAttrSurfaceAlignment               = 30   !!  /**< Alignment requirement for surfaces */
        enumerator :: cudaDevAttrConcurrentKernels              = 31   !!  /**< Device can possibly execute multiple kernels concurrently */
        enumerator :: cudaDevAttrEccEnabled                     = 32   !!  /**< Device has ECC support enabled */
        enumerator :: cudaDevAttrPciBusId                       = 33   !!  /**< PCI bus ID of the device */
        enumerator :: cudaDevAttrPciDeviceId                    = 34   !!  /**< PCI device ID of the device */
        enumerator :: cudaDevAttrTccDriver                      = 35   !!  /**< Device is using TCC driver model */
        enumerator :: cudaDevAttrMemoryClockRate                = 36   !!  /**< Peak memory clock frequency in kilohertz */
        enumerator :: cudaDevAttrGlobalMemoryBusWidth           = 37   !!  /**< Global memory bus width in bits */
        enumerator :: cudaDevAttrL2CacheSize                    = 38   !!  /**< Size of L2 cache in bytes */
        enumerator :: cudaDevAttrMaxThreadsPerMultiProcessor    = 39   !!  /**< Maximum resident threads per multiprocessor */
        enumerator :: cudaDevAttrAsyncEngineCount               = 40   !!  /**< Number of asynchronous engines */
        enumerator :: cudaDevAttrUnifiedAddressing              = 41   !!  /**< Device shares a unified address space with the host */
        enumerator :: cudaDevAttrMaxTexture1DLayeredWidth       = 42   !!  /**< Maximum 1D layered texture width */
        enumerator :: cudaDevAttrMaxTexture1DLayeredLayers      = 43   !!  /**< Maximum layers in a 1D layered texture */
        enumerator :: cudaDevAttrMaxTexture2DGatherWidth        = 45   !!  /**< Maximum 2D texture width if cudaArrayTextureGather is set */
        enumerator :: cudaDevAttrMaxTexture2DGatherHeight       = 46   !!  /**< Maximum 2D texture height if cudaArrayTextureGather is set */
        enumerator :: cudaDevAttrMaxTexture3DWidthAlt           = 47   !!  /**< Alternate maximum 3D texture width */
        enumerator :: cudaDevAttrMaxTexture3DHeightAlt          = 48   !!  /**< Alternate maximum 3D texture height */
        enumerator :: cudaDevAttrMaxTexture3DDepthAlt           = 49   !!  /**< Alternate maximum 3D texture depth */
        enumerator :: cudaDevAttrPciDomainId                    = 50   !!  /**< PCI domain ID of the device */
        enumerator :: cudaDevAttrTexturePitchAlignment          = 51   !!  /**< Pitch alignment requirement for textures */
        enumerator :: cudaDevAttrMaxTextureCubemapWidth         = 52   !!  /**< Maximum cubemap texture width/height */
        enumerator :: cudaDevAttrMaxTextureCubemapLayeredWidth  = 53   !!  /**< Maximum cubemap layered texture width/height */
        enumerator :: cudaDevAttrMaxTextureCubemapLayeredLayers = 54   !!  /**< Maximum layers in a cubemap layered texture */
        enumerator :: cudaDevAttrMaxSurface1DWidth              = 55   !!  /**< Maximum 1D surface width */
        enumerator :: cudaDevAttrMaxSurface2DWidth              = 56   !!  /**< Maximum 2D surface width */
        enumerator :: cudaDevAttrMaxSurface2DHeight             = 57   !!  /**< Maximum 2D surface height */
        enumerator :: cudaDevAttrMaxSurface3DWidth              = 58   !!  /**< Maximum 3D surface width */
        enumerator :: cudaDevAttrMaxSurface3DHeight             = 59   !!  /**< Maximum 3D surface height */
        enumerator :: cudaDevAttrMaxSurface3DDepth              = 60   !!  /**< Maximum 3D surface depth */
        enumerator :: cudaDevAttrMaxSurface1DLayeredWidth       = 61   !!  /**< Maximum 1D layered surface width */
        enumerator :: cudaDevAttrMaxSurface1DLayeredLayers      = 62   !!  /**< Maximum layers in a 1D layered surface */
        enumerator :: cudaDevAttrMaxSurface2DLayeredWidth       = 63   !!  /**< Maximum 2D layered surface width */
        enumerator :: cudaDevAttrMaxSurface2DLayeredHeight      = 64   !!  /**< Maximum 2D layered surface height */
        enumerator :: cudaDevAttrMaxSurface2DLayeredLayers      = 65   !!  /**< Maximum layers in a 2D layered surface */
        enumerator :: cudaDevAttrMaxSurfaceCubemapWidth         = 66   !!  /**< Maximum cubemap surface width */
        enumerator :: cudaDevAttrMaxSurfaceCubemapLayeredWidth  = 67   !!  /**< Maximum cubemap layered surface width */
        enumerator :: cudaDevAttrMaxSurfaceCubemapLayeredLayers = 68   !!  /**< Maximum layers in a cubemap layered surface */
        enumerator :: cudaDevAttrMaxTexture1DLinearWidth        = 69   !!  /**< Maximum 1D linear texture width */
        enumerator :: cudaDevAttrMaxTexture2DLinearWidth        = 70   !!  /**< Maximum 2D linear texture width */
        enumerator :: cudaDevAttrMaxTexture2DLinearHeight       = 71   !!  /**< Maximum 2D linear texture height */
        enumerator :: cudaDevAttrMaxTexture2DLinearPitch        = 72   !!  /**< Maximum 2D linear texture pitch in bytes */
        enumerator :: cudaDevAttrMaxTexture2DMipmappedWidth     = 73   !!  /**< Maximum mipmapped 2D texture width */
        enumerator :: cudaDevAttrMaxTexture2DMipmappedHeight    = 74   !!  /**< Maximum mipmapped 2D texture height */
        enumerator :: cudaDevAttrComputeCapabilityMajor         = 75   !!  /**< Major compute capability version number */
        enumerator :: cudaDevAttrComputeCapabilityMinor         = 76   !!  /**< Minor compute capability version number */
        enumerator :: cudaDevAttrMaxTexture1DMipmappedWidth     = 77   !!  /**< Maximum mipmapped 1D texture width */
        enumerator :: cudaDevAttrStreamPrioritiesSupported      = 78   !!  /**< Device supports stream priorities */
        enumerator :: cudaDevAttrGlobalL1CacheSupported         = 79   !!  /**< Device supports caching globals in L1 */
        enumerator :: cudaDevAttrLocalL1CacheSupported          = 80   !!  /**< Device supports caching locals in L1 */
        enumerator :: cudaDevAttrMaxSharedMemoryPerMultiprocessor = 81   !!  /**< Maximum shared memory available per multiprocessor in bytes */
        enumerator :: cudaDevAttrMaxRegistersPerMultiprocessor  = 82   !!  /**< Maximum number of 32-bit registers available per multiprocessor */
        enumerator :: cudaDevAttrManagedMemory                  = 83   !!  /**< Device can allocate managed memory on this system */
        enumerator :: cudaDevAttrIsMultiGpuBoard                = 84   !!  /**< Device is on a multi-GPU board */
        enumerator :: cudaDevAttrMultiGpuBoardGroupID           = 85   !!  /**< Unique identifier for a group of devices on the same multi-GPU board */
        enumerator :: cudaDevAttrHostNativeAtomicSupported      = 86   !!  /**< Link between the device and the host supports native atomic operations */
        enumerator :: cudaDevAttrSingleToDoublePrecisionPerfRatio = 87   !!  /**< Ratio of single precision performance (in floating-point operations per second) to double precision performance */
        enumerator :: cudaDevAttrPageableMemoryAccess           = 88   !!  /**< Device supports coherently accessing pageable memory without calling cudaHostRegister on it */
        enumerator :: cudaDevAttrConcurrentManagedAccess        = 89   !!  /**< Device can coherently access managed memory concurrently with the CPU */
        enumerator :: cudaDevAttrComputePreemptionSupported     = 90   !!  /**< Device supports Compute Preemption */
        enumerator :: cudaDevAttrCanUseHostPointerForRegisteredMem = 91   !!  /**< Device can access host registered memory at the same virtual address as the CPU */
        enumerator :: cudaDevAttrReserved92                     = 92   !!
        enumerator :: cudaDevAttrReserved93                     = 93   !!
        enumerator :: cudaDevAttrReserved94                     = 94   !!
        enumerator :: cudaDevAttrCooperativeLaunch              = 95   !!  /**< Device supports launching cooperative kernels via ::cudaLaunchCooperativeKernel*/
        enumerator :: cudaDevAttrCooperativeMultiDeviceLaunch   = 96   !!  /**< Device can participate in cooperative kernels launched via ::cudaLaunchCooperativeKernelMultiDevice */
        enumerator :: cudaDevAttrMaxSharedMemoryPerBlockOptin   = 97   !!  /**< The maximum optin shared memory per block. This value may vary by chip. See ::cudaFuncSetAttribute */
        enumerator :: cudaDevAttrCanFlushRemoteWrites           = 98   !!  /**< Device supports flushing of outstanding remote writes. */
        enumerator :: cudaDevAttrHostRegisterSupported          = 99   !!  /**< Device supports host memory registration via ::cudaHostRegister. */
        enumerator :: cudaDevAttrPageableMemoryAccessUsesHostPageTables = 100   !!  /**< Device accesses pageable memory via the host's page tables. */
        enumerator :: cudaDevAttrDirectManagedMemAccessFromHost = 101 !! /**< Host can directly access managed memory on the device without migration. */
    end enum ! cudaDeviceAttr

    enum, bind(C) !::cudaDeviceP2PAttr  CUDA device P2P attributes
        enumerator :: cudaDevP2PAttrPerformanceRank              = 1!! A relative value indicating the performance of the link between two devices
        enumerator ::  cudaDevP2PAttrAccessSupported              = 2!! Peer access is enabled
        enumerator :: cudaDevP2PAttrNativeAtomicSupported        = 3!! Native atomic operation over the link supported
        enumerator :: cudaDevP2PAttrCudaArrayAccessSupported     = 4!! Accessing CUDA arrays over the link supported
    end enum

    !! CUDA cooperative group scope
    enum, bind(C) !! :: cudaCGScope
        enumerator ::  cudaCGScopeInvalid   = 0 !! Invalid cooperative group scope */
        enumerator ::  cudaCGScopeGrid      = 1 !! Scope represented by a grid_group */
        enumerator :: cudaCGScopeMultiGrid = 2 !! Scope represented by a multi_grid_group */
    end enum


    enum, bind(C) !:: cudaSurfaceBoundaryMode
        enumerator :: cudaBoundaryModeZero=0
        enumerator :: cudaBoundaryModeClamp=1
        enumerator :: cudaBoundaryModeTrap=2
    end enum ! cudaSurfaceBoundaryMode

    enum, bind(C) !:: cudaSurfaceFormatMode
        enumerator :: cudaFormatModeForced=0
        enumerator :: cudaFormatModeAuto=1
    end enum ! cudaSurfaceFormatMode

    enum, bind(C) !:: cudaTextureAddressMode
        enumerator :: cudaTextureAddressMode=1000
        enumerator :: cudaAddressModeWrap=0
        enumerator :: cudaAddressModeClamp=1
        enumerator :: cudaAddressModeMirror=2
        enumerator :: cudaAddressModeBorder = 3
    end enum ! cudaTextureAddressMode

    enum, bind(C) !:: cudaTextureFilterMode
        enumerator :: cudaTextureFilterMode=1000
        enumerator :: cudaFilterModePoint=0
        enumerator :: cudaFilterModeLinear=1
    end enum ! cudaTextureFilterMode

    enum, bind(C) !:: cudaTextureReadMode
         enumerator :: cudaTextureReadMode=-1
        enumerator :: cudaReadModeElementType=0
        enumerator :: cudaReadModeNormalizedFloat=1
    end enum ! cudaTextureReadMode


    type, bind(C) :: cudaChannelFormatDesc
        integer(c_int) :: x
        integer(c_int) :: y
        integer(c_int) :: z
        integer(c_int) :: w
        integer (KIND(cudaChannelFormatKindSigned)) :: f
    end type cudaChannelFormatDesc

    type, bind(C) :: textureReference
        integer(c_int)                            :: normalized
        integer (KIND(cudaFilterModePoint))       :: filterMode
        integer (KIND(cudaTextureAddressMode))    :: addressMode(3);
        type(cudaChannelFormatDesc)               :: channelDesc;
        integer(c_int)                            :: sRGB;
        integer(c_int)                            :: maxAnisotropy;  !! unsigned int
        integer (KIND(cudaTextureFilterMode))     :: mipmapFilterMode;
        real(c_float)                             :: mipmapLevelBias;
        real(c_float)                             :: minMipmapLevelClamp;
        real(c_float)                             :: maxMipmapLevelClamp;
        integer(c_int) ,dimension(15)             :: cudaReserved;
    end type textureReference

    type, bind(C) :: cudaTextureDesc
        integer (KIND(cudaTextureAddressMode)) :: addressMode(3);
        integer (KIND(cudaTextureFilterMode) ) :: filterMode;
        integer (KIND(cudaTextureReadMode))    :: readMode;
        integer(c_int)                         :: sRGB;
        real(c_float)                          :: borderColor(4);
        integer(c_int)                         :: normalizedCoords;
        integer(c_int)                         :: maxAnisotropy;
        integer (KIND(cudaTextureFilterMode))  :: mipmapFilterMode;
        real(c_float)                          :: mipmapLevelBias;
        real(c_float )                         :: minMipmapLevelClamp;
        real(c_float)                          :: maxMipmapLevelClamp;
    end type cudaTextureDesc


   type, bind(C) :: surfaceReference
        type (cudaChannelFormatDesc) :: channelDesc
    end type surfaceReference

    type, bind(C) :: cudaPitchedPtr
        type(c_ptr) :: ptr
        integer(c_size_t) :: pitch
        integer(c_size_t) :: xsize
        integer(c_size_t) :: ysize
    end type cudaPitchedPtr

    type, bind(C) :: cudaExtent
        integer(c_size_t) :: width
        integer(c_size_t) :: height
        integer(c_size_t) :: depth
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

   type, bind(C) :: cudaMemcpy3DPeersParms
        type (c_ptr) :: srcArray
        type (cudaPos) :: srcPos
        type (cudaPitchedPtr) :: srcPtr
        integer(c_int) :: srcDevice
        type (c_ptr) :: dstArray
        type (cudaPos) :: dstPos
        type (cudaPitchedPtr) :: dstPtr
        integer(c_int) :: dstDevice
        type (cudaExtent) :: extent
    end type cudaMemcpy3DPeersParms

    type, bind(C) :: cudaPointerAttributes
        integer (KIND(cudaMemoryTypeHost)) :: memoryType
        integer(c_int) :: device
        type (c_ptr) ::devicePointer
        type (c_ptr) :: hostPointer
        integer(c_int) :: isManaged
    end type cudaPointerAttributes


    type, bind(C) :: cudaFuncAttributes
        integer(c_size_t) :: sharedSizeBytes
        integer(c_size_t) :: constSizeBytes
        integer(c_size_t) :: localSizeBytes
        integer(c_int) :: maxThreadsPerBlock
        integer(c_int) :: numRegs
        integer(c_int) :: ptxVersion
        integer(c_int) :: binaryVersion
        integer(c_int) :: cacheModeCA
        integer(c_int) :: maxDynamicSharedSizeBytes
        integer(c_int) :: preferredShmemCarveout
    end type cudaFuncAttributes

    type, bind(C) :: cudaDeviceProp
        character(c_char) :: name(256)
        integer(c_size_t) :: totalGlobalMem            ! amount of global memory available on the device in bytes
        integer(c_size_t) :: sharedMemPerBlock         ! maximum amount of shared memory available to a thread block in bytes
        integer(c_int) :: regsPerBlock                 ! maximum number of 32-bit registers available to a thread block
        integer(c_int) :: warpSize                     ! warp size in threads
        integer(c_size_t) :: memPitch! maximum pitch in  bytes allowed by the memory copy functions that involve memory regions allocated through ::cudaMallocPitch()
        integer(c_int) :: maxThreadsPerBlock! maximum number of threads per block
        integer(c_int) :: maxThreadsDim(3)! maximum size of each dimension of a block
        integer(c_int) :: maxGridSize(3)! maximum size of each dimension of a grid
        integer(c_int) :: clockRate ! clock frequency in   kilohertz
        integer(c_size_t) :: totalConstMem
        integer(c_int) :: major
        integer(c_int) :: minor
        integer(c_size_t) :: textureAlignment
        integer(c_size_t) :: texturePitchAlignment
        integer(c_int) :: deviceOverlap
        integer(c_int) :: multiProcessorCount
        integer(c_int) :: kernelExecTimeoutEnabled
        integer(c_int) :: integrated
        integer(c_int) :: canMapHostMemory
        integer(c_int) :: computeMode
        integer(c_int) :: maxTexture1D
        integer(c_int) :: maxTexture1DMipmap;
        integer(c_int) :: maxTexture1DLinear;
        integer(c_int) :: maxTexture2D(2)
        integer(c_int) :: maxTexture2DMipmap(2);
        integer(c_int) :: maxTexture2DLinear(3);
        integer(c_int) :: maxTexture2DGather(2);
        integer(c_int) :: maxTexture3D(3)
        integer(c_int) :: maxTexture3DAlt(3);
        !LEGACY   integer(c_int) :: maxTexture2DArray(3)
        integer(c_int) :: maxTextureCubemap;
        integer(c_int) :: maxTexture1DLayered(2);
        integer(c_int) :: maxTexture2DLayered(3);
        integer(c_int) :: maxTextureCubemapLayered(2)
        integer(c_int) :: maxSurface1D;
        integer(c_int) :: maxSurface2D(2)
        integer(c_int) :: maxSurface3D(3)
        integer(c_int) :: maxSurface1DLayered(2)
        integer(c_int) :: maxSurface2DLayered(3)
        integer(c_size_t) :: surfaceAlignment
        integer(c_int) :: concurrentKernels
        integer(c_int) :: ECCEnabled
        integer(c_int) :: pciBusID
        integer(c_int) :: pciDeviceID
        integer(c_int) :: pciDomainID;
        integer(c_int) :: tccDriver;
        integer(c_int) :: asyncEngineCount;! is 1 when the   device can concurrently copy memory between host and device while executing   a kernel. It is 2 when the device can concurrently copy memory between host and device in both directions and execute a kernel at the same time. It is 0 if neither of these is supported.
        integer(c_int) :: unifiedAddressing;
        integer(c_int) :: memoryClockRate;!the peak memory clock frequency in kilohertz
        integer(c_int) :: memoryBusWidth;
        integer(c_int) :: l2CacheSize;
        integer(c_int) :: maxThreadsPerMultiProcessor;
        integer(c_int) :: streamPrioritiesSupported;
        integer(c_int) :: globalL1CacheSupported; ! is 1 if the device supports caching of globals in L1 cache, or 0 if it is not supported
        integer(c_int) :: localL1CacheSupported;
        integer(c_size_t) :: sharedMemPerMultiprocessor;
        integer(c_int) :: regsPerMultiprocessor;
        integer(c_int) :: managedMemory; !! in CUDA 8.0 changed from managedMemorySupported
        integer(c_int) :: isMultiGpuBoard;
        integer(c_int) :: multiGpuBoardGroupID;
        integer(c_int) :: singleToDoublePrecisionPerfRatio;
        integer(c_int) :: pageableMemoryAccess;
        integer(c_int) :: concurrentManagedAccess;
        integer(c_int) :: computePreemptionSupported;
        integer(c_int) :: canUseHostPointerForRegisteredMem;
        integer(c_int) :: cooperativeLaunch;
        integer(c_int) :: cooperativeMultiDeviceLaunch;
        integer(c_int) :: pageableMemoryAccessUsesHostPageTables;
        integer(c_int) :: directManagedMemAccessFromHost;
        !LEGACY   integer(c_int) :: cudaReserved__(22)
    end type cudaDeviceProp

    !! Vector_types.h
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


    !! CUDA launch parameters
    type, bind(C) ::  cudaLaunchParams
        type(c_ptr) ::func ! Device function symbol, void *
        type( dim3) :: gridDim
        type( dim3) :: blockDim
        type(c_ptr) :: args  !! void**
        integer(c_size_t) :: sharedMem
        type(cudaStream_t) stream
    end type cudaLaunchParams


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

    interface ! [['cudaError_t', None], 'cudaFuncGetAttributes', [['struct', 'cudaDevAttributes', '*', 'attr'], ['const', 'char', '*', 'func']]]
        function cudaGetDeviceAttributes(retvalue,attr,device) result( res ) bind(C, name="cudaGetDeviceAttributes")
            use, intrinsic :: ISO_C_BINDING
            import cudaSuccess
            import cudaDeviceAttr
            implicit none
            integer(c_int) :: retvalue
            integer (KIND(cudaDeviceAttr)) :: attr
            integer(c_int) :: device
            integer (KIND(cudaSuccess)) :: res
        end function cudaGetDeviceAttributes
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
