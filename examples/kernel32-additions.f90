module kernel32_additions
    use IFWINTY
    use, intrinsic :: ISO_C_BINDING, only : C_INT
    implicit none
    
    
    enum, bind(C) ! LOGICAL_PROCESSOR_RELATIONSHIP
        enumerator :: RelationProcessorCore
        enumerator :: RelationNumaNode
        enumerator :: RelationCache
        enumerator :: RelationProcessorPackage
        enumerator :: RelationGroup
        enumerator :: RelationAll = INT(Z'ffff',C_INT)
    end enum
    
    enum, bind(C) ! PROCESSOR_CACHE_TYPE
        enumerator :: CacheUnified
        enumerator :: CacheInstruction
        enumerator :: CacheData
        enumerator :: CacheTrace
    end enum
    
    type :: T_CACHE_DESCRIPTOR
        sequence
        integer(BYTE) :: Level
        integer(BYTE) :: Associativity
        integer(WORD) :: LineSize
        integer(DWORD) :: Size
        integer(C_INT) :: Type ! Enumerator PROCESSOR_CACHE_TYPE
    end type
    
    ! This has to be declared as a STRUCTURE rather than TYPE
    ! because of the union.
    structure /T_SYSTEM_LOGICAL_PROCESSOR_INFORMATION/
        integer(ULONG_PTR) :: ProcessorMask
        integer(C_INT) :: Relationship ! Enumerator LOGICAL_PROCESSOR_RELATIONSHIP
        union
            map
                integer(BYTE) :: Flags 
            end map
            map
                integer(DWORD) :: NodeNumber
            end map
            map
                type(T_CACHE_DESCRIPTOR) :: Cache
            end map
            map
                integer(ULONGLONG) :: Reserved(2)
            end map
        end union
    end structure
    
    
    ! Functions
    interface
    
    function GetLogicalProcessorInformation (Buffer, ReturnLength)
    !DEC$ ATTRIBUTES STDCALL, REFERENCE, DECORATE, ALIAS:"GetLogicalProcessorInformation" :: GetLogicalProcessorInformation
    import
    integer(BOOL) :: GetLogicalProcessorInformation
    type(T_SYSTEM_LOGICAL_PROCESSOR_INFORMATION), DIMENSION(*), intent(OUT) :: Buffer
    integer(DWORD), intent(INOUT) :: ReturnLength
    end function GetLogicalProcessorInformation
    
    function GetLogicalProcessorInformationEx (RelationshipType, Buffer, ReturnLength)
    !DEC$ ATTRIBUTES STDCALL, REFERENCE, DECORATE, ALIAS:"GetLogicalProcessorInformationEx" :: GetLogicalProcessorInformationEx
    import
    integer(BOOL) :: GetLogicalProcessorInformationEx
    integer(C_INT), intent(in), value :: RelationshipType ! Enumerator LOGICAL_PROCESSOR_RELATIONSHIP
    type(T_SYSTEM_LOGICAL_PROCESSOR_INFORMATION), DIMENSION(*), intent(OUT) :: Buffer
    integer(DWORD), intent(INOUT) :: ReturnLength
    end function GetLogicalProcessorInformationEx
    
    end interface
    
    end module kernel32_additions