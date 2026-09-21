!@descr: versioned binary persistence and transaction primitives for the flex per-state weight files
!!
!! One committed `flex_weights_state_NNN.bin` per delivered flex_pca state, each holding that
!! state's per-particle weight over the FULL physical particle layout of the project (rows outside
!! the run's selection are zero) and a flag marking the particles hard-labelled to it, so a
!! workflow that selects one state has everything about it in one file. The files of one delivery
!! share generation, layout digest and `nstates`; column s of a run, label value s and `vol_flex`
!! state s in the same project's out segment describe the same state.
!!
!! The format follows the canonical sigma2 store (simple_sigma2_state_file): 16-byte magic + 32
!! int64 header words with a header checksum, fixed section offsets, an immutable committed file,
!! a candidate that is validated and atomically published, and worker range files for a future
!! distributed writer. As in the sigma2 store since 2026-09-21 there are no per-row checksums:
!! the particle sections are written contiguously; only the small per-state scalar section
!! carries a checksum.
!!
!! Layout: header(512) | scalars(real64, NFIELDS+ncomp) | weights(real32, nptcls) | flags(int32, nptcls)
module simple_flex_weights_file
use, intrinsic :: iso_fortran_env, only: int8, int32, int64, real32, real64
use simple_fileio, only: get_fpath
use simple_string, only: string
use simple_syslib, only: file_exists, simple_sync_file, simple_sync_dir, simple_atomic_replace
implicit none
private

public :: flex_weights_header
public :: flex_weights_init_header, flex_weights_read_header, flex_weights_create_candidate
public :: flex_weights_write_rows, flex_weights_read_rows
public :: flex_weights_write_scalars, flex_weights_read_scalars
public :: flex_weights_write_local_range, flex_weights_read_local_range
public :: flex_weights_validate_file, flex_weights_publish
public :: FLEX_WEIGHTS_FBODY, FLEX_WEIGHTS_EXT
public :: FLEX_WEIGHTS_KIND_PARTITION, FLEX_WEIGHTS_KIND_KERNEL
public :: FLEX_WEIGHTS_PROV_FLEX_PCA, FLEX_WEIGHTS_PROV_MERGED, FLEX_WEIGHTS_PROV_EXTERNAL
public :: FLEX_WEIGHTS_CANDIDATE, FLEX_WEIGHTS_COMMITTED
public :: FLEX_WEIGHTS_SCALAR_MASS, FLEX_WEIGHTS_SCALAR_NEFF, FLEX_WEIGHTS_SCALAR_POP
public :: FLEX_WEIGHTS_SCALAR_BW, FLEX_WEIGHTS_SCALAR_NFIELDS

character(len=*),  parameter :: FLEX_WEIGHTS_FBODY = 'flex_weights_state_'
character(len=*),  parameter :: FLEX_WEIGHTS_EXT   = '.bin'
character(len=16), parameter :: STATE_MAGIC = 'SIMPLE_FLEXW_V02'
character(len=16), parameter :: RANGE_MAGIC = 'SIMPLE_FXW_RNG_2'

integer(int32), parameter :: FLEX_WEIGHTS_VERSION = 2_int32
!> across the files of one delivery every assigned row sums to one: responsibilities, regions, hard labels
integer(int32), parameter :: FLEX_WEIGHTS_KIND_PARTITION = 1_int32
!> compact-support kernel coefficients: each in [0,1], rows need not sum to one across files
integer(int32), parameter :: FLEX_WEIGHTS_KIND_KERNEL    = 2_int32
integer(int32), parameter :: FLEX_WEIGHTS_PROV_FLEX_PCA = 1_int32
integer(int32), parameter :: FLEX_WEIGHTS_PROV_MERGED   = 2_int32
integer(int32), parameter :: FLEX_WEIGHTS_PROV_EXTERNAL = 3_int32
integer(int32), parameter :: FLEX_WEIGHTS_CANDIDATE = 1_int32
integer(int32), parameter :: FLEX_WEIGHTS_COMMITTED = 2_int32
!> scalar section fields; the latent target occupies fields NFIELDS+1 .. NFIELDS+ncomp
integer, parameter :: FLEX_WEIGHTS_SCALAR_MASS    = 1
integer, parameter :: FLEX_WEIGHTS_SCALAR_NEFF    = 2
integer, parameter :: FLEX_WEIGHTS_SCALAR_POP     = 3
integer, parameter :: FLEX_WEIGHTS_SCALAR_BW      = 4
integer, parameter :: FLEX_WEIGHTS_SCALAR_NFIELDS = 4

integer(int64), parameter :: STATE_HEADER_BYTES = 512_int64
integer(int64), parameter :: RANGE_HEADER_BYTES = 256_int64
integer, parameter :: STATE_NWORDS = 32
integer, parameter :: RANGE_NWORDS = 16
integer(int64), parameter :: FNV_OFFSET = int(z'CBF29CE484222325', int64)
integer(int64), parameter :: FNV_PRIME  = int(z'00000100000001B3', int64)

integer, parameter :: W_VERSION=1, W_HEADER_BYTES=2, W_REAL_BYTES=3, W_NPTCLS=4, W_NSTATES=5
integer, parameter :: W_STATE_INDEX=6, W_NCOMP=7, W_KIND=8, W_PROVENANCE=9, W_GENERATION=10
integer, parameter :: W_STATE=11, W_LAYOUT_DIGEST=12, W_BOX=13, W_SMPD_BITS=14, W_BOX_CROP=15
integer, parameter :: W_SMPD_CROP_BITS=16, W_SCALAR_OFFSET=17, W_WEIGHT_OFFSET=18, W_FLAG_OFFSET=19
integer, parameter :: W_FILE_BYTES=20, W_SCALAR_CHECKSUM=21, W_HEADER_CHECKSUM=32

integer, parameter :: RW_VERSION=1, RW_HEADER_BYTES=2, RW_STATE_INDEX=3, RW_FIRST=4, RW_LAST=5
integer, parameter :: RW_GENERATION=6, RW_LAYOUT_DIGEST=7, RW_WEIGHT_OFFSET=8, RW_FLAG_OFFSET=9
integer, parameter :: RW_FILE_BYTES=10, RW_HEADER_CHECKSUM=16

type :: flex_weights_header
    integer(int32) :: version         = FLEX_WEIGHTS_VERSION
    integer(int32) :: real_bytes      = 0_int32
    integer(int32) :: nptcls          = 0_int32
    integer(int32) :: nstates         = 0_int32   !< states in the delivery this file belongs to
    integer(int32) :: state_index     = 0_int32   !< which of them this file is
    integer(int32) :: ncomp           = 0_int32
    integer(int32) :: weight_kind     = 0_int32
    integer(int32) :: provenance      = 0_int32
    integer(int32) :: state           = FLEX_WEIGHTS_CANDIDATE
    integer(int32) :: box             = 0_int32
    integer(int32) :: box_crop        = 0_int32
    integer(int64) :: generation      = 0_int64
    integer(int64) :: layout_digest   = 0_int64
    integer(int64) :: scalar_offset   = 0_int64
    integer(int64) :: weight_offset   = 0_int64
    integer(int64) :: flag_offset     = 0_int64
    integer(int64) :: file_bytes      = 0_int64
    integer(int64) :: scalar_checksum = 0_int64
    real(real64)   :: smpd            = 0.0_real64
    real(real64)   :: smpd_crop       = 0.0_real64
end type flex_weights_header

contains

    subroutine flex_weights_init_header(header, nptcls, nstates, state_index, ncomp, box, smpd, box_crop, &
        &smpd_crop, weight_kind, provenance, generation, layout_digest)
        type(flex_weights_header), intent(out) :: header
        integer,                   intent(in)  :: nptcls, nstates, state_index, ncomp, box, box_crop
        real,                      intent(in)  :: smpd, smpd_crop
        integer(int32),            intent(in)  :: weight_kind, provenance
        integer(int64),            intent(in)  :: generation, layout_digest
        header%version       = FLEX_WEIGHTS_VERSION
        header%real_bytes    = int(storage_size(0.0_real32)/8, int32)
        header%nptcls        = int(nptcls, int32)
        header%nstates       = int(nstates, int32)
        header%state_index   = int(state_index, int32)
        header%ncomp         = int(ncomp, int32)
        header%box           = int(box, int32)
        header%box_crop      = int(box_crop, int32)
        header%smpd          = real(smpd, real64)
        header%smpd_crop     = real(smpd_crop, real64)
        header%weight_kind   = weight_kind
        header%provenance    = provenance
        header%generation    = generation
        header%layout_digest = layout_digest
        header%state         = FLEX_WEIGHTS_CANDIDATE
        call set_state_layout(header)
    end subroutine flex_weights_init_header

    subroutine set_state_layout(header)
        type(flex_weights_header), intent(inout) :: header
        integer(int64) :: scalar_bytes, weight_bytes, flag_bytes
        scalar_bytes = int(FLEX_WEIGHTS_SCALAR_NFIELDS + header%ncomp, int64) * 8_int64
        weight_bytes = int(header%nptcls, int64) * int(header%real_bytes, int64)
        flag_bytes   = int(header%nptcls, int64) * 4_int64
        header%scalar_offset = STATE_HEADER_BYTES + 1_int64
        header%weight_offset = header%scalar_offset + scalar_bytes
        header%flag_offset   = header%weight_offset + weight_bytes
        header%file_bytes    = header%flag_offset + flag_bytes - 1_int64
    end subroutine set_state_layout

    subroutine flex_weights_read_header(path, header, status, message)
        character(len=*),          intent(in)  :: path
        type(flex_weights_header), intent(out) :: header
        integer,                   intent(out) :: status
        character(len=*),          intent(out) :: message
        character(len=16) :: magic
        integer(int64) :: words(STATE_NWORDS)
        integer :: funit, io_stat
        header = flex_weights_header()
        status = 0
        message = ''
        if( .not. file_exists(path) )then
            status = 1; message = 'flex weights file does not exist:'//trim(path); return
        endif
        open(newunit=funit, file=trim(path), access='stream', form='unformatted', &
            &status='old', action='read', iostat=io_stat)
        if( io_stat /= 0 )then
            status = io_stat; message = 'cannot open flex weights header'; return
        endif
        read(funit, pos=1, iostat=io_stat) magic
        if( io_stat == 0 ) read(funit, pos=17, iostat=io_stat) words
        close(funit)
        if( io_stat /= 0 )then
            status = io_stat; message = 'cannot read flex weights header'; return
        endif
        if( magic /= STATE_MAGIC )then
            status = 1; message = 'invalid flex weights magic'; return
        endif
        if( checksum_words(words(:W_HEADER_CHECKSUM-1)) /= words(W_HEADER_CHECKSUM) )then
            status = 1; message = 'invalid flex weights header checksum'; return
        endif
        if( words(W_HEADER_BYTES) /= STATE_HEADER_BYTES )then
            status = 1; message = 'invalid flex weights header size'; return
        endif
        call unpack_state_header(words, header)
        call validate_header_values(header, status, message)
    end subroutine flex_weights_read_header

    subroutine pack_state_header(header, words)
        type(flex_weights_header), intent(in)  :: header
        integer(int64),            intent(out) :: words(STATE_NWORDS)
        words = 0_int64
        words(W_VERSION)         = int(header%version, int64)
        words(W_HEADER_BYTES)    = STATE_HEADER_BYTES
        words(W_REAL_BYTES)      = int(header%real_bytes, int64)
        words(W_NPTCLS)          = int(header%nptcls, int64)
        words(W_NSTATES)         = int(header%nstates, int64)
        words(W_STATE_INDEX)     = int(header%state_index, int64)
        words(W_NCOMP)           = int(header%ncomp, int64)
        words(W_KIND)            = int(header%weight_kind, int64)
        words(W_PROVENANCE)      = int(header%provenance, int64)
        words(W_GENERATION)      = header%generation
        words(W_STATE)           = int(header%state, int64)
        words(W_LAYOUT_DIGEST)   = header%layout_digest
        words(W_BOX)             = int(header%box, int64)
        words(W_SMPD_BITS)       = transfer(header%smpd, words(W_SMPD_BITS))
        words(W_BOX_CROP)        = int(header%box_crop, int64)
        words(W_SMPD_CROP_BITS)  = transfer(header%smpd_crop, words(W_SMPD_CROP_BITS))
        words(W_SCALAR_OFFSET)   = header%scalar_offset
        words(W_WEIGHT_OFFSET)   = header%weight_offset
        words(W_FLAG_OFFSET)     = header%flag_offset
        words(W_FILE_BYTES)      = header%file_bytes
        words(W_SCALAR_CHECKSUM) = header%scalar_checksum
        words(W_HEADER_CHECKSUM) = checksum_words(words(:W_HEADER_CHECKSUM-1))
    end subroutine pack_state_header

    subroutine unpack_state_header(words, header)
        integer(int64),            intent(in)  :: words(STATE_NWORDS)
        type(flex_weights_header), intent(out) :: header
        header%version         = int(words(W_VERSION), int32)
        header%real_bytes      = int(words(W_REAL_BYTES), int32)
        header%nptcls          = int(words(W_NPTCLS), int32)
        header%nstates         = int(words(W_NSTATES), int32)
        header%state_index     = int(words(W_STATE_INDEX), int32)
        header%ncomp           = int(words(W_NCOMP), int32)
        header%weight_kind     = int(words(W_KIND), int32)
        header%provenance      = int(words(W_PROVENANCE), int32)
        header%generation      = words(W_GENERATION)
        header%state           = int(words(W_STATE), int32)
        header%layout_digest   = words(W_LAYOUT_DIGEST)
        header%box             = int(words(W_BOX), int32)
        header%smpd            = transfer(words(W_SMPD_BITS), header%smpd)
        header%box_crop        = int(words(W_BOX_CROP), int32)
        header%smpd_crop       = transfer(words(W_SMPD_CROP_BITS), header%smpd_crop)
        header%scalar_offset   = words(W_SCALAR_OFFSET)
        header%weight_offset   = words(W_WEIGHT_OFFSET)
        header%flag_offset     = words(W_FLAG_OFFSET)
        header%file_bytes      = words(W_FILE_BYTES)
        header%scalar_checksum = words(W_SCALAR_CHECKSUM)
    end subroutine unpack_state_header

    subroutine validate_header_values(header, status, message)
        type(flex_weights_header), intent(in)  :: header
        integer,                   intent(out) :: status
        character(len=*),          intent(out) :: message
        type(flex_weights_header) :: expected
        status = 0
        message = ''
        if( header%version /= FLEX_WEIGHTS_VERSION )then
            status = 1; message = 'unsupported flex weights version'; return
        endif
        if( header%real_bytes /= storage_size(0.0_real32)/8 )then
            status = 1; message = 'unsupported flex weights scalar kind'; return
        endif
        if( header%nptcls < 1 .or. header%nstates < 1 .or. header%ncomp < 0 .or. header%box < 1 )then
            status = 1; message = 'invalid flex weights dimensions'; return
        endif
        if( header%state_index < 1 .or. header%state_index > header%nstates )then
            status = 1; message = 'invalid flex weights state index'; return
        endif
        if( header%smpd <= 0.0_real64 )then
            status = 1; message = 'invalid flex weights sampling distance'; return
        endif
        if( header%box_crop < 0 .or. header%smpd_crop < 0.0_real64 )then
            status = 1; message = 'invalid flex weights working lattice'; return
        endif
        if( header%weight_kind /= FLEX_WEIGHTS_KIND_PARTITION .and. &
            &header%weight_kind /= FLEX_WEIGHTS_KIND_KERNEL )then
            status = 1; message = 'invalid flex weights kind'; return
        endif
        select case(header%provenance)
            case(FLEX_WEIGHTS_PROV_FLEX_PCA, FLEX_WEIGHTS_PROV_MERGED, FLEX_WEIGHTS_PROV_EXTERNAL)
            case default
                status = 1; message = 'invalid flex weights provenance'; return
        end select
        if( header%generation < 1_int64 .or. header%layout_digest == 0_int64 )then
            status = 1; message = 'invalid flex weights generation or layout identity'; return
        endif
        if( header%state /= FLEX_WEIGHTS_CANDIDATE .and. header%state /= FLEX_WEIGHTS_COMMITTED )then
            status = 1; message = 'invalid flex weights transaction state'; return
        endif
        expected = header
        call set_state_layout(expected)
        if( header%scalar_offset /= expected%scalar_offset .or. &
            &header%weight_offset /= expected%weight_offset .or. &
            &header%flag_offset /= expected%flag_offset .or. &
            &header%file_bytes /= expected%file_bytes )then
            status = 1; message = 'invalid flex weights section layout'; return
        endif
    end subroutine validate_header_values

    !> Create a sized candidate. With a committed source of identical identity its sections are
    !! copied byte for byte, so an update that rewrites only some rows keeps the others.
    subroutine flex_weights_create_candidate(path, header, status, message, source_path)
        character(len=*),          intent(in)    :: path
        type(flex_weights_header), intent(inout) :: header
        integer,                   intent(out)   :: status
        character(len=*),          intent(out)   :: message
        character(len=*), optional, intent(in)   :: source_path
        type(flex_weights_header) :: source_header
        integer(int8), allocatable :: buffer(:)
        integer(int8) :: zero
        integer(int64) :: pos, remain
        integer :: src_unit, dst_unit, io_stat, ncopy
        status = 0
        message = ''
        header%state = FLEX_WEIGHTS_CANDIDATE
        call set_state_layout(header)
        call validate_header_values(header, status, message)
        if( status /= 0 ) return
        if( present(source_path) )then
            if( len_trim(source_path) > 0 .and. file_exists(source_path) )then
                call flex_weights_validate_file(source_path, status, message, deep=.true.)
                if( status /= 0 ) return
                call flex_weights_read_header(source_path, source_header, status, message)
                if( status /= 0 ) return
                call validate_copy_compatibility(source_header, header, status, message)
                if( status /= 0 ) return
                open(newunit=src_unit, file=trim(source_path), access='stream', form='unformatted', &
                    &status='old', action='read', iostat=io_stat)
                if( io_stat /= 0 )then
                    status = io_stat; message = 'cannot open committed flex weights'; return
                endif
                open(newunit=dst_unit, file=trim(path), access='stream', form='unformatted', &
                    &status='replace', action='readwrite', iostat=io_stat)
                if( io_stat /= 0 )then
                    close(src_unit)
                    status = io_stat; message = 'cannot create flex weights candidate'; return
                endif
                allocate(buffer(1048576))
                pos = 1_int64
                remain = source_header%file_bytes
                do while( remain > 0_int64 )
                    ncopy = int(min(remain, int(size(buffer),int64)))
                    read(src_unit, pos=pos, iostat=io_stat) buffer(:ncopy)
                    if( io_stat /= 0 ) exit
                    write(dst_unit, pos=pos, iostat=io_stat) buffer(:ncopy)
                    if( io_stat /= 0 ) exit
                    pos = pos + int(ncopy,int64)
                    remain = remain - int(ncopy,int64)
                enddo
                deallocate(buffer)
                close(src_unit)
                if( io_stat /= 0 )then
                    close(dst_unit)
                    status = io_stat; message = 'cannot copy committed flex weights'; return
                endif
                header%scalar_checksum = source_header%scalar_checksum
                call write_header_unit(dst_unit, header, io_stat)
                flush(dst_unit)
                close(dst_unit)
                if( io_stat /= 0 )then
                    status = io_stat; message = 'cannot finalize copied flex weights candidate'; return
                endif
                call simple_sync_file(path, io_stat)
                if( io_stat /= 0 )then
                    status = io_stat; message = 'cannot sync copied flex weights candidate'; return
                endif
                return
            endif
        endif
        header%scalar_checksum = 0_int64
        open(newunit=dst_unit, file=trim(path), access='stream', form='unformatted', &
            &status='replace', action='readwrite', iostat=io_stat)
        if( io_stat /= 0 )then
            status = io_stat; message = 'cannot create flex weights candidate'; return
        endif
        call write_header_unit(dst_unit, header, io_stat)
        zero = 0_int8
        if( io_stat == 0 ) write(dst_unit, pos=header%file_bytes, iostat=io_stat) zero
        flush(dst_unit)
        close(dst_unit)
        if( io_stat /= 0 )then
            status = io_stat; message = 'cannot size flex weights candidate'; return
        endif
        call simple_sync_file(path, io_stat)
        if( io_stat /= 0 )then
            status = io_stat; message = 'cannot sync flex weights candidate'
        endif
    end subroutine flex_weights_create_candidate

    subroutine validate_copy_compatibility(source, target, status, message)
        type(flex_weights_header), intent(in)  :: source, target
        integer,                   intent(out) :: status
        character(len=*),          intent(out) :: message
        status = 0
        message = ''
        if( source%state /= FLEX_WEIGHTS_COMMITTED )then
            status = 1; message = 'flex weights copy source is not committed'; return
        endif
        if( source%nptcls /= target%nptcls .or. source%nstates /= target%nstates .or. &
            &source%state_index /= target%state_index .or. source%ncomp /= target%ncomp .or. &
            &source%box /= target%box .or. source%layout_digest /= target%layout_digest .or. &
            &abs(source%smpd-target%smpd) > 1.0e-9_real64 )then
            status = 1; message = 'committed flex weights are incompatible with candidate'
        endif
    end subroutine validate_copy_compatibility

    subroutine write_header_unit(funit, header, io_stat)
        integer,                   intent(in)  :: funit
        type(flex_weights_header), intent(in)  :: header
        integer,                   intent(out) :: io_stat
        integer(int64) :: words(STATE_NWORDS)
        call pack_state_header(header, words)
        write(funit, pos=1, iostat=io_stat) STATE_MAGIC
        if( io_stat == 0 ) write(funit, pos=17, iostat=io_stat) words
    end subroutine write_header_unit

    !> Contiguous write of the rows first_row .. first_row+size(weights)-1 into a candidate
    subroutine flex_weights_write_rows(path, first_row, weights, flags, status, message)
        character(len=*), intent(in)  :: path
        integer,          intent(in)  :: first_row
        real(real32),     intent(in)  :: weights(:)
        integer(int32),   intent(in)  :: flags(:)
        integer,          intent(out) :: status
        character(len=*), intent(out) :: message
        type(flex_weights_header) :: header
        integer(int64) :: weight_pos, flag_pos
        integer :: funit, io_stat, nrows
        call flex_weights_read_header(path, header, status, message)
        if( status /= 0 ) return
        if( header%state /= FLEX_WEIGHTS_CANDIDATE )then
            status = 1; message = 'row writes require a flex weights candidate'; return
        endif
        nrows = size(weights)
        if( size(flags) /= nrows .or. first_row < 1 .or. nrows < 1 .or. first_row+nrows-1 > header%nptcls )then
            status = 1; message = 'invalid flex weights row write range'; return
        endif
        weight_pos = header%weight_offset + int(first_row-1,int64)*int(header%real_bytes,int64)
        flag_pos   = header%flag_offset   + int(first_row-1,int64)*4_int64
        open(newunit=funit, file=trim(path), access='stream', form='unformatted', &
            &status='old', action='readwrite', iostat=io_stat)
        if( io_stat /= 0 )then
            status = io_stat; message = 'cannot open flex weights candidate for row write'; return
        endif
        write(funit, pos=weight_pos, iostat=io_stat) weights
        if( io_stat == 0 ) write(funit, pos=flag_pos, iostat=io_stat) flags
        flush(funit)
        close(funit)
        if( io_stat /= 0 )then
            status = io_stat; message = 'cannot write flex weights rows'; return
        endif
        call simple_sync_file(path, io_stat)
        status = io_stat
        if( status /= 0 ) message = 'cannot sync flex weights rows'
    end subroutine flex_weights_write_rows

    subroutine flex_weights_read_rows(path, first_row, last_row, weights, flags, status, message)
        character(len=*), intent(in)  :: path
        integer,          intent(in)  :: first_row, last_row
        real(real32),   allocatable, intent(out) :: weights(:)
        integer(int32), allocatable, intent(out) :: flags(:)
        integer,          intent(out) :: status
        character(len=*), intent(out) :: message
        type(flex_weights_header) :: header
        integer(int64) :: weight_pos, flag_pos
        integer :: funit, io_stat, nrows
        call flex_weights_read_header(path, header, status, message)
        if( status /= 0 ) return
        if( first_row < 1 .or. last_row < first_row .or. last_row > header%nptcls )then
            status = 1; message = 'invalid flex weights row read range'; return
        endif
        nrows = last_row-first_row+1
        allocate(weights(nrows), flags(nrows))
        weight_pos = header%weight_offset + int(first_row-1,int64)*int(header%real_bytes,int64)
        flag_pos   = header%flag_offset   + int(first_row-1,int64)*4_int64
        open(newunit=funit, file=trim(path), access='stream', form='unformatted', &
            &status='old', action='read', iostat=io_stat)
        if( io_stat /= 0 )then
            status = io_stat; message = 'cannot open flex weights rows'; return
        endif
        read(funit, pos=weight_pos, iostat=io_stat) weights
        if( io_stat == 0 ) read(funit, pos=flag_pos, iostat=io_stat) flags
        close(funit)
        if( io_stat /= 0 )then
            status = io_stat; message = 'cannot read flex weights rows'; return
        endif
        status = 0
        message = ''
    end subroutine flex_weights_read_rows

    !> Scalar section: FLEX_WEIGHTS_SCALAR_NFIELDS + ncomp real64, checksummed in the header
    subroutine flex_weights_write_scalars(path, scalars, status, message)
        character(len=*), intent(in)  :: path
        real(real64),     intent(in)  :: scalars(:)
        integer,          intent(out) :: status
        character(len=*), intent(out) :: message
        type(flex_weights_header) :: header
        integer :: funit, io_stat
        call flex_weights_read_header(path, header, status, message)
        if( status /= 0 ) return
        if( header%state /= FLEX_WEIGHTS_CANDIDATE )then
            status = 1; message = 'scalar writes require a flex weights candidate'; return
        endif
        if( size(scalars) /= FLEX_WEIGHTS_SCALAR_NFIELDS+int(header%ncomp) )then
            status = 1; message = 'invalid flex weights scalar section size'; return
        endif
        header%scalar_checksum = checksum_reals64(scalars)
        open(newunit=funit, file=trim(path), access='stream', form='unformatted', &
            &status='old', action='readwrite', iostat=io_stat)
        if( io_stat /= 0 )then
            status = io_stat; message = 'cannot open flex weights scalar section'; return
        endif
        write(funit, pos=header%scalar_offset, iostat=io_stat) scalars
        if( io_stat == 0 ) call write_header_unit(funit, header, io_stat)
        flush(funit)
        close(funit)
        if( io_stat /= 0 )then
            status = io_stat; message = 'cannot write flex weights scalar section'; return
        endif
        call simple_sync_file(path, io_stat)
        status = io_stat
        if( status /= 0 ) message = 'cannot sync flex weights scalar section'
    end subroutine flex_weights_write_scalars

    subroutine flex_weights_read_scalars(path, scalars, status, message)
        character(len=*), intent(in)  :: path
        real(real64), allocatable, intent(out) :: scalars(:)
        integer,          intent(out) :: status
        character(len=*), intent(out) :: message
        type(flex_weights_header) :: header
        integer :: funit, io_stat
        call flex_weights_read_header(path, header, status, message)
        if( status /= 0 ) return
        allocate(scalars(FLEX_WEIGHTS_SCALAR_NFIELDS+int(header%ncomp)))
        open(newunit=funit, file=trim(path), access='stream', form='unformatted', &
            &status='old', action='read', iostat=io_stat)
        if( io_stat /= 0 )then
            status = io_stat; message = 'cannot open flex weights scalar section'; return
        endif
        read(funit, pos=header%scalar_offset, iostat=io_stat) scalars
        close(funit)
        if( io_stat /= 0 )then
            status = io_stat; message = 'cannot read flex weights scalar section'; return
        endif
        status = 0
        message = ''
    end subroutine flex_weights_read_scalars

    !> Shallow: header + file size. Deep: also the scalar section checksum.
    subroutine flex_weights_validate_file(path, status, message, deep)
        character(len=*), intent(in)  :: path
        integer,          intent(out) :: status
        character(len=*), intent(out) :: message
        logical, optional, intent(in) :: deep
        type(flex_weights_header) :: header
        real(real64), allocatable :: scalars(:)
        integer(int64) :: file_bytes
        logical :: l_deep
        call flex_weights_read_header(path, header, status, message)
        if( status /= 0 ) return
        inquire(file=trim(path), size=file_bytes, iostat=status)
        if( status /= 0 )then
            message = 'cannot inspect flex weights file size'; return
        endif
        if( file_bytes /= header%file_bytes )then
            status = 1; message = 'flex weights file size does not match header'; return
        endif
        l_deep = .false.
        if( present(deep) ) l_deep = deep
        if( .not. l_deep ) return
        call flex_weights_read_scalars(path, scalars, status, message)
        if( status /= 0 ) return
        if( checksum_reals64(scalars) /= header%scalar_checksum )then
            status = 1; message = 'flex weights scalar section checksum mismatch'; return
        endif
        status = 0
        message = ''
    end subroutine flex_weights_validate_file

    !> Mark the candidate committed and atomically rename it over the committed path
    subroutine flex_weights_publish(candidate_path, committed_path, status, message)
        character(len=*), intent(in)  :: candidate_path, committed_path
        integer,          intent(out) :: status
        character(len=*), intent(out) :: message
        type(flex_weights_header) :: header
        type(string) :: parent
        integer :: funit, io_stat
        call flex_weights_validate_file(candidate_path, status, message, deep=.true.)
        if( status /= 0 ) return
        call flex_weights_read_header(candidate_path, header, status, message)
        if( status /= 0 ) return
        if( header%state /= FLEX_WEIGHTS_CANDIDATE )then
            status = 1; message = 'only a flex weights candidate can be published'; return
        endif
        header%state = FLEX_WEIGHTS_COMMITTED
        open(newunit=funit, file=trim(candidate_path), access='stream', form='unformatted', &
            &status='old', action='readwrite', iostat=io_stat)
        if( io_stat /= 0 )then
            status = io_stat; message = 'cannot open flex weights candidate for publication'; return
        endif
        call write_header_unit(funit, header, io_stat)
        flush(funit)
        close(funit)
        if( io_stat /= 0 )then
            status = io_stat; message = 'cannot mark flex weights candidate committed'; return
        endif
        call simple_sync_file(candidate_path, io_stat)
        if( io_stat /= 0 )then
            status = io_stat; message = 'cannot sync committed flex weights candidate'; return
        endif
        call simple_atomic_replace(candidate_path, committed_path, io_stat)
        if( io_stat /= 0 )then
            status = io_stat; message = 'cannot atomically publish flex weights'; return
        endif
        parent = get_fpath(string(committed_path))
        call simple_sync_dir(parent, io_stat)
        call parent%kill
        status = io_stat
        if( status /= 0 ) message = 'cannot sync flex weights directory'
    end subroutine flex_weights_publish

    !> One worker's exclusive row range of one state for the update that commits `generation`
    subroutine flex_weights_write_local_range(path, generation, layout_digest, state_index, first_row, &
        &weights, flags, status, message)
        character(len=*), intent(in)  :: path
        integer(int64),   intent(in)  :: generation, layout_digest
        integer,          intent(in)  :: state_index, first_row
        real(real32),     intent(in)  :: weights(:)
        integer(int32),   intent(in)  :: flags(:)
        integer,          intent(out) :: status
        character(len=*), intent(out) :: message
        integer(int64) :: words(RANGE_NWORDS), weight_offset, flag_offset, file_bytes
        integer :: funit, io_stat, nrows, last_row
        status = 0
        message = ''
        nrows    = size(weights)
        last_row = first_row + nrows - 1
        if( first_row < 1 .or. nrows < 1 .or. state_index < 1 .or. size(flags) /= nrows )then
            status = 1; message = 'invalid local flex weights range'; return
        endif
        if( file_exists(path) )then
            status = 1; message = 'local flex weights range file already exists'; return
        endif
        weight_offset = RANGE_HEADER_BYTES + 1_int64
        flag_offset   = weight_offset + int(nrows,int64)*4_int64
        file_bytes    = flag_offset + int(nrows,int64)*4_int64 - 1_int64
        words = 0_int64
        words(RW_VERSION)         = int(FLEX_WEIGHTS_VERSION,int64)
        words(RW_HEADER_BYTES)    = RANGE_HEADER_BYTES
        words(RW_STATE_INDEX)     = int(state_index,int64)
        words(RW_FIRST)           = int(first_row,int64)
        words(RW_LAST)            = int(last_row,int64)
        words(RW_GENERATION)      = generation
        words(RW_LAYOUT_DIGEST)   = layout_digest
        words(RW_WEIGHT_OFFSET)   = weight_offset
        words(RW_FLAG_OFFSET)     = flag_offset
        words(RW_FILE_BYTES)      = file_bytes
        words(RW_HEADER_CHECKSUM) = checksum_words(words(:RW_HEADER_CHECKSUM-1))
        open(newunit=funit, file=trim(path), access='stream', form='unformatted', &
            &status='new', action='readwrite', iostat=io_stat)
        if( io_stat /= 0 )then
            status = io_stat; message = 'cannot create local flex weights range file'; return
        endif
        write(funit, pos=1, iostat=io_stat) RANGE_MAGIC
        if( io_stat == 0 ) write(funit, pos=17, iostat=io_stat) words
        if( io_stat == 0 ) write(funit, pos=weight_offset, iostat=io_stat) weights
        if( io_stat == 0 ) write(funit, pos=flag_offset, iostat=io_stat) flags
        flush(funit)
        close(funit)
        if( io_stat /= 0 )then
            status = io_stat; message = 'cannot write local flex weights range file'; return
        endif
        call simple_sync_file(path, io_stat)
        status = io_stat
        if( status /= 0 ) message = 'cannot sync local flex weights range file'
    end subroutine flex_weights_write_local_range

    subroutine flex_weights_read_local_range(path, generation, layout_digest, state_index, first_row, last_row, &
        &weights, flags, status, message)
        character(len=*), intent(in)  :: path
        integer(int64),   intent(out) :: generation, layout_digest
        integer,          intent(out) :: state_index, first_row, last_row
        real(real32),   allocatable, intent(out) :: weights(:)
        integer(int32), allocatable, intent(out) :: flags(:)
        integer,          intent(out) :: status
        character(len=*), intent(out) :: message
        character(len=16) :: magic
        integer(int64) :: words(RANGE_NWORDS), actual_bytes, expected_flag_offset, expected_file_bytes
        integer :: funit, io_stat, nrows
        status = 0
        message = ''
        generation = 0_int64; layout_digest = 0_int64
        state_index = 0; first_row = 0; last_row = -1
        if( .not. file_exists(path) )then
            status = 1; message = 'local flex weights range file does not exist'; return
        endif
        open(newunit=funit, file=trim(path), access='stream', form='unformatted', &
            &status='old', action='read', iostat=io_stat)
        if( io_stat /= 0 )then
            status = io_stat; message = 'cannot open local flex weights range file'; return
        endif
        read(funit, pos=1, iostat=io_stat) magic
        if( io_stat == 0 ) read(funit, pos=17, iostat=io_stat) words
        if( io_stat /= 0 )then
            close(funit)
            status = io_stat; message = 'cannot read local flex weights range header'; return
        endif
        if( magic /= RANGE_MAGIC .or. words(RW_VERSION) /= FLEX_WEIGHTS_VERSION .or. &
            &words(RW_HEADER_BYTES) /= RANGE_HEADER_BYTES )then
            close(funit); status = 1; message = 'invalid local flex weights range header'; return
        endif
        if( checksum_words(words(:RW_HEADER_CHECKSUM-1)) /= words(RW_HEADER_CHECKSUM) )then
            close(funit); status = 1; message = 'invalid local flex weights range header checksum'; return
        endif
        state_index = int(words(RW_STATE_INDEX))
        first_row   = int(words(RW_FIRST)); last_row = int(words(RW_LAST))
        generation  = words(RW_GENERATION); layout_digest = words(RW_LAYOUT_DIGEST)
        nrows = last_row-first_row+1
        inquire(file=trim(path), size=actual_bytes, iostat=io_stat)
        expected_flag_offset = RANGE_HEADER_BYTES + 1_int64 + int(nrows,int64)*4_int64
        expected_file_bytes  = expected_flag_offset + int(nrows,int64)*4_int64 - 1_int64
        if( io_stat /= 0 .or. actual_bytes /= words(RW_FILE_BYTES) .or. nrows < 1 .or. first_row < 1 .or. &
            &state_index < 1 .or. words(RW_WEIGHT_OFFSET) /= RANGE_HEADER_BYTES+1_int64 .or. &
            &words(RW_FLAG_OFFSET) /= expected_flag_offset .or. &
            &words(RW_FILE_BYTES) /= expected_file_bytes .or. generation < 1_int64 .or. &
            &layout_digest == 0_int64 )then
            close(funit); status = 1; message = 'invalid local flex weights range layout'; return
        endif
        allocate(weights(nrows), flags(nrows))
        read(funit, pos=words(RW_WEIGHT_OFFSET), iostat=io_stat) weights
        if( io_stat == 0 ) read(funit, pos=words(RW_FLAG_OFFSET), iostat=io_stat) flags
        close(funit)
        if( io_stat /= 0 )then
            status = io_stat; message = 'cannot read local flex weights range data'; return
        endif
    end subroutine flex_weights_read_local_range

    ! ---- FNV-1a helpers, the same hash the sigma2 store uses (kept private there) ----

    pure integer(int64) function checksum_reals64(values) result(hash)
        real(real64), intent(in) :: values(:)
        integer(int64) :: bits
        integer :: i
        hash = FNV_OFFSET
        do i = 1, size(values)
            bits = transfer(values(i), bits)
            hash = fnv_integer(hash, bits, 8)
        enddo
        if( hash == 0_int64 ) hash = 1_int64
    end function checksum_reals64

    pure integer(int64) function checksum_words(words) result(hash)
        integer(int64), intent(in) :: words(:)
        integer :: i
        hash = FNV_OFFSET
        do i = 1, size(words)
            hash = fnv_integer(hash, words(i), 8)
        enddo
        if( hash == 0_int64 ) hash = 1_int64
    end function checksum_words

    pure integer(int64) function fnv_integer(seed, value, nbytes) result(hash)
        integer(int64), intent(in) :: seed, value
        integer,        intent(in) :: nbytes
        integer :: i
        hash = seed
        do i = 0, nbytes-1
            hash = ieor(hash, ibits(value, 8*i, 8))
            hash = hash * FNV_PRIME
        enddo
    end function fnv_integer

end module simple_flex_weights_file
