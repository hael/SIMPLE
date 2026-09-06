!@descr: versioned binary persistence and transaction primitives for canonical sigma2 state
module simple_sigma2_state_file
use, intrinsic :: iso_fortran_env, only: int8, int32, int64, real32, real64
use simple_fileio, only: get_fpath
use simple_string, only: string
use simple_syslib, only: file_exists, simple_sync_file, simple_sync_dir, simple_atomic_replace
implicit none
private

public :: sigma2_state_header
public :: sigma2_state_init_header, sigma2_state_read_header, sigma2_state_create_candidate
public :: sigma2_state_write_particles, sigma2_state_read_particles
public :: sigma2_state_write_groups, sigma2_state_read_groups
public :: sigma2_state_write_local_range, sigma2_state_read_local_range
public :: sigma2_state_refresh_integrity, sigma2_state_validate_file, sigma2_state_publish
public :: sigma2_state_digest_begin, sigma2_state_digest_text, sigma2_state_digest_integer
public :: SIGMA2_STATE_FNAME, SIGMA2_STATE_NEXT_FNAME
public :: SIGMA2_GROUP_GLOBAL, SIGMA2_GROUP_STACK
public :: SIGMA2_PROV_PSPEC, SIGMA2_PROV_RESIDUAL, SIGMA2_PROV_LEGACY_PARTS, SIGMA2_PROV_STAR_SEED
public :: SIGMA2_RANGE_LOCAL, SIGMA2_RANGE_DIRECT
public :: SIGMA2_STATE_CANDIDATE, SIGMA2_STATE_COMMITTED

character(len=*), parameter :: SIGMA2_STATE_FNAME      = 'sigma2_state.bin'
character(len=*), parameter :: SIGMA2_STATE_NEXT_FNAME = 'sigma2_state.next'
character(len=16), parameter :: STATE_MAGIC = 'SIMPLE_SIGMA2_V1'
character(len=16), parameter :: RANGE_MAGIC = 'SIMPLE_S2_RANGE1'

integer(int32), parameter :: SIGMA2_STATE_VERSION = 1_int32
integer(int32), parameter :: SIGMA2_GROUP_GLOBAL  = 1_int32
integer(int32), parameter :: SIGMA2_GROUP_STACK   = 2_int32
integer(int32), parameter :: SIGMA2_PROV_PSPEC        = 1_int32
integer(int32), parameter :: SIGMA2_PROV_RESIDUAL     = 2_int32
integer(int32), parameter :: SIGMA2_PROV_LEGACY_PARTS = 3_int32
integer(int32), parameter :: SIGMA2_PROV_STAR_SEED    = 4_int32
integer(int32), parameter :: SIGMA2_RANGE_LOCAL  = 1_int32
integer(int32), parameter :: SIGMA2_RANGE_DIRECT = 2_int32
integer(int32), parameter :: SIGMA2_STATE_CANDIDATE = 1_int32
integer(int32), parameter :: SIGMA2_STATE_COMMITTED = 2_int32

integer(int64), parameter :: STATE_HEADER_BYTES = 512_int64
integer(int64), parameter :: RANGE_HEADER_BYTES = 256_int64
integer, parameter :: STATE_NWORDS = 32
integer, parameter :: RANGE_NWORDS = 16
integer(int64), parameter :: FNV_OFFSET = int(z'CBF29CE484222325', int64)
integer(int64), parameter :: FNV_PRIME  = int(z'00000100000001B3', int64)

integer, parameter :: W_VERSION=1, W_HEADER_BYTES=2, W_REAL_BYTES=3, W_KFROM=4, W_KTO=5
integer, parameter :: W_NPTCLS=6, W_BOX=7, W_NGROUPS=8, W_GROUPING=9, W_GENERATION=10
integer, parameter :: W_PROVENANCE=11, W_STATE=12, W_RANGE_IO=13, W_LAYOUT_DIGEST=14
integer, parameter :: W_GROUPED_OFFSET=15, W_PARTICLE_OFFSET=16, W_INTEGRITY_OFFSET=17
integer, parameter :: W_FILE_BYTES=18, W_GROUP_CHECKSUM=19, W_PARTICLE_CHECKSUM=20
integer, parameter :: W_SMPD_BITS=21, W_HEADER_CHECKSUM=32

integer, parameter :: RW_VERSION=1, RW_HEADER_BYTES=2, RW_REAL_BYTES=3, RW_KFROM=4, RW_KTO=5
integer, parameter :: RW_FIRST=6, RW_LAST=7, RW_GENERATION=8, RW_LAYOUT_DIGEST=9
integer, parameter :: RW_DATA_OFFSET=10, RW_INTEGRITY_OFFSET=11, RW_FILE_BYTES=12
integer, parameter :: RW_HEADER_CHECKSUM=16

type :: sigma2_state_header
    integer(int32) :: version       = SIGMA2_STATE_VERSION
    integer(int32) :: real_bytes    = 0_int32
    integer(int32) :: kfrom         = 0_int32
    integer(int32) :: kto           = -1_int32
    integer(int32) :: nptcls        = 0_int32
    integer(int32) :: box           = 0_int32
    integer(int32) :: ngroups       = 0_int32
    integer(int32) :: grouping      = 0_int32
    integer(int32) :: provenance    = 0_int32
    integer(int32) :: state         = SIGMA2_STATE_CANDIDATE
    integer(int32) :: range_io      = SIGMA2_RANGE_LOCAL
    integer(int64) :: generation    = 0_int64
    integer(int64) :: layout_digest = 0_int64
    integer(int64) :: grouped_offset   = 0_int64
    integer(int64) :: particle_offset  = 0_int64
    integer(int64) :: integrity_offset = 0_int64
    integer(int64) :: file_bytes       = 0_int64
    integer(int64) :: group_checksum    = 0_int64
    integer(int64) :: particle_checksum = 0_int64
    real(real64)   :: smpd = 0.0_real64
end type sigma2_state_header

contains

    subroutine sigma2_state_init_header(header, kfrom, kto, nptcls, box, smpd, ngroups, grouping, &
        &generation, layout_digest, provenance)
        type(sigma2_state_header), intent(out) :: header
        integer,                   intent(in)  :: kfrom, kto, nptcls, box, ngroups
        real,                      intent(in)  :: smpd
        integer(int32),            intent(in)  :: grouping, provenance
        integer(int64),            intent(in)  :: generation, layout_digest
        header%version       = SIGMA2_STATE_VERSION
        header%real_bytes    = int(storage_size(0.0_real32)/8, int32)
        header%kfrom         = int(kfrom, int32)
        header%kto           = int(kto, int32)
        header%nptcls        = int(nptcls, int32)
        header%box           = int(box, int32)
        header%ngroups       = int(ngroups, int32)
        header%grouping      = grouping
        header%generation    = generation
        header%layout_digest = layout_digest
        header%provenance    = provenance
        header%state         = SIGMA2_STATE_CANDIDATE
        header%range_io      = SIGMA2_RANGE_LOCAL
        header%smpd          = real(smpd, real64)
        call set_state_layout(header)
    end subroutine sigma2_state_init_header

    subroutine set_state_layout(header)
        type(sigma2_state_header), intent(inout) :: header
        integer(int64) :: nshell, grouped_bytes, particle_bytes
        nshell = int(header%kto, int64) - int(header%kfrom, int64) + 1_int64
        grouped_bytes  = nshell * 2_int64 * int(header%ngroups, int64) * int(header%real_bytes, int64)
        particle_bytes = nshell * int(header%nptcls, int64) * int(header%real_bytes, int64)
        header%grouped_offset   = STATE_HEADER_BYTES + 1_int64
        header%particle_offset  = header%grouped_offset + grouped_bytes
        header%integrity_offset = header%particle_offset + particle_bytes
        header%file_bytes       = header%integrity_offset + int(header%nptcls, int64)*8_int64 - 1_int64
    end subroutine set_state_layout

    subroutine sigma2_state_read_header(path, header, status, message)
        character(len=*),          intent(in)  :: path
        type(sigma2_state_header), intent(out) :: header
        integer,                   intent(out) :: status
        character(len=*),          intent(out) :: message
        character(len=16) :: magic
        integer(int64) :: words(STATE_NWORDS)
        integer :: funit, io_stat
        header = sigma2_state_header()
        status = 0
        message = ''
        if( .not. file_exists(path) )then
            status = 1; message = 'sigma2 state file does not exist'; return
        endif
        open(newunit=funit, file=trim(path), access='stream', form='unformatted', &
            &status='old', action='read', iostat=io_stat)
        if( io_stat /= 0 )then
            status = io_stat; message = 'cannot open sigma2 state header'; return
        endif
        read(funit, pos=1, iostat=io_stat) magic
        if( io_stat == 0 ) read(funit, pos=17, iostat=io_stat) words
        close(funit)
        if( io_stat /= 0 )then
            status = io_stat; message = 'cannot read sigma2 state header'; return
        endif
        if( magic /= STATE_MAGIC )then
            status = 1; message = 'invalid sigma2 state magic'; return
        endif
        if( checksum_words(words(:W_HEADER_CHECKSUM-1)) /= words(W_HEADER_CHECKSUM) )then
            status = 1; message = 'invalid sigma2 state header checksum'; return
        endif
        if( words(W_HEADER_BYTES) /= STATE_HEADER_BYTES )then
            status = 1; message = 'invalid sigma2 state header size'; return
        endif
        call unpack_state_header(words, header)
        call validate_header_values(header, status, message)
    end subroutine sigma2_state_read_header

    subroutine pack_state_header(header, words)
        type(sigma2_state_header), intent(in)  :: header
        integer(int64),            intent(out) :: words(STATE_NWORDS)
        words = 0_int64
        words(W_VERSION)       = int(header%version, int64)
        words(W_HEADER_BYTES)  = STATE_HEADER_BYTES
        words(W_REAL_BYTES)    = int(header%real_bytes, int64)
        words(W_KFROM)         = int(header%kfrom, int64)
        words(W_KTO)           = int(header%kto, int64)
        words(W_NPTCLS)        = int(header%nptcls, int64)
        words(W_BOX)           = int(header%box, int64)
        words(W_NGROUPS)       = int(header%ngroups, int64)
        words(W_GROUPING)      = int(header%grouping, int64)
        words(W_GENERATION)    = header%generation
        words(W_PROVENANCE)    = int(header%provenance, int64)
        words(W_STATE)         = int(header%state, int64)
        words(W_RANGE_IO)      = int(header%range_io, int64)
        words(W_LAYOUT_DIGEST) = header%layout_digest
        words(W_GROUPED_OFFSET)   = header%grouped_offset
        words(W_PARTICLE_OFFSET)  = header%particle_offset
        words(W_INTEGRITY_OFFSET) = header%integrity_offset
        words(W_FILE_BYTES)       = header%file_bytes
        words(W_GROUP_CHECKSUM)    = header%group_checksum
        words(W_PARTICLE_CHECKSUM) = header%particle_checksum
        words(W_SMPD_BITS) = transfer(header%smpd, words(W_SMPD_BITS))
        words(W_HEADER_CHECKSUM) = checksum_words(words(:W_HEADER_CHECKSUM-1))
    end subroutine pack_state_header

    subroutine unpack_state_header(words, header)
        integer(int64),            intent(in)  :: words(STATE_NWORDS)
        type(sigma2_state_header), intent(out) :: header
        header%version       = int(words(W_VERSION), int32)
        header%real_bytes    = int(words(W_REAL_BYTES), int32)
        header%kfrom         = int(words(W_KFROM), int32)
        header%kto           = int(words(W_KTO), int32)
        header%nptcls        = int(words(W_NPTCLS), int32)
        header%box           = int(words(W_BOX), int32)
        header%ngroups       = int(words(W_NGROUPS), int32)
        header%grouping      = int(words(W_GROUPING), int32)
        header%generation    = words(W_GENERATION)
        header%provenance    = int(words(W_PROVENANCE), int32)
        header%state         = int(words(W_STATE), int32)
        header%range_io      = int(words(W_RANGE_IO), int32)
        header%layout_digest = words(W_LAYOUT_DIGEST)
        header%grouped_offset   = words(W_GROUPED_OFFSET)
        header%particle_offset  = words(W_PARTICLE_OFFSET)
        header%integrity_offset = words(W_INTEGRITY_OFFSET)
        header%file_bytes       = words(W_FILE_BYTES)
        header%group_checksum    = words(W_GROUP_CHECKSUM)
        header%particle_checksum = words(W_PARTICLE_CHECKSUM)
        header%smpd = transfer(words(W_SMPD_BITS), header%smpd)
    end subroutine unpack_state_header

    subroutine validate_header_values(header, status, message)
        type(sigma2_state_header), intent(in)  :: header
        integer,                   intent(out) :: status
        character(len=*),          intent(out) :: message
        type(sigma2_state_header) :: expected
        status = 0
        message = ''
        if( header%version /= SIGMA2_STATE_VERSION )then
            status = 1; message = 'unsupported sigma2 state version'; return
        endif
        if( header%real_bytes /= storage_size(0.0_real32)/8 )then
            status = 1; message = 'unsupported sigma2 state scalar kind'; return
        endif
        if( header%kfrom < 0 .or. header%kto < header%kfrom )then
            status = 1; message = 'invalid sigma2 state shell bounds'; return
        endif
        if( header%nptcls < 1 .or. header%box < 1 .or. header%ngroups < 1 )then
            status = 1; message = 'invalid sigma2 state dimensions'; return
        endif
        if( header%smpd <= 0.0_real64 )then
            status = 1; message = 'invalid sigma2 state sampling distance'; return
        endif
        if( header%grouping /= SIGMA2_GROUP_GLOBAL .and. header%grouping /= SIGMA2_GROUP_STACK )then
            status = 1; message = 'invalid sigma2 state grouping policy'; return
        endif
        select case(header%provenance)
            case(SIGMA2_PROV_PSPEC, SIGMA2_PROV_RESIDUAL, SIGMA2_PROV_LEGACY_PARTS, SIGMA2_PROV_STAR_SEED)
            case default
                status = 1; message = 'invalid sigma2 state provenance'; return
        end select
        if( header%generation < 1_int64 .or. header%layout_digest == 0_int64 )then
            status = 1; message = 'invalid sigma2 state generation or layout identity'; return
        endif
        if( header%range_io /= SIGMA2_RANGE_LOCAL .and. header%range_io /= SIGMA2_RANGE_DIRECT )then
            status = 1; message = 'invalid sigma2 state range I/O mode'; return
        endif
        if( header%state /= SIGMA2_STATE_CANDIDATE .and. header%state /= SIGMA2_STATE_COMMITTED )then
            status = 1; message = 'invalid sigma2 state transaction state'; return
        endif
        expected = header
        call set_state_layout(expected)
        if( header%grouped_offset /= expected%grouped_offset .or. &
            &header%particle_offset /= expected%particle_offset .or. &
            &header%integrity_offset /= expected%integrity_offset .or. &
            &header%file_bytes /= expected%file_bytes )then
            status = 1; message = 'invalid sigma2 state section layout'; return
        endif
    end subroutine validate_header_values

    subroutine sigma2_state_create_candidate(path, header, status, message, source_path)
        character(len=*),          intent(in)    :: path
        type(sigma2_state_header), intent(inout) :: header
        integer,                   intent(out)   :: status
        character(len=*),          intent(out)   :: message
        character(len=*), optional, intent(in)   :: source_path
        type(sigma2_state_header) :: source_header
        integer(int8), allocatable :: buffer(:)
        integer(int8) :: zero
        integer(int64) :: pos, remain
        integer :: src_unit, dst_unit, io_stat, ncopy
        status = 0
        message = ''
        header%state    = SIGMA2_STATE_CANDIDATE
        header%range_io = SIGMA2_RANGE_LOCAL
        call set_state_layout(header)
        call validate_header_values(header, status, message)
        if( status /= 0 ) return
        if( present(source_path) )then
            if( len_trim(source_path) > 0 .and. file_exists(source_path) )then
                call sigma2_state_validate_file(source_path, status, message, deep=.true.)
                if( status /= 0 ) return
                call sigma2_state_read_header(source_path, source_header, status, message)
                if( status /= 0 ) return
                call validate_copy_compatibility(source_header, header, status, message)
                if( status /= 0 ) return
                open(newunit=src_unit, file=trim(source_path), access='stream', form='unformatted', &
                    &status='old', action='read', iostat=io_stat)
                if( io_stat /= 0 )then
                    status = io_stat; message = 'cannot open committed sigma2 state'; return
                endif
                open(newunit=dst_unit, file=trim(path), access='stream', form='unformatted', &
                    &status='replace', action='readwrite', iostat=io_stat)
                if( io_stat /= 0 )then
                    close(src_unit)
                    status = io_stat; message = 'cannot create sigma2 candidate'; return
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
                    status = io_stat; message = 'cannot copy committed sigma2 state'; return
                endif
                header%group_checksum    = source_header%group_checksum
                header%particle_checksum = source_header%particle_checksum
                call write_header_unit(dst_unit, header, io_stat)
                flush(dst_unit)
                close(dst_unit)
                if( io_stat /= 0 )then
                    status = io_stat; message = 'cannot finalize copied sigma2 candidate'; return
                endif
                call simple_sync_file(path, io_stat)
                if( io_stat /= 0 )then
                    status = io_stat; message = 'cannot sync copied sigma2 candidate'; return
                endif
                return
            endif
        endif
        header%group_checksum    = 0_int64
        header%particle_checksum = 0_int64
        open(newunit=dst_unit, file=trim(path), access='stream', form='unformatted', &
            &status='replace', action='readwrite', iostat=io_stat)
        if( io_stat /= 0 )then
            status = io_stat; message = 'cannot create sigma2 candidate'; return
        endif
        call write_header_unit(dst_unit, header, io_stat)
        zero = 0_int8
        if( io_stat == 0 ) write(dst_unit, pos=header%file_bytes, iostat=io_stat) zero
        flush(dst_unit)
        close(dst_unit)
        if( io_stat /= 0 )then
            status = io_stat; message = 'cannot size sigma2 candidate'; return
        endif
        call simple_sync_file(path, io_stat)
        if( io_stat /= 0 )then
            status = io_stat; message = 'cannot sync sigma2 candidate'
        endif
    end subroutine sigma2_state_create_candidate

    subroutine validate_copy_compatibility(source, target, status, message)
        type(sigma2_state_header), intent(in)  :: source, target
        integer,                   intent(out) :: status
        character(len=*),          intent(out) :: message
        status = 0
        message = ''
        if( source%state /= SIGMA2_STATE_COMMITTED )then
            status = 1; message = 'sigma2 copy source is not committed'; return
        endif
        if( source%kfrom /= target%kfrom .or. source%kto /= target%kto .or. &
            &source%nptcls /= target%nptcls .or. source%box /= target%box .or. &
            &source%ngroups /= target%ngroups .or. source%grouping /= target%grouping .or. &
            &source%layout_digest /= target%layout_digest .or. &
            &abs(source%smpd-target%smpd) > 1.0e-9_real64 )then
            status = 1; message = 'committed sigma2 state is incompatible with candidate'
        endif
    end subroutine validate_copy_compatibility

    subroutine write_header_unit(funit, header, io_stat)
        integer,                   intent(in)  :: funit
        type(sigma2_state_header), intent(in)  :: header
        integer,                   intent(out) :: io_stat
        integer(int64) :: words(STATE_NWORDS)
        call pack_state_header(header, words)
        write(funit, pos=1, iostat=io_stat) STATE_MAGIC
        if( io_stat == 0 ) write(funit, pos=17, iostat=io_stat) words
    end subroutine write_header_unit

    subroutine sigma2_state_write_particles(path, first_row, spectra, status, message)
        character(len=*), intent(in) :: path
        integer,          intent(in) :: first_row
        real(real32),     intent(in) :: spectra(:,:)
        integer,          intent(out) :: status
        character(len=*), intent(out) :: message
        type(sigma2_state_header) :: header
        integer(int64) :: data_pos, checksum_pos, checksum
        integer :: funit, io_stat, i, row, nshell
        call sigma2_state_read_header(path, header, status, message)
        if( status /= 0 ) return
        if( header%state /= SIGMA2_STATE_CANDIDATE )then
            status = 1; message = 'particle writes require a sigma2 candidate'; return
        endif
        nshell = int(header%kto-header%kfrom+1)
        if( size(spectra,1) /= nshell .or. first_row < 1 .or. &
            &first_row+size(spectra,2)-1 > header%nptcls )then
            status = 1; message = 'invalid sigma2 particle write range'; return
        endif
        open(newunit=funit, file=trim(path), access='stream', form='unformatted', &
            &status='old', action='readwrite', iostat=io_stat)
        if( io_stat /= 0 )then
            status = io_stat; message = 'cannot open sigma2 candidate for particle write'; return
        endif
        do i = 1, size(spectra,2)
            row = first_row + i - 1
            data_pos = header%particle_offset + int(row-1,int64)*int(nshell*header%real_bytes,int64)
            checksum_pos = header%integrity_offset + int(row-1,int64)*8_int64
            checksum = checksum_reals(spectra(:,i))
            write(funit, pos=data_pos, iostat=io_stat) spectra(:,i)
            if( io_stat == 0 ) write(funit, pos=checksum_pos, iostat=io_stat) checksum
            if( io_stat /= 0 ) exit
        enddo
        flush(funit)
        close(funit)
        if( io_stat /= 0 )then
            status = io_stat; message = 'cannot write sigma2 particle range'; return
        endif
        call simple_sync_file(path, io_stat)
        status = io_stat
        if( status /= 0 ) message = 'cannot sync sigma2 particle range'
    end subroutine sigma2_state_write_particles

    subroutine sigma2_state_read_particles(path, first_row, last_row, spectra, status, message, checksums)
        character(len=*), intent(in) :: path
        integer,          intent(in) :: first_row, last_row
        real(real32), allocatable, intent(out) :: spectra(:,:)
        integer,          intent(out) :: status
        character(len=*), intent(out) :: message
        integer(int64), allocatable, optional, intent(out) :: checksums(:)
        type(sigma2_state_header) :: header
        integer(int64) :: data_pos, checksum_pos
        integer(int64), allocatable :: sums(:)
        integer :: funit, io_stat, nrows, nshell
        call sigma2_state_read_header(path, header, status, message)
        if( status /= 0 ) return
        if( first_row < 1 .or. last_row < first_row .or. last_row > header%nptcls )then
            status = 1; message = 'invalid sigma2 particle read range'; return
        endif
        nshell = int(header%kto-header%kfrom+1)
        nrows  = last_row-first_row+1
        allocate(spectra(nshell,nrows), sums(nrows))
        data_pos = header%particle_offset + int(first_row-1,int64)*int(nshell*header%real_bytes,int64)
        checksum_pos = header%integrity_offset + int(first_row-1,int64)*8_int64
        open(newunit=funit, file=trim(path), access='stream', form='unformatted', &
            &status='old', action='read', iostat=io_stat)
        if( io_stat /= 0 )then
            status = io_stat; message = 'cannot open sigma2 particle range'; return
        endif
        read(funit, pos=data_pos, iostat=io_stat) spectra
        if( io_stat == 0 ) read(funit, pos=checksum_pos, iostat=io_stat) sums
        close(funit)
        if( io_stat /= 0 )then
            status = io_stat; message = 'cannot read sigma2 particle range'; return
        endif
        if( present(checksums) ) call move_alloc(sums, checksums)
        status = 0
        message = ''
    end subroutine sigma2_state_read_particles

    subroutine sigma2_state_write_groups(path, groups, status, message)
        character(len=*), intent(in) :: path
        real(real32),     intent(in) :: groups(:,:,:)
        integer,          intent(out) :: status
        character(len=*), intent(out) :: message
        type(sigma2_state_header) :: header
        integer :: funit, io_stat, nshell
        call sigma2_state_read_header(path, header, status, message)
        if( status /= 0 ) return
        if( header%state /= SIGMA2_STATE_CANDIDATE )then
            status = 1; message = 'group writes require a sigma2 candidate'; return
        endif
        nshell = int(header%kto-header%kfrom+1)
        if( any(shape(groups) /= [nshell,2,int(header%ngroups)]) )then
            status = 1; message = 'invalid sigma2 grouped section shape'; return
        endif
        header%group_checksum = checksum_reals(reshape(groups, [size(groups)]))
        open(newunit=funit, file=trim(path), access='stream', form='unformatted', &
            &status='old', action='readwrite', iostat=io_stat)
        if( io_stat /= 0 )then
            status = io_stat; message = 'cannot open sigma2 grouped section'; return
        endif
        write(funit, pos=header%grouped_offset, iostat=io_stat) groups
        if( io_stat == 0 ) call write_header_unit(funit, header, io_stat)
        flush(funit)
        close(funit)
        if( io_stat /= 0 )then
            status = io_stat; message = 'cannot write sigma2 grouped section'; return
        endif
        call simple_sync_file(path, io_stat)
        status = io_stat
        if( status /= 0 ) message = 'cannot sync sigma2 grouped section'
    end subroutine sigma2_state_write_groups

    subroutine sigma2_state_read_groups(path, groups, status, message)
        character(len=*), intent(in) :: path
        real(real32), allocatable, intent(out) :: groups(:,:,:)
        integer,          intent(out) :: status
        character(len=*), intent(out) :: message
        type(sigma2_state_header) :: header
        integer :: funit, io_stat, nshell
        call sigma2_state_read_header(path, header, status, message)
        if( status /= 0 ) return
        nshell = int(header%kto-header%kfrom+1)
        allocate(groups(nshell,2,int(header%ngroups)))
        open(newunit=funit, file=trim(path), access='stream', form='unformatted', &
            &status='old', action='read', iostat=io_stat)
        if( io_stat /= 0 )then
            status = io_stat; message = 'cannot open sigma2 grouped section'; return
        endif
        read(funit, pos=header%grouped_offset, iostat=io_stat) groups
        close(funit)
        if( io_stat /= 0 )then
            status = io_stat; message = 'cannot read sigma2 grouped section'; return
        endif
        status = 0
        message = ''
    end subroutine sigma2_state_read_groups

    subroutine sigma2_state_refresh_integrity(path, status, message)
        character(len=*), intent(in) :: path
        integer,          intent(out) :: status
        character(len=*), intent(out) :: message
        type(sigma2_state_header) :: header
        integer(int64), allocatable :: checksums(:)
        integer :: funit, io_stat
        call sigma2_state_read_header(path, header, status, message)
        if( status /= 0 ) return
        allocate(checksums(header%nptcls))
        open(newunit=funit, file=trim(path), access='stream', form='unformatted', &
            &status='old', action='readwrite', iostat=io_stat)
        if( io_stat /= 0 )then
            deallocate(checksums)
            status = io_stat; message = 'cannot open sigma2 integrity metadata'; return
        endif
        read(funit, pos=header%integrity_offset, iostat=io_stat) checksums
        if( io_stat == 0 )then
            header%particle_checksum = checksum_words(checksums)
            call write_header_unit(funit, header, io_stat)
        endif
        flush(funit)
        close(funit)
        deallocate(checksums)
        if( io_stat /= 0 )then
            status = io_stat; message = 'cannot refresh sigma2 integrity metadata'; return
        endif
        call simple_sync_file(path, io_stat)
        status = io_stat
        if( status /= 0 ) message = 'cannot sync sigma2 integrity metadata'
    end subroutine sigma2_state_refresh_integrity

    subroutine sigma2_state_validate_file(path, status, message, deep)
        character(len=*), intent(in) :: path
        integer,          intent(out) :: status
        character(len=*), intent(out) :: message
        logical, optional, intent(in) :: deep
        type(sigma2_state_header) :: header
        real(real32), allocatable :: groups(:,:,:), spectra(:,:)
        integer(int64), allocatable :: checksums(:)
        integer(int64) :: file_bytes
        integer :: i
        logical :: l_deep
        call sigma2_state_read_header(path, header, status, message)
        if( status /= 0 ) return
        inquire(file=trim(path), size=file_bytes, iostat=status)
        if( status /= 0 )then
            message = 'cannot inspect sigma2 state file size'; return
        endif
        if( file_bytes /= header%file_bytes )then
            status = 1; message = 'sigma2 state file size does not match header'; return
        endif
        l_deep = .false.
        if( present(deep) ) l_deep = deep
        if( .not. l_deep ) return
        call sigma2_state_read_groups(path, groups, status, message)
        if( status /= 0 ) return
        if( checksum_reals(reshape(groups,[size(groups)])) /= header%group_checksum )then
            status = 1; message = 'sigma2 grouped checksum mismatch'; return
        endif
        call sigma2_state_read_particles(path, 1, int(header%nptcls), spectra, status, message, checksums)
        if( status /= 0 ) return
        do i = 1, int(header%nptcls)
            if( checksums(i) /= 0_int64 )then
                if( checksum_reals(spectra(:,i)) /= checksums(i) )then
                    status = 1; message = 'sigma2 particle checksum mismatch'; return
                endif
            endif
        enddo
        if( checksum_words(checksums) /= header%particle_checksum )then
            status = 1; message = 'sigma2 integrity-section checksum mismatch'; return
        endif
        status = 0
        message = ''
    end subroutine sigma2_state_validate_file

    subroutine sigma2_state_publish(candidate_path, committed_path, status, message)
        character(len=*), intent(in) :: candidate_path, committed_path
        integer,          intent(out) :: status
        character(len=*), intent(out) :: message
        type(sigma2_state_header) :: header
        type(string) :: parent
        integer :: funit, io_stat
        call sigma2_state_validate_file(candidate_path, status, message, deep=.true.)
        if( status /= 0 ) return
        call sigma2_state_read_header(candidate_path, header, status, message)
        if( status /= 0 ) return
        if( header%state /= SIGMA2_STATE_CANDIDATE )then
            status = 1; message = 'only a sigma2 candidate can be published'; return
        endif
        header%state = SIGMA2_STATE_COMMITTED
        open(newunit=funit, file=trim(candidate_path), access='stream', form='unformatted', &
            &status='old', action='readwrite', iostat=io_stat)
        if( io_stat /= 0 )then
            status = io_stat; message = 'cannot open sigma2 candidate for publication'; return
        endif
        call write_header_unit(funit, header, io_stat)
        flush(funit)
        close(funit)
        if( io_stat /= 0 )then
            status = io_stat; message = 'cannot mark sigma2 candidate committed'; return
        endif
        call simple_sync_file(candidate_path, io_stat)
        if( io_stat /= 0 )then
            status = io_stat; message = 'cannot sync committed sigma2 candidate'; return
        endif
        call simple_atomic_replace(candidate_path, committed_path, io_stat)
        if( io_stat /= 0 )then
            status = io_stat; message = 'cannot atomically publish sigma2 state'; return
        endif
        parent = get_fpath(string(committed_path))
        call simple_sync_dir(parent, io_stat)
        call parent%kill
        status = io_stat
        if( status /= 0 ) message = 'cannot sync sigma2 state directory'
    end subroutine sigma2_state_publish

    subroutine sigma2_state_write_local_range(path, generation, layout_digest, first_row, spectra, &
        &kfrom, kto, status, message)
        character(len=*), intent(in) :: path
        integer(int64),   intent(in) :: generation, layout_digest
        integer,          intent(in) :: first_row, kfrom, kto
        real(real32),     intent(in) :: spectra(:,:)
        integer,          intent(out) :: status
        character(len=*), intent(out) :: message
        integer(int64) :: words(RANGE_NWORDS), data_offset, integrity_offset, file_bytes
        integer(int64), allocatable :: checksums(:)
        integer :: funit, io_stat, i, last_row
        status = 0
        message = ''
        last_row = first_row + size(spectra,2) - 1
        if( first_row < 1 .or. kto < kfrom .or. size(spectra,1) /= kto-kfrom+1 .or. &
            &size(spectra,2) < 1 )then
            status = 1; message = 'invalid local sigma2 range'; return
        endif
        if( file_exists(path) )then
            status = 1; message = 'local sigma2 range file already exists'; return
        endif
        data_offset = RANGE_HEADER_BYTES + 1_int64
        integrity_offset = data_offset + int(size(spectra),int64)*4_int64
        file_bytes = integrity_offset + int(size(spectra,2),int64)*8_int64 - 1_int64
        words = 0_int64
        words(RW_VERSION)       = int(SIGMA2_STATE_VERSION,int64)
        words(RW_HEADER_BYTES)  = RANGE_HEADER_BYTES
        words(RW_REAL_BYTES)    = int(storage_size(0.0_real32)/8,int64)
        words(RW_KFROM)         = int(kfrom,int64)
        words(RW_KTO)           = int(kto,int64)
        words(RW_FIRST)         = int(first_row,int64)
        words(RW_LAST)          = int(last_row,int64)
        words(RW_GENERATION)    = generation
        words(RW_LAYOUT_DIGEST) = layout_digest
        words(RW_DATA_OFFSET)      = data_offset
        words(RW_INTEGRITY_OFFSET) = integrity_offset
        words(RW_FILE_BYTES)       = file_bytes
        words(RW_HEADER_CHECKSUM)  = checksum_words(words(:RW_HEADER_CHECKSUM-1))
        allocate(checksums(size(spectra,2)))
        do i = 1, size(spectra,2)
            checksums(i) = checksum_reals(spectra(:,i))
        enddo
        open(newunit=funit, file=trim(path), access='stream', form='unformatted', &
            &status='new', action='readwrite', iostat=io_stat)
        if( io_stat /= 0 )then
            deallocate(checksums)
            status = io_stat; message = 'cannot create local sigma2 range file'; return
        endif
        write(funit, pos=1, iostat=io_stat) RANGE_MAGIC
        if( io_stat == 0 ) write(funit, pos=17, iostat=io_stat) words
        if( io_stat == 0 ) write(funit, pos=data_offset, iostat=io_stat) spectra
        if( io_stat == 0 ) write(funit, pos=integrity_offset, iostat=io_stat) checksums
        flush(funit)
        close(funit)
        deallocate(checksums)
        if( io_stat /= 0 )then
            status = io_stat; message = 'cannot write local sigma2 range file'; return
        endif
        call simple_sync_file(path, io_stat)
        status = io_stat
        if( status /= 0 ) message = 'cannot sync local sigma2 range file'
    end subroutine sigma2_state_write_local_range

    subroutine sigma2_state_read_local_range(path, generation, layout_digest, first_row, last_row, &
        &kfrom, kto, spectra, status, message)
        character(len=*), intent(in) :: path
        integer(int64),   intent(out) :: generation, layout_digest
        integer,          intent(out) :: first_row, last_row, kfrom, kto
        real(real32), allocatable, intent(out) :: spectra(:,:)
        integer,          intent(out) :: status
        character(len=*), intent(out) :: message
        character(len=16) :: magic
        integer(int64) :: words(RANGE_NWORDS), actual_bytes
        integer(int64) :: expected_integrity_offset, expected_file_bytes
        integer(int64), allocatable :: checksums(:)
        integer :: funit, io_stat, i, nshell, nrows
        status = 0
        message = ''
        generation = 0_int64; layout_digest = 0_int64
        first_row = 0; last_row = -1; kfrom = 0; kto = -1
        if( .not. file_exists(path) )then
            status = 1; message = 'local sigma2 range file does not exist'; return
        endif
        open(newunit=funit, file=trim(path), access='stream', form='unformatted', &
            &status='old', action='read', iostat=io_stat)
        if( io_stat /= 0 )then
            status = io_stat; message = 'cannot open local sigma2 range file'; return
        endif
        read(funit, pos=1, iostat=io_stat) magic
        if( io_stat == 0 ) read(funit, pos=17, iostat=io_stat) words
        if( io_stat /= 0 )then
            close(funit)
            status = io_stat; message = 'cannot read local sigma2 range header'; return
        endif
        if( magic /= RANGE_MAGIC .or. words(RW_VERSION) /= SIGMA2_STATE_VERSION .or. &
            &words(RW_HEADER_BYTES) /= RANGE_HEADER_BYTES .or. words(RW_REAL_BYTES) /= 4_int64 )then
            close(funit); status = 1; message = 'invalid local sigma2 range header'; return
        endif
        if( checksum_words(words(:RW_HEADER_CHECKSUM-1)) /= words(RW_HEADER_CHECKSUM) )then
            close(funit); status = 1; message = 'invalid local sigma2 range header checksum'; return
        endif
        kfrom = int(words(RW_KFROM)); kto = int(words(RW_KTO))
        first_row = int(words(RW_FIRST)); last_row = int(words(RW_LAST))
        generation = words(RW_GENERATION); layout_digest = words(RW_LAYOUT_DIGEST)
        nshell = kto-kfrom+1; nrows = last_row-first_row+1
        inquire(file=trim(path), size=actual_bytes, iostat=io_stat)
        expected_integrity_offset = RANGE_HEADER_BYTES + 1_int64 + &
            &int(nshell,int64)*int(nrows,int64)*4_int64
        expected_file_bytes = expected_integrity_offset + int(nrows,int64)*8_int64 - 1_int64
        if( io_stat /= 0 .or. actual_bytes /= words(RW_FILE_BYTES) .or. nshell < 1 .or. nrows < 1 .or. &
            &words(RW_DATA_OFFSET) /= RANGE_HEADER_BYTES+1_int64 .or. &
            &words(RW_INTEGRITY_OFFSET) /= expected_integrity_offset .or. &
            &words(RW_FILE_BYTES) /= expected_file_bytes .or. generation < 1_int64 .or. &
            &layout_digest == 0_int64 )then
            close(funit); status = 1; message = 'invalid local sigma2 range layout'; return
        endif
        allocate(spectra(nshell,nrows), checksums(nrows))
        read(funit, pos=words(RW_DATA_OFFSET), iostat=io_stat) spectra
        if( io_stat == 0 ) read(funit, pos=words(RW_INTEGRITY_OFFSET), iostat=io_stat) checksums
        close(funit)
        if( io_stat /= 0 )then
            status = io_stat; message = 'cannot read local sigma2 range data'; return
        endif
        do i = 1, nrows
            if( checksum_reals(spectra(:,i)) /= checksums(i) )then
                status = 1; message = 'local sigma2 range checksum mismatch'; return
            endif
        enddo
        deallocate(checksums)
    end subroutine sigma2_state_read_local_range

    pure integer(int64) function checksum_reals(values) result(hash)
        real(real32), intent(in) :: values(:)
        integer(int32) :: bits
        integer :: i
        hash = FNV_OFFSET
        do i = 1, size(values)
            bits = transfer(values(i), bits)
            hash = fnv_integer(hash, int(bits,int64), 4)
        enddo
        if( hash == 0_int64 ) hash = 1_int64
    end function checksum_reals

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

    pure integer(int64) function sigma2_state_digest_begin() result(hash)
        hash = FNV_OFFSET
    end function sigma2_state_digest_begin

    pure subroutine sigma2_state_digest_text(hash, text)
        integer(int64),   intent(inout) :: hash
        character(len=*), intent(in)    :: text
        integer :: i
        do i = 1, len_trim(text)
            hash = ieor(hash, int(iachar(text(i:i)),int64))
            hash = hash * FNV_PRIME
        enddo
        hash = ieor(hash, int(10,int64))
        hash = hash * FNV_PRIME
    end subroutine sigma2_state_digest_text

    pure subroutine sigma2_state_digest_integer(hash, value)
        integer(int64), intent(inout) :: hash
        integer,        intent(in)    :: value
        hash = fnv_integer(hash, int(value,int64), 8)
    end subroutine sigma2_state_digest_integer

end module simple_sigma2_state_file
