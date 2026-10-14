!@descr: versioned binary persistence and transaction primitives for canonical sigma2 state
module simple_sigma2_state_file
use, intrinsic :: iso_fortran_env, only: int8, int32, int64, real32, real64
use simple_fileio, only: get_fpath, fopen, fclose
use simple_string, only: string
use simple_syslib, only: file_exists, simple_sync_file, simple_sync_dir, simple_atomic_replace
implicit none
private

public :: sigma2_state_header
public :: sigma2_state_init_header, sigma2_state_read_header, sigma2_state_create_candidate
public :: sigma2_state_write_particles, sigma2_state_read_particles, sigma2_state_read_particles_into
public :: sigma2_state_open_candidate, sigma2_state_write_particles_unit, sigma2_state_close_candidate
public :: sigma2_state_open_store, sigma2_state_read_particles_unit
public :: sigma2_state_open_local_range, sigma2_state_read_local_range_rows
public :: sigma2_state_write_global, sigma2_state_read_global
public :: sigma2_state_write_local_range, sigma2_state_read_local_range
public :: sigma2_state_validate_file, sigma2_state_publish
public :: sigma2_state_digest_begin, sigma2_state_digest_text, sigma2_state_digest_integer
public :: sigma2_state_digest_file, sigma2_state_artifact_signature
public :: SIGMA2_STATE_FNAME, SIGMA2_STATE_NEXT_FNAME
public :: SIGMA2_PROV_PSPEC, SIGMA2_PROV_RESIDUAL
public :: SIGMA2_STATE_CANDIDATE, SIGMA2_STATE_COMMITTED

character(len=*), parameter :: SIGMA2_STATE_FNAME      = 'sigma2_state.bin'
character(len=*), parameter :: SIGMA2_STATE_NEXT_FNAME = 'sigma2_state.next'
character(len=16), parameter :: STATE_MAGIC = 'SIMPLE_SIGMA2_V2'
character(len=16), parameter :: RANGE_MAGIC = 'SIMPLE_S2_RANGE2'

! version 2 (release 4): no reserved particle-integrity section; version 1 files are refused.
! One global spectrum per even/odd half: the grouping word is the global policy (1) with one group,
! and the range-mode word is the local mode (1); the retired stack policy (2), STAR-seed
! provenance (3) and direct range mode (2) are refused.
integer(int32), parameter :: SIGMA2_STATE_VERSION = 2_int32
integer(int32), parameter :: SIGMA2_GROUP_GLOBAL  = 1_int32
integer(int32), parameter :: SIGMA2_PROV_PSPEC     = 1_int32
integer(int32), parameter :: SIGMA2_PROV_RESIDUAL  = 2_int32
integer(int32), parameter :: SIGMA2_RANGE_LOCAL  = 1_int32
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
integer, parameter :: W_GLOBAL_OFFSET=15, W_PARTICLE_OFFSET=16, W_FILE_BYTES=17
! The retired checksum words (18 and 32) stay in the version-two layout: zero in new files,
! ignored on read.
integer, parameter :: W_RETIRED_CHECKSUM=18, W_SMPD_BITS=19, W_RETIRED_HEADER_CHECKSUM=32

integer, parameter :: RW_VERSION=1, RW_HEADER_BYTES=2, RW_REAL_BYTES=3, RW_KFROM=4, RW_KTO=5
integer, parameter :: RW_FIRST=6, RW_LAST=7, RW_GENERATION=8, RW_LAYOUT_DIGEST=9
integer, parameter :: RW_DATA_OFFSET=10, RW_FILE_BYTES=11
integer, parameter :: RW_RETIRED_CHECKSUM=16 !< zero in new range files, ignored on read

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
    integer(int64) :: global_offset   = 0_int64
    integer(int64) :: particle_offset  = 0_int64
    integer(int64) :: file_bytes       = 0_int64
    real(real64)   :: smpd = 0.0_real64
end type sigma2_state_header

contains

    subroutine sigma2_state_init_header(header, kfrom, kto, nptcls, box, smpd, generation, layout_digest, provenance)
        type(sigma2_state_header), intent(out) :: header
        integer,                   intent(in)  :: kfrom, kto, nptcls, box
        real,                      intent(in)  :: smpd
        integer(int32),            intent(in)  :: provenance
        integer(int64),            intent(in)  :: generation, layout_digest
        header%version       = SIGMA2_STATE_VERSION
        header%real_bytes    = int(storage_size(0.0_real32)/8, int32)
        header%kfrom         = int(kfrom, int32)
        header%kto           = int(kto, int32)
        header%nptcls        = int(nptcls, int32)
        header%box           = int(box, int32)
        header%ngroups       = 1_int32
        header%grouping      = SIGMA2_GROUP_GLOBAL
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
        integer(int64) :: nshell, global_bytes, particle_bytes
        nshell = int(header%kto, int64) - int(header%kfrom, int64) + 1_int64
        global_bytes  = nshell * 2_int64 * int(header%ngroups, int64) * int(header%real_bytes, int64)
        particle_bytes = nshell * int(header%nptcls, int64) * int(header%real_bytes, int64)
        header%global_offset   = STATE_HEADER_BYTES + 1_int64
        header%particle_offset  = header%global_offset + global_bytes
        header%file_bytes       = header%particle_offset + particle_bytes - 1_int64
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
            status = 1; message = 'sigma2 state file does not exist:'//trim(path); return
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
        words(W_GLOBAL_OFFSET)   = header%global_offset
        words(W_PARTICLE_OFFSET)  = header%particle_offset
        words(W_FILE_BYTES)       = header%file_bytes
        words(W_RETIRED_CHECKSUM) = 0_int64
        words(W_SMPD_BITS) = transfer(header%smpd, words(W_SMPD_BITS))
        words(W_RETIRED_HEADER_CHECKSUM) = 0_int64
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
        header%global_offset   = words(W_GLOBAL_OFFSET)
        header%particle_offset  = words(W_PARTICLE_OFFSET)
        header%file_bytes       = words(W_FILE_BYTES)
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
        if( header%nptcls < 1 .or. header%box < 1 )then
            status = 1; message = 'invalid sigma2 state dimensions'; return
        endif
        if( header%smpd <= 0.0_real64 )then
            status = 1; message = 'invalid sigma2 state sampling distance'; return
        endif
        if( header%grouping /= SIGMA2_GROUP_GLOBAL .or. header%ngroups /= 1 )then
            status = 1; message = 'sigma2 state is not one global spectrum per half (a retired grouping policy)'; return
        endif
        select case(header%provenance)
            case(SIGMA2_PROV_PSPEC, SIGMA2_PROV_RESIDUAL)
            case default
                status = 1; message = 'invalid sigma2 state provenance'; return
        end select
        if( header%generation < 1_int64 .or. header%layout_digest == 0_int64 )then
            status = 1; message = 'invalid sigma2 state generation or layout identity'; return
        endif
        if( header%range_io /= SIGMA2_RANGE_LOCAL )then
            status = 1; message = 'invalid sigma2 state range I/O mode'; return
        endif
        if( header%state /= SIGMA2_STATE_CANDIDATE .and. header%state /= SIGMA2_STATE_COMMITTED )then
            status = 1; message = 'invalid sigma2 state transaction state'; return
        endif
        expected = header
        call set_state_layout(expected)
        if( header%global_offset /= expected%global_offset .or. &
            &header%particle_offset /= expected%particle_offset .or. &
            &header%file_bytes /= expected%file_bytes )then
            status = 1; message = 'invalid sigma2 state section layout'; return
        endif
    end subroutine validate_header_values

    !> A fresh candidate at path of the given identity: header written, file at its full size, synced
    subroutine sigma2_state_create_candidate(path, header, status, message)
        character(len=*),          intent(in)    :: path
        type(sigma2_state_header), intent(inout) :: header
        integer,                   intent(out)   :: status
        character(len=*),          intent(out)   :: message
        integer :: funit, io_stat
        call sigma2_state_open_candidate(path, header, funit, status, message)
        if( status /= 0 ) return
        flush(funit)
        close(funit)
        call simple_sync_file(path, io_stat)
        if( io_stat /= 0 )then
            status = io_stat; message = 'cannot sync sigma2 candidate'
        endif
    end subroutine sigma2_state_create_candidate

    !> Opens a fresh candidate at path for one writer: the header is written
    !! and the file has its full size. Records (sigma2_state_write_particles_unit) and the global
    !! spectra (sigma2_state_close_candidate, which closes the unit and syncs once) follow through funit.
    subroutine sigma2_state_open_candidate(path, header, funit, status, message)
        character(len=*),          intent(in)    :: path
        type(sigma2_state_header), intent(inout) :: header
        integer,                   intent(out)   :: funit
        integer,                   intent(out)   :: status
        character(len=*),          intent(out)   :: message
        integer(int8) :: zero
        integer :: io_stat
        status  = 0
        message = ''
        funit   = -1
        header%state           = SIGMA2_STATE_CANDIDATE
        header%range_io        = SIGMA2_RANGE_LOCAL
        call set_state_layout(header)
        call validate_header_values(header, status, message)
        if( status /= 0 ) return
        open(newunit=funit, file=trim(path), access='stream', form='unformatted', &
            &status='replace', action='readwrite', iostat=io_stat)
        if( io_stat /= 0 )then
            funit = -1
            status = io_stat; message = 'cannot create sigma2 candidate'; return
        endif
        call write_header_unit(funit, header, io_stat)
        zero = 0_int8
        if( io_stat == 0 ) write(funit, pos=header%file_bytes, iostat=io_stat) zero
        if( io_stat /= 0 )then
            close(funit)
            funit = -1
            status = io_stat; message = 'cannot size sigma2 candidate'
        endif
    end subroutine sigma2_state_open_candidate

    !> writes the records spectra(shells, rows) from row first_row through the open candidate unit
    subroutine sigma2_state_write_particles_unit(funit, header, first_row, spectra, status, message)
        integer,                   intent(in)  :: funit
        type(sigma2_state_header), intent(in)  :: header
        integer,                   intent(in)  :: first_row
        real(real32),              intent(in)  :: spectra(:,:)
        integer,                   intent(out) :: status
        character(len=*),          intent(out) :: message
        integer(int64) :: data_pos
        integer :: nshell
        status  = 0
        message = ''
        nshell = int(header%kto-header%kfrom+1)
        if( size(spectra,1) /= nshell .or. first_row < 1 .or. &
            &first_row+size(spectra,2)-1 > header%nptcls )then
            status = 1; message = 'invalid sigma2 particle write range'; return
        endif
        data_pos = header%particle_offset + int(first_row-1,int64)*int(nshell*header%real_bytes,int64)
        write(funit, pos=data_pos, iostat=status) spectra
        if( status /= 0 ) message = 'cannot write sigma2 particle range'
    end subroutine sigma2_state_write_particles_unit

    !> writes the two global spectra (shells, half: even, odd) through the open
    !! candidate unit, closes it and syncs the candidate once
    subroutine sigma2_state_close_candidate(path, funit, header, spectra, status, message)
        character(len=*),          intent(in)    :: path
        integer,                   intent(in)    :: funit
        type(sigma2_state_header), intent(inout) :: header
        real(real32),              intent(in)    :: spectra(:,:)
        integer,                   intent(out)   :: status
        character(len=*),          intent(out)   :: message
        integer :: io_stat
        status  = 0
        message = ''
        if( any(shape(spectra) /= [int(header%kto-header%kfrom+1),2]) )then
            close(funit)
            status = 1; message = 'invalid sigma2 global section shape'; return
        endif
        write(funit, pos=header%global_offset, iostat=io_stat) spectra
        flush(funit)
        close(funit)
        if( io_stat /= 0 )then
            status = io_stat; message = 'cannot write sigma2 global section'; return
        endif
        call simple_sync_file(path, io_stat)
        status = io_stat
        if( status /= 0 ) message = 'cannot sync sigma2 candidate'
    end subroutine sigma2_state_close_candidate

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
        integer :: funit, io_stat
        call sigma2_state_read_header(path, header, status, message)
        if( status /= 0 ) return
        if( header%state /= SIGMA2_STATE_CANDIDATE )then
            status = 1; message = 'particle writes require a sigma2 candidate'; return
        endif
        open(newunit=funit, file=trim(path), access='stream', form='unformatted', &
            &status='old', action='readwrite', iostat=io_stat)
        if( io_stat /= 0 )then
            status = io_stat; message = 'cannot open sigma2 candidate for particle write'; return
        endif
        call sigma2_state_write_particles_unit(funit, header, first_row, spectra, status, message)
        flush(funit)
        close(funit)
        if( status /= 0 ) return
        call simple_sync_file(path, io_stat)
        status = io_stat
        if( status /= 0 ) message = 'cannot sync sigma2 particle range'
    end subroutine sigma2_state_write_particles

    !> the records of rows first_row:last_row, into an array the reader allocates (shells, rows)
    subroutine sigma2_state_read_particles(path, first_row, last_row, spectra, status, message)
        character(len=*), intent(in) :: path
        integer,          intent(in) :: first_row, last_row
        real(real32), allocatable, intent(out) :: spectra(:,:)
        integer,          intent(out) :: status
        character(len=*), intent(out) :: message
        type(sigma2_state_header) :: header
        call sigma2_state_read_header(path, header, status, message)
        if( status /= 0 ) return
        if( first_row < 1 .or. last_row < first_row .or. last_row > header%nptcls )then
            status = 1; message = 'invalid sigma2 particle read range'; return
        endif
        allocate(spectra(int(header%kto-header%kfrom+1),last_row-first_row+1))
        call sigma2_state_read_particles_into(path, first_row, spectra, status, message)
    end subroutine sigma2_state_read_particles

    !> the records of rows first_row:first_row+size(spectra,2)-1 straight into the caller's array
    !! (shells, rows), whose first extent must be the store's shell count; any lower bounds
    subroutine sigma2_state_read_particles_into(path, first_row, spectra, status, message)
        character(len=*), intent(in)  :: path
        integer,          intent(in)  :: first_row
        real(real32),     intent(out) :: spectra(:,:)
        integer,          intent(out) :: status
        character(len=*), intent(out) :: message
        type(sigma2_state_header) :: header
        integer :: funit
        call sigma2_state_open_store(path, header, funit, status, message)
        if( status /= 0 ) return
        call sigma2_state_read_particles_unit(funit, header, first_row, spectra, status, message)
        close(funit)
    end subroutine sigma2_state_read_particles_into

    !> opens the store at path for reading through one handle; its header is returned
    subroutine sigma2_state_open_store(path, header, funit, status, message)
        character(len=*),          intent(in)  :: path
        type(sigma2_state_header), intent(out) :: header
        integer,                   intent(out) :: funit
        integer,                   intent(out) :: status
        character(len=*),          intent(out) :: message
        integer :: io_stat
        funit = -1
        call sigma2_state_read_header(path, header, status, message)
        if( status /= 0 ) return
        open(newunit=funit, file=trim(path), access='stream', form='unformatted', &
            &status='old', action='read', iostat=io_stat)
        if( io_stat /= 0 )then
            funit = -1
            status = io_stat; message = 'cannot open sigma2 particle range'
        endif
    end subroutine sigma2_state_open_store

    !> the records of rows first_row:first_row+size(spectra,2)-1 through an open store unit, into the
    !! caller's array (shells, rows), whose first extent must be the store's shell count
    subroutine sigma2_state_read_particles_unit(funit, header, first_row, spectra, status, message)
        integer,                   intent(in)  :: funit
        type(sigma2_state_header), intent(in)  :: header
        integer,                   intent(in)  :: first_row
        real(real32),              intent(out) :: spectra(:,:)
        integer,                   intent(out) :: status
        character(len=*),          intent(out) :: message
        integer(int64) :: data_pos
        integer :: nrows, nshell, last_row
        nshell   = int(header%kto-header%kfrom+1)
        nrows    = size(spectra,2)
        last_row = first_row + nrows - 1
        if( first_row < 1 .or. nrows < 1 .or. last_row > header%nptcls )then
            status = 1; message = 'invalid sigma2 particle read range'; return
        endif
        if( size(spectra,1) /= nshell )then
            status = 1; message = 'sigma2 particle records have another shell count'; return
        endif
        data_pos = header%particle_offset + int(first_row-1,int64)*int(nshell*header%real_bytes,int64)
        read(funit, pos=data_pos, iostat=status) spectra
        if( status /= 0 )then
            message = 'cannot read sigma2 particle range'; return
        endif
        message = ''
    end subroutine sigma2_state_read_particles_unit

    subroutine sigma2_state_write_global(path, spectra, status, message)
        character(len=*), intent(in) :: path
        real(real32),     intent(in) :: spectra(:,:) !< (shells, half): even, odd
        integer,          intent(out) :: status
        character(len=*), intent(out) :: message
        type(sigma2_state_header) :: header
        integer :: funit, io_stat, nshell
        call sigma2_state_read_header(path, header, status, message)
        if( status /= 0 ) return
        if( header%state /= SIGMA2_STATE_CANDIDATE )then
            status = 1; message = 'global spectra writes require a sigma2 candidate'; return
        endif
        nshell = int(header%kto-header%kfrom+1)
        if( any(shape(spectra) /= [nshell,2]) )then
            status = 1; message = 'invalid sigma2 global section shape'; return
        endif
        open(newunit=funit, file=trim(path), access='stream', form='unformatted', &
            &status='old', action='readwrite', iostat=io_stat)
        if( io_stat /= 0 )then
            status = io_stat; message = 'cannot open sigma2 global section'; return
        endif
        write(funit, pos=header%global_offset, iostat=io_stat) spectra
        flush(funit)
        close(funit)
        if( io_stat /= 0 )then
            status = io_stat; message = 'cannot write sigma2 global section'; return
        endif
        call simple_sync_file(path, io_stat)
        status = io_stat
        if( status /= 0 ) message = 'cannot sync sigma2 global section'
    end subroutine sigma2_state_write_global

    subroutine sigma2_state_read_global(path, spectra, status, message)
        character(len=*), intent(in) :: path
        real(real32), allocatable, intent(out) :: spectra(:,:) !< (shells, half): even, odd
        integer,          intent(out) :: status
        character(len=*), intent(out) :: message
        type(sigma2_state_header) :: header
        integer :: funit, io_stat, nshell
        call sigma2_state_read_header(path, header, status, message)
        if( status /= 0 ) return
        nshell = int(header%kto-header%kfrom+1)
        allocate(spectra(nshell,2))
        open(newunit=funit, file=trim(path), access='stream', form='unformatted', &
            &status='old', action='read', iostat=io_stat)
        if( io_stat /= 0 )then
            status = io_stat; message = 'cannot open sigma2 global section'; return
        endif
        read(funit, pos=header%global_offset, iostat=io_stat) spectra
        close(funit)
        if( io_stat /= 0 )then
            status = io_stat; message = 'cannot read sigma2 global section'; return
        endif
        status = 0
        message = ''
    end subroutine sigma2_state_read_global

    subroutine sigma2_state_validate_file(path, status, message)
        character(len=*), intent(in) :: path
        integer,          intent(out) :: status
        character(len=*), intent(out) :: message
        type(sigma2_state_header) :: header
        integer(int64) :: file_bytes
        call sigma2_state_read_header(path, header, status, message)
        if( status /= 0 ) return
        inquire(file=trim(path), size=file_bytes, iostat=status)
        if( status /= 0 )then
            message = 'cannot inspect sigma2 state file size'; return
        endif
        if( file_bytes /= header%file_bytes )then
            status = 1; message = 'sigma2 state file size does not match header'; return
        endif
    end subroutine sigma2_state_validate_file

    subroutine sigma2_state_publish(candidate_path, committed_path, status, message)
        character(len=*), intent(in) :: candidate_path, committed_path
        integer,          intent(out) :: status
        character(len=*), intent(out) :: message
        type(sigma2_state_header) :: header
        type(string) :: parent
        integer :: funit, io_stat
        call sigma2_state_validate_file(candidate_path, status, message)
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
        integer(int64) :: words(RANGE_NWORDS), data_offset, file_bytes
        integer :: funit, io_stat, last_row
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
        file_bytes  = data_offset + int(size(spectra),int64)*4_int64 - 1_int64
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
        words(RW_DATA_OFFSET)     = data_offset
        words(RW_FILE_BYTES)      = file_bytes
        words(RW_RETIRED_CHECKSUM) = 0_int64
        open(newunit=funit, file=trim(path), access='stream', form='unformatted', &
            &status='new', action='readwrite', iostat=io_stat)
        if( io_stat /= 0 )then
            status = io_stat; message = 'cannot create local sigma2 range file'; return
        endif
        write(funit, pos=1, iostat=io_stat) RANGE_MAGIC
        if( io_stat == 0 ) write(funit, pos=17, iostat=io_stat) words
        if( io_stat == 0 ) write(funit, pos=data_offset, iostat=io_stat) spectra
        flush(funit)
        close(funit)
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
        integer :: funit
        call sigma2_state_open_local_range(path, funit, generation, layout_digest, first_row, last_row, &
            &kfrom, kto, status, message)
        if( status /= 0 ) return
        allocate(spectra(kto-kfrom+1,last_row-first_row+1))
        call sigma2_state_read_local_range_rows(funit, first_row, first_row, spectra, status, message)
        close(funit)
    end subroutine sigma2_state_read_local_range

    !> opens the local range at path and validates its header; the unit stays open for
    !! sigma2_state_read_local_range_rows and the caller closes it
    subroutine sigma2_state_open_local_range(path, funit, generation, layout_digest, first_row, last_row, &
        &kfrom, kto, status, message)
        character(len=*), intent(in)  :: path
        integer,          intent(out) :: funit
        integer(int64),   intent(out) :: generation, layout_digest
        integer,          intent(out) :: first_row, last_row, kfrom, kto
        integer,          intent(out) :: status
        character(len=*), intent(out) :: message
        character(len=16) :: magic
        integer(int64) :: words(RANGE_NWORDS), actual_bytes
        integer(int64) :: expected_file_bytes
        integer :: io_stat, nshell, nrows
        status = 0
        message = ''
        funit = -1
        generation = 0_int64; layout_digest = 0_int64
        first_row = 0; last_row = -1; kfrom = 0; kto = -1
        if( .not. file_exists(path) )then
            status = 1; message = 'local sigma2 range file does not exist'; return
        endif
        open(newunit=funit, file=trim(path), access='stream', form='unformatted', &
            &status='old', action='read', iostat=io_stat)
        if( io_stat /= 0 )then
            funit = -1
            status = io_stat; message = 'cannot open local sigma2 range file'; return
        endif
        read(funit, pos=1, iostat=io_stat) magic
        if( io_stat == 0 ) read(funit, pos=17, iostat=io_stat) words
        if( io_stat /= 0 )then
            call close_range
            status = io_stat; message = 'cannot read local sigma2 range header'; return
        endif
        if( magic /= RANGE_MAGIC .or. words(RW_VERSION) /= SIGMA2_STATE_VERSION .or. &
            &words(RW_HEADER_BYTES) /= RANGE_HEADER_BYTES .or. words(RW_REAL_BYTES) /= 4_int64 )then
            call close_range
            status = 1; message = 'invalid local sigma2 range header'; return
        endif
        kfrom = int(words(RW_KFROM)); kto = int(words(RW_KTO))
        first_row = int(words(RW_FIRST)); last_row = int(words(RW_LAST))
        generation = words(RW_GENERATION); layout_digest = words(RW_LAYOUT_DIGEST)
        nshell = kto-kfrom+1; nrows = last_row-first_row+1
        inquire(file=trim(path), size=actual_bytes, iostat=io_stat)
        expected_file_bytes = RANGE_HEADER_BYTES + int(nshell,int64)*int(nrows,int64)*4_int64
        if( io_stat /= 0 .or. actual_bytes /= words(RW_FILE_BYTES) .or. nshell < 1 .or. nrows < 1 .or. &
            &words(RW_DATA_OFFSET) /= RANGE_HEADER_BYTES+1_int64 .or. &
            &words(RW_FILE_BYTES) /= expected_file_bytes .or. generation < 1_int64 .or. &
            &layout_digest == 0_int64 )then
            call close_range
            status = 1; message = 'invalid local sigma2 range layout'; return
        endif

      contains

        subroutine close_range
            close(funit)
            funit = -1
        end subroutine close_range

    end subroutine sigma2_state_open_local_range

    !> rows first_row:first_row+size(spectra,2)-1 of the open local range whose first row is range_first,
    !! into the caller's array (shells, rows); the range's bounds were validated when it was opened
    subroutine sigma2_state_read_local_range_rows(funit, range_first, first_row, spectra, status, message)
        integer,          intent(in)  :: funit, range_first, first_row
        real(real32),     intent(out) :: spectra(:,:)
        integer,          intent(out) :: status
        character(len=*), intent(out) :: message
        integer(int64) :: data_pos
        data_pos = RANGE_HEADER_BYTES + 1_int64 + &
            &int(first_row-range_first,int64)*int(size(spectra,1),int64)*4_int64
        read(funit, pos=data_pos, iostat=status) spectra
        if( status /= 0 )then
            message = 'cannot read local sigma2 range data'; return
        endif
        message = ''
    end subroutine sigma2_state_read_local_range_rows

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

    !> Generation of a structurally valid committed sigma2 artifact, 0 when invalid.
    !! The manifest does not compare header or payload checksums.
    integer(int64) function sigma2_state_artifact_signature( fname ) result( signature )
        class(string), intent(in) :: fname
        type(sigma2_state_header) :: header
        character(len=256) :: message
        integer :: status
        signature = 0_int64
        call sigma2_state_validate_file(fname%to_char(), status, message)
        if( status /= 0 ) return
        call sigma2_state_read_header(fname%to_char(), header, status, message)
        if( status /= 0 .or. header%state /= SIGMA2_STATE_COMMITTED ) return
        signature = header%generation
    end function sigma2_state_artifact_signature

    !> FNV-1a 64 over the bytes of a file (8-byte words, then the tail), 0 when unreadable
    integer(int64) function sigma2_state_digest_file( fname ) result( digest )
        class(string), intent(in) :: fname
        integer, parameter :: CHUNK = 65536
        integer(int64), allocatable :: words(:)
        integer(int64) :: fsize, pos, nwords, i
        integer(int8)  :: tail(8)
        integer :: funit, io_stat, ntail, j
        digest = 0_int64
        if( .not. file_exists(fname) ) return
        inquire(file=fname%to_char(), size=fsize)
        if( fsize < 0 ) return
        call fopen(funit, file=fname, status='OLD', action='READ', access='STREAM', iostat=io_stat)
        if( io_stat /= 0 ) return
        digest = sigma2_state_digest_begin()
        allocate(words(CHUNK))
        pos    = 1
        nwords = fsize / 8_int64
        do while( nwords > 0 )
            i = min(int(CHUNK,int64), nwords)
            read(funit, pos=pos, iostat=io_stat) words(1:i)
            if( io_stat /= 0 )then
                digest = 0_int64
                exit
            endif
            do j = 1, int(i)
                call digest_word(digest, words(j))
            enddo
            pos    = pos + 8_int64*i
            nwords = nwords - i
        enddo
        ntail = int(fsize - 8_int64*(fsize/8_int64))
        if( digest /= 0_int64 .and. ntail > 0 )then
            read(funit, pos=pos, iostat=io_stat) tail(1:ntail)
            if( io_stat /= 0 )then
                digest = 0_int64
            else
                do j = 1, ntail
                    call sigma2_state_digest_integer(digest, int(tail(j)))
                enddo
            endif
        endif
        call fclose(funit)
        if( digest /= 0_int64 ) call sigma2_state_digest_integer(digest, int(mod(fsize, 2147483647_int64)))
        if( digest == 0_int64 .and. io_stat == 0 ) digest = 1_int64
    end function sigma2_state_digest_file

    subroutine digest_word( digest, word )
        integer(int64), intent(inout) :: digest
        integer(int64), intent(in)    :: word
        call sigma2_state_digest_integer(digest, int(ibits(word, 0, 31)))
        call sigma2_state_digest_integer(digest, int(ibits(word, 31, 31)))
        call sigma2_state_digest_integer(digest, int(ibits(word, 62, 2)))
    end subroutine digest_word

end module simple_sigma2_state_file
