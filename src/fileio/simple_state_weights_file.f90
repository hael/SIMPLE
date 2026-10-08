!@descr: file codec of a state weight set: one binary file per state and a text manifest published last
! A state file holds one state's weight over the full particle layout and a hard-label flag per row:
!   magic(16) | header words(int64 x 8) | weights(real32 x nptcls) | flags(int8 x nptcls)
! The manifest (state_weights.txt) names the files of one generation with their byte sizes, checksums
! and per-state populations, and the set it was derived from (a work project's single-state set). Policy (selection, thresholds, publication order) is in simple_state_weight_set.
module simple_state_weights_file
use, intrinsic :: iso_fortran_env, only: int8, int64, real32, real64
use simple_fileio,            only: fopen, fclose, get_fpath
use simple_string,            only: string
use simple_string_utils,      only: int2str_pad
use simple_syslib,            only: file_exists, simple_sync_file, simple_sync_dir, simple_atomic_replace
use simple_sigma2_state_file, only: sigma2_state_digest_file
implicit none
private

public :: state_weights_manifest, state_weights_file_header
public :: state_weights_fname, state_weights_kind_name, state_weights_kind_id
public :: state_weights_write_state_file, state_weights_read_state_file, state_weights_file_checksum
public :: state_weights_write_manifest, state_weights_read_manifest
public :: STATE_WEIGHTS_MANIFEST_FNAME, STATE_WEIGHTS_KIND_PARTITION, STATE_WEIGHTS_KIND_KERNEL

character(len=*),  parameter :: STATE_WEIGHTS_MANIFEST_FNAME       = 'state_weights.txt'
character(len=*),  parameter :: STATE_WEIGHTS_FBODY          = 'state_weights_g'
character(len=16), parameter :: FILE_MAGIC                   = 'SIMPLE_STWGT_V01'
character(len=*),  parameter :: MANIFEST_FORMAT              = 'simple_state_weights'
integer,           parameter :: MANIFEST_VERSION             = 2
integer,           parameter :: NHEADER_WORDS                = 8
integer(int64),    parameter :: HEADER_BYTES                 = 16_int64 + 8_int64 * NHEADER_WORDS
!> rows are mixture responsibilities: every weighted row sums to one across the states
integer,           parameter :: STATE_WEIGHTS_KIND_PARTITION = 1
!> rows are kernel coefficients in [0,1] that need not sum to one
integer,           parameter :: STATE_WEIGHTS_KIND_KERNEL    = 2

!> The header of one state file; the manifest's identity is repeated so a mixed delivery is detected
type :: state_weights_file_header
    integer(int64) :: version       = 1_int64
    integer(int64) :: nptcls        = 0_int64
    integer(int64) :: nstates       = 0_int64
    integer(int64) :: state         = 0_int64
    integer(int64) :: generation    = 0_int64
    integer(int64) :: layout_digest = 0_int64
    integer(int64) :: kind          = 0_int64
    integer(int64) :: file_bytes    = 0_int64
end type state_weights_file_header

!> The manifest record of one published generation
type :: state_weights_manifest
    integer(int64)              :: generation    = 0_int64
    integer                     :: kind          = 0
    character(len=64)           :: producer      = ''
    integer                     :: nstates       = 0
    integer                     :: nptcls        = 0
    integer(int64)              :: layout_digest = 0_int64
    ! provenance: the set this one was derived from (one state of it, rows remapped); zero for none
    integer(int64)              :: parent_generation    = 0_int64
    integer(int64)              :: parent_layout_digest = 0_int64
    integer                     :: parent_state         = 0
    type(string),   allocatable :: fnames(:)        !< state file names, relative to the manifest's directory
    integer(int64), allocatable :: nbytes(:)
    integer(int64), allocatable :: checksums(:)
    real(real64),   allocatable :: mass(:)          !< sum of the weights
    real(real64),   allocatable :: ess(:)           !< effective sample size, mass**2 / sum of squared weights
    integer,        allocatable :: pop(:)           !< hard population (flagged rows)
end type state_weights_manifest

contains

    !> state_weights_gGGGGG_sNNN.bin
    function state_weights_fname( generation, state ) result( fname )
        integer(int64), intent(in) :: generation
        integer,        intent(in) :: state
        type(string) :: fname
        fname = STATE_WEIGHTS_FBODY//int2str_pad(int(generation),6)//'_s'//int2str_pad(state,3)//'.bin'
    end function state_weights_fname

    function state_weights_kind_name( kind ) result( name )
        integer, intent(in) :: kind
        character(len=:), allocatable :: name
        select case(kind)
            case(STATE_WEIGHTS_KIND_PARTITION)
                name = 'partition'
            case(STATE_WEIGHTS_KIND_KERNEL)
                name = 'kernel'
            case DEFAULT
                name = 'unknown'
        end select
    end function state_weights_kind_name

    integer function state_weights_kind_id( name ) result( kind )
        character(len=*), intent(in) :: name
        select case(trim(name))
            case('partition')
                kind = STATE_WEIGHTS_KIND_PARTITION
            case('kernel')
                kind = STATE_WEIGHTS_KIND_KERNEL
            case DEFAULT
                kind = 0
        end select
    end function state_weights_kind_id

    !> Write and flush one state file; status /= 0 on any failure
    subroutine state_weights_write_state_file( path, header, weights, flags, status, message )
        class(string),                   intent(in)    :: path
        type(state_weights_file_header), intent(inout) :: header
        real(real32),                    intent(in)    :: weights(:)
        integer(int8),                   intent(in)    :: flags(:)
        integer,                         intent(out)   :: status
        character(len=*),                intent(out)   :: message
        integer :: funit
        message = ''
        if( size(weights) /= header%nptcls .or. size(flags) /= header%nptcls )then
            status = 1; message = 'state weight file rows do not match its header'; return
        endif
        header%file_bytes = HEADER_BYTES + 5_int64 * header%nptcls
        open(newunit=funit, file=path%to_char(), access='stream', form='unformatted', status='replace', &
            &action='write', iostat=status)
        if( status /= 0 )then
            message = 'cannot create state weight file '//path%to_char(); return
        endif
        write(funit, iostat=status) FILE_MAGIC, header_words(header), weights, flags
        close(funit)
        if( status /= 0 )then
            message = 'cannot write state weight file '//path%to_char(); return
        endif
        call simple_sync_file(path, status)
        if( status /= 0 ) message = 'cannot flush state weight file '//path%to_char()
    end subroutine state_weights_write_state_file

    !> Read one state file; the header must be consistent with the file's size
    subroutine state_weights_read_state_file( path, header, weights, flags, status, message )
        class(string),                   intent(in)  :: path
        type(state_weights_file_header), intent(out) :: header
        real(real32),  allocatable,      intent(out) :: weights(:)
        integer(int8), allocatable,      intent(out) :: flags(:)
        integer,                         intent(out) :: status
        character(len=*),                intent(out) :: message
        character(len=16) :: magic
        integer(int64)    :: words(NHEADER_WORDS), fsize
        integer           :: funit
        message = ''
        status  = 1
        if( .not. file_exists(path) )then
            message = 'state weight file does not exist: '//path%to_char(); return
        endif
        inquire(file=path%to_char(), size=fsize)
        open(newunit=funit, file=path%to_char(), access='stream', form='unformatted', status='old', &
            &action='read', iostat=status)
        if( status /= 0 )then
            message = 'cannot open state weight file '//path%to_char(); return
        endif
        read(funit, iostat=status) magic, words
        if( status /= 0 .or. magic /= FILE_MAGIC )then
            close(funit)
            status = 1; message = 'not a state weight file: '//path%to_char(); return
        endif
        header = header_from_words(words)
        if( header%nptcls < 1 .or. header%file_bytes /= HEADER_BYTES + 5_int64 * header%nptcls .or. &
            &fsize /= header%file_bytes )then
            close(funit)
            status = 1; message = 'state weight file is truncated or inconsistent: '//path%to_char(); return
        endif
        allocate(weights(header%nptcls), flags(header%nptcls))
        read(funit, iostat=status) weights, flags
        close(funit)
        if( status /= 0 ) message = 'cannot read state weight rows of '//path%to_char()
    end subroutine state_weights_read_state_file

    !> Whole-file checksum (the FNV-1a digest of the sigma2 store); 0 when the file cannot be read
    integer(int64) function state_weights_file_checksum( path ) result( checksum )
        class(string), intent(in) :: path
        checksum = sigma2_state_digest_file(path)
    end function state_weights_file_checksum

    !> Publish a manifest: write it under a temporary name, flush, rename it over the final name
    subroutine state_weights_write_manifest( path, manifest, status, message )
        class(string),                intent(in)  :: path
        type(state_weights_manifest), intent(in)  :: manifest
        integer,                      intent(out) :: status
        character(len=*),             intent(out) :: message
        type(string) :: tmp, dir
        integer :: funit, s
        message = ''
        tmp = path%to_char()//'.tmp'
        call fopen(funit, file=tmp, status='REPLACE', action='WRITE', iostat=status)
        if( status /= 0 )then
            message = 'cannot create the state weight manifest'; return
        endif
        write(funit,'(A,1X,I0)',iostat=status) MANIFEST_FORMAT, MANIFEST_VERSION
        if( status == 0 ) write(funit,'(A,1X,I0)',iostat=status) 'generation', manifest%generation
        if( status == 0 ) write(funit,'(A,1X,A)', iostat=status) 'kind', state_weights_kind_name(manifest%kind)
        if( status == 0 ) write(funit,'(A,1X,A)', iostat=status) 'producer', trim(manifest%producer)
        if( status == 0 ) write(funit,'(A,1X,I0)',iostat=status) 'nstates', manifest%nstates
        if( status == 0 ) write(funit,'(A,1X,I0)',iostat=status) 'nptcls', manifest%nptcls
        if( status == 0 ) write(funit,'(A,1X,I0)',iostat=status) 'layout_digest', manifest%layout_digest
        if( status == 0 ) write(funit,'(A,3(1X,I0))',iostat=status) 'provenance', manifest%parent_generation, &
            &manifest%parent_layout_digest, manifest%parent_state
        do s = 1, manifest%nstates
            if( status /= 0 ) exit
            write(funit,'(A,1X,I0,1X,A,1X,I0,1X,I0,1X,ES24.16E3,1X,ES24.16E3,1X,I0)',iostat=status) 'state', s, &
                &manifest%fnames(s)%to_char(), manifest%nbytes(s), manifest%checksums(s), manifest%mass(s), &
                &manifest%ess(s), manifest%pop(s)
        enddo
        call fclose(funit)
        if( status /= 0 )then
            message = 'cannot write the state weight manifest'; return
        endif
        call simple_sync_file(tmp, status)
        if( status == 0 ) call simple_atomic_replace(tmp, path, status)
        if( status /= 0 )then
            message = 'cannot publish the state weight manifest'; return
        endif
        dir = get_fpath(path)
        call simple_sync_dir(dir, status)
        if( status /= 0 ) message = 'cannot flush the directory of the state weight manifest'
    end subroutine state_weights_write_manifest

    !> Read and check the syntax of a manifest; status /= 0 when it is missing, of another version or malformed
    subroutine state_weights_read_manifest( path, manifest, status, message )
        class(string),                intent(in)  :: path
        type(state_weights_manifest), intent(out) :: manifest
        integer,                      intent(out) :: status
        character(len=*),             intent(out) :: message
        character(len=1024) :: line
        character(len=64)   :: key, word
        character(len=512)  :: fname
        integer :: funit, version, s, idx
        message = ''
        status  = 1
        if( .not. file_exists(path) )then
            message = 'state weight manifest does not exist: '//path%to_char(); return
        endif
        call fopen(funit, file=path, status='OLD', action='READ', iostat=status)
        if( status /= 0 )then
            message = 'cannot open the state weight manifest'; return
        endif
        read(funit,'(A)',iostat=status) line
        if( status == 0 ) read(line,*,iostat=status) word, version
        if( status /= 0 .or. trim(word) /= MANIFEST_FORMAT .or. version /= MANIFEST_VERSION )then
            call fclose(funit)
            status = 1; message = 'unreadable state weight manifest (format or version)'; return
        endif
        call read_keyed(funit, 'generation', line, status)
        if( status == 0 ) read(line,*,iostat=status) manifest%generation
        if( status == 0 ) call read_keyed(funit, 'kind', line, status)
        if( status == 0 ) read(line,*,iostat=status) word
        if( status == 0 ) manifest%kind = state_weights_kind_id(word)
        if( status == 0 ) call read_keyed(funit, 'producer', line, status)
        if( status == 0 ) read(line,*,iostat=status) manifest%producer
        if( status == 0 ) call read_keyed(funit, 'nstates', line, status)
        if( status == 0 ) read(line,*,iostat=status) manifest%nstates
        if( status == 0 ) call read_keyed(funit, 'nptcls', line, status)
        if( status == 0 ) read(line,*,iostat=status) manifest%nptcls
        if( status == 0 ) call read_keyed(funit, 'layout_digest', line, status)
        if( status == 0 ) read(line,*,iostat=status) manifest%layout_digest
        if( status == 0 ) call read_keyed(funit, 'provenance', line, status)
        if( status == 0 ) read(line,*,iostat=status) manifest%parent_generation, manifest%parent_layout_digest, &
            &manifest%parent_state
        if( status /= 0 .or. manifest%kind == 0 .or. manifest%nstates < 1 .or. manifest%nptcls < 1 .or. &
            &manifest%generation < 1_int64 )then
            call fclose(funit)
            status = 1; message = 'malformed state weight manifest header'; return
        endif
        allocate(manifest%fnames(manifest%nstates), manifest%nbytes(manifest%nstates), &
            &manifest%checksums(manifest%nstates), manifest%mass(manifest%nstates), &
            &manifest%ess(manifest%nstates), manifest%pop(manifest%nstates))
        do s = 1, manifest%nstates
            call read_keyed(funit, 'state', line, status)
            if( status == 0 ) read(line,*,iostat=status) idx, fname, manifest%nbytes(s), manifest%checksums(s), &
                &manifest%mass(s), manifest%ess(s), manifest%pop(s)
            if( status /= 0 .or. idx /= s )then
                call fclose(funit)
                status = 1; message = 'malformed state line in the state weight manifest'; return
            endif
            manifest%fnames(s) = trim(fname)
        enddo
        call fclose(funit)
        status = 0
    contains

        !> next line must start with key; returns the rest of the line
        subroutine read_keyed( funit, key_expected, rest, io_stat )
            integer,          intent(in)  :: funit
            character(len=*), intent(in)  :: key_expected
            character(len=*), intent(out) :: rest
            integer,          intent(out) :: io_stat
            character(len=1024) :: full
            integer :: ipos
            rest = ''
            read(funit,'(A)',iostat=io_stat) full
            if( io_stat /= 0 ) return
            full = adjustl(full)
            ipos = index(full, ' ')
            key  = full(:ipos-1)
            if( trim(key) /= key_expected )then
                io_stat = 1
                return
            endif
            rest = full(ipos+1:)
        end subroutine read_keyed

    end subroutine state_weights_read_manifest

    ! ---- header words ----

    function header_words( header ) result( words )
        type(state_weights_file_header), intent(in) :: header
        integer(int64) :: words(NHEADER_WORDS)
        words = [header%version, header%nptcls, header%nstates, header%state, header%generation, &
            &header%layout_digest, header%kind, header%file_bytes]
    end function header_words

    function header_from_words( words ) result( header )
        integer(int64), intent(in) :: words(NHEADER_WORDS)
        type(state_weights_file_header) :: header
        header%version       = words(1)
        header%nptcls        = words(2)
        header%nstates       = words(3)
        header%state         = words(4)
        header%generation    = words(5)
        header%layout_digest = words(6)
        header%kind          = words(7)
        header%file_bytes    = words(8)
    end function header_from_words

end module simple_state_weights_file
