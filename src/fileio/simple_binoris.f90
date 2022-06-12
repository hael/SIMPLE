! for manging orientation data using binary files
module simple_binoris
!$ use omp_lib
use, intrinsic :: ISO_C_BINDING
use simple_defs
use simple_strings
use simple_error
use simple_syslib
use simple_fileio
use simple_defs_ori
use simple_map_reduce, only: split_nobjs_even
implicit none

public :: binoris, binoris_seginfo
private
#include "simple_local_flags.inc"

integer(kind(ENUM_ORISEG)), parameter :: MAX_N_SEGMENTS     = 20
integer(kind=8),            parameter :: N_VARS_HEAD_SEG    = 5
integer(kind=8),            parameter :: N_BYTES_HEADER     = MAX_N_SEGMENTS * N_VARS_HEAD_SEG * 8 ! because dp integer
integer,                    parameter :: THREAD_NSTRINGS    = 2000
integer,                    parameter :: PTCL_BYTES_PER_REC = N_PTCL_ORIPARAMS * 4

type binoris_seginfo
    integer(kind=8) :: fromto(2)          = 0
    integer(kind=8) :: n_bytes_per_record = 0
    integer(kind=8) :: n_records          = 0
    integer(kind=8) :: first_data_byte    = 0
end type binoris_seginfo

type binoris
    private
    type(binoris_seginfo)         :: header(MAX_N_SEGMENTS)
    character(len=:), allocatable :: fname
    integer                       :: n_segments  = 0
    integer                       :: funit       = 0
    logical                       :: l_open      = .false.
  contains
    ! I/O
    procedure          :: open
    procedure, private :: clear_segments
    procedure          :: close
    procedure, private :: read_header
    procedure          :: write_header
    procedure, private :: write_segment_inside_1
    procedure, private :: write_segment_inside_2
    generic            :: write_segment_inside => write_segment_inside_1, write_segment_inside_2
    procedure, private :: byte_manager4seg_inside_1
    procedure, private :: byte_manager4seg_inside_2
    procedure, private :: write_segment_1
    procedure, private :: write_segment_2
    generic            :: write_segment => write_segment_1, write_segment_2
    procedure, private :: add_segment_1
    procedure, private :: add_segment_2
    generic            :: add_segment => add_segment_1, add_segment_2
    procedure          :: update_byte_ranges
    procedure          :: read_first_segment_record
    procedure, private :: read_segment_1, read_segment_2
    generic            :: read_segment => read_segment_1, read_segment_2
    procedure          :: read_record
    ! getters
    procedure          :: get_segments_info
    procedure          :: get_n_segments
    procedure          :: get_fromto
    procedure          :: get_n_records
    procedure          :: get_n_bytes_per_record
    procedure          :: get_n_bytes_tot
    procedure          :: get_first_data_byte
    procedure          :: is_opened
end type binoris

contains

    ! I/O

    subroutine open( self, fname, del_if_exists )
        class(binoris),    intent(inout) :: self          !< instance
        character(len=*),  intent(in)    :: fname         !< filename
        logical, optional, intent(in)    :: del_if_exists !< If the file already exists on disk, replace
        integer(kind=8) :: filesz
        integer(kind(ENUM_ORISEG)) :: isegment
        if( present(del_if_exists) )then
            if( del_if_exists )then
                call del_file(trim(fname))
            endif
        endif
        if( .not. file_exists(trim(fname)) )then
            call self%clear_segments
            call open_local
            return
        endif
        call open_local
        self%fname = trim(fname)
        ! check size
        filesz = funit_size(self%funit)
        if( filesz == -1 )then
            THROW_HARD('file_size cannot be inquired')
        else if( filesz >= N_BYTES_HEADER )then
            ! ok
        else
            THROW_HARD('size of file: '//trim(fname)//' too small to contain a header')
        endif
        ! clear segments before reading header
        call self%clear_segments
        call self%read_header
        self%n_segments = 0 ! for counting # segments
        do isegment=1,MAX_N_SEGMENTS
            ! update # segments counter
            if( self%header(isegment)%n_bytes_per_record > 0 .and. self%header(isegment)%n_records > 0&
                &.and. self%header(isegment)%first_data_byte > 0 ) self%n_segments = isegment
        end do

        contains

            subroutine open_local
                integer :: io_stat, tmpunit
                if( .not. self%l_open )then
                    call fopen(tmpunit, trim(fname), access='STREAM', action='READWRITE',&
                        &status='UNKNOWN', form='UNFORMATTED', iostat=io_stat)
                    call fileiochk('binoris ; open_local '//trim(fname), io_stat)
                    self%funit  = tmpunit
                    self%l_open = .true.
                endif
            end subroutine open_local

    end subroutine open

    subroutine clear_segments( self )
        class(binoris), intent(inout) :: self
        if( self%n_segments <= 0 ) return
        ! clear header
        self%header(:)%fromto(1)          = 0
        self%header(:)%fromto(2)          = 0
        self%header(:)%n_bytes_per_record = 0
        self%header(:)%n_records          = 0
        self%header(:)%first_data_byte    = 0
        ! zero # segments
        self%n_segments = 0
    end subroutine clear_segments

    subroutine close( self )
        class(binoris), intent(inout) :: self !< instance
        integer :: io_stat
        if( self%l_open )then
            io_stat = fsync(fnum(self%funit))
            if( io_stat /= 0 ) write(logfhandle,'(a)')'FSYNC ERROR FILE UNIT: '//int2str(self%funit)
            call fclose(self%funit)
            self%l_open = .false.
        end if
    end subroutine close

    subroutine read_header( self )
        class(binoris), intent(inout) :: self  !< instance
        integer :: io_status
        if( .not. self%l_open ) THROW_HARD('file needs to be open')
        read(unit=self%funit,pos=1,iostat=io_status) self%header
        if( io_status .ne. 0 ) call fileiochk('binoris ::read_header, ERROR reading header bytes ', io_status)
    end subroutine read_header

    subroutine write_header( self )
        class(binoris), intent(inout) :: self  !< instance
        integer :: io_status
        if( .not. self%l_open ) THROW_HARD('file needs to be open')
        write(unit=self%funit,pos=1,iostat=io_status) self%header
        if( io_status .ne. 0 ) call fileiochk('binoris :: write_header, ERROR writing header bytes ', io_status)
    end subroutine write_header

    subroutine write_segment_inside_1( self, isegment, os, fromto )
        use simple_oris, only: oris
        class(binoris),             intent(inout) :: self
        integer(kind(ENUM_ORISEG)), intent(in)    :: isegment
        class(oris), optional,      intent(inout) :: os ! indexed from 1 to nptcls
        integer,     optional,      intent(in)    :: fromto(2)
        character(len=:), allocatable :: str_dyn
        integer :: i, nspaces, noris
        integer(kind(ENUM_ORISEG)) :: iseg
        integer(kind=8) :: end_part1, start_part3, end_part3
        integer(kind=8) :: ibytes, first_data_byte
        real            :: ptcl_record(N_PTCL_ORIPARAMS)
        character(len=1), allocatable ::bytearr_part3(:)
        noris = os%get_noris()
        if( noris == 0 ) return
        if( present(fromto) )then
            if( fromto(1) < 1 .or. fromto(2) > noris )then
                write(logfhandle,*) 'filename',trim(self%fname)
                write(logfhandle,*) 'noris : ', noris
                write(logfhandle,*) 'fromto: ', fromto
                THROW_HARD('fromto out of range')
            endif
        endif
        if( .not. self%l_open ) THROW_HARD('file needs to be open: '//trim(self%fname))
        ! ranges and raw bytes
        call self%byte_manager4seg_inside_1(isegment, end_part1, start_part3, end_part3, bytearr_part3)
        ! add segment to stack, this sets all the information needed for allocation
        call self%add_segment_1(isegment, os, fromto)
        ! error checks
        call self%byte_manager4seg_inside_2(isegment, end_part1, start_part3, end_part3, bytearr_part3)
        ! WRITE FILE
        ! write orientation data (2nd part)
        ibytes = self%header(isegment)%first_data_byte
        if( is_particle_seg(isegment) )then
            do i=self%header(isegment)%fromto(1),self%header(isegment)%fromto(2)
                call os%ori2prec(i, ptcl_record)
                write(unit=self%funit,pos=ibytes) ptcl_record
                ibytes = ibytes + self%header(isegment)%n_bytes_per_record
            end do
        else
            do i=self%header(isegment)%fromto(1),self%header(isegment)%fromto(2)
                str_dyn = os%ori2str(i)
                nspaces = self%header(isegment)%n_bytes_per_record - len_trim(str_dyn)
                if( nspaces > 0 )then
                    write(unit=self%funit,pos=ibytes) trim(str_dyn)//spaces(nspaces)
                else
                    write(unit=self%funit,pos=ibytes) trim(str_dyn)
                endif
                ibytes = ibytes + self%header(isegment)%n_bytes_per_record
            end do
        endif
        ! 3d part
        if( allocated(bytearr_part3) )then
            ! find next nonzero first_data_byte
            do iseg=isegment + 1,self%n_segments
                if( self%header(iseg)%first_data_byte > 0 )then
                    first_data_byte = self%header(iseg)%first_data_byte
                    exit
                endif
            end do
            if( first_data_byte > 0 )then
                write(unit=self%funit,pos=first_data_byte) bytearr_part3
            else
                write(logfhandle,*) 'filename',trim(self%fname)
                write(logfhandle,*) 'first_data_byte: ', first_data_byte
                THROW_HARD('first_data_byte must be > 0 for non-empty 3d segment')
            endif
        endif
        ! write updated header
        write(unit=self%funit,pos=1) self%header
        ! so no need to update header in file after this operation
    end subroutine write_segment_inside_1

    subroutine write_segment_inside_2( self, isegment, os_strings, fromto, strlen_max )
        use simple_oris, only: oris
        class(binoris),             intent(inout) :: self
        integer(kind(ENUM_ORISEG)), intent(in)    :: isegment
        type(str4arr),              intent(in)    :: os_strings(:)
        integer,                    intent(in)    :: fromto(2), strlen_max
        character(len=:), allocatable :: str_dyn
        integer(kind(ENUM_ORISEG)) :: iseg
        integer         :: i, nspaces, noris
        integer(kind=8) :: end_part1, start_part3, end_part3
        integer(kind=8) :: ibytes, first_data_byte
        character(len=1), allocatable :: bytearr_part3(:)
        noris = size(os_strings)
        if( noris == 0 ) return
        if( .not. self%l_open ) THROW_HARD('file needs to be open: '//trim(self%fname))
        ! ranges and raw bytes
        call self%byte_manager4seg_inside_1(isegment, end_part1, start_part3, end_part3, bytearr_part3)
        ! add segment to stack, this sets all the information needed for allocation
        call self%add_segment_2(isegment, fromto, strlen_max)
        ! error checks
        call self%byte_manager4seg_inside_2(isegment, end_part1, start_part3, end_part3, bytearr_part3)
        ! WRITE FILE
        ! write orientation data (2nd part)
        ibytes = self%header(isegment)%first_data_byte
        do i=self%header(isegment)%fromto(1),self%header(isegment)%fromto(2)
            str_dyn = os_strings(i)%str
            nspaces = self%header(isegment)%n_bytes_per_record - len_trim(str_dyn)
            if( nspaces > 0 )then
                write(unit=self%funit,pos=ibytes) trim(str_dyn)//spaces(nspaces)
            else
                write(unit=self%funit,pos=ibytes) trim(str_dyn)
            endif
            ibytes = ibytes + self%header(isegment)%n_bytes_per_record
        end do
        ! 3d part
        if( allocated(bytearr_part3) )then
            ! find next nonzero first_data_byte
            do iseg=isegment + 1,self%n_segments
                if( self%header(iseg)%first_data_byte > 0 )then
                    first_data_byte = self%header(iseg)%first_data_byte
                    exit
                endif
            end do
            if( first_data_byte > 0 )then
                write(unit=self%funit,pos=first_data_byte) bytearr_part3
            else
                write(logfhandle,*) 'filename',trim(self%fname)
                write(logfhandle,*) 'first_data_byte: ', first_data_byte
                THROW_HARD('first_data_byte must be > 0 for non-empty 3d segment')
            endif
        endif
        ! write updated header
        write(unit=self%funit,pos=1) self%header
        ! so no need to update header in file after this operation
    end subroutine write_segment_inside_2

    subroutine byte_manager4seg_inside_1( self, isegment, end_part1, start_part3, end_part3, bytearr_part3 )
        class(binoris),                intent(inout) :: self
        integer(kind(ENUM_ORISEG)),    intent(in)    :: isegment
        integer(kind=8),               intent(out)   :: end_part1, start_part3, end_part3
        character(len=1), allocatable, intent(out)   :: bytearr_part3(:)
        integer(kind(ENUM_ORISEG))  :: iseg
        ! READ RAW BYTES OF PART 1 OF FILE
        ! figure out byte range
        end_part1 = N_BYTES_HEADER
        if( isegment > 1 )then
            do iseg=1,isegment-1
                if( self%header(iseg)%n_records > 0 .and. self%header(iseg)%n_bytes_per_record > 0 )then
                    end_part1 = end_part1 + self%header(iseg)%n_records * self%header(iseg)%n_bytes_per_record
                endif
            end do
        endif
        ! READ RAW BYTES OF PART 3 OF FILE
        ! figure out byte range
        if( isegment < self%n_segments )then
            start_part3 = N_BYTES_HEADER
            do iseg=1,isegment
                if( self%header(iseg)%n_records > 0 .and. self%header(iseg)%n_bytes_per_record > 0 )then
                    start_part3 = start_part3 + self%header(iseg)%n_records * self%header(iseg)%n_bytes_per_record
                endif
            end do
            start_part3 = start_part3 + 1
            end_part3 = self%get_n_bytes_tot()
            ! allocate & read 3d part
            if( allocated(bytearr_part3) ) deallocate(bytearr_part3)
            allocate( bytearr_part3(start_part3:end_part3) )
            read(unit=self%funit,pos=start_part3) bytearr_part3
        endif
    end subroutine byte_manager4seg_inside_1

    subroutine byte_manager4seg_inside_2( self, isegment, end_part1, start_part3, end_part3, bytearr_part3 )
        class(binoris),                intent(inout) :: self
        integer(kind(ENUM_ORISEG)),       intent(in)    :: isegment
        integer(kind=8),               intent(in)    :: end_part1, start_part3, end_part3
        character(len=1), allocatable, intent(in)    :: bytearr_part3(:)
        integer(kind=8) :: n_bytes_part3_orig, n_bytes_part3
        integer(kind(ENUM_ORISEG)) :: iseg
        ! update byte ranges in header
        call self%update_byte_ranges
        ! validate byte ranges
        if( self%header(isegment)%first_data_byte - 1 /= end_part1 )then
            write(logfhandle,*) 'filename',trim(self%fname)
            write(logfhandle,*) 'first data byte of segment: ', self%header(isegment)%first_data_byte
            write(logfhandle,*) 'end of part 1 (bytes)     : ', end_part1
            THROW_HARD('end of part 1 of file does not match first data byte of segment')
        endif
        if( allocated(bytearr_part3) )then
            ! calculate byte size of second part of file, given the updated header
            n_bytes_part3 = 0
            do iseg=isegment + 1,self%n_segments
                if( self%header(iseg)%n_records > 0 .and. self%header(iseg)%n_bytes_per_record > 0 )then
                    n_bytes_part3 = n_bytes_part3 + self%header(iseg)%n_bytes_per_record * self%header(iseg)%n_records
                endif
            end do
            ! compare with original
            n_bytes_part3_orig = end_part3 - start_part3 + 1
            if( n_bytes_part3_orig /= n_bytes_part3 )then
                write(logfhandle,*) 'filename',trim(self%fname)
                write(logfhandle,*) '# bytes of part 3 in original: ', n_bytes_part3_orig
                write(logfhandle,*) '# bytes of part 3 in updated : ', n_bytes_part3
                THROW_HARD('byte sizes of part3 in original and updated do not match')
            endif
        endif
    end subroutine byte_manager4seg_inside_2

    subroutine write_segment_1( self, isegment, os, fromto )
        use simple_oris,   only: oris
        class(binoris),          intent(inout) :: self
        integer(kind(ENUM_ORISEG)), intent(in)    :: isegment
        class(oris),             intent(in) :: os ! indexed from 1 to nptcls
        integer, optional,       intent(in)    :: fromto(2)
        character(len=:), allocatable :: str_dyn
        real    :: ptcl_record(N_PTCL_ORIPARAMS)
        integer :: i, nspaces, noris
        integer(kind=8) :: ibytes
        noris = os%get_noris()
        if( noris == 0 ) return
        if( present(fromto) )then
            if( fromto(1) < 1 .or. fromto(2) > noris )then
                write(logfhandle,*) 'filename',trim(self%fname)
                write(logfhandle,*) 'noris : ', noris
                write(logfhandle,*) 'fromto: ', fromto
                THROW_HARD('fromto out of range')
            endif
        endif
        if( .not. self%l_open ) THROW_HARD('file needs to be open: '//trim(self%fname))
        ! add segment to stack, this sets all the information needed for allocation
        call self%add_segment_1(isegment, os, fromto)
        ! update byte ranges in header
        call self%update_byte_ranges
        ! write orientation data
        ibytes = self%header(isegment)%first_data_byte
        if( is_particle_seg(isegment) )then ! is ptcl2D or ptcl3D segment, see simple_sp_project
            do i=self%header(isegment)%fromto(1),self%header(isegment)%fromto(2)
                call os%ori2prec(i, ptcl_record)
                write(unit=self%funit,pos=ibytes) ptcl_record
                ibytes = ibytes + self%header(isegment)%n_bytes_per_record
            end do
        else
            do i=self%header(isegment)%fromto(1),self%header(isegment)%fromto(2)
                str_dyn = os%ori2str(i)
                nspaces = self%header(isegment)%n_bytes_per_record - len_trim(str_dyn)
                if( nspaces > 0 )then
                    write(unit=self%funit,pos=ibytes) str_dyn//spaces(nspaces)
                else
                    write(unit=self%funit,pos=ibytes) str_dyn
                endif
                ibytes = ibytes + self%header(isegment)%n_bytes_per_record
            end do
        endif
    end subroutine write_segment_1

    subroutine write_segment_2( self, isegment, fromto, strlen_max, sarr )
        class(binoris),             intent(inout) :: self
        integer(kind(ENUM_ORISEG)), intent(in)    :: isegment
        integer,                    intent(in)    :: fromto(2), strlen_max
        type(str4arr),              intent(inout) :: sarr(:) ! indexed from 1 to nptcls
        integer :: i, nspaces
        integer(kind=8) :: ibytes
        if( .not. self%l_open ) THROW_HARD('file needs to be open: '//trim(self%fname))
        ! add segment to stack
        call self%add_segment_2(isegment, fromto, strlen_max)
        ! update byte ranges in header
        call self%update_byte_ranges
        ! write orientation data
        ibytes = self%header(isegment)%first_data_byte
        do i=self%header(isegment)%fromto(1),self%header(isegment)%fromto(2)
            nspaces = self%header(isegment)%n_bytes_per_record - len_trim(sarr(i)%str)
            if( nspaces > 0 )then
                write(unit=self%funit,pos=ibytes) sarr(i)%str//spaces(nspaces)
            else
                write(unit=self%funit,pos=ibytes) sarr(i)%str
            endif
            ibytes = ibytes + self%header(isegment)%n_bytes_per_record
        end do
    end subroutine write_segment_2

    subroutine add_segment_1( self, isegment, os, fromto )
        use simple_oris, only: oris
        class(binoris),             intent(inout) :: self
        integer(kind(ENUM_ORISEG)), intent(in)    :: isegment
        class(oris),                intent(in)    :: os
        integer, optional,          intent(in)    :: fromto(2)
        integer :: strlen_max
        ! sanity check isegment
        if( isegment < 1 .or. isegment > MAX_N_SEGMENTS )then
            write(logfhandle,*) 'filename',trim(self%fname)
            write(logfhandle,*) 'isegment: ', isegment
            THROW_HARD('isegment out of range')
        endif
        if( isegment > self%n_segments ) self%n_segments = isegment
        ! set range in segment
        if( present(fromto) )then
            self%header(isegment)%fromto = fromto
        else
            self%header(isegment)%fromto(1) = 1
            self%header(isegment)%fromto(2) = os%get_noris()
        endif
        if( is_particle_seg(isegment) )then ! is ptcl2D or ptcl3D segment, see simple_sp_project
            ! fixed byte length per record
            self%header(isegment)%n_bytes_per_record = PTCL_BYTES_PER_REC
        else
            ! set maximum trimmed string lenght to n_bytes_per_record
            strlen_max = os%max_ori_strlen_trim()
            self%header(isegment)%n_bytes_per_record = strlen_max
        endif
        ! set n_records
        self%header(isegment)%n_records = self%header(isegment)%fromto(2) - self%header(isegment)%fromto(1) + 1
        if( self%header(isegment)%n_records < 1 ) THROW_HARD('input oritab (os) empty')
    end subroutine add_segment_1

    subroutine add_segment_2( self, isegment, fromto, strlen_max )
        class(binoris),             intent(inout) :: self
        integer(kind(ENUM_ORISEG)), intent(in)    :: isegment
        integer,                    intent(in)    :: fromto(2)
        integer,                    intent(in)    :: strlen_max
        ! sanity check isegment
        if( isegment < 1 .or. isegment > MAX_N_SEGMENTS )then
            write(logfhandle,*) 'filename',trim(self%fname)
            write(logfhandle,*) 'isegment: ', isegment
            THROW_HARD('isegment out of range')
        endif
        if( isegment > self%n_segments ) self%n_segments = isegment
        ! set range in segment
        self%header(isegment)%fromto = fromto
        ! set maximum trimmed string lenght to n_bytes_per_record
        self%header(isegment)%n_bytes_per_record = strlen_max
        ! set n_records
        self%header(isegment)%n_records = self%header(isegment)%fromto(2) - self%header(isegment)%fromto(1) + 1
    end subroutine add_segment_2

    subroutine update_byte_ranges( self )
        class(binoris), intent(inout) :: self
        integer(kind=8)         :: n_bytes_tot
        integer(kind(ENUM_ORISEG)) :: isegment
        n_bytes_tot = N_BYTES_HEADER
        if( self%n_segments <= 0 ) return
        do isegment=1,self%n_segments
            if( self%header(isegment)%n_records > 0 .and. self%header(isegment)%n_bytes_per_record > 0 )then
                self%header(isegment)%first_data_byte = n_bytes_tot + 1
                n_bytes_tot = n_bytes_tot + self%header(isegment)%n_bytes_per_record * self%header(isegment)%n_records
            endif
        end do
    end subroutine update_byte_ranges

    subroutine read_first_segment_record( self, isegment, o )
        use simple_ori, only: ori
        class(binoris),             intent(inout) :: self
        integer(kind(ENUM_ORISEG)), intent(in)    :: isegment
        class(ori),                 intent(inout) :: o
        character(len=self%header(isegment)%n_bytes_per_record) :: str_os_line ! string with static lenght (set to max(strlen))
        if( .not. self%l_open ) THROW_HARD('file needs to be open: '//trim(self%fname))
        if( isegment < 1 .or. isegment > self%n_segments )then
            write(logfhandle,*) 'filename',trim(self%fname)
            write(logfhandle,*) 'isegment: ', isegment
            write(logfhandle,*) 'n_segments: ', self%n_segments
            THROW_HARD('isegment out of range')
        endif
        if( is_particle_seg(isegment) ) THROW_HARD('Not intended for ptcl2D/3D segments')
        if( self%header(isegment)%n_records > 0 .and. self%header(isegment)%n_bytes_per_record > 0 )then
            read(unit=self%funit,pos=self%header(isegment)%first_data_byte) str_os_line
            call o%str2ori(str_os_line, is_particle_seg(isegment))
        else
            ! empty segment, nothing to do
        endif
    end subroutine read_first_segment_record

    subroutine read_segment_1( self, isegment, os, fromto, only_ctfparams_state_eo, wthreads )
        use simple_oris, only: oris
        class(binoris),             intent(inout) :: self
        integer(kind(ENUM_ORISEG)), intent(in)    :: isegment
        class(oris),                intent(inout) :: os
        integer,          optional, intent(in)    :: fromto(2)
        logical,          optional, intent(in)    :: only_ctfparams_state_eo, wthreads
        integer,       allocatable :: batches(:,:)
        character(len=self%header(isegment)%n_bytes_per_record), allocatable :: str_os_line(:) ! string with static length (set to max(strlen))
        character(len=:), allocatable ::  tmp_string
        integer(kind=8) :: ibytes, ipos
        real            :: ptcl_record(N_PTCL_ORIPARAMS)
        integer         :: fromto_here(2), i, irec, n, nl, nbatches, ibatch, nbatch, nthr
        logical         :: present_fromto, oonly_ctfparams_state_eo
        if( .not. self%l_open ) THROW_HARD('file needs to be open: '//trim(self%fname))
        if( isegment < 1 .or. isegment > self%n_segments )then
            write(logfhandle,*) 'filename',trim(self%fname)
            write(logfhandle,*) 'isegment: ', isegment
            THROW_HARD('isegment out of range')
        endif
        fromto_here = self%header(isegment)%fromto
        if( present(fromto) ) fromto_here = fromto
        if( fromto_here(1)<self%header(isegment)%fromto(1) .or. fromto_here(1)>self%header(isegment)%fromto(2) )then
            write(logfhandle,*) 'filename',trim(self%fname)
            THROW_HARD('Invalid fromto(1) index, out of range')
        endif
        if( fromto_here(2)<self%header(isegment)%fromto(1) .or. fromto_here(2)>self%header(isegment)%fromto(2) )then
            write(logfhandle,*) 'filename',trim(self%fname)
            THROW_HARD('Invalid fromto(2) index, out of range')
        endif
        nthr = 1
        if( present(wthreads) )then
            if( wthreads ) nthr = omp_get_max_threads()
        endif
        oonly_ctfparams_state_eo = .false.
        if( present(only_ctfparams_state_eo) ) oonly_ctfparams_state_eo = only_ctfparams_state_eo
        if( self%header(isegment)%n_records > 0 .and. self%header(isegment)%n_bytes_per_record > 0 )then
            if( is_particle_seg(isegment) )then ! ptcl2D/3D segment, see simple_sp_project
                if( .not.os%is_particle() ) THROW_HARD('os needs to be of particle kind')
                ibytes = self%header(isegment)%first_data_byte
                ibytes = ibytes + (fromto_here(1)-self%header(isegment)%fromto(1))*self%header(isegment)%n_bytes_per_record
                do i=fromto_here(1),fromto_here(2)
                    read(unit=self%funit,pos=ibytes) ptcl_record
                    call os%prec2ori(i, ptcl_record)
                    ibytes = ibytes + self%header(isegment)%n_bytes_per_record
                end do
            else
                ibytes = self%header(isegment)%first_data_byte
                ! read orientation data as strings
                nl       = self%header(isegment)%fromto(2) - self%header(isegment)%fromto(1) + 1
                n        = nthr * THREAD_NSTRINGS
                nbatches = ceiling(real(nl)/real(n))
                batches  = split_nobjs_even(nl,nbatches)
                do ibatch = 1,nbatches
                    nbatch = batches(ibatch,2)-batches(ibatch,1)+1
                    ! read
                    allocate(character(len=nbatch*self%header(isegment)%n_bytes_per_record) :: tmp_string)
                    read(unit=self%funit,pos=ibytes) tmp_string
                    ibytes = ibytes + nbatch*self%header(isegment)%n_bytes_per_record
                    ! parse
                    !omp parallel do default(shared) private(i,irec,ipos) schedule(static) proc_bind(close)
                    do i = 1,nbatch
                        irec  = batches(ibatch,1) + i - 1
                        ipos  = (i-1) * self%header(isegment)%n_bytes_per_record + 1
                        if( oonly_ctfparams_state_eo )then
                            if( present_fromto )then
                                call os%str2ori_ctfparams_state_eo(i,&
                                    &tmp_string(ipos:ipos+self%header(isegment)%n_bytes_per_record-1))
                            else
                                call os%str2ori_ctfparams_state_eo(irec,&
                                    &tmp_string(ipos:ipos+self%header(isegment)%n_bytes_per_record-1))
                            endif
                        else
                            if( present_fromto )then
                                call os%str2ori(i, tmp_string(ipos:ipos+self%header(isegment)%n_bytes_per_record-1))
                            else
                                call os%str2ori(irec, tmp_string(ipos:ipos+self%header(isegment)%n_bytes_per_record-1))
                            endif
                        endif
                    enddo
                    !omp end parallel do
                    deallocate(tmp_string)
                enddo
            endif
        else
            ! empty segment, nothing to do
        endif
    end subroutine read_segment_1

    subroutine read_segment_2( self, isegment, sarr )
        class(binoris),             intent(inout) :: self
        integer(kind(ENUM_ORISEG)), intent(in)    :: isegment
        type(str4arr),              intent(inout) :: sarr(:)
        integer :: i
        integer(kind=8) :: ibytes
        if( .not. self%l_open ) THROW_HARD('file needs to be open: '//trim(self%fname))
        if( isegment < 1 .or. isegment > self%n_segments )then
            write(logfhandle,*) 'filename',trim(self%fname)
            write(logfhandle,*) 'isegment: ', isegment
            THROW_HARD('isegment out of range')
        endif
        if( self%header(isegment)%n_records > 0 .and. self%header(isegment)%n_bytes_per_record > 0 )then
            ! read orientation data into array of allocatable strings
            ibytes = self%header(isegment)%first_data_byte
            do i=self%header(isegment)%fromto(1),self%header(isegment)%fromto(2)
                allocate(character(len=self%header(isegment)%n_bytes_per_record) :: sarr(i)%str)
                read(unit=self%funit,pos=ibytes) sarr(i)%str
                ibytes = ibytes + self%header(isegment)%n_bytes_per_record
            end do
        else
            ! empty segment, nothing to do
        endif
    end subroutine read_segment_2

    subroutine read_record( self, isegment, ibytes, str )
        use simple_ori, only: ori
        class(binoris),                intent(inout) :: self
        integer(kind(ENUM_ORISEG)),    intent(in)    :: isegment
        integer(kind=8),               intent(inout) :: ibytes
        character(len=:), allocatable, intent(inout) :: str
        type(ori) :: o
        real      :: ptcl_record(N_PTCL_ORIPARAMS)
        if( allocated(str) ) deallocate(str)
        if( is_particle_seg(isegment) )then ! ptcl2D/3D segment, see simple_sp_project
            call o%new(is_ptcl=.true.)
            read(unit=self%funit,pos=ibytes) ptcl_record
            call o%prec2ori(ptcl_record)
            str = o%ori2str()
        else
            allocate(character(len=self%header(isegment)%n_bytes_per_record) :: str)
            read(unit=self%funit,pos=ibytes) str
        endif
        ibytes = ibytes + self%header(isegment)%n_bytes_per_record
    end subroutine read_record

    ! getters

    subroutine get_segments_info( self, seginds, info_struct )
        class(binoris),                     intent(in)  :: self
        integer,               allocatable, intent(out) :: seginds(:)
        type(binoris_seginfo), allocatable, intent(out) :: info_struct(:)
        integer(kind(ENUM_ORISEG)) :: isegment, cnt
        if( allocated(seginds)     ) deallocate(seginds)
        if( allocated(info_struct) ) deallocate(info_struct)
        if( self%n_segments <= 0 ) return
        ! count # populated segments
        cnt = 0
        do isegment=1,MAX_N_SEGMENTS
            if(   self%header(isegment)%n_records > 0&
            .and. self%header(isegment)%n_bytes_per_record > 0 )then
                cnt = cnt + 1
            endif
        end do
        ! allocate arrays
        allocate(info_struct(cnt), seginds(cnt))
        cnt = 0
        do isegment=1,MAX_N_SEGMENTS
            if(   self%header(isegment)%n_records > 0&
            .and. self%header(isegment)%n_bytes_per_record > 0 )then
                ! segment is populateed
                cnt = cnt + 1
                seginds(cnt)                        = isegment
                info_struct(cnt)%fromto             = self%header(isegment)%fromto
                info_struct(cnt)%n_records          = self%header(isegment)%n_records
                info_struct(cnt)%n_bytes_per_record = self%header(isegment)%n_bytes_per_record
                info_struct(cnt)%first_data_byte    = self%header(isegment)%first_data_byte
            endif
        end do
    end subroutine get_segments_info

    pure integer function get_n_segments( self )
        class(binoris), intent(in) :: self
        get_n_segments = self%n_segments
    end function get_n_segments

    pure function get_fromto( self, isegment ) result( fromto )
        class(binoris),             intent(in) :: self
        integer(kind(ENUM_ORISEG)), intent(in) :: isegment
        integer :: fromto(2)
        fromto = self%header(isegment)%fromto
    end function get_fromto

    pure integer function get_n_records( self, isegment )
        class(binoris),             intent(in) :: self
        integer(kind(ENUM_ORISEG)), intent(in) :: isegment
        get_n_records = self%header(isegment)%n_records
    end function get_n_records

    pure integer function get_n_bytes_per_record( self, isegment )
        class(binoris),             intent(in) :: self
        integer(kind(ENUM_ORISEG)), intent(in) :: isegment
        get_n_bytes_per_record = self%header(isegment)%n_bytes_per_record
    end function get_n_bytes_per_record

    pure integer(kind=8) function get_n_bytes_tot( self )
        class(binoris), intent(in) :: self
        integer(kind(ENUM_ORISEG)) :: isegment
        get_n_bytes_tot = N_BYTES_HEADER
        if( self%n_segments <= 0 ) return
        do isegment=1,self%n_segments
            if( self%header(isegment)%n_records > 0 .and. self%header(isegment)%n_bytes_per_record > 0 )then
                get_n_bytes_tot = get_n_bytes_tot + self%header(isegment)%n_bytes_per_record * self%header(isegment)%n_records
            endif
        end do
    end function get_n_bytes_tot

    pure integer(kind=8) function get_first_data_byte( self, isegment )
        class(binoris),          intent(in) :: self
        integer(kind(ENUM_ORISEG)), intent(in) :: isegment
        get_first_data_byte = self%header(isegment)%first_data_byte
    end function get_first_data_byte

    logical function is_opened( self )
        class(binoris), intent(in) :: self
        is_opened = self%l_open
    end function is_opened

    logical function is_particle_seg( isegment )
        integer, intent(in) :: isegment
        is_particle_seg = .false.
        ! In simple_sp_project:
        ! type(oris) :: os_ptcl2D ! per-particle 2D os, segment 3
        ! type(oris) :: os_ptcl3D ! per-particle 3D os, segment 6
        if( isegment == 3 .or. isegment == 6 ) is_particle_seg = .true.
    end function is_particle_seg

end module simple_binoris
