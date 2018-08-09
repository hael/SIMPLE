! for manging orientation data using binary files
module simple_binoris
use, intrinsic :: ISO_C_BINDING
use simple_defs     ! use all in there
use simple_strings, only: spaces
use simple_error,   only: allocchk
use simple_fileio,  only: fopen, fclose, fileiochk, funit_size
use simple_syslib,  only: del_file, file_exists
implicit none

public :: binoris
private

integer(kind=8), parameter :: MAX_N_SEGEMENTS = 20
integer(kind=8), parameter :: N_VARS_HEAD_SEG = 5
integer(kind=8), parameter :: N_BYTES_HEADER  = MAX_N_SEGEMENTS * N_VARS_HEAD_SEG * 8 ! because dp integer

#include "simple_local_flags.inc"

type file_header_segment
    integer(kind=8)   :: fromto(2)          = 0
    integer(kind=8)   :: n_bytes_per_record = 0
    integer(kind=8)   :: n_records          = 0
    integer(kind=8)   :: first_data_byte    = 0
end type file_header_segment

type binoris
    private
    type(file_header_segment) :: header(MAX_N_SEGEMENTS)
    integer                   :: n_segments  = 0
    integer                   :: funit       = 0
    logical                   :: l_open      = .false.
  contains
    ! I/O
    procedure          :: open
    procedure, private :: clear_segments
    procedure          :: close
    procedure, private :: read_header
    procedure          :: write_header
    procedure          :: print_header
    procedure, private :: write_segment_1
    procedure, private :: write_segment_2
    generic            :: write_segment => write_segment_1, write_segment_2
    procedure, private :: write_segment_inside_1
    procedure, private :: write_segment_inside_2
    generic            :: write_segment_inside => write_segment_inside_1, write_segment_inside_2
    procedure, private :: byte_manager4seg_inside_1
    procedure, private :: byte_manager4seg_inside_2
    procedure, private :: add_segment_1
    procedure, private :: add_segment_2
    procedure, private :: update_byte_ranges
    procedure          :: read_first_segment_record
    procedure          :: read_segment_1
    procedure          :: read_segment_2
    generic            :: read_segment => read_segment_1, read_segment_2
    ! getters
    procedure          :: get_n_segments
    procedure          :: get_fromto
    procedure          :: get_n_records
    procedure          :: get_n_bytes_per_record
    procedure          :: get_n_bytes_tot
    procedure          :: is_opened
end type binoris

contains

    ! I/O

    subroutine open( self, fname, del_if_exists )
        class(binoris),    intent(inout) :: self          !< instance
        character(len=*),  intent(in)    :: fname         !< filename
        logical, optional, intent(in)    :: del_if_exists !< If the file already exists on disk, replace
        integer(kind=8)    :: filesz
        integer            :: isegment
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
        ! check size
        filesz = funit_size(self%funit)
        if( filesz == -1 )then
            stop 'file_size cannot be inquired; binoris :: open'
        else if( filesz >= N_BYTES_HEADER )then
            ! ok
        else
            write(*,*) 'file: ', trim(fname)
            stop 'file_size too small to contain a header; binoris :: open'
        endif
        ! clear segments before reading header
        call self%clear_segments
        call self%read_header
        self%n_segments = 0 ! for counting # segments
        do isegment=1,MAX_N_SEGEMENTS
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
            call fclose(self%funit,io_stat,errmsg='binoris ; close ')
            self%l_open = .false.
        end if
    end subroutine close

    subroutine read_header( self )
        class(binoris), intent(inout) :: self  !< instance
        integer :: io_status
        if( .not. self%l_open ) stop 'file needs to be open; binoris :: read_header'
        read(unit=self%funit,pos=1,iostat=io_status) self%header
        if( io_status .ne. 0 ) call fileiochk('binoris ::read_header, ERROR reading header bytes ', io_status)
    end subroutine read_header

    subroutine write_header( self )
        class(binoris), intent(inout) :: self  !< instance
        integer :: io_status
        if( .not. self%l_open ) stop 'file needs to be open; binoris :: write_header'
        write(unit=self%funit,pos=1,iostat=io_status) self%header
        if( io_status .ne. 0 ) call fileiochk('binoris :: write_header, ERROR writing header bytes ', io_status)
        DebugPrint  'wrote: ', sizeof(self%header), ' header bytes'
    end subroutine write_header

    subroutine print_header( self )
        class(binoris), intent(in) :: self
        integer :: isegment
        do isegment=1,MAX_N_SEGEMENTS
            write(*,*) '*****  HEADER, segment: ', isegment
            write(*,*) 'fromto(1)         : ', self%header(isegment)%fromto(1)
            write(*,*) 'fromto(2)         : ', self%header(isegment)%fromto(2)
            write(*,*) 'n_bytes_per_record: ', self%header(isegment)%n_bytes_per_record
            write(*,*) 'n_records         : ', self%header(isegment)%n_records
            write(*,*) 'first_data_byte   : ', self%header(isegment)%first_data_byte
        end do
    end subroutine print_header

    subroutine write_segment_inside_1( self, isegment, os, fromto )
        use simple_oris,   only: oris
        class(binoris),        intent(inout) :: self
        integer,               intent(in)    :: isegment
        class(oris), optional, intent(inout) :: os ! indexed from 1 to nptcls
        integer,     optional, intent(in)    :: fromto(2)
        character(len=:), allocatable :: str_dyn
        integer         :: i, nspaces, noris, iseg
        integer(kind=8) :: end_part1, start_part3, end_part3
        integer(kind=8) :: ibytes, first_data_byte
        character(len=1), allocatable ::bytearr_part3(:)
        noris = os%get_noris()
        if( noris == 0 ) return
        if( present(fromto) )then
            if( fromto(1) < 1 .or. fromto(2) > noris )then
                write(*,*) 'noris : ', noris
                write(*,*) 'fromto: ', fromto
                stop 'fromto out of range; binoris :: write_segment_inside'
            endif
        endif
        if( .not. self%l_open ) stop 'file needs to be open; binoris :: write_segment_inside'
        ! ranges and raw bytes
        call self%byte_manager4seg_inside_1(isegment, end_part1, start_part3, end_part3, bytearr_part3)
        ! add segment to stack, this sets all the information needed for allocation
        call self%add_segment_1(isegment, os, fromto)
        ! error checks
        call self%byte_manager4seg_inside_2(isegment, end_part1, start_part3, end_part3, bytearr_part3)
        ! WRITE FILE
        ! write orientation data (2nd part)
        ibytes = self%header(isegment)%first_data_byte
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
                write(*,*) 'first_data_byte: ', first_data_byte
                write(*,*) 'ERROR! first_data_byte must be > 0 for non-empty 3d segment'
                stop 'simple_binoris :: write_segment_inside'
            endif
        endif
        ! write updated header
        write(unit=self%funit,pos=1) self%header
        ! so no need to update header in file after this operation
    end subroutine write_segment_inside_1

    subroutine write_segment_inside_2( self, isegment, os_strings, fromto, strlen_max )
        use simple_oris,   only: oris
        class(binoris), intent(inout) :: self
        integer,        intent(in)    :: isegment
        type(str4arr),  intent(in)    :: os_strings(:)
        integer,        intent(in)    :: fromto(2), strlen_max
        character(len=:), allocatable :: str_dyn
        integer         :: i, nspaces, noris, iseg
        integer(kind=8) :: end_part1, start_part3, end_part3
        integer(kind=8) :: ibytes, first_data_byte
        character(len=1), allocatable :: bytearr_part3(:)
        noris = size(os_strings)
        if( noris == 0 ) return
        if( .not. self%l_open ) stop 'file needs to be open; binoris :: write_segment_inside'
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
                write(unit=self%funit,pos=ibytes) str_dyn//spaces(nspaces)
            else
                write(unit=self%funit,pos=ibytes) str_dyn
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
                write(*,*) 'first_data_byte: ', first_data_byte
                write(*,*) 'ERROR! first_data_byte must be > 0 for non-empty 3d segment'
                stop 'simple_binoris :: write_segment_inside'
            endif
        endif
        ! write updated header
        write(unit=self%funit,pos=1) self%header
        ! so no need to update header in file after this operation
    end subroutine write_segment_inside_2

    subroutine byte_manager4seg_inside_1( self, isegment, end_part1, start_part3, end_part3, bytearr_part3 )
        class(binoris),                intent(inout) :: self
        integer,                       intent(in)    :: isegment
        integer(kind=8),               intent(out)   :: end_part1, start_part3, end_part3
        character(len=1), allocatable, intent(out)   :: bytearr_part3(:)
        integer :: iseg
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
        integer,                       intent(in)    :: isegment
        integer(kind=8),               intent(in)    :: end_part1, start_part3, end_part3
        character(len=1), allocatable, intent(in)    :: bytearr_part3(:)
        integer(kind=8) :: n_bytes_part3_orig, n_bytes_part3
        integer         :: iseg
        ! update byte ranges in header
        call self%update_byte_ranges
        ! validate byte ranges
        if( self%header(isegment)%first_data_byte - 1 /= end_part1 )then
            write(*,*) 'first data byte of segment: ', self%header(isegment)%first_data_byte
            write(*,*) 'end of part 1 (bytes)     : ', end_part1
            write(*,*) 'ERROR! end of part 1 of file does not match first data byte of segment'
            stop 'simple_binoris :: byte_manager4seg_inside_2'
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
                write(*,*) '# bytes of part 3 in original: ', n_bytes_part3_orig
                write(*,*) '# bytes of part 3 in updated : ', n_bytes_part3
                write(*,*) 'ERROR! byte sizes of part3 in original and updated do not match'
                stop 'simple_binoris :: byte_manager4seg_inside_2'
            endif
        endif
    end subroutine byte_manager4seg_inside_2

    subroutine write_segment_1( self, isegment, os, fromto )
        use simple_oris,   only: oris
        class(binoris),    intent(inout) :: self
        integer,           intent(in)    :: isegment
        class(oris),       intent(inout) :: os ! indexed from 1 to nptcls
        integer, optional, intent(in)    :: fromto(2)
        character(len=:), allocatable :: str_dyn
        integer :: i, nspaces, noris
        integer(kind=8) :: ibytes
        noris = os%get_noris()
        if( noris == 0 ) return
        if( present(fromto) )then
            if( fromto(1) < 1 .or. fromto(2) > noris )then
                write(*,*) 'noris : ', noris
                write(*,*) 'fromto: ', fromto
                stop 'fromto out of range; binoris :: write_segment_1'
            endif
        endif
        if( .not. self%l_open ) stop 'file needs to be open; binoris :: write_segment_1'
        ! add segment to stack, this sets all the information needed for allocation
        call self%add_segment_1(isegment, os, fromto)
        ! update byte ranges in header
        call self%update_byte_ranges
        ! write orientation data
        ibytes = self%header(isegment)%first_data_byte
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
    end subroutine write_segment_1

    subroutine write_segment_2( self, isegment, fromto, strlen_max, sarr )
        class(binoris), intent(inout) :: self
        integer,        intent(in)    :: isegment
        integer,        intent(in)    :: fromto(2), strlen_max
        type(str4arr),  intent(inout) :: sarr(:) ! indexed from 1 to nptcls
        integer :: i, nspaces
        integer(kind=8) :: ibytes
        if( .not. self%l_open ) stop 'file needs to be open; binoris :: write_segment_2'
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
                DebugPrint 'wrote: ', sizeof(sarr(i)%str//spaces(nspaces)),&
                    &' segment: ', isegment, ' bytes, starting @: ', ibytes
            else
                write(unit=self%funit,pos=ibytes) sarr(i)%str
                DebugPrint 'wrote: ', sizeof(sarr(i)%str),&
                    &' segment: ', isegment, ' bytes, starting @: ', ibytes
            endif
            ibytes = ibytes + self%header(isegment)%n_bytes_per_record
        end do
    end subroutine write_segment_2

    subroutine add_segment_1( self, isegment, os, fromto )
        use simple_oris,   only: oris
        class(binoris),    intent(inout) :: self
        integer,           intent(in)    :: isegment
        class(oris),       intent(inout) :: os
        integer, optional, intent(in)    :: fromto(2)
        integer :: strlen_max
        ! sanity check isegment
        if( isegment < 1 .or. isegment > MAX_N_SEGEMENTS )then
            write(*,*) 'isegment: ', isegment
            stop 'ERROR, isegment out of range; binoris :: add_segment_1'
        endif
        if( isegment > self%n_segments ) self%n_segments = isegment
        ! set range in segment
        if( present(fromto) )then
            self%header(isegment)%fromto = fromto
        else
            self%header(isegment)%fromto(1) = 1
            self%header(isegment)%fromto(2) = os%get_noris()
        endif
        ! set maximum trimmed string lenght to n_bytes_per_record
        strlen_max = os%max_ori_strlen_trim()
        self%header(isegment)%n_bytes_per_record = strlen_max
        ! set n_records
        self%header(isegment)%n_records = self%header(isegment)%fromto(2) - self%header(isegment)%fromto(1) + 1
        if( self%header(isegment)%n_records < 1 ) stop 'ERROR, input oritab (os) empty; binoris :: add_segment_1'
    end subroutine add_segment_1

    subroutine add_segment_2( self, isegment, fromto, strlen_max )
        class(binoris), intent(inout) :: self
        integer,        intent(in)    :: isegment
        integer,        intent(in)    :: fromto(2)
        integer,        intent(in)    :: strlen_max
        ! sanity check isegment
        if( isegment < 1 .or. isegment > MAX_N_SEGEMENTS )then
            write(*,*) 'isegment: ', isegment
            stop 'ERROR, isegment out of range; binoris :: add_segment_2'
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
        integer(kind=8) :: n_bytes_tot
        integer :: isegment
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
        class(binoris), intent(inout) :: self
        integer,        intent(in)    :: isegment
        class(ori),     intent(inout) :: o
        character(len=self%header(isegment)%n_bytes_per_record) :: str_os_line ! string with static lenght (set to max(strlen))
        if( .not. self%l_open ) stop 'file needs to be open; binoris :: read_first_segment_record'
        if( isegment < 1 .or. isegment > self%n_segments )then
            write(*,*) 'isegment: ', isegment
            write(*,*) 'n_segments: ', self%n_segments
            stop 'isegment out of bound; binoris :: read_first_segment_record'
        endif
        if( self%header(isegment)%n_records > 0 .and. self%header(isegment)%n_bytes_per_record > 0 )then
            read(unit=self%funit,pos=self%header(isegment)%first_data_byte) str_os_line
            call o%str2ori(str_os_line)
        else
            ! empty segment, nothing to do
        endif
    end subroutine read_first_segment_record

    subroutine read_segment_1( self, isegment, os, fromto, only_ctfparams_state_eo )
        use simple_oris,   only: oris
        class(binoris),    intent(inout) :: self
        integer,           intent(in)    :: isegment
        class(oris),       intent(inout) :: os
        integer, optional, intent(in) :: fromto(2)
        logical, optional, intent(in) :: only_ctfparams_state_eo
        character(len=self%header(isegment)%n_bytes_per_record) :: str_os_line ! string with static lenght (set to max(strlen))
        integer :: i, irec
        integer(kind=8) :: ibytes
        logical :: present_fromto, oonly_ctfparams_state_eo
        if( .not. self%l_open ) stop 'file needs to be open; binoris :: read_segment_1'
        if( isegment < 1 .or. isegment > self%n_segments )then
            write(*,*) 'isegment: ', isegment
            stop 'isegment out of bound; binoris :: read_segment_1'
        endif
        present_fromto = present(fromto)
        if( present_fromto )then
            if( .not. all(fromto .eq. self%header(isegment)%fromto) )&
                &stop 'passed dummy fromto not consistent with self%header; binoris :: read_segment_1'
        endif
        oonly_ctfparams_state_eo = .false.
        if( present(only_ctfparams_state_eo) ) oonly_ctfparams_state_eo = only_ctfparams_state_eo
        if( self%header(isegment)%n_records > 0 .and. self%header(isegment)%n_bytes_per_record > 0 )then
            ! read orientation data
            ibytes = self%header(isegment)%first_data_byte
            irec   = 0
            do i=self%header(isegment)%fromto(1),self%header(isegment)%fromto(2)
                irec = irec + 1
                read(unit=self%funit,pos=ibytes) str_os_line
                if( oonly_ctfparams_state_eo )then
                    if( present_fromto )then
                        call os%str2ori_ctfparams_state_eo(i, str_os_line)
                    else
                        call os%str2ori_ctfparams_state_eo(irec, str_os_line)
                    endif
                else
                    if( present_fromto )then
                        call os%str2ori(i, str_os_line)
                    else
                        call os%str2ori(irec, str_os_line)
                    endif
                endif
                ibytes = ibytes + self%header(isegment)%n_bytes_per_record
            end do
        else
            ! empty segment, nothing to do
        endif
    end subroutine read_segment_1

    subroutine read_segment_2( self, isegment, sarr )
        class(binoris), intent(inout) :: self
        integer,        intent(in)    :: isegment
        type(str4arr),  intent(inout) :: sarr(:)
        integer :: i
        integer(kind=8) :: ibytes
        if( .not. self%l_open ) stop 'file needs to be open; binoris :: read_segment_2'
        if( isegment < 1 .or. isegment > self%n_segments )then
            write(*,*) 'isegment: ', isegment
            stop 'isegment out of bound; binoris :: read_segment_2'
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

    ! getters

    pure integer function get_n_segments( self )
        class(binoris), intent(in) :: self
        get_n_segments = self%n_segments
    end function get_n_segments

    pure function get_fromto( self, isegment ) result( fromto )
        class(binoris), intent(in) :: self
        integer,        intent(in) :: isegment
        integer :: fromto(2)
        fromto = self%header(isegment)%fromto
    end function get_fromto

    pure integer function get_n_records( self, isegment )
        class(binoris), intent(in) :: self
        integer,        intent(in) :: isegment
        get_n_records = self%header(isegment)%n_records
    end function get_n_records

    pure integer function get_n_bytes_per_record( self, isegment )
        class(binoris), intent(in) :: self
        integer,        intent(in) :: isegment
        get_n_bytes_per_record = self%header(isegment)%n_bytes_per_record
    end function get_n_bytes_per_record

    pure integer(kind=8) function get_n_bytes_tot( self )
        class(binoris), intent(in) :: self
        integer :: isegment
        get_n_bytes_tot = N_BYTES_HEADER
        if( self%n_segments <= 0 ) return
        do isegment=1,self%n_segments
            if( self%header(isegment)%n_records > 0 .and. self%header(isegment)%n_bytes_per_record > 0 )then
                get_n_bytes_tot = get_n_bytes_tot + self%header(isegment)%n_bytes_per_record * self%header(isegment)%n_records
            endif
        end do
    end function get_n_bytes_tot

    logical function is_opened( self )
        class(binoris), intent(in) :: self
        is_opened = self%l_open
    end function is_opened

end module simple_binoris
