! for manging orientation data using binary files
module simple_binoris
use simple_defs   ! use all in there
use simple_fileio ! use all in there
use simple_oris,  only: oris
#include "simple_lib.f08"
implicit none

public :: binoris
private

integer(kind=8), parameter :: MAX_N_SEGEMENTS = 20
integer(kind=8), parameter :: N_VARS_HEAD_SEG = 5
integer(kind=8), parameter :: N_BYTES_HEADER  = MAX_N_SEGEMENTS * N_VARS_HEAD_SEG * 8 ! because dp integer
logical,         parameter :: DEBUG = .false.

type file_header_segment
    integer(kind=8)   :: fromto(2)          = 0
    integer(kind=8)   :: n_bytes_per_record = 0
    integer(kind=8)   :: n_records          = 0
    integer(kind=8)   :: first_data_byte    = 0
end type file_header_segment

type binoris
    private
    type(file_header_segment) :: header(MAX_N_SEGEMENTS)
    integer                   :: n_segments = 0
    integer                   :: funit      = 0
    logical                   :: l_open     = .false.
  contains
    ! I/O
    procedure          :: open
    procedure, private :: clear_segments
    procedure          :: close
    procedure          :: read_header
    procedure          :: write_header
    procedure          :: print_header
    procedure, private :: write_segment_1
    procedure, private :: write_segment_2
    generic            :: write_segment => write_segment_1, write_segment_2
    procedure, private :: add_segment_1
    procedure, private :: add_segment_2

    procedure, private :: update_byte_ranges
    procedure          :: read_segment_1
    procedure          :: read_segment_2
    generic            :: read_segment => read_segment_1, read_segment_2
    procedure          :: read_segment_ctfparams_state_eo
    ! getters
    procedure          :: get_n_segments
    procedure          :: get_fromto
    procedure          :: get_n_records
end type binoris

contains

    ! I/O

    subroutine open( self, fname, del_if_exists )
        class(binoris),    intent(inout) :: self          !< instance
        character(len=*),  intent(in)    :: fname         !< filename
        logical, optional, intent(in)    :: del_if_exists !< If the file already exists on disk, replace
        integer(kind=8)    :: filesz
        integer            :: io_status, isegment
        character(len=512) :: io_message
        ! deletion logics
        if( present(del_if_exists) )then
            if( del_if_exists )then
                call del_file(trim(fname))
            endif
        endif
        ! existence logics
        if( .not. file_exists(trim(fname)) )then
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
                &.and. self%header(isegment)%first_data_byte > 0 ) self%n_segments = isegment ! to allow empty in-between segments
        end do

        contains

            subroutine open_local
                integer :: io_stat, tmpunit
                if( .not. self%l_open )then
                    call fopen(tmpunit, trim(fname), access='STREAM', action='READWRITE',&
                        &status='UNKNOWN', form='UNFORMATTED', iostat=io_stat)
                    call fileio_errmsg('binoris ; open_local '//trim(fname), io_stat)
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
        if( io_status .ne. 0 ) call fileio_errmsg('binoris ::read_header, ERROR reading header bytes ', io_status)
    end subroutine read_header

    subroutine write_header( self )
        class(binoris), intent(inout) :: self  !< instance
        integer :: io_status
        if( .not. self%l_open ) stop 'file needs to be open; binoris :: write_header'
        write(unit=self%funit,pos=1,iostat=io_status) self%header
        if( io_status .ne. 0 ) call fileio_errmsg('binoris :: write_header, ERROR writing header bytes ', io_status)
        if( DEBUG ) print *, 'wrote: ', sizeof(self%header), ' header bytes'
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

    subroutine write_segment_1( self, isegment, os, fromto )
        class(binoris),    intent(inout) :: self
        integer,           intent(in)    :: isegment
        class(oris),       intent(inout) :: os ! indexed from 1 to nptcls
        integer, optional, intent(in)    :: fromto(2)
        character(len=:), allocatable :: str_dyn
        integer :: i, ibytes, io_status, nspaces, noris
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
            nspaces = self%header(isegment)%n_bytes_per_record - len(str_dyn)
            if( nspaces > 0 )then
                write(unit=self%funit,pos=ibytes) str_dyn//spaces(nspaces)
                if( DEBUG ) print *, 'wrote: ', sizeof(str_dyn//spaces(nspaces)),&
                    &'bytes, segment: ', isegment, ' bytes, starting @: ', ibytes
            else
                write(unit=self%funit,pos=ibytes) str_dyn
                if( DEBUG ) print *, 'wrote: ', sizeof(str_dyn),&
                    &'bytes, segment: ', isegment, ' bytes, starting @: ', ibytes
            endif
            ibytes = ibytes + self%header(isegment)%n_bytes_per_record
        end do
    end subroutine write_segment_1

    subroutine write_segment_2( self, isegment, fromto, strlen_max, sarr )
        class(binoris), intent(inout) :: self
        integer,        intent(in)    :: isegment
        integer,        intent(in)    :: fromto(2), strlen_max
        type(str4arr),  intent(inout) :: sarr(:) ! indexed from 1 to nptcls
        integer :: i, ibytes, io_status, nspaces
        if( .not. self%l_open ) stop 'file needs to be open; binoris :: write_segment_2'
        ! add segment to stack
        call self%add_segment_2(isegment, fromto, strlen_max)
        ! update byte ranges in header
        call self%update_byte_ranges
        ! write orientation data
        ibytes = self%header(isegment)%first_data_byte
        do i=self%header(isegment)%fromto(1),self%header(isegment)%fromto(2)
            nspaces = self%header(isegment)%n_bytes_per_record - len(sarr(i)%str)
            if( nspaces > 0 )then
                write(unit=self%funit,pos=ibytes) sarr(i)%str//spaces(nspaces)
                if( DEBUG ) print *, 'wrote: ', sizeof(sarr(i)%str//spaces(nspaces)),&
                    &' segment: ', isegment, ' bytes, starting @: ', ibytes
            else
                write(unit=self%funit,pos=ibytes) sarr(i)%str
                if( DEBUG ) print *, 'wrote: ', sizeof(sarr(i)%str),&
                    &' segment: ', isegment, ' bytes, starting @: ', ibytes
            endif
            ibytes = ibytes + self%header(isegment)%n_bytes_per_record
        end do
    end subroutine write_segment_2

    subroutine add_segment_1( self, isegment, os, fromto )
        class(binoris),    intent(inout) :: self
        integer,           intent(in)    :: isegment
        class(oris),       intent(inout) :: os
        integer, optional, intent(in)    :: fromto(2)
        integer :: strlen_max
        ! sanity check isegment
        if( isegment < 1 .or. isegment > MAX_N_SEGEMENTS ) stop 'ERROR, isegment out of range; binoris :: add_segment_1'
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
        if( isegment < 1 .or. isegment > MAX_N_SEGEMENTS ) stop 'ERROR, isegment out of range; binoris :: add_segment_2'
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

    subroutine read_segment_1( self, isegment, os, fromto )
        class(binoris), intent(inout) :: self
        integer,        intent(in)    :: isegment
        class(oris),    intent(inout) :: os
        integer, optional, intent(in) :: fromto(2)
        character(len=self%header(isegment)%n_bytes_per_record) :: str_os_line ! string with static lenght (set to max(strlen))
        integer :: i, ibytes, irec
        logical :: present_fromto
        if( .not. self%l_open ) stop 'file needs to be open; binoris :: read_segment_1'
        if( isegment < 1 .or. isegment > self%n_segments ) stop 'isegment out of bound; binoris :: read_segment_1'
        present_fromto = present(fromto)
        if( present_fromto )then
            if( .not. all(fromto .eq. self%header(isegment)%fromto) )&
                &stop 'passed dummy fromto not consistent with self%header; binoris :: read_segment_1'
        endif
        if( self%header(isegment)%n_records > 0 .and. self%header(isegment)%n_bytes_per_record > 0 )then
            ! read orientation data
            ibytes = self%header(isegment)%first_data_byte
            irec   = 0
            do i=self%header(isegment)%fromto(1),self%header(isegment)%fromto(2)
                irec = irec + 1
                read(unit=self%funit,pos=ibytes) str_os_line
                if( present(fromto) )then
                    call os%str2ori(i, str_os_line)
                else
                    call os%str2ori(irec, str_os_line)
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
        integer :: i, ibytes
        if( .not. self%l_open ) stop 'file needs to be open; binoris :: read_segment_2'
        if( isegment < 1 .or. isegment > self%n_segments ) stop 'isegment out of bound; binoris :: read_segment_2'
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

    subroutine read_segment_ctfparams_state_eo( self, isegment, os )
        use simple_ori, only: ori
        class(binoris), intent(inout) :: self
        integer,        intent(in)    :: isegment
        class(oris),    intent(inout) :: os
        integer, parameter :: NFLAGS = 11
        character(len=32)  :: flags(NFLAGS)
        type(ori)          :: o
        character(len=self%header(isegment)%n_bytes_per_record) :: str_os_line ! string with static lenght (set to max(strlen))
        integer :: i, j, ibytes
        if( .not. self%l_open ) stop 'file needs to be open; binoris :: read_segment_ctfparams_state_eo'
        if( isegment < 1 .or. isegment > self%n_segments ) stop 'isegment out of bound; binoris :: read_segment_ctfparams_state_eo'
        if( self%header(isegment)%n_records > 0 .and. self%header(isegment)%n_bytes_per_record > 0 )then
            ! set flags for ctfparams, state & eo
            flags(1)  = 'smpd'
            flags(2)  = 'kv'
            flags(3)  = 'cs'
            flags(4)  = 'fraca'
            flags(5)  = 'dfx'
            flags(6)  = 'dfy'
            flags(7)  = 'angast'
            flags(8)  = 'bfac'
            flags(9)  = 'state'
            flags(10) = 'eo'
            flags(11) = 'phshift'
            ! read orientation data
            ibytes = self%header(isegment)%first_data_byte
            do i=self%header(isegment)%fromto(1),self%header(isegment)%fromto(2)
                read(unit=self%funit,pos=ibytes) str_os_line
                call o%str2ori(str_os_line)
                do j=1,NFLAGS
                    if( o%isthere(trim(flags(j))) ) call os%set(i, trim(flags(j)), o%get(trim(flags(j))))
                end do
                ibytes = ibytes + self%header(isegment)%n_bytes_per_record
            end do
        else
            ! empty segment, nothing to do
        endif
    end subroutine read_segment_ctfparams_state_eo

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

end module simple_binoris
