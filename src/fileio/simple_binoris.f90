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
integer(kind=8), parameter :: N_BYTES_HEADER  = MAX_N_SEGEMENTS * 5 * 8
! max(# segments) * # kind=8 variables in header segment * bytes per variable + string lenght of descriptor

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
    procedure          :: close
    procedure          :: write_header
    procedure          :: print_header
    procedure          :: write_segment
    procedure          :: read_segment
    procedure          :: read_segment_ctfparams_state_eo
    ! getters
    procedure          :: get_n_segments
    procedure          :: get_fromto
    procedure          :: get_n_records
    ! helper routines
    procedure, private :: clear_segments
    procedure, private :: add_segment
    procedure, private :: update_byte_ranges
end type binoris

contains

    ! I/O

    subroutine open( self, fname, del_if_exists )
        class(binoris),    intent(inout) :: self          !< instance
        character(len=*),  intent(in)    :: fname         !< filename
        logical, optional, intent(in)    :: del_if_exists !< If the file already exists on disk, replace
        integer(kind=1)    :: byte_array(N_BYTES_HEADER)
        integer(kind=8)    :: filesz
        integer            :: io_status, isegment, first_byte, sz
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
            stop 'file_size too small to contain a header; binoris :: open'
        endif
        ! clear segments before reading header
        call self%clear_segments
        ! read header
        read(unit=self%funit,pos=1,iostat=io_status,iomsg=io_message) byte_array
        if( io_status .ne. 0 ) call fileio_errmsg('binoris :: open, ERROR reading header records '//trim(io_message), io_status)
        self%n_segments = 0 ! for counting # segments
        first_byte = 1
        do isegment=1,MAX_N_SEGEMENTS
            sz = sizeof(self%header(isegment)%fromto(1))
            self%header(isegment)%fromto(1)          = transfer(byte_array(first_byte:first_byte + sz - 1),  self%header(isegment)%fromto(1))
            first_byte = first_byte + sz
            self%header(isegment)%fromto(2)          = transfer(byte_array(first_byte:first_byte + sz - 1),  self%header(isegment)%fromto(2))
            first_byte = first_byte + sz
            self%header(isegment)%n_bytes_per_record = transfer(byte_array(first_byte:first_byte + sz - 1),  self%header(isegment)%n_bytes_per_record)
            first_byte = first_byte + sz
            self%header(isegment)%n_records          = transfer(byte_array(first_byte:first_byte + sz - 1),  self%header(isegment)%n_records)
            first_byte = first_byte + sz
            self%header(isegment)%first_data_byte    = transfer(byte_array(first_byte:first_byte + sz - 1),  self%header(isegment)%first_data_byte)
            first_byte = first_byte + sz
            ! update # segments counter
            if( self%header(isegment)%n_bytes_per_record > 0 .and. self%header(isegment)%n_records > 0&
                &.and. self%header(isegment)%first_data_byte > 0 ) self%n_segments = isegment ! to allow empty in-between segments
        end do

        contains

            subroutine open_local
                integer :: io_stat, tmpunit
                if( .not. self%l_open )then
                    call fopen(tmpunit,fname,access='STREAM', action='READWRITE', status='UNKNOWN', iostat=io_stat)
                    call fileio_errmsg('binoris ; open_local '//trim(fname), io_stat)
                    self%funit  = tmpunit
                    self%l_open = .true.
                endif
            end subroutine open_local

    end subroutine open

    subroutine close( self )
        class(binoris), intent(inout) :: self !< instance
        integer :: io_stat
        if( self%l_open )then
            call fclose(self%funit,io_stat,errmsg='binoris ; close ')
            self%l_open = .false.
        end if
    end subroutine close

    subroutine write_header( self )
        class(binoris), intent(inout) :: self  !< instance
        integer(kind=1) :: byte_array(N_BYTES_HEADER)
        integer :: first_byte, sz
        integer :: io_status, isegment
        if( .not. self%l_open ) stop 'file needs to be open; binoris :: write_header'
        ! transfer header to byte array
        first_byte = 1
        do isegment=1,MAX_N_SEGEMENTS
            sz = sizeof(self%header(isegment)%fromto(1))
            byte_array(first_byte:first_byte + sz - 1) = transfer(self%header(isegment)%fromto(1),          byte_array(first_byte:first_byte + sz - 1))
            first_byte = first_byte + sz
            byte_array(first_byte:first_byte + sz - 1) = transfer(self%header(isegment)%fromto(2),          byte_array(first_byte:first_byte + sz - 1))
            first_byte = first_byte + sz
            byte_array(first_byte:first_byte + sz - 1) = transfer(self%header(isegment)%n_bytes_per_record, byte_array(first_byte:first_byte + sz - 1))
            first_byte = first_byte + sz
            byte_array(first_byte:first_byte + sz - 1) = transfer(self%header(isegment)%n_records,          byte_array(first_byte:first_byte + sz - 1))
            first_byte = first_byte + sz
            byte_array(first_byte:first_byte + sz - 1) = transfer(self%header(isegment)%first_data_byte,    byte_array(first_byte:first_byte + sz - 1))
            first_byte = first_byte + sz
        end do
        ! write header
        write(unit=self%funit,pos=1,iostat=io_status) byte_array
        if( io_status .ne. 0 ) call fileio_errmsg('binoris :: write_header, ERROR writing header bytes ', io_status)
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

    subroutine write_segment( self, isegment, os, fromto )
        class(binoris),    intent(inout) :: self
        integer,           intent(in)    :: isegment
        class(oris),       intent(inout) :: os
        integer, optional, intent(in)    :: fromto(2)
        integer(kind=1),  allocatable :: byte_array(:,:)
        character(len=:), allocatable :: str_os_line, str_dyn
        integer :: i, irec, io_status
        if( .not. self%l_open ) stop 'file needs to be open; binoris :: write_segment'
        ! add segment to stack, this sets all the information needed for allocation
        call self%add_segment(isegment, os, fromto)
        ! update byte ranges in header
        call self%update_byte_ranges
        ! allocate byte array
        allocate(byte_array(self%header(isegment)%n_records,self%header(isegment)%n_bytes_per_record))
        ! allocate string with static lenght (set to max(strlen))
        allocate(character(len=self%header(isegment)%n_bytes_per_record) :: str_os_line)
        ! transfer orientation data to raw byte array via string representation
        irec = 0
        do i=self%header(isegment)%fromto(1),self%header(isegment)%fromto(2)
            irec               = irec + 1
            str_dyn            = os%ori2str(i)
            str_os_line        = str_dyn ! string of lenght that matches record, since different oris will have different strlen
            byte_array(irec,:) = transfer(str_os_line, byte_array(irec,:))
        end do
        ! write raw bytes
        write(unit=self%funit,pos=self%header(isegment)%first_data_byte,iostat=io_status) byte_array
        if( io_status .ne. 0 ) call fileio_errmsg('binoris :: write_segment, ERROR when writing 2D byte array to disk', io_status)
    end subroutine write_segment

    subroutine read_segment( self, isegment, os )
        class(binoris), intent(inout) :: self
        integer,        intent(in)    :: isegment
        class(oris),    intent(inout) :: os
        integer(kind=1) :: byte_array(self%header(isegment)%n_records,self%header(isegment)%n_bytes_per_record)
        character(len=self%header(isegment)%n_bytes_per_record) :: str_os_line ! string with static lenght (set to max(strlen))
        character(len=512) :: io_message
        integer :: i, irec, io_status
        if( .not. self%l_open ) stop 'file needs to be open; binoris :: read_segment'
        if( isegment < 1 .or. isegment > self%n_segments ) stop 'isegment out of bound; binoris :: write_segment'
        if( self%header(isegment)%n_records > 0 .and. self%header(isegment)%n_bytes_per_record > 0 )then
            ! read raw byte array
            read(unit=self%funit,pos=self%header(isegment)%first_data_byte,iostat=io_status,iomsg=io_message) byte_array
            if( io_status .ne. 0 ) call fileio_errmsg('binoris :: open, ERROR when reading 2D byte array from disk '//trim(io_message), io_status)
            ! transfer raw bytes to oris via string representation
            irec = 0
            do i=self%header(isegment)%fromto(1),self%header(isegment)%fromto(2)
                irec        = irec + 1
                str_os_line = transfer(byte_array(irec,:), str_os_line)
                call os%str2ori(irec, str_os_line) ! irec because of sp_project implementation
            end do
        else
            ! empty segment, nothing to do
        endif
    end subroutine read_segment

    subroutine read_segment_ctfparams_state_eo( self, isegment, os )
        use simple_ori, only: ori
        class(binoris), intent(inout) :: self
        integer,        intent(in)    :: isegment
        class(oris),    intent(inout) :: os
        integer, parameter :: NFLAGS = 11
        character(len=32)  :: flags(NFLAGS)
        type(ori)          :: o
        integer(kind=1)    :: byte_array(self%header(isegment)%n_records,self%header(isegment)%n_bytes_per_record)
        character(len=self%header(isegment)%n_bytes_per_record) :: str_os_line ! string with static lenght (set to max(strlen))
        character(len=512) :: io_message
        integer :: i, j, irec, io_status
        if( .not. self%l_open ) stop 'file needs to be open; binoris :: read_segment'
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
        ! read raw byte array
        read(unit=self%funit,pos=self%header(isegment)%first_data_byte,iostat=io_status,iomsg=io_message) byte_array
        if( io_status .ne. 0 ) call fileio_errmsg('binoris :: open, ERROR when reading 2D byte array from disk '//trim(io_message), io_status)
        ! transfer raw bytes to oris via string representation
        irec = 0
        do i=self%header(isegment)%fromto(1),self%header(isegment)%fromto(2)
            irec        = irec + 1
            str_os_line = transfer(byte_array(irec,:), str_os_line)
            call o%str2ori(str_os_line)
            do j=1,NFLAGS
                if( o%isthere(trim(flags(j))) ) call os%set(i, trim(flags(j)), o%get(trim(flags(j))))
            end do
        end do
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

    ! private routines

    subroutine clear_segments( self )
        class(binoris), intent(inout) :: self
        if( self%n_segments <= 0 ) return
        ! clear header
        self%header(:)%fromto(1)          = 0
        self%header(:)%fromto(2)          = 0
        self%header(:)%n_bytes_per_record = 0
        self%header(:)%n_records          = 0
        self%header(:)%first_data_byte    = 0
        ! clear the rest
        self%n_segments = 0
    end subroutine clear_segments

    subroutine add_segment( self, isegment, os, fromto )
        class(binoris),    intent(inout) :: self
        integer,           intent(in)    :: isegment
        class(oris),       intent(inout) :: os
        integer, optional, intent(in)    :: fromto(2)
        integer :: strlen_max, n_oris
        ! sanity check isegment
        if( isegment < 1 .or. isegment > MAX_N_SEGEMENTS ) stop 'ERROR, isegment out of range; binoris :: add_segment'
        if( isegment > self%n_segments ) self%n_segments = isegment
        ! set range in segment
        n_oris = os%get_noris()
        if( present(fromto) )then
            if( fromto(1) < 1 .or. fromto(2) > n_oris ) stop 'fromto out of range; binoris :: add_segment'
            self%header(isegment)%fromto = fromto
        else
            self%header(isegment)%fromto(1) = 1
            self%header(isegment)%fromto(2) = n_oris
        endif
        ! set maximum trimmed string lenght to n_bytes_per_record
        strlen_max = os%max_ori_strlen_trim()
        self%header(isegment)%n_bytes_per_record = strlen_max
        ! set n_records
        self%header(isegment)%n_records = self%header(isegment)%fromto(2) - self%header(isegment)%fromto(1) + 1
        if( self%header(isegment)%n_records < 1 ) stop 'ERROR, input oritab (os) empty; binoris :: add_segment'
    end subroutine add_segment

    subroutine update_byte_ranges( self )
        class(binoris), intent(inout) :: self
        integer(kind=8) :: n_bytes_tot
        integer :: isegment
        n_bytes_tot = N_BYTES_HEADER
        if( self%n_segments <= 0 ) return
        do isegment=1,self%n_segments
            self%header(isegment)%first_data_byte = n_bytes_tot + 1
            n_bytes_tot = n_bytes_tot + self%header(isegment)%n_bytes_per_record * self%header(isegment)%n_records
        end do
    end subroutine update_byte_ranges

end module simple_binoris
