! for manging orientation data using binary files
module simple_binoris
use, intrinsic :: iso_c_binding
use simple_defs  ! use all in there
use simple_ori,  only: ori
use simple_oris, only: oris
use simple_fileio
implicit none

public :: binoris
private

type binoris
    private 
    ! in header
    integer                        :: n_bytes_header  = 0
    integer                        :: n_hash_vals     = 0 
    integer                        :: first_data_byte = 0
    integer                        :: n_records       = 0
    integer                        :: n_peaks         = 0
    integer                        :: fromto(2)
    character(len=32), allocatable :: hash_keys(:)
    ! derived
    integer                        :: n_bytes_hash_keys  = 0
    integer                        :: n_reals_per_record = 0
    ! header byte array
    integer(kind=1),   allocatable :: byte_array_header(:)
    ! record
    real(kind=4),      allocatable :: record(:)
    ! for on-line use
    integer                        :: funit  = 0
    logical                        :: l_open = .false.
    logical                        :: exists = .false.
  contains
    ! constructors
    procedure          :: new
    ! checkers/setters/getters
    procedure, private :: same_dims
    generic            :: operator(.eqdims.) => same_dims
    procedure          :: set_fromto
    procedure          :: get_n_records
    procedure          :: get_fromto
    procedure          :: get_n_peaks
    ! I/O
    procedure          :: open
    procedure, private :: open_local
    procedure          :: close
    procedure          :: print_header
    procedure          :: print_hash_keys
    procedure          :: write_header
    procedure, private :: write_record_1
    procedure, private :: write_record_2
    generic            :: write_record => write_record_1, write_record_2
    procedure, private :: read_record_1
    procedure, private :: read_record_2
    generic            :: read_record => read_record_1, read_record_2
    procedure          :: read_ctfparams_and_state
    ! byte indexing
    procedure          :: first_byte
    procedure          :: last_byte
    ! header byte array conversions for mixed format read/write
    procedure, private :: header2byte_array
    procedure, private :: byte_array2header
    ! destructor
    procedure          :: kill
end type binoris

integer, parameter :: NOPEAKFLAGS = 9
character(len=32)  :: o_peak_flags(NOPEAKFLAGS)

contains

    ! constructors

    subroutine new( self, a, fromto, os_peak )
        use simple_syslib, only: alloc_errchk
        class(binoris),        intent(inout) :: self
        class(oris),           intent(inout) :: a
        integer,     optional, intent(in)    :: fromto(2)
        class(oris), optional, intent(in)    :: os_peak
        type(ori) :: o
        ! destruct possibly pre-existing
        call self%kill
        ! set n_hash_vals
        o = a%get_ori(1)
        self%n_hash_vals = o%hash_size()
        ! set hash keys
        self%hash_keys = o%hash_keys()               ! Intel warn: code for reallocating lhs
        if( size(self%hash_keys) /= self%n_hash_vals )&
        &stop 'ERROR, n_hash_vals /= n_keys; binoris :: new_1'
        ! set range
        if( present(fromto) )then
            self%fromto = fromto
        else
            self%fromto(1) = 1
            self%fromto(2) = a%get_noris()
        endif
        ! set n_records
        self%n_records = self%fromto(2) - self%fromto(1) + 1
        if( self%n_records < 1 ) stop 'ERROR, input oritab (a) empty; binoris :: new_1'
        ! set n_peaks
        self%n_peaks = 0
        if( present(os_peak) )then
            self%n_peaks = os_peak%get_noris()
        endif
        ! set derived
        self%n_bytes_hash_keys  = 32 * self%n_hash_vals
        self%n_bytes_header     = 7  * 4 + self%n_bytes_hash_keys
        self%first_data_byte    = self%n_bytes_header + 1
        self%n_reals_per_record = self%n_hash_vals + self%n_peaks * NOPEAKFLAGS
        ! allocate
        allocate( self%byte_array_header(self%n_bytes_header),&
                 &self%record(self%n_reals_per_record), stat=alloc_stat )
        if(alloc_stat/=0)call alloc_errchk( 'In: binoris :: new_1', alloc_stat )
        ! set
        call self%header2byte_array
        self%record = 0.0
        call set_o_peak_flags ! class variable
        ! flag existence
        self%exists = .true.
    end subroutine new

    ! checkers/setters/getters

    logical function same_dims( self1, self2 )
        class(binoris), intent(in) :: self1, self2 !< instances
        same_dims = all([self1%n_bytes_header == self2%n_bytes_header,&
                        &self1%n_hash_vals    == self2%n_hash_vals,&
                        &self1%n_peaks        == self2%n_peaks])
    end function same_dims

    subroutine set_fromto( self, fromto )
        class(binoris), intent(inout) :: self      !< instance
        integer,        intent(in)    :: fromto(2) !< range
        if( fromto(1) > 1 .or. fromto(2) < fromto(1) )&
        &stop 'unallowed fromto range; binoris :: set_fromto'
        self%fromto = fromto 
        ! set n_records
        self%n_records = self%fromto(2) - self%fromto(1) + 1
    end subroutine set_fromto

    integer function get_n_records( self )
        class(binoris), intent(in) :: self !< instance
        get_n_records = self%n_records
    end function get_n_records

    function get_fromto( self ) result( fromto )
        class(binoris), intent(in) :: self !< instance
        integer :: fromto(2)
        fromto = self%fromto
    end function get_fromto

    integer function get_n_peaks( self )
        class(binoris), intent(in) :: self !< instance
        get_n_peaks = self%n_peaks
    end function get_n_peaks

    ! I/O supporting ori + oris (peaks)

    subroutine open( self, fname, del_if_exists )
        use simple_syslib, only: alloc_errchk
        class(binoris),    intent(inout) :: self          !< instance
        character(len=*),  intent(in)    :: fname         !< filename
        logical, optional, intent(in)    :: del_if_exists !< If the file already exists on disk, replace 
        integer(dp)            :: filesz
        integer            :: io_status
        character(len=512) :: io_message
        integer(kind=1)    :: bytes(8)
        ! deletion logics
        if( present(del_if_exists) )then
            if( del_if_exists )then
                call del_file(trim(fname))
            endif
        endif
        ! existence logics
        if( .not. file_exists(trim(fname)) )then
            if( self%exists )then
                call self%open_local(trim(fname))
                return
            else
                stop 'cannot open non-existing file using non-existing object; binoris :: open'
            endif
        endif
        call self%kill        
        ! open file
        call self%open_local(trim(fname))
        ! check size
        filesz = funit_size(self%funit)
        if( filesz == -1 )then 
            stop 'file_size cannot be inquired; binoris :: new_2'
        else if( filesz >= 28 )then
            ! ok
        else
            stop 'file_size too small to contain a header; binoris :: new_2'
        endif
        ! read first two header records (n_bytes_header & n_hash_vals)
        read(unit=self%funit,pos=1,iostat=io_status,iomsg=io_message) bytes
        if( io_status .ne. 0 )then
            write(*,'(a,i0,2a)') '**error(binoris::new_2): error ', io_status,&
            &' when reading first two header record from disk: ', trim(io_message)
            stop 'I/O error; new_2; simple_binoris'
        endif
        ! allocate header byte array and hash_keys
        self%n_bytes_header = transfer(bytes(1:4), self%n_bytes_header)
        self%n_hash_vals    = transfer(bytes(5:8), self%n_hash_vals)
        allocate( self%byte_array_header(self%n_bytes_header),&
            &self%hash_keys(self%n_hash_vals), stat=alloc_stat )
        if(alloc_stat/=0)call alloc_errchk( 'In: binoris :: new_2, 1', alloc_stat )
        read(unit=self%funit,pos=1,iostat=io_status,iomsg=io_message) self%byte_array_header
        if( io_status .ne. 0 )then
            write(*,'(a,i0,2a)') '**error(binoris::new_2): error ', io_status,&
            &' when reading header bytes from disk: ', trim(io_message)
            stop 'I/O error; new_2; simple_binoris'
        endif
        call self%byte_array2header
        ! set derived
        self%n_bytes_hash_keys  = 32 * self%n_hash_vals
        self%n_reals_per_record = self%n_hash_vals + self%n_peaks * NOPEAKFLAGS
        ! allocate
        allocate( self%record(self%n_reals_per_record), stat=alloc_stat )
        if(alloc_stat/=0)call alloc_errchk( 'In: binoris :: new_2, 2', alloc_stat )
        ! set
        self%record = 0.0
        call set_o_peak_flags ! class variable
        ! flag existence
        self%exists = .true.
    end subroutine open

    subroutine open_local( self, fname, rwaction )
        class(binoris),             intent(inout) :: self     !< instance
        character(len=*),           intent(in)    :: fname    !< filename
        character(len=*), optional, intent(in)    :: rwaction !< read/write flag
        character(len=9) :: rw_str
        character(len=7) :: stat_str
        integer          :: io_stat, tmpunit
        if( .not. self%l_open )then
            if( present(rwaction) )then
                rw_str = trim(rwaction)
            else
                rw_str = 'READWRITE'
            endif
            stat_str = 'UNKNOWN'
            call fopen(tmpunit,fname,access='STREAM',action=rw_str,status=stat_str, iostat=io_stat)
            call fileio_errmsg('binoris ; open_local '// trim(fname), io_stat)
            self%funit  = tmpunit
            self%l_open = .true.
        endif
    end subroutine open_local

    subroutine close( self )
        class(binoris), intent(inout) :: self !< instance
        integer          :: io_stat
        if( self%l_open )then
            call fclose(self%funit,io_stat,errmsg='binoris ; close ')
            self%l_open =.false.
        end if
    end subroutine close

    subroutine print_header( self )
        class(binoris), intent(in) :: self
        write(*,*) '*****  HEADER VALUES  *****'
        write(*,*) 'n_bytes_header    : ', self%n_bytes_header
        write(*,*) 'n_hash_vals       : ', self%n_hash_vals
        write(*,*) 'first_data_byte   : ', self%first_data_byte
        write(*,*) 'n_records         : ', self%n_records
        write(*,*) 'n_peaks           : ', self%n_peaks
        write(*,*) 'fromto(1)         : ', self%fromto(1)
        write(*,*) 'fromto(2)         : ', self%fromto(2)
        write(*,*) 'n_bytes_hash_keys : ', self%n_bytes_hash_keys
        write(*,*) 'n_reals_per_record: ', self%n_reals_per_record
    end subroutine print_header

    subroutine print_hash_keys( self )
        class(binoris), intent(in) :: self
        integer :: i
        write(*,*) '*****  HASH KEYS  *****'
        do i=1,size(self%hash_keys)
            write(*,*) 'key ', i, ' is ', trim(self%hash_keys(i))
        end do
    end subroutine print_hash_keys

    subroutine write_header( self )
        class(binoris), intent(inout) :: self  !< instance
        integer :: io_status
        if( .not. self%l_open ) stop 'file needs to be open; binoris :: write_header'
        call self%header2byte_array
        write(unit=self%funit,pos=1,iostat=io_status) self%byte_array_header
        if( io_status .ne. 0 )&
            call fileio_errmsg('simple_binoris::write_header: error when writing header bytes to disk', io_status) 
    end subroutine write_header

    subroutine write_record_1( self, i, self2 )
        class(binoris), intent(inout) :: self
        integer,        intent(in)    :: i
        class(binoris), intent(in)    :: self2
        integer :: io_status
        if( .not. self%l_open         ) stop 'file needs to be open; binoris :: write_record_1'
        if( .not. (self.eqdims.self2) ) stop 'filehandlers have different dims; binoris :: write_record_1'
        if( i < self%fromto(1) .or. i > self%fromto(2) ) stop 'index i out of bound; binoris :: write_record_1'
        write(unit=self%funit,pos=self%first_byte(i),iostat=io_status) self2%record
        if( io_status .ne. 0 )&
             call fileio_errmsg('simple_binoris::write_record_1: error when writing record bytes to disk', io_status) 
    end subroutine write_record_1

    subroutine write_record_2( self, i, a, os_peak )
        class(binoris),        intent(inout) :: self
        integer,               intent(in)    :: i
        class(oris),           intent(inout) :: a
        class(oris), optional, intent(inout) :: os_peak
        integer                   :: io_status, iflag, ipeak, cnt, sz_vals
        type(ori)                 :: o
        real(kind=4), allocatable :: vals(:)
        if( .not. self%l_open ) stop 'file needs to be open; binoris :: write_record_2'
        if( i < self%fromto(1) .or. i > self%fromto(2) ) stop 'index i out of bound; binoris :: write_record_2'
        ! transfer hash data to self%record
        o       = a%get_ori(i)
        vals    = o%hash_vals()                  ! Intel warn: code for reallocating lhs
        sz_vals = size(vals)
        if( self%n_hash_vals /= sz_vals )then
            print *, 'self%n_hash_vals: ', self%n_hash_vals
            print *, 'sz_vals         : ', sz_vals
            stop 'nonconforming hash size; binoris :: write_record_2'
        endif
        self%record = 0.0
        self%record(:self%n_hash_vals) = vals
        if( present(os_peak) )then
            ! transfer os_peak data to self%record
            if( self%n_peaks /= os_peak%get_noris() ) stop 'nonconforming os_peak size; binoris :: write_record_2'        
            cnt = self%n_hash_vals                
            do ipeak=1,self%n_peaks
                do iflag=1,NOPEAKFLAGS
                    cnt = cnt + 1
                    if( .not. os_peak%isthere(trim(o_peak_flags(iflag))) )then
                        write(*,'(a)') 'WARNING! The '//trim(o_peak_flags(iflag))//' is missing from os_peak'
                        write(*,'(a)') 'In: simple_binors; write_record_2'
                    endif
                    self%record(cnt) = os_peak%get(ipeak, trim(o_peak_flags(iflag)))
                end do
            end do
        endif
        write(unit=self%funit,pos=self%first_byte(i),iostat=io_status) self%record
        if( io_status .ne. 0 )&
            call fileio_errmsg('binoris::write_record_2: error when writing record bytes to disk', io_status) 
    end subroutine write_record_2

    subroutine read_record_1( self, i )
        class(binoris), intent(inout) :: self
        integer,        intent(in)    :: i
        integer :: io_status
        if( .not. self%l_open ) stop 'file needs to be open; binoris :: read_record_1'
        if( i < self%fromto(1) .or. i > self%fromto(2) ) stop 'index i out of bound; binoris :: read_record_1'
        read(unit=self%funit,pos=self%first_byte(i),iostat=io_status) self%record
        if( io_status .ne. 0 )&
            call fileio_errmsg("binoris::read_record_1  error when reading record bytes from disk", io_status)
    end subroutine read_record_1

    subroutine read_record_2( self, i, a, os_peak, nst )
        class(binoris),        intent(inout) :: self
        integer,               intent(in)    :: i
        class(oris),           intent(inout) :: a
        class(oris), optional, intent(inout) :: os_peak
        integer,     optional, intent(out)   :: nst
        integer   :: iflag, ipeak, cnt, j, state
        real      :: euls(3)
        call self%read_record_1(i)
        ! transfer hash data
        euls = 0.
        do j=1,self%n_hash_vals
            select case(trim(self%hash_keys(j)))
                case('e1')
                    euls(1) = self%record(j)
                case('e2')
                    euls(2) = self%record(j)
                case('e3')
                    euls(3) = self%record(j)
                case DEFAULT
                    call a%set(i, trim(self%hash_keys(j)), self%record(j))
            end select            
        end do
        call a%set_euler(i, euls)
        if( present(nst) )then
            state = int(a%get(i, 'state'))
            nst = max(1,max(state,nst))
        endif
        if( present(os_peak) )then
            ! transfer os_peak data
            cnt = self%n_hash_vals
            call os_peak%new(self%n_peaks)
            do ipeak=1,self%n_peaks
                euls = 0.
                do iflag=1,NOPEAKFLAGS
                    cnt = cnt + 1
                    select case(trim(o_peak_flags(iflag)))
                        case('e1')
                            euls(1) = self%record(cnt)
                        case('e2')
                            euls(2) = self%record(cnt)
                        case('e3')
                            euls(3) = self%record(cnt)
                        case DEFAULT
                            call os_peak%set(ipeak, trim(o_peak_flags(iflag)), self%record(cnt))
                    end select
                end do
                call os_peak%set_euler(ipeak, euls)
            end do
        endif
    end subroutine read_record_2

    subroutine read_ctfparams_and_state( self, i, a )
        class(binoris), intent(inout) :: self
        integer,        intent(in)    :: i
        class(oris),    intent(inout) :: a
        integer :: j
        call self%read_record_1(i)
        ! transfer hash data
        do j=1,self%n_hash_vals
            select case(trim(self%hash_keys(j)))
                case('smpd','kv','cs','fraca','phaseplate','dfx','dfy','angast','bfac','state')
                    call a%set(i, trim(self%hash_keys(j)), self%record(j))
            end select            
        end do
    end subroutine read_ctfparams_and_state

    ! byte indexing

    integer function first_byte( self, irec )
        class(binoris), intent(in) :: self
        integer,        intent(in) :: irec
        integer :: ind
        ind = irec - self%fromto(1) + 1
        first_byte = self%first_data_byte + (ind - 1) * self%n_reals_per_record * 4
    end function first_byte

    integer function last_byte( self )
        class(binoris), intent(in) :: self
        integer :: ind
        ind = self%fromto(2) - self%fromto(1) + 1
        last_byte = self%n_bytes_header + ind * self%n_reals_per_record * 4
    end function last_byte

    ! header byte array conversions for mixed format read/write

    subroutine header2byte_array( self )
        class(binoris), intent(inout) :: self
        self%byte_array_header(1:4)   = transfer(self%n_bytes_header,  self%byte_array_header(1:4))
        self%byte_array_header(5:8)   = transfer(self%n_hash_vals,     self%byte_array_header(5:8))
        self%byte_array_header(9:12)  = transfer(self%n_records,       self%byte_array_header(9:12))
        self%byte_array_header(13:16) = transfer(self%first_data_byte, self%byte_array_header(13:16))
        self%byte_array_header(17:20) = transfer(self%n_peaks,         self%byte_array_header(17:20))
        self%byte_array_header(21:24) = transfer(self%fromto(1),       self%byte_array_header(21:24))
        self%byte_array_header(25:28) = transfer(self%fromto(2),       self%byte_array_header(25:28))
        self%byte_array_header(29:self%n_bytes_header)&
        &= transfer(self%hash_keys, self%byte_array_header(29:self%n_bytes_header))
    end subroutine header2byte_array

    subroutine byte_array2header( self )
        class(binoris), intent(inout) :: self
        self%n_bytes_header  = transfer(self%byte_array_header(1:4),   self%n_bytes_header)
        self%n_hash_vals     = transfer(self%byte_array_header(5:8),   self%n_hash_vals)
        self%n_records       = transfer(self%byte_array_header(9:12),  self%n_records)
        self%first_data_byte = transfer(self%byte_array_header(13:16), self%first_data_byte)
        self%n_peaks         = transfer(self%byte_array_header(17:20), self%n_peaks)
        self%fromto(1)       = transfer(self%byte_array_header(21:24), self%fromto(1))
        self%fromto(2)       = transfer(self%byte_array_header(25:28), self%fromto(2))
        ! Intel warn: code for reallocating lhs
        self%hash_keys       = transfer(self%byte_array_header(29:self%n_bytes_header), self%hash_keys)
    end subroutine byte_array2header

    ! destructor

    subroutine kill( self )
        class(binoris), intent(inout) :: self
        call self%close
        if( self%exists )then
            deallocate(self%hash_keys, self%byte_array_header, self%record)
            self%exists = .false.
        endif
    end subroutine kill

    ! local subs

    subroutine set_o_peak_flags
        o_peak_flags(1) = 'e1'
        o_peak_flags(2) = 'e2'
        o_peak_flags(3) = 'e3'
        o_peak_flags(4) = 'x'
        o_peak_flags(5) = 'y'
        o_peak_flags(6) = 'state'
        o_peak_flags(7) = 'proj'
        o_peak_flags(8) = 'corr'
        o_peak_flags(9) = 'ow'
    end subroutine set_o_peak_flags

end module simple_binoris
