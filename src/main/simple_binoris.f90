! for manging orientation data using binary files
module simple_binoris
use, intrinsic :: iso_c_binding
use simple_ori,  only: ori
use simple_oris, only: oris
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
	! in record
	real,              allocatable :: hash_vals(:)
	real,              allocatable :: peak_array(:,:)
	! byte arrays
	integer(kind=1),   allocatable :: byte_array_header(:)
	real(kind=4),      allocatable :: record(:)
	! for on-line use
	integer                        :: funit  = 0
	logical                        :: exists = .false.
  contains
  	procedure, private :: new_1
  	procedure, private :: new_2
  	generic            :: new => new_1, new_2
  	procedure, private :: header2byte_array
  	procedure, private :: byte_array2header
  	procedure          :: open
  	procedure          :: close
  	procedure          :: write_header
  	procedure          :: write_record
  	
  	procedure          :: kill
end type binoris

enum, bind(c)
	enumerator :: E1=1, E2, E3, X, Y, STATE
	enumerator :: PROJ, CORR, OW
end enum

contains

	subroutine new_1( self, a, fromto, os_peak )
		use simple_jiffys, only: alloc_err
		class(binoris),        intent(inout) :: self
		class(oris),           intent(inout) :: a
		integer,     optional, intent(in)    :: fromto(2)
		class(oris), optional, intent(in)    :: os_peak
		type(ori) :: o
		integer   :: alloc_stat
		! destruct possibly pre-existing
		call self%kill
		! set n_hash_vals
		o = a%get_ori(1)
		self%n_hash_vals = o%hash_size()
		! set hash keys
		self%hash_keys = o%hash_keys()
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
		self%n_reals_per_record = self%n_hash_vals + self%n_peaks * OW
		! allocate
		allocate( self%peak_array(self%n_peaks,OW),&
			     &self%byte_array_header(self%n_bytes_header),&
			     &self%record(self%n_reals_per_record), stat=alloc_stat )
		call alloc_err( 'In: binoris :: new_1', alloc_stat )
		! set
		self%peak_array = 0.0
		call self%header2byte_array
        self%record     = 0.0
        ! flag existence
        self%exists     = .true.
	end subroutine new_1

	subroutine new_2( self, fname )
		use simple_filehandling, only: file_exists
		class(binoris),   intent(inout) :: self  !< instance
		character(len=*), intent(in)    :: fname !< filename
		integer            :: file_size, io_status, alloc_stat
		character(len=512) :: io_message
		integer(kind=1)    :: bytes(8)
		! destruct possibly pre-existing
		call self%kill
		! check file existence
		if( .not. file_exists(fname) ) stop 'input file does not exists; binoris :: new_2'
		! open file
		call self%open(fname)
		! check size
		inquire(unit=self%funit, size=file_size)
		if( file_size == -1 )then 
			stop 'file_size cannot be inquired; binoris :: new_2'
		else if( file_size >= 28 )then
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
		call alloc_err( 'In: binoris :: new_2, 1', alloc_stat )
		read(unit=self%funit,pos=1,iostat=io_status,iomsg=io_message) self%byte_array_header
        if( io_status .ne. 0 )then
            write(*,'(a,i0,2a)') '**error(binoris::new_2): error ', io_status,&
            &' when reading header bytes from disk: ', trim(io_message)
            stop 'I/O error; new_2; simple_binoris'
        endif
        call self%close
        call self%byte_array2header
        ! set derived
		self%n_bytes_hash_keys  = 32 * self%n_hash_vals
		self%n_reals_per_record = self%n_hash_vals + self%n_peaks * OW
		! allocate
		allocate( self%peak_array(self%n_peaks,OW),&
			     &self%record(self%n_reals_per_record), stat=alloc_stat )
		call alloc_err( 'In: binoris :: new_2, 2', alloc_stat )
		! set
		self%peak_array = 0.
        self%record     = 0.0
        ! flag existence
        self%exists     = .true.
    end subroutine new_2

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
        self%hash_keys       = transfer(self%byte_array_header(29:self%n_bytes_header), self%hash_keys)
	end subroutine byte_array2header

    subroutine open( self, fname, rwaction )
		use simple_filehandling, only: get_fileunit, is_open
		class(binoris),             intent(inout) :: self     !< instance
		character(len=*),           intent(in)    :: fname    !< filename
        character(len=*), optional, intent(in)    :: rwaction !< read/write flag
        character(len=9) :: rw_str
        character(len=7) :: stat_str
        if( .not. is_open(self%funit) )then
	        if( present(rwaction) )then
	            rw_str = trim(rwaction)
	        else
	            rw_str = 'READWRITE'
	        endif
	        stat_str = 'UNKNOWN'
	        self%funit = get_fileunit()
	        open(unit=self%funit,access='STREAM',file=fname,action=rw_str,status=stat_str)
	    endif
	end subroutine open

	subroutine close( self )
		use simple_filehandling, only: is_open
		class(binoris), intent(in) :: self !< instance
		if( is_open(self%funit) ) close(self%funit)
	end subroutine close

	subroutine write_header( self, fname )
		class(binoris),   intent(inout) :: self  !< instance
		character(len=*), intent(in)    :: fname !< filename
		integer :: io_status
		! assuming file open
		call self%header2byte_array
		write(unit=self%funit,pos=1,iostat=io_status) self%byte_array_header
        if( io_status .ne. 0 )then
            write(*,'(a,i0,a)') '**error(binoris::write_header): error ', io_status, ' when writing header bytes to disk'
            stop 'I/O error; write_header; simple_binoris'
        endif
	end subroutine write_header

	subroutine write_record( self, a, i, os_peak )
		class(binoris),        intent(inout) :: self
		class(oris),           intent(inout) :: a
		integer,               intent(in)    :: i
		class(oris), optional, intent(in)    :: os_peak
		integer :: first_byte, io_status

		! assuming file open

		! transfer data to self%record


		first_byte = self%first_data_byte + (i - 1) * self%n_reals_per_record * 4
		write(unit=self%funit,pos=first_byte,iostat=io_status) self%record
        if( io_status .ne. 0 )then
            write(*,'(a,i0,a)') '**error(binoris::write_record): error ', io_status, ' when writing record bytes to disk'
            stop 'I/O error; write_record; simple_binoris'
        endif
	end subroutine write_record


	subroutine kill( self )
		class(binoris), intent(inout) :: self
		call self%close
		if( self%exists )then
			deallocate(self%hash_keys, self%hash_vals, self%peak_array,&
				&self%byte_array_header, self%record)
			self%exists = .false.
		endif
	end subroutine kill

end module simple_binoris
