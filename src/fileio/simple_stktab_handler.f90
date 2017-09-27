module simple_stktab_handler
use simple_defs
implicit none

public :: stktab_handler
private

type :: mic_stk
	character(len=:), allocatable :: stkname
	integer :: fromp=0, top=0
end type mic_stk

type :: ptcl_entry
	type(mic_stk), pointer :: msp => null()
end type ptcl_entry

type :: stktab_handler
	private
	integer                       :: nmics   = 0
	integer                       :: nptcls  = 0
	integer                       :: ldim(3) = [0,0,0]
	type(ptcl_entry), allocatable :: pens(:)
	type(mic_stk),    allocatable :: mics(:)
	logical                       :: exists=.false.
  contains
  	! constructor
    procedure :: new
    ! getters
    procedure :: get_nmics
    procedure :: get_nptcls
    procedure :: get_ldim
    ! reader
    procedure :: get_stkname_and_ind
    ! destructor
    procedure :: kill
end type

contains

	subroutine new( self, filetabname )
		use simple_fileio,  only: read_filetable
		use simple_imghead, only: find_ldim_nptcls 
		class(stktab_handler), target, intent(inout) :: self
		character(len=*),              intent(in)    :: filetabname
		character(len=STDLEN), allocatable :: micnames(:) 
		integer :: imic, ldim(3), pind_cnt, istart, istop, nptcls, iptcl
		! take care of micrograph stacks
		call read_filetable(filetabname, micnames)
		self%nmics = size(micnames)
		allocate(self%mics(self%nmics))
		istart = 1
		istop  = 0
		do imic=1,self%nmics
			allocate(self%mics(imic)%stkname, source=trim(micnames(imic)))
			call find_ldim_nptcls(trim(micnames(imic)), ldim, nptcls)
			ldim(3) = 1
			if( imic == 1 )then
				self%ldim = ldim
			else
				if( .not. all(self%ldim == ldim) )then
					print *, 'micrograph stack #:   ', imic
					print *, 'ldim in object:       ', self%ldim
					print *, 'ldim read from stack: ', ldim
					stop 'inconsistent logical dimensions; stktab_handler :: new'
				endif
			endif
			istop  = istop + nptcls
			self%mics(imic)%fromp = istart
			self%mics(imic)%top   = istop
			istart = istart + nptcls
		end do
		self%nptcls = self%mics(self%nmics)%top
		! take care of particle entries
		allocate(self%pens(self%nptcls))
		pind_cnt = 0
		do imic=1,self%nmics
			do iptcl=self%mics(imic)%fromp,self%mics(imic)%top
				pind_cnt = pind_cnt + 1
				self%pens(pind_cnt)%msp => self%mics(imic)
			end do
		end do
		self%exists = .true.
	end subroutine new

	integer function get_nmics( self )
		class(stktab_handler), intent(in) :: self
		get_nmics = self%nmics
	end function get_nmics

	integer function get_nptcls( self )
		class(stktab_handler), intent(in) :: self
		get_nptcls = self%nptcls
	end function get_nptcls

	function get_ldim( self ) result( ldim )
		class(stktab_handler), intent(in) :: self
		integer :: ldim(3)
		ldim = self%ldim
	end function get_ldim

	subroutine get_stkname_and_ind( self, iptcl, stkname, ind )
		class(stktab_handler),         intent(in)  :: self
		integer,                       intent(in)  :: iptcl
		character(len=:), allocatable, intent(out) :: stkname
		integer,                       intent(out) :: ind
		integer :: ldim(3), mic_stk_ind
		real    :: smpd
		! first sanity check
		if( iptcl < 1 .or. iptcl > self%nptcls )then
			print *, 'iptcl:       ', iptcl
			print *, 'self%nptcls: ', self%nptcls
			stop 'iptcl index out of overall range; simple_stktab_handler :: get_stkname_and_ind'
		endif
		! second sanity check
		if( iptcl < self%pens(iptcl)%msp%fromp .or. iptcl > self%pens(iptcl)%msp%top )then
			print *, 'iptcl:             ', iptcl
			print *, 'prange for micstk: ', self%pens(iptcl)%msp%fromp, self%pens(iptcl)%msp%top
			stop 'iptcl index out of micstk range; simple_stktab_handler :: get_stkname_and_ind'
		endif
		! calculate index in mic stk
		ind = iptcl - self%pens(iptcl)%msp%fromp + 1
		! read image
		if( allocated(stkname) ) deallocate(stkname)
		allocate(stkname, source=self%pens(iptcl)%msp%stkname)
	end subroutine get_stkname_and_ind

	subroutine kill( self )
		class(stktab_handler), intent(inout) :: self
		integer :: imic, iptcl
		if( self%exists )then
			do imic=1,self%nmics
				deallocate(self%mics(imic)%stkname)
			end do
			do iptcl=1,self%nptcls
				self%pens(iptcl)%msp => null()
			end do
			self%nmics  = 0
			self%nptcls = 0
			self%ldim   = 0
			deallocate(self%pens, self%mics)
			self%exists = .false.
		endif
	end subroutine kill

end module simple_stktab_handler
