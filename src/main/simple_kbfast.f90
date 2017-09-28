module simple_kbfast
#include "simple_lib.f08"
implicit none

public :: kbfast, test_kbfast
private

type :: kbfast
	private
	real, allocatable :: grid(:)
	real, allocatable :: vals(:)
	integer           :: ngrid = 0 
	logical           :: exists=.false.
  contains
    procedure :: new
    procedure :: lookup
    procedure :: kill
end type kbfast

contains
	
	subroutine new( self, kinterpol_obj, ngrid, which )
		use simple_kbinterpol, only: kbinterpol
		class(kbfast),     intent(inout) :: self
		class(kbinterpol), intent(in)    :: kinterpol_obj
		integer,           intent(in)    :: ngrid
		character(len=*),  intent(in)    :: which
		real    :: Whalf, stepsz, x
		integer :: cnt
		logical :: l_calc_apod
		self%ngrid = ngrid
		allocate(self%grid(ngrid), self%vals(ngrid), stat=alloc_stat)
		allocchk('In: new; simple_kbfast')
		select case(which)
			case('apod')
				l_calc_apod = .true. 
			case('instr')
				l_calc_apod = .false.
			case DEFAULT
				stop 'unsupported which flag (apod or instr); new; simple_kbfast'
		end select
		Whalf  = kinterpol_obj%get_winsz()
		stepsz = (2.0 * Whalf) / real(ngrid)
		x   = -Whalf
		cnt = 0
		do while(x <= Whalf)
			cnt = cnt + 1
			self%grid(cnt) = x
			if( l_calc_apod )then
				self%vals(cnt) = kinterpol_obj%apod(x)
			else
				self%vals(cnt) = kinterpol_obj%instr(x)
			endif
			x = x + stepsz
		end do
		self%exists = .true.
	end subroutine new

	real function lookup( self, x )
		use simple_math, only: locate
		class(kbfast), intent(in) :: self
		real,          intent(in) :: x
		integer :: ind
		real    :: x0, x1, y0, y1
    	ind = locate(self%grid, self%ngrid, x)
    	! ind means x is between self%vals(ind) and self%vals(ind + 1)
    	! ind = 0 or ind == self%ngrid means out of bound
    	if( ind == 0 .or. ind == self%ngrid )then
    		print *, 'WARNING! out of bounds; simple_kbfast :: lookup'
    		if( ind == 0 ) ind = 1
    		lookup = self%vals(ind)
    		return
    	endif
    	x0 = self%grid(ind)
    	x1 = self%grid(ind + 1)
    	y0 = self%vals(ind)
    	y1 = self%vals(ind + 1)
    	lookup = (y0 * (x1 - x) + y1 * (x - x0)) / (x1 - x0)
	end function lookup

	subroutine kill( self )
		class(kbfast), intent(inout) :: self
		if( self%exists )then
			deallocate(self%grid, self%vals)
			self%exists = .false.
		endif
	end subroutine kill

	subroutine test_kbfast
		use simple_kbinterpol, only: kbinterpol
		use simple_rnd,        only: ran3
		type(kbinterpol) :: kbobj
		type(kbfast)     :: kbfastobj
		integer, parameter :: NGRID_OPT = 100000, NTESTS=1000
		integer :: i
		real    :: x, val, diff, val_true, diffsum, maxdiff, mindiff
		call kbobj%new(1.5, 2.0)
		diffsum=0
		maxdiff=1e-20
		mindiff=1e20
		call kbfastobj%new(kbobj, NGRID_OPT, 'apod')
		do i=1,NTESTS
			x        = ran3() * 3.0 -1.5
			val      = kbfastobj%lookup(x)
			val_true = kbobj%apod(x)
			diff     = abs(val_true - val)
			diffsum  = diffsum + diff
			if (diff > maxdiff) maxdiff = diff
			if (diff < mindiff) mindiff = diff
		end do
		print *,'avg/max,min ', diffsum/real(NTESTS), maxdiff,mindiff
		call kbfastobj%kill
	end subroutine test_kbfast

end module simple_kbfast
