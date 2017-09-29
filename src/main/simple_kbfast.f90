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

integer, parameter :: NGRID_DEFAULT = 50000
! ngrid = 2500   => avgerr = 3e-2 maxerr = 0.1  speedup = 3.0
! ngrid = 5000   => avgerr = 2e-2 maxerr = 5e-2 speedup = 3.0
! ngrid = 10000  => avgerr = 9e-3 maxerr = 3e-2 speedup = 3.0
! ngrid = 50000  => avgerr = 2e-3 maxerr = 5e-3 speedup = 2.6
! ngrid = 100000 => avgerr = 8e-4 maxerr = 3e-3 speedup = 2.5
logical, parameter :: WARN=.true.

contains
    
    subroutine new( self, kbexact, which, ngrid )
        use simple_kbinterpol, only: kbinterpol
        class(kbfast),     intent(inout) :: self
        class(kbinterpol), intent(in)    :: kbexact
        character(len=*),  intent(in)    :: which
        integer, optional, intent(in)    :: ngrid
        real    :: Whalf, stepsz, x
        integer :: cnt
        logical :: l_calc_apod
        call self%kill
        if( present(ngrid) )then
        	self%ngrid = ngrid
        else
        	self%ngrid = NGRID_DEFAULT
        endif
        select case(which)
            case('apod')
                l_calc_apod = .true. 
            case('instr')
                l_calc_apod = .false.
            case DEFAULT
                stop 'unsupported which flag (apod or instr); new; simple_kbfast'
        end select
        allocate(self%grid(self%ngrid), self%vals(self%ngrid), source=0., stat=alloc_stat)
        allocchk('In: new; simple_kbfast')
        Whalf  = kbexact%get_winsz()
        stepsz = (2.0 * Whalf) / real(self%ngrid)
        x      = -Whalf
        cnt    = 0
        do while(x <= Whalf)
            cnt = cnt + 1
            if( cnt > self%ngrid ) exit
            self%grid(cnt) = x
            if( l_calc_apod )then
                self%vals(cnt) = kbexact%apod(x)
            else
                self%vals(cnt) = kbexact%instr(x)
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
        real    :: dist
        ind = locate(self%grid, self%ngrid, x)
        if( ind < self%ngrid )then
            lookup = (self%vals(ind) + self%vals(ind + 1)) / 2.0
        else
            lookup = self%vals(ind)
        endif
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
        use simple_rnd,        only: ran3, ran3arr
        use simple_timer,      only: tic, toc, timer_int_kind
        type(kbinterpol)        :: kbobj
        type(kbfast)            :: kbfastobj
        integer, parameter      :: NTESTS=1000, NTESTS_TIME=1000000
        integer(timer_int_kind) :: tslow, tfast
        real    :: ranarr(NTESTS_TIME)
        integer :: i, j, k
        real    :: x, val, diff, val_true, diffsum, maxdiff, mindiff
        call kbobj%new(1.5, 2.0)
        diffsum = 0
        maxdiff = 1e-20
        mindiff = 1e20
        call kbfastobj%new(kbobj, 'apod')
        do i=1,NTESTS
            x        = ran3() * 3.0 - 1.5
            val      = kbfastobj%lookup(x)
            val_true = kbobj%apod(x)
            diff     = abs(val_true - val)
            diffsum  = diffsum + diff
            if( diff > maxdiff ) maxdiff = diff
            if( diff < mindiff ) mindiff = diff
        end do
        print *,'avg/max,min ', diffsum/real(NTESTS), maxdiff,mindiff
        call ran3arr(ranarr)
        ranarr = ranarr * 3.0 - 1.5
        tslow=tic()
        do i=1,NTESTS_TIME
            do j=1,NTESTS_TIME
                do k=1,NTESTS_TIME
                    val_true = kbobj%apod(x)
                end do
            end do
        end do
        print *, 'toc(tslow): ', toc(tslow)
        tfast=tic()
        do i=1,NTESTS_TIME
            do j=1,NTESTS_TIME
                do k=1,NTESTS_TIME
                    val = kbfastobj%lookup(ranarr(k))
                end do
            end do
        end do
        print *, 'toc(tfast): ', toc(tfast)
        call kbfastobj%kill
    end subroutine test_kbfast

end module simple_kbfast
