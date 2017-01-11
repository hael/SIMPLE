!==Class simple_comlin_symsrch
!
! The simple_comlin_sym singleton provides the common line hybrid discrete/continuous symmetry
! search. The code is distributed with the hope 
! that it will be useful, but _WITHOUT_ _ANY_ _WARRANTY_. Redistribution or modification 
! is regulated by the GNU General Public License. *Author:* Hans Elmlund, 2009-06-11.
module simple_comlin_sym
use simple_defs
use simple_build,  only: build
use simple_params, only: params
use simple_jiffys, only: progress
use simple_ori,    only: ori
use simple_oris,   only: oris
use simple_sym,    only: sym
use simple_comlin_symsrch
implicit none

public :: comlin_sym_init, comlin_sym_axis
private

class(build), pointer  :: bp=>null()      !< pointer to builder
class(params), pointer :: pp=>null()      !< pointer to params
type(sym)              :: ico             !< symmetry operator for sampling
type(oris)             :: resoris         !< for results storage and ranking
real                   :: optlims(5,2)=0. !< Euler limits; only used to setup optimizer
real                   :: hp=0.           !< fixed high-pass limit
real                   :: lp=0.           !< fixed low-pass limit
logical                :: debug=.false.

real, parameter    :: ico_alpha = 37.38  !< face-vertex angle, leads to ~4% oversampling, from Baldwin & Penczek JSB2007
real, parameter    :: delta1 = ico_alpha
real, parameter    :: delta2 = ico_alpha
real, parameter    :: delta3 = 36.       !< 2pi/5/2; icosahedral inplane symmetry

contains

    !>  \brief  is a constructor
    subroutine comlin_sym_init( b, p )
        use simple_build,  only: build
        use simple_params, only: params
        class(build), intent(in), target  :: b
        class(params), intent(in), target :: p
        integer                           :: n
        call ico%new('ico')
        n = ico%get_nsym()
        call resoris%new(n)
        hp = p%hp
        lp = p%lp
        bp => b
        pp => p
    end subroutine

    !>  \brief  is for finding the symmetry axis given an aligned set of images
    subroutine comlin_sym_axis( orientation_best, mode, doprint )
        use simple_math,       only: hpsort
        use simple_jiffys,     only: alloc_err
        class(ori), intent(inout)    :: orientation_best
        character(len=*), intent(in) :: mode
        logical, intent(in)          :: doprint
        type(ori)            :: o, o_zero
        real, allocatable    :: corrs(:)
        integer, allocatable :: corrs_inds(:)
        integer              :: i, ntot, ncorr, alloc_stat
        real                 :: corr, corr_best!, prev_corr
        real                 :: opto(6)
        if( doprint ) write(*,'(A)') '>>> SYMMETRY AXIS SEARCH'
        corr_best = -1.
        ntot      = ico%get_nsym() 
        ncorr     = ntot
        allocate( corrs(ncorr), corrs_inds(ncorr), stat=alloc_stat)
        call alloc_err("In: comlin_sym_axis; simple_comlin_sym", alloc_stat )
        call o%new
        call o_zero%new
        ! First pass Simplex search
        do i=1,ncorr
            o = ico%apply(o_zero, i)         ! o is the vertex at play
            if( doprint ) call progress(i, ncorr)
            call update_lims(o, [delta1, delta2, delta3])
            if( mode .eq. 'sym' )then
                call comlin_symsrch_init(bp, pp, 'simplex', optlims, o, ico_alpha, 'sym')
                opto(:4) = comlin_symsrch_minimize()
            else if( mode .eq. 'pair' )then
                call comlin_symsrch_init(bp, pp, 'simplex', optlims, o, ico_alpha, 'pair')
                opto = comlin_pairsrch_minimize()
                call o%set('x', opto(5))
                call o%set('y', opto(6))
            else
                stop 'unsupported mode in simple_comlin_sym::comlin_sym_axis'
            endif
            call o%set_euler(opto(2:4))
            corr = opto(1)
            corrs(i) = corr
            call o%set('corr', corr)
            call resoris%set_ori(i, o)
        end do
        do i=1, ncorr
            corrs_inds(i) = i
        end do
        call hpsort(ncorr, corrs, corrs_inds)
        orientation_best = resoris%get_ori(corrs_inds(ncorr))
        if( doprint )then
           write(*,'(A)') '>>> FOUND SYMMETRY AXIS ORIENTATION:'
           call orientation_best%print
        endif
        deallocate(corrs, corrs_inds)
   end subroutine

   !< Update limits when visiting new vertex
   subroutine update_lims( o, lims)
        type(ori), intent(in)    :: o
        real, intent(in)         :: lims(3)          
        real                     :: euls(3)          
        euls = o%get_euler()
        optlims(1:3,1) = euls-lims
        optlims(1:3,2) = euls+lims
        optlims(4:5,1) = -pp%trs
        optlims(4:5,2) = pp%trs
   end subroutine
    
end module simple_comlin_sym