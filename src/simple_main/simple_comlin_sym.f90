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
use simple_jiffys, only: progress, alloc_err
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
integer                :: fromv=0, tov=0  !< vertices used for search
logical                :: debug=.false.

real, parameter    :: ico_alpha = 37.38  !< face-vertex angle, leads to ~4% oversampling, from Baldwin & Penczek JSB2007
real, parameter    :: delta1 = ico_alpha
real, parameter    :: delta2 = ico_alpha
real, parameter    :: delta3 = 36.       !< 2pi/5/2; icosahedral inplane symmetry

contains

    !>  \brief  is a constructor
    subroutine comlin_sym_init( b, p )
        use simple_map_reduce, only: split_nobjs_even
        use simple_build,      only: build
        class(build),  target, intent(in) :: b
        class(params), target, intent(in) :: p
        integer, allocatable :: parts(:,:)
        integer              :: i, n
        call ico%new('ico')
        n = ico%get_nsym()
        call resoris%new(n)
        hp = p%hp
        lp = p%lp
        bp => b
        pp => p
        ! automated partioning
        if( p%l_distr_exec .and. p%nparts > 1 )then
            if( p%nparts > n )then
                fromv = 0
                tov   = 0
            else
                parts = split_nobjs_even(n, p%nparts)
                fromv = parts(p%part,1)
                tov   = parts(p%part,2)
                deallocate(parts)
            endif
        else
            fromv = 1
            tov   = n
        endif
    end subroutine comlin_sym_init

    !>  \brief  is for finding the symmetry axis given an aligned set of images
    subroutine comlin_sym_axis( p, orientation_best, mode, doprint )
        use simple_math,    only: hpsort
        use simple_strings, only: int2str_pad
        class(params),    intent(in)    :: p
        class(ori),       intent(inout) :: orientation_best
        character(len=*), intent(in)    :: mode
        logical,          intent(in)    :: doprint
        type(ori)            :: o, o_zero
        real,    allocatable :: corrs(:)
        integer, allocatable :: corr_inds(:)
        integer              :: i, ncorr, alloc_stat
        real                 :: corr
        real                 :: opto(6)
        if( doprint ) write(*,'(A)') '>>> SYMMETRY AXIS SEARCH'
        ncorr = ico%get_nsym() 
        allocate(corrs(ncorr), corr_inds(ncorr), stat=alloc_stat)
        call alloc_err("In: comlin_sym_axis; simple_comlin_sym", alloc_stat )
        call o%new
        call o_zero%new
        ! Simplex search
        do i = fromv, tov   ! vertices range from partitioning
            ! init with vertex
            o = ico%apply(o_zero, i)
            if( doprint ) call progress(i, ncorr)
            call update_lims(o, [delta1, delta2, delta3])
            ! calculate
            if( mode .eq. 'sym' )then
                call comlin_symsrch_init(bp, pp, 'simplex', optlims, o, ico_alpha, 'sym')
                opto(:4) = comlin_symsrch_minimize()
            else if( mode .eq. 'pair' )then
                call comlin_symsrch_init(bp, pp, 'simplex', optlims, o, ico_alpha, 'pair')
                opto = comlin_pairsrch_minimize()
                call o%set_shift(opto(5:6))
            else
                stop 'unsupported mode in simple_comlin_sym::comlin_sym_axis'
            endif
            ! stash
            call o%set_euler(opto(2:4))
            corr     = opto(1)
            corrs(i) = corr
            call o%set('corr', corr)
            call resoris%set_ori(i, o)
        end do
        ! output & sorting
        if(p%l_distr_exec .and. mode.eq.'sym')then
            call resoris%write(p%outfile, [fromv,tov])
        else
            corr_inds = (/ (i,i=1,ncorr) /)
            call hpsort(ncorr, corrs, corr_inds)
            orientation_best = resoris%get_ori(corr_inds(ncorr))
            if( doprint )then
               write(*,'(A)') '>>> FOUND SYMMETRY AXIS ORIENTATION:'
               call orientation_best%print
            endif
            deallocate(corrs, corr_inds)
        endif
   end subroutine comlin_sym_axis

   !< Update limits when visiting new vertex
   subroutine update_lims( o, lims)
        type(ori), intent(in) :: o
        real,      intent(in) :: lims(3)          
        real :: euls(3)          
        euls = o%get_euler()
        optlims(1:3,1) = euls-lims
        optlims(1:3,2) = euls+lims
        optlims(4:5,1) = -pp%trs
        optlims(4:5,2) = pp%trs
   end subroutine update_lims
    
end module simple_comlin_sym