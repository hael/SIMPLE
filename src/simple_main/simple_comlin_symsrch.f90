!==Class simple_comlin_symsrch
!
! The simple_comlin_symsrch singleton provides the  common line continuous symmetry
! search. The code is distributed with the hope 
! that it will be useful, but _WITHOUT_ _ANY_ _WARRANTY_. Redistribution or modification 
! is regulated by the GNU General Public License. *Author:* Hans Elmlund, 2009-06-11.
module simple_comlin_symsrch
use simple_defs
use simple_build,       only: build
use simple_opt_factory, only: opt_factory
use simple_opt_spec,    only: opt_spec
use simple_ori,         only: ori
use simple_oris,        only: oris
use simple_optimizer,   only: optimizer
implicit none

public :: comlin_symsrch_init, comlin_symsrch_minimize, comlin_symsrch_get_nevals, comlin_pairsrch_minimize
private

type(opt_factory)         :: ofac            !< optimizer factory
type(opt_spec)            :: ospec           !< optimizer specification object
class(build), pointer     :: bp=>null()      !< pointer to builder
class(optimizer), pointer :: nlopt=>null()   !< pointer to nonlinear optimizer
type(oris)                :: a_copy          !< safe keeping for cost function
type(ori)                 :: optori          !< axis ori being optimised
type(ori)                 :: vertex          !< current vertex (initial position)
integer, pointer          :: iptcl=>null()   !< ptcl index
integer, pointer          :: jptcl=>null()   !< ptcl index
real, pointer             :: dynlim=>null()  !< dynamic lowpass limit
integer                   :: nptcls=0        !< nr of particles
real                      :: eullims(3,2)=0. !< Euler limits; only used to setup optimizer
real                      :: hp=0.           !< fixed high-pass limit
real                      :: lp=0.           !< fixed low-pass limit
real                      :: euldistmax=1.   !< max Euler distance
real                      :: trs=0.          !< half shift range

contains

    !>  \brief  is a constructor
    subroutine comlin_symsrch_init( b, p, opt_str, lims, o_vertex, dmax, mode )
        use simple_build,        only: build
        use simple_params,       only: params
        class(build), target,  intent(in)    :: b
        class(params), target, intent(in)    :: p
        character(len=*),      intent(in)    :: opt_str
        real,                  intent(in)    :: lims(5,2)
        type(ori),             intent(inout) :: o_vertex    ! inital position (is an ico vertex)
        real,                  intent(in)    :: dmax        ! Max angular distance from inital position
        character(len=*),      intent(in)    :: mode
        character(len=8) :: str_opt= 'simplex'
        ! set constants & pointers:
        nptcls     = p%nptcls
        hp         = p%hp
        lp         = p%lp
        euldistmax = dmax
        bp         => b
        iptcl      => p%iptcl
        jptcl      => p%jptcl
        dynlim     => p%lp_dyn
        trs        = p%trs
        eullims    = lims(:3,:2)
        ! OPTIMISER INIT
        a_copy = bp%a
        vertex = o_vertex
        optori = o_vertex
        ! make optimizer spec
        if( opt_str.ne.'simplex' .and. opt_str.ne.'de' .and. opt_str.ne.'oasis')then
            write(*,*)'Unsupported minimizer in simple_comlin_symsrch; comlin_symsrch_init'
            stop
        else
            str_opt = opt_str
        endif
        ! Optimizer init
        if( mode .eq. 'sym' )then
            call ospec%specify(str_opt, 3, ftol=1e-4, gtol=1e-4, limits=lims(:3,:), nrestarts=3)
            call ospec%set_costfun(comlin_symsrch_cost)
        elseif( mode .eq. 'pair')then
            call ospec%specify(str_opt, 5, ftol=1e-4, gtol=1e-4, limits=lims, nrestarts=1, maxits=30)
            call ospec%set_costfun(comlin_pairsrch_cost)
        else
            stop 'Unsupported mode; simple_comlin_symsrch::comlin_symsrch_init'
        endif
        call ofac%new(ospec, nlopt)
    end subroutine comlin_symsrch_init

    function comlin_symsrch_minimize() result( cxy )
        real :: cxy(4)
        ospec%x = vertex%get_euler()
        ospec%nevals = 0
        call nlopt%minimize(ospec, cxy(1))
        cxy(1)  = -cxy(1) ! correlation 
        cxy(2:) = ospec%x ! euler set
    end function comlin_symsrch_minimize
    
    function comlin_pairsrch_minimize() result( cxy )
        real :: cxy(6)
        ospec%x(1:3) = vertex%get_euler()
        ospec%x(4)   = 0.
        ospec%x(5)   = 0.
        ospec%nevals = 0
        call nlopt%minimize(ospec, cxy(1))
        cxy(1)  = -cxy(1) ! correlation 
        cxy(2:) = ospec%x ! ori set
    end function comlin_pairsrch_minimize
    
    function comlin_symsrch_cost( vec, D) result( cost )
        integer, intent(in)  :: D
        real,    intent(in)  :: vec(D)
        real    :: cost
        integer :: i
        cost = 0.
        call optori%set_euler(vec)
        if(any(vec(1:3)-eullims(:,1) < 0.))then     ! lower limts
            cost = 1.
        elseif(any(vec(1:3)-eullims(:,2) > 0.))then ! upper limits
            cost = 1.
        elseif((optori.euldist.vertex)>euldistmax)then
            cost = 1.
        else
            call bp%a%rot(optori)
            call bp%se%apply2all(bp%a)
            cost = -bp%clins%corr(dynlim)
            bp%a = a_copy
        endif
    end function comlin_symsrch_cost
    
    function comlin_pairsrch_cost( vec, D) result( cost )
        integer, intent(in)  :: D
        real,    intent(in)  :: vec(D)
        real :: cost
        cost = 0.
        call optori%set_euler(vec(1:3))
        if(any(vec(1:3)-eullims(:,1) < 0.))then     ! lower limts
            cost = 1.
        elseif(any(vec(1:3)-eullims(:,2) > 0.))then ! upper limits
            cost = 1.
        elseif((optori.euldist.vertex) > euldistmax)then
            cost = 1.
        elseif( vec(4) < -trs .or. vec(4) > trs )then 
            cost = 1.
        elseif( vec(5) < -trs .or. vec(5) > trs )then 
            cost = 1.
        else
            call bp%a%set_euler(jptcl, vec(1:3))
            call bp%a%set_shift(jptcl, vec(4:5))
            cost = -bp%clins%pcorr(iptcl, jptcl, dynlim) 
        endif
    end function comlin_pairsrch_cost

    function comlin_symsrch_get_nevals() result( nevals )
        integer :: nevals
        nevals = ospec%nevals
    end function comlin_symsrch_get_nevals

    ! OLD SYMMETRY SRCH ROUTINE

    !>  \brief  is for finding the symmetry axis give an aligned set of images
    subroutine comlin_srch_symaxis( orientation_best, doprint )
        use simple_ori,  only: ori
        use simple_oris, only: oris
        class(ori), intent(inout) :: orientation_best
        logical,    intent(in)    :: doprint
        integer    :: i, j, k
        type(ori)  :: orientation
        type(oris) :: a_copy
        real       :: corr, corr_best
        integer    :: ntot, cnt, lims(3,2)
        if( doprint ) write(*,'(A)') '>>> GLOBAL SYMMETRY AXIS SEARCH'
        a_copy    = bp%a
        corr_best = -1.
        ntot      = 24624
        cnt       = 0
        lims(1,1) = 0
        lims(1,2) = 359
        lims(2,1) = 0
        lims(2,2) = 180
        lims(3,1) = 0
        lims(3,2) = 359
        call orientation%new
        do i=lims(1,1),lims(1,2),10
            call orientation%e1set(real(i))
            do j=lims(2,1),lims(2,2),10
                call orientation%e2set(real(j))
                do k=lims(3,1),lims(3,2),10
                    cnt = cnt+1
                    if( doprint ) call progress(cnt, ntot)
                    call orientation%e3set(real(k))
                    bp%a = a_copy
                    call bp%a%rot(orientation)
                    call bp%se%apply2all(bp%a)
                    corr = bp%clins%corr(dynlim)
                    call orientation%set('corr',corr)
                    if( corr > corr_best )then
                        corr_best = corr
                        orientation_best = orientation
                    endif
                end do
            end do
        end do
        if( doprint )then
            write(*,'(A)') '>>> FOUND FIRST SYMMETRY AXIS ORIENTATION'
            call orientation_best%print
            write(*,'(A)') '>>> REFINED SYMMETRY AXIS SEARCH'
        endif
        lims(1,1) = max(nint(orientation_best%e1get()-10.),lims(1,1))
        lims(1,2) = min(nint(orientation_best%e1get()+10.),lims(1,2))
        lims(2,1) = max(nint(orientation_best%e2get()-10.),lims(2,1))
        lims(2,2) = min(nint(orientation_best%e2get()+10.),lims(2,2))
        lims(3,1) = max(nint(orientation_best%e3get()-10.),lims(3,1))
        lims(3,2) = min(nint(orientation_best%e3get()+10.),lims(3,2))
        ntot      = (lims(1,2)-lims(1,1)+1)*(lims(2,2)-lims(2,1)+1)*(lims(3,2)-lims(3,1)+1)
        cnt       = 0 
        do i=lims(1,1),lims(1,2)
            call orientation%e1set(real(i))
            do j=lims(2,1),lims(2,2)
                call orientation%e2set(real(j))
                do k=lims(3,1),lims(3,2)
                    cnt = cnt+1
                    if( doprint ) call progress(cnt, ntot)
                    call orientation%e3set(real(k))
                    bp%a = a_copy
                    call bp%a%rot(orientation)
                    call bp%se%apply2all(bp%a)
                    corr = bp%clins%corr(dynlim)
                    call orientation%set('corr',corr)
                    if( corr > corr_best )then
                        corr_best = corr
                        orientation_best = orientation
                    endif
                end do
            end do
        end do
        if( doprint )then
            write(*,'(A)') '>>> FOUND REFINED SYMMETRY AXIS ORIENTATION'
            call orientation_best%print
        endif
        bp%a = a_copy
        call a_copy%kill
    end subroutine comlin_srch_symaxis
    
end module simple_comlin_symsrch
