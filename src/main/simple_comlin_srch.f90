! for common-lines-based search
module simple_comlin_srch
include 'simple_lib.f08'
use simple_singletons
use simple_optimizer,   only: optimizer
use simple_opt_factory, only: opt_factory
use simple_opt_spec,    only: opt_spec
use simple_ori,         only: ori
use simple_oris,        only: oris
use simple_sym,         only: sym
implicit none

public :: comlin_srch_init, comlin_srch_get_nproj, comlin_srch_get_nbest,&
comlin_coarsesrch_symaxis, comlin_singlesrch_symaxis, comlin_srch_pair
private
#include "simple_local_flags.inc"
integer, parameter :: NPROJC1    = 400
integer, parameter :: NPROJ      = 200
integer, parameter :: NBEST      = 30
integer, parameter :: NBEST_PAIR = 10
integer, parameter :: ANGSTEP    = 10

type(opt_factory)         :: ofac             !< optimizer factory
type(opt_spec)            :: ospec            !< optimizer specification object
class(optimizer), pointer :: nlopt=>null()    !< pointer to nonlinear optimizer
type(oris)                :: a_copy           !< safe keeping for cost function
type(oris)                :: resoris          !< for results storage and ranking
type(oris)                :: espace           !< projection directions
type(ori)                 :: oref             !< starting projection direction
integer                   :: nptcls=0         !< nr of particles
integer                   :: nproj_sym=0      !< number of reference projections within the asu
real                      :: eullims(3,2)=0.  !< Euler limits; only used to setup optimizer
real                      :: optlims(5,2)=0.  !< Euler limits; only used to setup optimizer
real                      :: hp=0.            !< fixed high-pass limit
real                      :: lp=0.            !< fixed low-pass limit
real                      :: trs=0.           !< half shift range

contains

    !> common-line mode search constructor
    subroutine comlin_srch_init( opt_str, mode )
        character(len=*),      intent(in) :: opt_str !< 'simplex', 'de' or 'oasis' search type
        character(len=*),      intent(in) :: mode    !< 'sym' or 'pair' mode
        character(len=8) :: str_opt= 'simplex'
        nptcls    =  p%nptcls
        hp        =  p%hp
        lp        =  p%lp
        trs       =  p%trs
        a_copy    =  b%a
        nproj_sym = comlin_srch_get_nproj( pgrp=trim(p%pgrp) )
        write(*,'(A,I3)') '>>> NUMBER OF REFERENCE PROJECTIONS IN ASU: ', nproj_sym
        call resoris%new(nproj_sym)
        call espace%new(nproj_sym)
        ! make optimizer spec
        if( opt_str.ne.'simplex' .and. opt_str.ne.'de' .and. opt_str.ne.'oasis' )then
            call simple_stop ('Unsupported minimizer in simple_comlin_srch; comlin_srch_init')
        else
            str_opt = opt_str
        endif
        ! set limits & create projection diection search space
        eullims        = 0.
        eullims(:,2)   = [360., 180., 360.]
        if( mode .eq. 'sym' )then
            eullims(2,2) = 90.
            call espace%spiral(2, eullims)
        else
            call espace%spiral
        endif
        optlims(:3,:) = eullims
        optlims(4,:)  = [-p%trs,p%trs]
        optlims(5,:)  = [-p%trs,p%trs]
        ! Optimizer init
        if( mode .eq. 'sym' )then
            call ospec%specify(str_opt, 3, ftol=1e-4, gtol=1e-4, limits=eullims, nrestarts=3)
            call ospec%set_costfun(comlin_srch_cost)
        elseif( mode .eq. 'pair' )then
            call ospec%specify(str_opt, 5, ftol=1e-4, gtol=1e-4, limits=optlims, nrestarts=1, maxits=30)
            call ospec%set_costfun(comlin_pairsrch_cost)
        else
            call simple_stop ('Unsupported mode; simple_comlin_srch::comlin_srch_init')
        endif
        if (associated(nlopt)) then
            call nlopt%kill
            deallocate(nlopt)
        end if
        call ofac%new(ospec, nlopt)
    end subroutine comlin_srch_init

    function comlin_srch_get_nproj( pgrp ) result( nprojout )
        character(len=*), optional, intent(in) :: pgrp
        type(sym) :: se
        integer   :: nprojout
        if( present(pgrp) )then
            call se%new(trim(pgrp))
            nprojout = ceiling( real(NPROJC1) / real(se%get_nsym()) )
        else
            nprojout = NPROJ
        endif
    end function comlin_srch_get_nproj

    function comlin_srch_get_nbest( ) result( nbestout )
        integer :: nbestout
        nbestout = NBEST
    end function comlin_srch_get_nbest

    !> common-line mode search symmetry axis
    subroutine comlin_coarsesrch_symaxis( fromto, resoris )
        integer,     intent(in)    :: fromto(2)
        class(oris), intent(inout) :: resoris
        type(ori)            :: orientation, orientation_best
        real                 :: corr, corr_best
        integer              :: iproj, inpl, cnt, ntot
        if( fromto(1) < 1 .or. fromto(2) > nproj_sym )then
            stop 'range out of bound; simple_comlin_srch :: comlin_coarsesrch_symaxis'
        endif
        ntot = fromto(2) - fromto(1) + 1
        call resoris%new(nproj_sym)
        ! grid search using the spiral geometry & ANGSTEP degree in-plane resolution
        write(*,'(A)') '>>> GLOBAL GRID SYMMETRY AXIS SEARCH'
        cnt = 0
        do iproj=fromto(1),fromto(2)
            cnt = cnt + 1
            call progress(cnt, ntot)
            orientation = espace%get_ori(iproj)
            corr_best   = -1.
            do inpl=0,359,ANGSTEP
                call orientation%e3set(real(inpl))
                b%a = a_copy
                call b%a%rot(orientation)
                call b%se%apply2all(b%a)
                corr = b%clins%corr()
                call orientation%print_ori
                call orientation%set('corr',corr)
                if( corr > corr_best )then
                    corr_best = corr
                    orientation_best = orientation
                endif
            end do
            ! set local in-plane optimum for iproj
            call resoris%set_ori(iproj, orientation_best)
            call resoris%print_(iproj)
        end do
    end subroutine comlin_coarsesrch_symaxis

    !> continuous single common-line search symmetry axis
    subroutine comlin_singlesrch_symaxis( orientation )
        class(ori), intent(inout) :: orientation
        real :: cxy(4)
        cxy = comlin_srch_minimize(orientation)
        call orientation%set_euler(cxy(2:))
        call orientation%set('corr', cxy(1))
    end subroutine comlin_singlesrch_symaxis

    !> Calculate similarity between pairs
    function comlin_srch_pair() result( similarity )
        integer, allocatable :: order(:)
        type(ori) :: orientation, orientation_best
        real      :: corr, corr_best, similarity, cxy(6)
        integer   :: iproj, inpl, iloc
        ! grid search using the spiral geometry & ANGSTEP degree in-plane resolution
        do iproj = 1, nproj
            orientation = espace%get_ori(iproj)
            corr_best   = -1.
            do inpl=0,359,ANGSTEP
                call orientation%e3set(real(inpl))
                call b%a%set_ori(p%jptcl, orientation)
                corr = b%clins%pcorr(p%iptcl, p%jptcl)
                call orientation%set('corr',corr)
                if( corr > corr_best )then
                    corr_best = corr
                    orientation_best = orientation
                endif
            end do
            ! set local in-plane optimum for iproj
            call resoris%set_ori(iproj,orientation_best)
        end do
        ! refine the NBEST_PAIR local optima
        ! continuous refinement with rotational origin optimization
        ! first, order the local optima according to correlation
        order = resoris%order_corr()
        ! refine the NBEST_PAIR solutions
        do iloc=1,NBEST_PAIR
            orientation = resoris%get_ori(order(iloc))
            cxy = comlin_pairsrch_minimize(orientation)
            call resoris%set(order(iloc), 'corr', cxy(1))
        end do
        ! re-order
        order = resoris%order_corr()
        ! we are only interested in a similarity value, return best
        similarity = resoris%get(order(1), 'corr')
    end function comlin_srch_pair

    !> common-line mode search minimisation function
    function comlin_srch_minimize( otarget ) result( cxy )
        class(ori)        :: otarget  !< orientation target
        real              :: cxy(4)   !< output correlation and spec coordinates
        class(*), pointer :: fun_self => null()
        oref = otarget
        call oref%swape1e3
        ospec%x = otarget%get_euler()
        ospec%nevals = 0
        call nlopt%minimize(ospec, fun_self, cxy(1))
        cxy(1)  = -cxy(1) ! correlation
        cxy(2:) = ospec%x ! euler set
    end function comlin_srch_minimize

    !> Pair search minimisation
    function comlin_pairsrch_minimize( otarget ) result( cxy )
        class(ori)        :: otarget
        real              :: cxy(6)
        class(*), pointer :: fun_self => null()
        ospec%x(1:3) = otarget%get_euler()
        ospec%x(4)   = 0.
        ospec%x(5)   = 0.
        ospec%nevals = 0
        call nlopt%minimize(ospec, fun_self, cxy(1))
        cxy(1)  = -cxy(1) ! correlation
        cxy(2:) = ospec%x ! ori set
    end function comlin_pairsrch_minimize

    !> Common-line mode search cost function
    function comlin_srch_cost( fun_self, vec, D ) result( cost )
        class(*), intent(inout) :: fun_self      !< self-pointer for cost function
        integer,  intent(in)    :: D             !< size of parameter vector
        real,     intent(in)    :: vec(D)        !<  parameter vector
        type(ori)               :: o, oswap
        real                    :: cost, euldist
        cost = 0.
        if(any(vec(1:3)-eullims(:,1) < 0.))then     ! lower limits
            cost = 1.
        elseif(any(vec(1:3)-eullims(:,2) > 0.))then ! upper limits
            cost = 1.
        else
            call o%new
            call o%set_euler(vec)
            oswap = o
            call oswap%swape1e3
            euldist = oswap.euldist.oref
            if( euldist > PI/3. )then
                cost = 1.
                return
            endif
            !
            call b%a%rot(o)
            call b%se%apply2all(b%a)
            cost = -b%clins%corr()
            b%a = a_copy ! puts back unmodified oris
        endif
    end function comlin_srch_cost

     !> Pair search cost function
    function comlin_pairsrch_cost( fun_self, vec, D ) result( cost )
        class(*), intent(inout) :: fun_self !< self-pointer for cost function
        integer,  intent(in)    :: D        !< size of parameter vector
        real,     intent(in)    :: vec(D)   !< parameter vector
        real                    :: cost
        cost = 0.
        if(any(vec(1:3)-eullims(:,1) < 0.))then     ! lower limits
            cost = 1.
        elseif(any(vec(1:3)-eullims(:,2) > 0.))then ! upper limits
            cost = 1.
        elseif( vec(4) < -trs .or. vec(4) > trs )then
            cost = 1.
        elseif( vec(5) < -trs .or. vec(5) > trs )then
            cost = 1.
        else
            call b%a%set_euler(p%jptcl, vec(1:3))
            call b%a%set_shift(p%jptcl, vec(4:5))
            cost = -b%clins%pcorr(p%iptcl, p%jptcl)
        endif
    end function comlin_pairsrch_cost

end module simple_comlin_srch
