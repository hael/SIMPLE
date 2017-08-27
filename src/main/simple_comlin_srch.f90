! for common-lines-based search
module simple_comlin_srch
use simple_defs
use simple_build,       only: build
use simple_params,      only: params
use simple_optimizer,   only: optimizer
use simple_opt_factory, only: opt_factory
use simple_opt_spec,    only: opt_spec
use simple_syslib,      only: alloc_errchk, simple_stop
use simple_jiffys,      only: progress
use simple_ori,         only: ori
use simple_oris,        only: oris
use simple_sym,         only: sym
implicit none

public :: comlin_srch_init, comlin_srch_get_nproj, comlin_srch_get_nbest, comlin_srch_write_resoris,&
comlin_srch_symaxis, comlin_srch_pair
private
#include "simple_local_flags.inc"
integer, parameter :: NPROJ      = 200
integer, parameter :: NBEST      = 20
integer, parameter :: NBEST_PAIR = 10
integer, parameter :: ANGRES     = 10

class(build),     pointer :: bp=>null()       !< pointer to builder
class(params),    pointer :: pp=>null()       !< pointer to params
integer,          pointer :: iptcl=>null()    !< ptcl index
integer,          pointer :: jptcl=>null()    !< ptcl index
type(opt_factory)         :: ofac             !< optimizer factory
type(opt_spec)            :: ospec            !< optimizer specification object
class(optimizer), pointer :: nlopt=>null()    !< pointer to nonlinear optimizer
type(oris)                :: a_copy           !< safe keeping for cost function
type(oris)                :: resoris          !< for results storage and ranking
type(oris)                :: espace           !< projection directions
type(ori)                 :: oref             !< starting projection direction
integer                   :: nptcls=0         !< nr of particles
real                      :: eullims(3,2)=0.  !< Euler limits; only used to setup optimizer
real                      :: optlims(5,2)=0.  !< Euler limits; only used to setup optimizer
real                      :: hp=0.            !< fixed high-pass limit
real                      :: lp=0.            !< fixed low-pass limit
real                      :: trs=0.           !< half shift range

contains

    !> common-line mode search constructor
    subroutine comlin_srch_init( b, p, opt_str, mode )
        class(build),  target, intent(in) :: b           !< build object
        class(params), target, intent(in) :: p           !< Parameters
        character(len=*),      intent(in) :: opt_str     !< 'simplex', 'de' or 'oasis' search type
        character(len=*),      intent(in) :: mode        !< 'sym' or 'pair' mode
        character(len=8) :: str_opt= 'simplex'
        bp     => b
        pp     => p
        iptcl  => p%iptcl
        jptcl  => p%jptcl
        nptcls =  p%nptcls
        hp     =  p%hp
        lp     =  p%lp
        trs    =  p%trs
        a_copy = bp%a
        call resoris%new(NPROJ)
        call espace%new(NPROJ)
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
            eullims(2,2)  = 90.
            call espace%spiral(2,eullims)
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
        call ofac%new(ospec, nlopt)
    end subroutine comlin_srch_init

    function comlin_srch_get_nproj( ) result( nprojout )
        integer :: nprojout
        nprojout = NPROJ
    end function comlin_srch_get_nproj

    function comlin_srch_get_nbest( ) result( nbestout )
        integer :: nbestout
        nbestout = NBEST
    end function comlin_srch_get_nbest

    !> Write results of common-line mode search to file
    subroutine comlin_srch_write_resoris( fname, fromto )
        use simple_binoris_io, only: binwrite_oritab
        character(len=*) :: fname
        integer          :: fromto(2)
        call binwrite_oritab(fname, resoris, fromto)
    end subroutine comlin_srch_write_resoris

    !> common-line mode search symmetry axis
    subroutine comlin_srch_symaxis( orientation_best, fromto )
        class(ori)           :: orientation_best
        integer, optional    :: fromto(2)
        type(ori)            :: orientation
        real                 :: corr, corr_best, cxy(4)
        integer              :: iproj, inpl, iloc, ffromto(2), cnt, ntot
        integer, allocatable :: order(:)
        logical              :: distr_exec
        ! flag distributed/shmem exec & set range
        distr_exec = present(fromto)
        if( distr_exec )then
            ffromto = fromto
        else
            ffromto(1) = 1
            ffromto(2) = NPROJ
        endif
        if( ffromto(1) < 1 .or. ffromto(2) > NPROJ )then
            stop 'range out of bound; simple_comlin_srch :: comlin_srch_symaxis'
        endif
        ntot = ffromto(2) - ffromto(1) + 1
        ! grid search using the spiral geometry & ANGRES degree in-plane resolution
        write(*,'(A)') '>>> GLOBAL GRID SYMMETRY AXIS SEARCH'
        cnt = 0
        do iproj=ffromto(1),ffromto(2)
            cnt = cnt + 1
            call progress(cnt, ntot)
            orientation = espace%get_ori(iproj)
            corr_best   = -1.
            do inpl=0,359,ANGRES
                call orientation%e3set(real(inpl))
                bp%a = a_copy
                call bp%a%rot(orientation)
                call bp%se%apply2all(bp%a)
                corr = bp%clins%corr()
                call orientation%set('corr',corr)
                if( corr > corr_best )then
                    corr_best = corr
                    orientation_best = orientation
                endif
            end do
            ! set local in-plane optimum for iproj
            call resoris%set_ori(iproj,orientation_best)
        end do
        if( distr_exec )then
            ! refine all local optima
            do iloc=ffromto(1),ffromto(2)
                orientation = resoris%get_ori(iloc)
                cxy = comlin_srch_minimize(orientation)
                call resoris%set_euler(iloc, cxy(2:))
                call resoris%set(iloc, 'corr', cxy(1))
            end do
            ! order the local optima according to correlation
            order = resoris%order_corr()
            ! return best
            orientation_best = resoris%get_ori(order(1))
        else
            write(*,'(A)') '>>> CONTINOUS SYMMETRY AXIS REFINEMENT'
            ! refine the NBEST local optima
            ! order the local optima according to correlation
            order = resoris%order_corr()
            ! refine the NBEST solutions
            do iloc=1,NBEST
                orientation = resoris%get_ori(order(iloc))
                cxy = comlin_srch_minimize(orientation)
                call resoris%set_euler(order(iloc), cxy(2:))
                call resoris%set(order(iloc), 'corr', cxy(1))
            end do
            ! re-order
            order = resoris%order_corr()
            ! return best
            orientation_best = resoris%get_ori(order(1))
            write(*,'(A)') '>>> FOUND REFINED SYMMETRY AXIS ORIENTATION'
            call orientation_best%print_ori()
        endif
    end subroutine comlin_srch_symaxis

    !> Calculate similarity between pairs
    function comlin_srch_pair() result( similarity )
        integer, allocatable :: order(:)
        type(ori) :: orientation, orientation_best
        real      :: corr, corr_best, similarity, cxy(6)
        integer   :: iproj, inpl, iloc
        ! grid search using the spiral geometry & ANGRES degree in-plane resolution
        do iproj=1,NPROJ
            orientation = espace%get_ori(iproj)
            corr_best   = -1.
            do inpl=0,359,ANGRES
                call orientation%e3set(real(inpl))
                call bp%a%set_ori(jptcl, orientation)
                corr = bp%clins%pcorr(iptcl, jptcl) 
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
        class(ori) :: otarget  !< orientation target
        real       :: cxy(4)   !< output correlation and spec coordinates
        oref    = otarget
        call oref%swape1e3
        ospec%x = otarget%get_euler()
        ospec%nevals = 0
        call nlopt%minimize(ospec, cxy(1))
        cxy(1)  = -cxy(1) ! correlation 
        cxy(2:) = ospec%x ! euler set
    end function comlin_srch_minimize

    !> Pair search minimisation
    function comlin_pairsrch_minimize( otarget ) result( cxy )
        class(ori) :: otarget
        real       :: cxy(6)
        ospec%x(1:3) = otarget%get_euler()
        ospec%x(4)   = 0.
        ospec%x(5)   = 0.
        ospec%nevals = 0
        call nlopt%minimize(ospec, cxy(1))
        cxy(1)  = -cxy(1) ! correlation 
        cxy(2:) = ospec%x ! ori set
    end function comlin_pairsrch_minimize

    !> Common-line mode search cost function
    function comlin_srch_cost( vec, D ) result( cost )
        integer, intent(in)  :: D     !< size of parameter vector
        real,    intent(in)  :: vec(D)!<  parameter vector
        type(ori) :: o, oswap
        real      :: cost, euldist
        cost = 0.
        if(any(vec(1:3)-eullims(:,1) < 0.))then     ! lower limits
            cost = 1.
        elseif(any(vec(1:3)-eullims(:,2) > 0.))then ! upper limits
            cost = 1.
        else
            call o%set_euler(vec)
            !
            oswap = o
            call oswap%swape1e3
            euldist = oswap.euldist.oref
            if( euldist > PI/3. )then
                cost = 1.
                return
            endif
            !
            call bp%a%rot(o)
            call bp%se%apply2all(bp%a)
            cost = -bp%clins%corr()
            bp%a = a_copy ! puts back unmodified oris
        endif
    end function comlin_srch_cost
    
     !> Pair search cost function
    function comlin_pairsrch_cost( vec, D ) result( cost )
        integer, intent(in)  :: D   !< size of parameter vector
        real,    intent(in)  :: vec(D) !< parameter vector
        real :: cost
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
            call bp%a%set_euler(jptcl, vec(1:3))
            call bp%a%set_shift(jptcl, vec(4:5))
            cost = -bp%clins%pcorr(iptcl, jptcl) 
        endif
    end function comlin_pairsrch_cost

end module simple_comlin_srch
