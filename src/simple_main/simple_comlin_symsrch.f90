!==Class simple_comlin_symsrch
!
! The simple_comlin_symsrch singleton provides the  common line continuous symmetry
! search. The code is distributed with the hope 
! that it will be useful, but _WITHOUT_ _ANY_ _WARRANTY_. Redistribution or modification 
! is regulated by the GNU General Public License. *Author:* Hans Elmlund, 2009-06-11.
module simple_comlin_symsrch
use simple_defs
use simple_comlin_corr  ! use all in there
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
integer, pointer          :: ptcl=>null()    !< ptcl index
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
        class(build), intent(in), target  :: b
        class(params), intent(in), target :: p
        type(ori)                         :: o_vertex    ! inital position (is an ico vertex)
        real, intent(in)                  :: lims(5,2)
        real, intent(in)                  :: dmax        ! Max angular distance from inital position
        character(len=*), intent(in)      :: opt_str, mode
        character(len=8)                  :: str_opt= 'simplex'
        ! set constants & pointers:
        nptcls     = p%nptcls
        hp         = p%hp
        lp         = p%lp
        euldistmax = dmax
        bp         => b
        ptcl       => p%ptcl
        iptcl      => p%iptcl
        jptcl      => p%jptcl
        dynlim     => p%lp_dyn
        trs        = p%trs
        eullims    = lims(:3,:2)
        ! make comlin_corr functionality:
        call comlin_corr_init( bp , ptcl, dynlim )
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
    end subroutine

    function comlin_symsrch_minimize() result( cxy )
        real :: cxy(4)
        ospec%x = vertex%get_euler()
        ospec%nevals = 0
        call nlopt%minimize(ospec, cxy(1))
        cxy(1)  = -cxy(1) ! correlation 
        cxy(2:) = ospec%x ! euler set
    end function
    
    function comlin_pairsrch_minimize() result( cxy )
        real :: cxy(6)
        ospec%x(1:3) = vertex%get_euler()
        ospec%x(4)   = 0.
        ospec%x(5)   = 0.
        ospec%nevals = 0
        call nlopt%minimize(ospec, cxy(1))
        cxy(1)  = -cxy(1) ! correlation 
        cxy(2:) = ospec%x ! ori set
    end function
    
    function comlin_symsrch_cost( vec, D) result( cost )
        integer, intent(in)  :: D
        real, intent(in)     :: vec(D)
        real                 :: cost
        integer              :: i
        cost = 0.
        call optori%set_euler(vec)
        if( (vec(3)<eullims(3,1)) .or. (vec(3)>eullims(3,2)) )then
            cost = 1.
        elseif( (optori.euldist.vertex)>euldistmax )then
            cost = 1.
        else
            call bp%a%rot(optori)
            call bp%se%apply2all(bp%a)
            do i=1,nptcls
                ptcl = i
                cost = cost+pcorr_comlin()
            end do
            cost = -cost/real(nptcls)
            bp%a = a_copy
        endif
    end function
    
    function comlin_pairsrch_cost( vec, D) result( cost )
        integer, intent(in)  :: D
        real, intent(in)     :: vec(D)
        real                 :: cost
        cost = 0.
        call optori%set_euler(vec(1:3))
        if( (vec(3)<eullims(3,1)) .or. (vec(3)>eullims(3,2)) )then
            cost = 1.
        elseif( (optori.euldist.vertex)>euldistmax )then
            cost = 1.
        elseif( vec(4) < -trs .or. vec(4) > trs )then 
            cost = 1.
        elseif( vec(5) < -trs .or. vec(5) > trs )then 
            cost = 1.
        else
            call bp%a%set_euler(iptcl, [0.,0.,0.])
            call bp%a%set(iptcl, 'x', 0.)
            call bp%a%set(iptcl, 'y', 0.)
            call bp%a%set_euler(jptcl, [vec(1),vec(2),vec(3)])
            call bp%a%set(jptcl, 'x', vec(4))
            call bp%a%set(jptcl, 'y', vec(5))
            cost = -pcorr_comlin(iptcl, jptcl) 
        endif
    end function

    function comlin_symsrch_get_nevals() result( nevals )
        integer :: nevals
        nevals = ospec%nevals
    end function
    
end module simple_comlin_symsrch
