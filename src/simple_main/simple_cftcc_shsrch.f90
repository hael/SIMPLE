!------------------------------------------------------------------------------!
! SIMPLE v2.5         Elmlund & Elmlund Lab          simplecryoem.com          !
!------------------------------------------------------------------------------!
!> Simple CTF optimisation module:  basic shift search
module simple_cftcc_shsrch
use simple_opt_factory,     only: opt_factory
use simple_opt_spec,        only: opt_spec
use simple_optimizer,       only: optimizer
use simple_cartft_corrcalc, only: cartft_corrcalc
use simple_image,           only: image
use simple_ori,             only: ori
use simple_oris,             only: oris
use simple_defs
implicit none

public :: cftcc_shsrch_init, cftcc_shsrch_set_state, cftcc_shsrch_minimize
private

type(opt_factory)               :: ofac              !< optimizer factory
type(opt_spec)                  :: ospec             !< optimizer specification object
class(optimizer), pointer       :: nlopt    =>null() !< pointer to nonlinear optimizer
class(cartft_corrcalc), pointer :: cftcc_ptr=>null() !< pointer to cftcc object
class(image), pointer           :: pimg     =>null() !< pointer to image
integer                         :: state=1           !< state to evaluate
contains

    subroutine cftcc_shsrch_init( cftcc, img, opt_str, trs, nrestarts )
        class(cartft_corrcalc), target, intent(in) :: cftcc
        class(image),           target, intent(in) :: img
        character(len=*),               intent(in) :: opt_str
        real,                           intent(in) :: trs
        integer,                        intent(in) :: nrestarts
        real :: lims(2,2)
        ! make optimizer spec
        lims(:,1) = -trs
        lims(:,2) =  trs
        call ospec%specify(opt_str, 2, ftol=1e-4, gtol=1e-4, limits=lims, nrestarts=nrestarts)
        ! set optimizer cost function
        call ospec%set_costfun(cftcc_shsrch_cost)
        ! generate optimizer object with the factory
        call ofac%new(ospec, nlopt)
        ! set pointers
        cftcc_ptr => cftcc
        pimg      => img
    end subroutine cftcc_shsrch_init
    
    subroutine cftcc_shsrch_set_state( state_in )
        integer, intent(in) :: state_in
        state = state_in
    end subroutine cftcc_shsrch_set_state
    
    function cftcc_shsrch_get_nevals() result( nevals )
        integer :: nevals
        nevals = ospec%nevals
    end function cftcc_shsrch_get_nevals
    
    function cftcc_shsrch_cost( vec, D ) result( cost )
        use simple_ori, only: ori
        use simple_math, only: rad2deg
        integer, intent(in) :: D
        real,    intent(in) :: vec(D)
        real      :: cost, shvec(3)
        integer   :: i
        ! enforce the barrier constraint for the shifts
        do i=1,2
            if( vec(i) < ospec%limits(i,1) .or. vec(i) > ospec%limits(i,2) )then
                cost = 1.
                return
            endif
        end do
        ! calculate cost
        shvec(1:2) = vec(1:2)
        shvec(3)   = 0.0
        cost = -cftcc_ptr%correlate(pimg, 1, shvec)
    end function cftcc_shsrch_cost
    
    subroutine cftcc_shsrch_minimize( o, os )
        use simple_math, only: rad2deg, arg
        class(ori),  intent(inout) :: o
        class(oris), intent(inout) :: os
        type(ori)         :: o_new
        real, allocatable :: vertices(:,:), costs(:)
        real    :: corr, cost, dist, prev_shvec(2), wcorr
        integer :: i
        prev_shvec = o%get_shift()
        ! project the input orientation
        call cftcc_ptr%project(o, 1)
        ! initialise optimiser
        ospec%x = 0.
        ! search
        call nlopt%minimize(ospec, cost)
        if( ospec%str_opt.eq.'simplex')then
            ! soft matching
            o_new = o
            call nlopt%get_vertices( ospec, vertices, costs )
            call os%new(3)
            do i = 1,3
                o_new = o
                call o_new%set('corr', -costs(i))
                call o_new%set_shift(prev_shvec-vertices(i,:))
                call os%set_ori(i, o_new)
            enddo
            call prep_soft_oris( wcorr, os )
            corr = wcorr
            call o%set('ow', os%get(1,'ow'))
            deallocate(vertices, costs)
        else
            ! hard matching
            corr = -cost
            call o%set('ow', 1.)
        endif
        call o%set('corr', corr)
        ! shifts must be obtained by vector addition
        ! the reference is rotated upon projection: no need to rotate the shift
        call o%set_shift(prev_shvec - ospec%x)
        ! distance
        dist = arg( ospec%x )
        call o%set('dist',dist)
        ! set the overlaps
        call o%set('mi_class', 1.)
        call o%set('mi_inpl',  1.)
        call o%set('mi_state', 1.)
        call o%set('mi_joint', 1.)
        call o%set('dist_inpl', 0.)
        ! all the other stuff
        call o%set('frac', 100.)
        call o%set('mirr', 0.)
        call o%set('sdev', 0.)
    end subroutine cftcc_shsrch_minimize

    subroutine prep_soft_oris( wcorr, os )
        type(oris), intent(inout) :: os
        real,       intent(out)   :: wcorr
        real              :: ws(3)
        real, allocatable :: corrs(:)
        integer           :: i
        if( os%get_noris() /= 3 )stop 'invalid number of orientations; simple_cftcc_shsrch::prep_soft_oris'
        ws    = 0.
        wcorr = 0.
        ! get unnormalised correlations
        corrs = os%get_all('corr')
        ! calculate normalised weights and weighted corr
        where( corrs > TINY ) ws = exp(corrs) ! ignore invalid corrs
        ws    = ws/sum(ws)
        wcorr = sum(ws*corrs) 
        ! update npeaks individual weights
        do i=1,3
            call os%set(i,'ow',ws(i))
        enddo
        deallocate(corrs)
    end subroutine prep_soft_oris

end module simple_cftcc_shsrch
