! Cartesian volume-volume alignment based on band-pass limited cross-correlation
module simple_vol_srch
include 'simple_lib.f08'
use simple_image,       only: image
use simple_opt_spec,    only: opt_spec
use simple_opt_simplex, only: opt_simplex
use simple_projector,   only: projector
implicit none

public :: vol_srch_init, vol_shsrch_minimize
private

logical, parameter :: INI_W_DISCR_SRCH = .false.

type(opt_spec)        :: ospec           !< optimizer specification object
type(opt_simplex)     :: nlopt           !< optimizer object
integer               :: nrestarts = 3   !< simplex restarts (randomized bounds)
real                  :: hp, lp          !< srch ctrl params
real                  :: lims(3,2)       !< variable bounds
class(image), pointer :: vref  => null() !< reference volume
class(image), pointer :: vtarg => null() !< target volume (subjected to shift)

contains

    subroutine vol_srch_init( vol_ref, vol_target, hp_in, lp_in, trs_in, nrestarts_in )
        use simple_ori, only: ori
        class(image), target, intent(inout) :: vol_ref, vol_target
        real,                 intent(in)    :: hp_in, lp_in, trs_in
        integer, optional, intent(in)       :: nrestarts_in
        integer :: ldim(3), ldim_pd(3), boxpd
        real    :: smpd, eul(3)
        hp  = hp_in
        lp  = lp_in
        nrestarts = 3
        if( present(nrestarts_in) ) nrestarts = nrestarts_in
        vref  => vol_ref
        vtarg => vol_target
        lims(1:3,1) = -trs_in
        lims(1:3,2) =  trs_in
        call ospec%specify('simplex', 3, ftol=1e-4,&
        &gtol=1e-4, limits=lims, nrestarts=nrestarts, maxits=100)
        call ospec%set_costfun(vol_shsrch_costfun)
        call nlopt%new(ospec)
    end subroutine vol_srch_init

    function vol_shsrch_costfun( fun_self, vec, D ) result( cost )
        class(*), intent(inout) :: fun_self
        integer,  intent(in)    :: D
        real,     intent(in)    :: vec(D)
        real :: cost
        if( any(vec < lims(:,1)) .or. any(vec > lims(:,2)) )then
            cost = 1.0
            return
        endif
        cost = - vref%corr_shifted(vtarg, vec(1:3), lp, hp)
    end function vol_shsrch_costfun

    function vol_shsrch_minimize( ) result( cxyz )
        real :: cost_init, cost, cxyz(4), cost_best
        real :: shvec(3), shvec_best(3), xsh, ysh, zsh
        class(*), pointer :: fun_self => null()
        shvec_best = 0.
        if( INI_W_DISCR_SRCH )then
            ! discrete search to start-off with (half-pixel resolution)
            cost_best = 1.
            xsh = lims(1,1)
            ysh = lims(2,1)
            zsh = lims(3,1)
            do while( xsh <= lims(1,2) )
                do while( ysh <= lims(2,2) )
                    do while( zsh <= lims(3,2) )
                        shvec = [xsh,ysh,zsh]
                        cost = vol_shsrch_costfun(fun_self, shvec, ospec%ndim)
                        if( cost < cost_best )then
                            shvec_best = shvec
                            cost_best  = cost
                        endif
                        zsh = zsh + 0.5
                    end do
                    ysh = ysh + 0.5
                end do
                xsh = xsh + 0.5
            end do
        endif
        ! refinement with simplex
        ospec%x   = shvec_best ! assumed that vol is shifted to previous centre
        cost_init = vol_shsrch_costfun(fun_self, ospec%x, ospec%ndim)
        call nlopt%minimize(ospec, fun_self, cost)
        if( cost < cost_init )then
            cxyz(1)  = -cost    ! correlation
            cxyz(2:) =  ospec%x ! shift
        else
            cxyz(1)  = -cost_init
            cxyz(2:) =  0.
        endif
    end function vol_shsrch_minimize

end module simple_vol_srch
