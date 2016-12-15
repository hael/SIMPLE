module simple_volpft_srch
use simple_opt_factory,     only: opt_factory
use simple_opt_spec,        only: opt_spec
use simple_optimizer,       only: optimizer
use simple_volpft_corrcalc, only: volpft_corrcalc
use simple_ori,             only: ori
implicit none

public :: volpft_srch_init, volpft_6dimsrch
private

type(opt_factory)               :: ofac                !< optimizer factory for Euler angles
type(opt_spec)                  :: ospec               !< optimizer specification object for Euler angles
class(optimizer), pointer       :: nlopt      =>null() !< pointer to nonlinear Euler optimizer
class(volpft_corrcalc), pointer :: vpftcc_ptr =>null() !< pointer to pftcc object
real, allocatable               :: corrs(:)            !< correlations
integer, allocatable            :: inds(:)             !< orientation indices
real                            :: trs                 !< shift half-range
integer, parameter              :: NRESTARTS=5         !< number of restarts

contains
    
    subroutine volpft_srch_init( vpftcc, trs_in )
        class(volpft_corrcalc), target, intent(in) :: vpftcc
        real,                           intent(in) :: trs_in
        real :: lims(6,2)
        trs = trs_in
        ! make optimizer specs
        lims      = 0.
        lims(1,2) = 359.99
        lims(2,2) = 180.
        lims(3,2) = 359.99
        lims(4,1) = -trs
        lims(4,2) =  trs
        lims(5,1) = -trs
        lims(5,2) =  trs
        lims(6,1) = -trs
        lims(6,2) =  trs
        call ospec%specify('simplex', 6, ftol=1e-4, gtol=1e-4, limits=lims, nrestarts=NRESTARTS)
        call ospec%set_costfun(vpftcc_cost)
        ! generate optimizer objects with the factory
        call ofac%new(ospec, nlopt)
        ! set pointer to corrcalc object
        vpftcc_ptr => vpftcc
    end subroutine volpft_srch_init
    
    subroutine volpft_6dimsrch( npeaks, corr_best, o_best )
        use simple_math, only: hpsort
        integer,   intent(in)  :: npeaks
        real,      intent(out) :: corr_best
        type(ori), intent(out) :: o_best
        integer   :: noris, ipeak, iori
        type(ori) :: e
        real      :: cost, cost_best, corr, shvec_best(3), rotmat(3,3)
        ! discrete search over Euler angles
        noris = vpftcc_ptr%get_nspace()
        if( allocated(corrs) ) deallocate(corrs)
        if( allocated(inds)  ) deallocate(inds)
        allocate( corrs(noris), inds(noris) )
        do iori=1,noris
            corrs(iori) = vpftcc_ptr%corr(iori)
            inds(iori)  = iori
        end do
        call hpsort(noris, corrs, inds)
        ! continuous search over all parameters
        call o_best%new
        cost_best = 1.
        do ipeak=noris,noris-npeaks+1,-1
            e = vpftcc_ptr%get_ori(inds(ipeak))
            ospec%x      = 0.
            ospec%x(1:3) = e%get_euler()
            call nlopt%minimize(ospec, cost)
            if( cost < cost_best )then
                cost_best = cost
                corr_best = -cost
                call o_best%set_euler(ospec%x(1:3))
                shvec_best = ospec%x(4:6)
            endif
        end do
        call o_best%set('x',shvec_best(1))
        call o_best%set('y',shvec_best(2))
        call o_best%set('z',shvec_best(3))
    end subroutine volpft_6dimsrch

    function vpftcc_cost( vec, D ) result( cost )
        integer, intent(in) :: D
        real,    intent(in) :: vec(D)
        real                :: cost
        type(ori)           :: e
        cost = 1.
        if( any(abs(vec(4:6)) > trs) ) return ! barrier constraint
        call e%new
        call e%set_euler(vec(1:3))
        cost = -vpftcc_ptr%corr(e, vec(4:6))
    end function vpftcc_cost

end module simple_volpft_srch
