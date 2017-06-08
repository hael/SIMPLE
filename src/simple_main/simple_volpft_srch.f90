module simple_volpft_srch
use simple_opt_spec,        only: opt_spec
use simple_pftcc_opt,       only: pftcc_opt
use simple_volpft_corrcalc, only: volpft_corrcalc
use simple_simplex_opt,     only: simplex_opt
implicit none

public :: volpft_srch_init, volpft_srch_minimize
private

logical, parameter :: DEBUG  = .false.
integer, parameter :: NPROJ  = 200
integer, parameter :: NBEST  = 20
integer, parameter :: ANGRES = 10

type(volpft_corrcalc)           :: vpftcc               !< corr calculator
type(opt_spec)                  :: ospec                !< optimizer specification object
type(simplex_opt)               :: nlopt                !< optimizer object
class(volpft_corrcalc), pointer :: vpftcc_ptr => null() !< vpftcc object
integer                         :: nrestarts = 3        !< simplex restarts (randomized bounds)
logical                         :: serial = .true.      !< controls shared-mem parallellization of corr calc

contains
    
    subroutine volpft_srch_init( vol_ref, vol_target, hp, lp, nrestarts_in )
        use simple_projector, only: projector
        class(projector),  intent(in) :: vol_ref, vol_target
        real,              intent(in) :: hp, lp
        integer, optional, intent(in) :: nrestarts_in
        real :: lims(3,2)
        ! make the corr calculator
        call vpftcc%new( vol_ref, vol_target, hp, lp )
        ! set nrestarts
        nrestarts = 3
        if( present(nrestarts_in) ) nrestarts = nrestarts_in
        ! make optimizer spec
        lims = 0.
        lims(1,2) = 359.99
        lims(2,2) = 180.
        lims(3,2) = 359.99
        call ospec%specify('simplex', 3, ftol=1e-4,&
        &gtol=1e-4, limits=lims, nrestarts=nrestarts, maxits=30)
        call ospec%set_costfun(volpft_srch_costfun)
        ! generate the simplex optimizer object 
        call nlopt%new(ospec)
    end subroutine volpft_srch_init

    function volpft_srch_costfun( vec, D ) result( cost )
        use simple_ori, only: ori
        integer, intent(in) :: D
        real,    intent(in) :: vec(D)
        type(ori) :: e
        real      :: cost
        call e%new
        call e%set_euler(vec)
        cost = -vpftcc%corr(e, serial)
    end function volpft_srch_costfun

    !> \brief  minimization of the cost function
    function volpft_srch_minimize( fromto ) result( orientation_best )
        use simple_oris, only: oris
        use simple_ori,  only: ori
        integer, optional, intent(in) :: fromto(2)
        integer                       :: ffromto(2), ntot, iproj, config_best(2), iloc, inpl
        real                          :: euls(3), corr_best, cost
        logical                       :: distr_exec
        type(oris)                    :: espace, resoris
        type(ori)                     :: orientation, orientation_best
        integer, allocatable          :: order(:)
        real,    allocatable          :: corrs(:,:)
        ! flag distributed/shmem exec & set range
        distr_exec = present(fromto)
        if( distr_exec )then
            ffromto = fromto
        else
            ffromto(1) = 1
            ffromto(2) = NPROJ
        endif
        if( ffromto(1) < 1 .or. ffromto(2) > NPROJ )then
            stop 'range out of bound; simple_volpft_srch :: volpft_srch_minimize'
        endif
        ntot = ffromto(2) - ffromto(1) + 1
        ! create
        call orientation%new
        call resoris%new(NPROJ)
        call espace%new(NPROJ)
        call espace%spiral
        allocate(corrs(ffromto(1):ffromto(2),0:359))
        ! grid search using the spiral geometry & ANGRES degree in-plane resolution
        serial = .true. ! since the below loop is parallel
        corrs  = -1.
        !$omp parallel do schedule(static) default(shared) private(iproj,euls,inpl) proc_bind(close)
        do iproj=ffromto(1),ffromto(2)
            euls = espace%get_euler(iproj)
            do inpl=0,359,ANGRES
                euls(3) = real(inpl)
                corrs(iproj,inpl) = vpftcc%corr(euls, serial)
            end do
        end do
        !$omp end parallel do
        ! identify the best candidates (serial code)
        do iproj=ffromto(1),ffromto(2)
            corr_best   = -1.
            config_best = 0
            do inpl=0,359,ANGRES
                if( corrs(iproj,inpl) > corr_best )then
                    corr_best = corrs(iproj,inpl) 
                    config_best(1) = iproj
                    config_best(2) = inpl
                endif
            end do
            ! set local in-plane optimum for iproj
            orientation = espace%get_ori(config_best(1))
            call orientation%e3set(real(config_best(2)))
            call orientation%set('corr', corr_best)
            call resoris%set_ori(iproj,orientation)
        end do
        serial = .false. ! since the simplex search (below) is not parallel, we parallelise the cost eval
        if( distr_exec )then
            ! refine all local optima
            do iloc=ffromto(1),ffromto(2)
                ospec%x = resoris%get_euler(iloc)
                call nlopt%minimize(ospec, cost)
                call resoris%set_euler(iloc, ospec%x)
                call resoris%set(iloc, 'corr', -cost)
            end do
        else
            ! refine the NBEST local optima
            ! order the local optima according to correlation
            order = resoris%order_corr()
            ! refine the NBEST solutions
            do iloc=1,NBEST
                ospec%x = resoris%get_euler(order(iloc))
                call nlopt%minimize(ospec, cost)
                call resoris%set_euler(order(iloc), ospec%x)
                call resoris%set(order(iloc), 'corr', -cost)
            end do
        endif
        ! order the local optima according to correlation
        order = resoris%order_corr()
        ! return best
        orientation_best = resoris%get_ori(order(1))
    end function volpft_srch_minimize

end module simple_volpft_srch
