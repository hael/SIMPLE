! volume alignment based on band-pass limited cross-correlation
module simple_volpft_srch
use simple_opt_spec,        only: opt_spec
use simple_pftcc_opt,       only: pftcc_opt
use simple_volpft_corrcalc, only: volpft_corrcalc
use simple_simplex_opt,     only: simplex_opt
use simple_ori,             only: ori
implicit none

public :: volpft_srch_init, volpft_srch_minimize_eul, volpft_srch_minimize_shift, volpft_srch_minimize_all
private

integer, parameter :: NPROJ  = 200
integer, parameter :: NBEST  = 20
integer, parameter :: ANGSTEP = 10

type(volpft_corrcalc) :: vpftcc               !< corr calculator
type(opt_spec)        :: ospec_eul            !< optimizer specification object, Euler angles
type(opt_spec)        :: ospec_shift          !< optimizer specification object, rotational origin shifts
type(opt_spec)        :: ospec_all            !< optimizer specification object, all df:s
type(simplex_opt)     :: nlopt_eul            !< optimizer object, Euler angles
type(simplex_opt)     :: nlopt_shift          !< optimizer object, rotational origin shifts
type(simplex_opt)     :: nlopt_all            !< optimizer object, all df:s
type(ori)             :: e_glob               !< ori of best solution found so far
real                  :: shift_scale =  1.0   !< shift scale factor
logical               :: shbarr      = .true. !< shift barrier constraint or not
integer               :: nrestarts   = 3      !< simplex restarts (randomized bounds)
logical               :: serial      = .true. !< controls shared-mem parallellization of corr calc

contains

    subroutine volpft_srch_init( vol_ref, vol_target, hp, lp, trs, shbarrier, nrestarts_in )
        use simple_projector, only: projector
        class(projector),           intent(in) :: vol_ref, vol_target
        real,                       intent(in) :: hp, lp, trs
        character(len=*), optional, intent(in) :: shbarrier
        integer,          optional, intent(in) :: nrestarts_in
        real :: lims_eul(3,2), lims_shift(3,2), lims_all(6,2)
        ! create the correlator
        call vpftcc%new( vol_ref, vol_target, hp, lp )
        ! flag the barrier constraint
        shbarr = .true.
        if( present(shbarrier) )then
            if( shbarrier .eq. 'no' ) shbarr = .false.
        endif
        ! set nrestarts
        nrestarts = 3
        if( present(nrestarts_in) ) nrestarts = nrestarts_in
        ! make optimizer specs
        lims_eul        = 0.
        lims_eul(1,2)   = 359.99
        lims_eul(2,2)   = 180.
        lims_eul(3,2)   = 359.99
        lims_shift(:,1) = -trs
        lims_shift(:,2) =  trs
        lims_all(:3,:)  = lims_eul
        lims_all(4:,:)  = lims_shift
        call ospec_eul%specify('simplex', 3, ftol=1e-4,&
        &gtol=1e-4, limits=lims_eul, nrestarts=nrestarts, maxits=30)
        call ospec_shift%specify('simplex', 3, ftol=1e-4,&
        &gtol=1e-4, limits=lims_shift, nrestarts=nrestarts, maxits=30)
        call ospec_all%specify('simplex', 6, ftol=1e-4,&
        &gtol=1e-4, limits=lims_all, nrestarts=nrestarts, maxits=30)
        ! point to costfuns
        call ospec_eul%set_costfun(volpft_srch_costfun_eul)
        call ospec_shift%set_costfun(volpft_srch_costfun_shift)
        call ospec_all%set_costfun(volpft_srch_costfun_all)
        ! create simplex optimiser objects
        call nlopt_eul%new(ospec_eul)
        call nlopt_shift%new(ospec_shift)
        call nlopt_all%new(ospec_all)
        ! create global ori
        call e_glob%new()
    end subroutine volpft_srch_init

    function volpft_srch_minimize_eul( fromto ) result( orientation_best )
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
        ! grid search using the spiral geometry & ANGSTEP degree in-plane resolution
        serial = .true. ! since the below loop is parallel
        corrs  = -1.
        !$omp parallel do schedule(static) default(shared) private(iproj,euls,inpl) proc_bind(close)
        do iproj=ffromto(1),ffromto(2)
            euls = espace%get_euler(iproj)
            do inpl=0,359,ANGSTEP
                euls(3) = real(inpl)
                corrs(iproj,inpl) = vpftcc%corr(euls, serial)
            end do
        end do
        !$omp end parallel do
        ! identify the best candidates (serial code)
        do iproj=ffromto(1),ffromto(2)
            corr_best   = -1.
            config_best = 0
            do inpl=0,359,ANGSTEP
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
                ospec_eul%x = resoris%get_euler(iloc)
                call nlopt_eul%minimize(ospec_eul, cost)
                call resoris%set_euler(iloc, ospec_eul%x)
                call resoris%set(iloc, 'corr', -cost)
            end do
        else
            ! refine the NBEST local optima
            ! order the local optima according to correlation
            order = resoris%order_corr()
            ! refine the NBEST solutions
            do iloc=1,NBEST
                ospec_eul%x = resoris%get_euler(order(iloc))
                call nlopt_eul%minimize(ospec_eul, cost)
                call resoris%set_euler(order(iloc), ospec_eul%x)
                call resoris%set(order(iloc), 'corr', -cost)
            end do
        endif
        ! order the local optima according to correlation
        order = resoris%order_corr()
        ! update global ori
        e_glob = resoris%get_ori(order(1))
        ! return best
        orientation_best = resoris%get_ori(order(1))
    end function volpft_srch_minimize_eul

    function volpft_srch_minimize_shift( ) result( orientation_best )
        type(ori) :: orientation_best
        real      :: cost_init, cost
        orientation_best = e_glob
        ospec_shift%x    = 0.
        serial           = .false. ! since we want to parallellize the cost calc
        cost_init        = volpft_srch_costfun_shift(ospec_shift%x, ospec_shift%ndim)
        call nlopt_shift%minimize(ospec_shift, cost)
        if( cost < cost_init )then
            ! set corr
            call orientation_best%set('corr', -cost)
            ! set shift
            call orientation_best%set('x', ospec_shift%x(1))
            call orientation_best%set('y', ospec_shift%x(2))
            call orientation_best%set('z', ospec_shift%x(3))
        else
            ! set corr
            call orientation_best%set('corr', -cost_init)
            ! set shift
            call orientation_best%set('x', 0.0)
            call orientation_best%set('y', 0.0)
            call orientation_best%set('z', 0.0)
        endif
        ! update global ori
        e_glob = orientation_best  
    end function volpft_srch_minimize_shift

    function volpft_srch_minimize_all( ) result( orientation_best )
        type(ori) :: orientation_best
        real      :: cost_init, cost
        orientation_best = e_glob
        ospec_all%x(1:3) = e_glob%get_euler()
        ospec_all%x(4)   = e_glob%get('x')
        ospec_all%x(5)   = e_glob%get('y')
        ospec_all%x(6)   = e_glob%get('z')
        ! determines & applies shift scaling
        shift_scale = (maxval(ospec_all%limits(1:3,2)) - minval(ospec_all%limits(1:3,1)))/&
                     &(maxval(ospec_all%limits(4:6,2)) - minval(ospec_all%limits(4:6,1)))
        ospec_all%limits(4:,:) = ospec_all%limits(4:,:) * shift_scale
        ospec_all%x(4:)        = ospec_all%x(4:)        * shift_scale
        ! initialize
        serial    = .false. ! since we want to parallellize the cost calc
        cost_init = volpft_srch_costfun_all(ospec_all%x, ospec_all%ndim)
        ! minimisation
        call nlopt_all%minimize(ospec_all, cost)
        ! solution and un-scaling
        if( cost < cost_init )then
            ! set corr
            call orientation_best%set('corr', -cost)
            ! set Euler
            call orientation_best%set_euler(ospec_all%x(1:3))
            ! unscale & set shift
            ospec_all%x(4:6) = ospec_all%x(4:6) / shift_scale
            call orientation_best%set('x', ospec_all%x(4))
            call orientation_best%set('y', ospec_all%x(5))
            call orientation_best%set('z', ospec_all%x(6))
            ! update global ori
            e_glob = orientation_best  
        endif
    end function volpft_srch_minimize_all

    function volpft_srch_costfun_eul( vec, D ) result( cost )
        use simple_ori, only: ori
        integer, intent(in) :: D
        real,    intent(in) :: vec(D)
        real :: cost
        cost = -vpftcc%corr(vec(:3), serial)
    end function volpft_srch_costfun_eul

    function volpft_srch_costfun_shift( vec, D ) result( cost )
        use simple_ori, only: ori
        integer, intent(in) :: D
        real,    intent(in) :: vec(D)
        real :: vec_here(3), cost
        vec_here = vec
        where( abs(vec) < 1.e-6 ) vec_here = 0.0
        if( shbarr )then
            if( any(vec_here(:) < ospec_shift%limits(:,1)) .or.&
               &any(vec_here(:) > ospec_shift%limits(:,2)) )then
                cost = 1.
                return
            endif
        endif
        cost = -vpftcc%corr(e_glob, vec_here(:3), serial)
    end function volpft_srch_costfun_shift

    function volpft_srch_costfun_all( vec, D ) result( cost )
        integer, intent(in) :: D
        real,    intent(in) :: vec(D)
        real :: vec_here(6), cost
        vec_here = vec
        if( shbarr )then
            if( any(vec_here(4:6) < ospec_all%limits(4:6,1)) .or.&
               &any(vec_here(4:6) > ospec_all%limits(4:6,2)) )then
                cost = 1.
                return
            endif
        endif
        ! unscale shift
        vec_here(4:6) = vec_here(4:6) / shift_scale
        ! zero small shifts
        if( abs(vec_here(4)) < 1.e-6 ) vec_here(4) = 0.0
        if( abs(vec_here(5)) < 1.e-6 ) vec_here(5) = 0.0
        if( abs(vec_here(6)) < 1.e-6 ) vec_here(6) = 0.0
        ! cost
        cost = -vpftcc%corr(vec_here(:3), vec_here(4:), serial)
    end function volpft_srch_costfun_all

end module simple_volpft_srch
