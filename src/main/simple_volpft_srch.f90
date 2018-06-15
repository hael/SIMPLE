! volume alignment based on band-pass limited cross-correlation
! the all df:s & shift routines need further testing
! need to introduce openMP:ed cost function in volpft_corrcalc for the serial optimisers
module simple_volpft_srch
include 'simple_lib.f08'
!$ use omp_lib
!$ use omp_lib_kinds
use simple_opt_spec,        only: opt_spec
use simple_volpft_corrcalc, only: volpft_corrcalc
use simple_optimizer,       only: optimizer
use simple_oris,            only: oris
use simple_ori,             only: ori, euler2m
implicit none

public :: volpft_srch_init, volpft_srch_minimize_eul, volpft_srch_minimize_shift, volpft_srch_minimize_all
private

logical, parameter :: DEBUG   = .false.
integer, parameter :: NPROJ   = 200
integer, parameter :: NBEST   = 20
integer, parameter :: ANGSTEP = 10

type opt4openMP
    type(opt_spec)            :: ospec                 !< optimizer specification object
    class(optimizer), pointer :: nlopt => null()       !< optimizer object
end type opt4openMP

type(volpft_corrcalc)         :: vpftcc                !< corr calculator
type(opt4openMP), allocatable :: opt_eul(:)            !< parallel optimisation, Euler angles
type(opt_spec)                :: ospec_shift           !< optimizer specification object, rotational origin shifts
type(opt_spec)                :: ospec_all             !< optimizer specification object, all df:s
class(optimizer), pointer     :: nlopt_shift => null() !< optimizer object, rotational origin shifts
class(optimizer), pointer     :: nlopt_all   => null() !< optimizer object, all df:s
type(ori)                     :: e_glob                !< ori of best solution found so far
real                          :: shift_scale =  1.0    !< shift scale factor
logical                       :: shbarr      = .true.  !< shift barrier constraint or not
integer                       :: nrestarts   = 3       !< simplex restarts (randomized bounds)

contains

    subroutine volpft_srch_init( vol_ref, vol_target, hp, lp, trs, shbarrier, nrestarts_in )
        use simple_projector,   only: projector
        use simple_opt_factory, only: opt_factory
        class(projector),           intent(in) :: vol_ref, vol_target
        real,                       intent(in) :: hp, lp, trs
        character(len=*), optional, intent(in) :: shbarrier
        integer,          optional, intent(in) :: nrestarts_in
        type(opt_factory) :: ofac
        real              :: lims_eul(3,2), lims_shift(3,2), lims_all(6,2)
        integer           :: ithr
        call volpft_srch_kill
        ! create the correlator
        call vpftcc%new( vol_ref, vol_target, hp, lp, KBALPHA )
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
        ! make parallel optimiser struct
        allocate(opt_eul(nthr_glob))
        do ithr=1,nthr_glob
            ! optimiser spec
            call opt_eul(ithr)%ospec%specify('simplex', 3, ftol=1e-4,&
            &gtol=1e-4, limits=lims_eul, nrestarts=nrestarts, maxits=30)
            ! point to costfun
            call opt_eul(ithr)%ospec%set_costfun(volpft_srch_costfun_eul)
            ! generate optimizer object with the factory
            call ofac%new(opt_eul(ithr)%ospec, opt_eul(ithr)%nlopt)
        end do
        ! serial optimizers
        call ospec_shift%specify('simplex', 3, ftol=1e-4,&
        &gtol=1e-4, limits=lims_shift, nrestarts=nrestarts, maxits=30)
        call ospec_all%specify('simplex', 6, ftol=1e-4,&
        &gtol=1e-4, limits=lims_all, nrestarts=nrestarts, maxits=30)
        call ospec_shift%set_costfun(volpft_srch_costfun_shift)
        call ospec_all%set_costfun(volpft_srch_costfun_all)
        call ofac%new(ospec_shift, nlopt_shift)
        call ofac%new(ospec_all, nlopt_all)
        ! create global ori
        call e_glob%new()
        if( DEBUG ) write(*,*) 'debug(volpft_srch); volpft_srch_init, DONE'
    end subroutine volpft_srch_init

    function volpft_srch_minimize_eul( fromto ) result( orientation_best )
        integer, optional, intent(in) :: fromto(2)
        real,    allocatable :: inpl_angs(:), rmats(:,:,:,:), corrs(:,:)
        integer, allocatable :: order(:)
        class(*), pointer    :: fun_self => null()
        type(oris) :: espace, cand_oris
        type(ori)  :: orientation, orientation_best
        integer    :: ffromto(2), ntot, iproj, iproj_best, inpl_best
        integer    :: iloc, inpl, ithr, istart, istop, n_inpls
        real       :: eul(3), corr_best, cost
        logical    :: distr_exec
        ! flag distributed/shmem exec & set range
        distr_exec = present(fromto)
        if( distr_exec )then
            ffromto = fromto
        else
            ffromto(1) = 1
            ffromto(2) = NPROJ
        endif
        if( ffromto(1) < 1 .or. ffromto(2) > NPROJ )then
            stop 'range out of bound; simple_volpft_srch :: volpft_srch_minimize_eul'
        endif
        ntot = ffromto(2) - ffromto(1) + 1
        ! create
        call orientation%new
        call cand_oris%new(NPROJ)
        call cand_oris%set_all2single('corr', -1.0) ! for later ordering
        call espace%new(NPROJ)
        call espace%spiral
        ! count # in-plane angles
        n_inpls = 0
        do inpl=0,359,ANGSTEP
            n_inpls = n_inpls + 1
        end do
        ! allocate
        allocate(inpl_angs(n_inpls), rmats(ffromto(1):ffromto(2),n_inpls,3,3),&
            &corrs(ffromto(1):ffromto(2),n_inpls))
        ! fill-in the in-plane angles
        n_inpls = 0
        do inpl=0,359,ANGSTEP
            n_inpls = n_inpls + 1
            inpl_angs(n_inpls) = real(inpl)
        end do
        ! fill-in the rotation matrices
        do iproj=ffromto(1),ffromto(2)
            eul = espace%get_euler(iproj)
            do inpl=1,n_inpls
                eul(3) = inpl_angs(inpl)
                rmats(iproj,inpl,:,:) = euler2m(eul)
            end do
        end do
        ! grid search using the spiral geometry & ANGSTEP degree in-plane resolution
        corrs  = -1.
        !omp parallel do schedule(static) default(shared) private(iproj,inpl) proc_bind(close) collapse(2)
        do iproj=ffromto(1),ffromto(2)
            do inpl=1,n_inpls
                corrs(iproj,inpl) = vpftcc%corr(rmats(iproj,inpl,:,:))
            end do
        end do
        !omp end parallel do
        ! identify the best candidates (serial code)
        do iproj=ffromto(1),ffromto(2)
            corr_best  = -1.
            iproj_best = 0
            inpl_best  = 0
            do inpl=1,n_inpls
                if( corrs(iproj,inpl) > corr_best )then
                    corr_best  = corrs(iproj,inpl)
                    iproj_best = iproj
                    inpl_best  = inpl
                endif
            end do
            ! set local in-plane optimum for iproj
            orientation = espace%get_ori(iproj_best)
            call orientation%e3set(inpl_angs(inpl_best))
            call orientation%set('corr', corr_best)
            call cand_oris%set_ori(iproj,orientation)
        end do
        ! refine local optima
        ! order the local optima according to correlation
        order = cand_oris%order_corr()
        ! determine range
        istart = 1
        istop  = min(ffromto(2) - ffromto(1) + 1,NBEST)
        !$omp parallel do schedule(static) default(shared) private(iloc,ithr) proc_bind(close)
        do iloc=istart,istop
            ithr = omp_get_thread_num() + 1
            opt_eul(ithr)%ospec%x = cand_oris%get_euler(order(iloc))
            call opt_eul(ithr)%nlopt%minimize(opt_eul(ithr)%ospec, fun_self, cost)
            call cand_oris%set_euler(order(iloc), opt_eul(ithr)%ospec%x)
            call cand_oris%set(order(iloc), 'corr', -cost)
        end do
        !$omp end parallel do
        ! order the local optima according to correlation
        order = cand_oris%order_corr()
        ! update global ori
        e_glob = cand_oris%get_ori(order(1))
        ! return best
        orientation_best = cand_oris%get_ori(order(1))
        ! destruct
        call espace%kill
        call cand_oris%kill
        call orientation%kill
    end function volpft_srch_minimize_eul

    function volpft_srch_minimize_shift( ) result( orientation_best )
        type(ori)         :: orientation_best
        real              :: cost_init, cost
        class(*), pointer :: fun_self => null()
        orientation_best = e_glob
        ospec_shift%x    = 0.
        cost_init        = volpft_srch_costfun_shift(fun_self, ospec_shift%x, ospec_shift%ndim)
        call nlopt_shift%minimize(ospec_shift, fun_self, cost)
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
        type(ori)         :: orientation_best
        real              :: cost_init, cost
        class(*), pointer :: fun_self => null()
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
        cost_init = volpft_srch_costfun_all(fun_self, ospec_all%x, ospec_all%ndim)
        ! minimisation
        call nlopt_all%minimize(ospec_all, fun_self, cost)
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

    function volpft_srch_costfun_eul( fun_self, vec, D ) result( cost )
        use simple_ori, only: ori
        class(*), intent(inout) :: fun_self
        integer,  intent(in)    :: D
        real,     intent(in)    :: vec(D)
        real                    :: cost, rmat(3,3)
        rmat = euler2m(vec(:3))
        cost = -vpftcc%corr(rmat)
    end function volpft_srch_costfun_eul

    function volpft_srch_costfun_shift( fun_self, vec, D ) result( cost )
        use simple_ori, only: ori
        class(*), intent(inout) :: fun_self
        integer,  intent(in)    :: D
        real,     intent(in)    :: vec(D)
        real                    :: vec_here(3), cost
        vec_here = vec
        where( abs(vec) < 1.e-6 ) vec_here = 0.0
        if( shbarr )then
            if( any(vec_here(:) < ospec_shift%limits(:,1)) .or.&
               &any(vec_here(:) > ospec_shift%limits(:,2)) )then
                cost = 1.
                return
            endif
        endif
        cost = -vpftcc%corr(e_glob%get_mat(), vec_here(:3))
    end function volpft_srch_costfun_shift

    function volpft_srch_costfun_all( fun_self, vec, D ) result( cost )
        class(*), intent(inout) :: fun_self
        integer,  intent(in)    :: D
        real,     intent(in)    :: vec(D)
        real                    :: vec_here(6), cost, rmat(3,3)
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
        rmat = euler2m(vec_here(:3))
        cost = -vpftcc%corr(rmat, vec_here(4:))
    end function volpft_srch_costfun_all

    subroutine volpft_srch_kill
        integer :: ithr
        if( allocated(opt_eul) )then
            do ithr=1,size(opt_eul)
                call opt_eul(ithr)%ospec%kill
                if( associated(opt_eul(ithr)%nlopt) )then
                    call opt_eul(ithr)%nlopt%kill
                    nullify(opt_eul(ithr)%nlopt)
                endif
            end do
            deallocate(opt_eul)
        endif
        call vpftcc%kill
        call ospec_shift%kill
        call ospec_all%kill
        if( associated(nlopt_shift) )then
            call nlopt_shift%kill
            nullify(nlopt_shift)
        endif
        if( associated(nlopt_all) )then
            call nlopt_all%kill
            nullify(nlopt_all)
        endif
    end subroutine volpft_srch_kill

end module simple_volpft_srch
