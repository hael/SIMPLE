! for fast rotational docking of volumes
module simple_volpft_srch
include 'simple_lib.f08'
!$ use omp_lib
!$ use omp_lib_kinds
use simple_opt_spec,        only: opt_spec
use simple_volpft_corrcalc, only: volpft_corrcalc
use simple_optimizer,       only: optimizer
implicit none

public :: volpft_srch_init, volpft_srch_set_shvec, volpft_srch_minimize, volpft_srch_refine, volpft_srch_kill
private

logical, parameter :: DEBUG   = .false.
integer, parameter :: NPROJ   = 500
integer, parameter :: NBEST   = 10
integer, parameter :: ANGSTEP = 7

type opt4openMP
    type(opt_spec)            :: ospec             !< optimizer specification object
    class(optimizer), pointer :: opt_obj => null() !< optimizer object
end type opt4openMP

type(volpft_corrcalc)         :: vpftcc                  !< corr calculator
type(opt4openMP), allocatable :: opt_eul(:)              !< parallel optimisation, Euler angles
type(ori)                     :: e_glob                  !< ori of best global solution
real                          :: shvec_glob(3) = 0.      !< global origin shift
real                          :: lims_eul(3,2)           !< opt boundaries
logical                       :: doshift       = .false. !< shifted cross-correlation search
integer                       :: nrestarts     = 5       !< simplex restarts (randomized bounds)

contains

    subroutine volpft_srch_init( vol_ref, vol_target, hp, lp, nrestarts_in )
        use simple_projector,   only: projector
        use simple_opt_factory, only: opt_factory
        class(projector),  intent(in) :: vol_ref, vol_target
        real,              intent(in) :: hp, lp
        integer, optional, intent(in) :: nrestarts_in
        type(opt_factory) :: ofac
        integer           :: ithr
        call volpft_srch_kill
        ! create the correlator
        call vpftcc%new( vol_ref, vol_target, hp, lp, KBALPHA )
        ! set nrestarts
        nrestarts = 5
        if( present(nrestarts_in) ) nrestarts = nrestarts_in
        ! make optimizer specs
        lims_eul      = 0.
        lims_eul(1,2) = 359.99
        lims_eul(2,2) = 180.
        lims_eul(3,2) = 359.99
        ! make parallel optimiser struct
        allocate(opt_eul(nthr_glob))
        do ithr=1,nthr_glob
            ! optimiser spec
            call opt_eul(ithr)%ospec%specify('simplex', 3, ftol=1e-4,&
                &gtol=1e-4, limits=lims_eul, nrestarts=nrestarts, maxits=100)
            ! point to costfun
            call opt_eul(ithr)%ospec%set_costfun(volpft_srch_costfun)
            ! generate optimizer object with the factory
            call ofac%new(opt_eul(ithr)%ospec, opt_eul(ithr)%opt_obj)
        end do
        ! create global ori
        call e_glob%new(is_ptcl=.false.)
        ! unflag doshift
        doshift    = .false.
        shvec_glob = 0.
        if( DEBUG ) write(logfhandle,*) 'debug(volpft_srch); volpft_srch_init, DONE'
    end subroutine volpft_srch_init

    subroutine volpft_srch_set_shvec( shvec )
        real, intent(in) :: shvec(3)
        doshift    = .true.
        shvec_glob = shvec
    end subroutine volpft_srch_set_shvec

    function volpft_srch_minimize() result( orientation_best )
        real,    allocatable :: inpl_angs(:), rmats(:,:,:,:), corrs(:,:)
        integer, allocatable :: order(:)
        class(*),    pointer :: fun_self => null()
        type(ori),    target :: dummyobject ! a dummy object for minimizer, could be any type
        type(oris) :: espace, cand_oris
        type(ori)  :: orientation, orientation_best
        integer    :: iproj, iproj_best, inpl_best
        integer    :: iloc, inpl, ithr, n_inpls
        real       :: eul(3), corr_best, cost, prev_cost
        ! dummy association
        fun_self => dummyobject
        ! create
        call orientation%new(is_ptcl=.false.)
        call cand_oris%new(NPROJ, is_ptcl=.false.)
        call cand_oris%set_all2single('corr', -1.0) ! for later ordering
        call espace%new(NPROJ, is_ptcl=.false.)
        call espace%spiral
        ! count # in-plane angles
        n_inpls = 0
        do inpl=0,359,ANGSTEP
            n_inpls = n_inpls + 1
        end do
        ! allocate
        allocate(inpl_angs(n_inpls), rmats(NPROJ,n_inpls,3,3), corrs(NPROJ,n_inpls))
        ! fill-in the in-plane angles
        n_inpls = 0
        do inpl=0,359,ANGSTEP
            n_inpls = n_inpls + 1
            inpl_angs(n_inpls) = real(inpl)
        end do
        ! fill-in the rotation matrices
        do iproj=1,NPROJ
            eul = espace%get_euler(iproj)
            do inpl=1,n_inpls
                eul(3) = inpl_angs(inpl)
                rmats(iproj,inpl,:,:) = euler2m(eul)
            end do
        end do
        ! grid search using the spiral geometry & ANGSTEP degree in-plane resolution
        corrs  = -1.
        !$omp parallel do schedule(static) default(shared) private(iproj,inpl) proc_bind(close) collapse(2)
        do iproj=1,NPROJ
            do inpl=1,n_inpls
                if( doshift )then
                    corrs(iproj,inpl) = vpftcc%corr(rmats(iproj,inpl,:,:), shvec_glob)
                else
                    corrs(iproj,inpl) = vpftcc%corr(rmats(iproj,inpl,:,:))
                endif
            end do
        end do
        !$omp end parallel do
        ! identify the best candidates (serial code)
        do iproj=1,NPROJ
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
            call espace%get_ori(iproj_best, orientation)
            call orientation%e3set(inpl_angs(inpl_best))
            call orientation%set('corr', corr_best)
            call cand_oris%set_ori(iproj,orientation)
        end do
        ! refine local optima
        ! order the local optima according to correlation
        order = cand_oris%order_corr()
        call cand_oris%get_ori(order(1), orientation_best)
        !$omp parallel do schedule(static) default(shared) private(iloc,ithr,prev_cost,cost) proc_bind(close)
        do iloc=1,NBEST
            ithr = omp_get_thread_num() + 1
            opt_eul(ithr)%ospec%x = cand_oris%get_euler(order(iloc))
            prev_cost = -cand_oris%get(order(iloc),'corr')
            call opt_eul(ithr)%opt_obj%minimize(opt_eul(ithr)%ospec, fun_self, cost)
            if( cost < prev_cost )then
              call cand_oris%set_euler(order(iloc), opt_eul(ithr)%ospec%x)
              call cand_oris%set(order(iloc), 'corr', -cost)
            endif
        end do
        !$omp end parallel do
        ! order the local optima according to correlation
        order = cand_oris%order_corr()
        ! update global ori
        call cand_oris%get_ori(order(1), e_glob)
        ! return best
        call cand_oris%get_ori(order(1), orientation_best)
        ! destruct
        call espace%kill
        call cand_oris%kill
        call orientation%kill
        call orientation_best%kill
    end function volpft_srch_minimize

    function volpft_srch_refine( orientation_start, angres ) result( orientation_best )
        class(ori),     intent(in) :: orientation_start
        real, optional, intent(in) :: angres
        class(*), pointer      :: fun_self => null()
        type(ori) :: orientation_best
        real      :: cost, prev_cost
        integer   :: ithr
        call orientation_best%new(is_ptcl=.false.)
        ithr = omp_get_thread_num() + 1
        opt_eul(ithr)%ospec%x = orientation_start%get_euler()
        if( present(angres) )then
            opt_eul(ithr)%ospec%limits(:,1) = opt_eul(ithr)%ospec%x(:) - angres
            opt_eul(ithr)%ospec%limits(:,2) = opt_eul(ithr)%ospec%x(:) + angres
        else
            opt_eul(ithr)%ospec%limits = lims_eul
        endif
        prev_cost = -vpftcc%corr(euler2m(opt_eul(ithr)%ospec%x(:)), shvec_glob)
        call opt_eul(ithr)%opt_obj%minimize(opt_eul(ithr)%ospec, fun_self, cost)
        if(cost < prev_cost) then
          call orientation_best%set_euler(opt_eul(ithr)%ospec%x)
          call orientation_best%set('corr', -cost)
          call orientation_best%print_ori
        else
          orientation_best = orientation_start
          call orientation_best%set('corr', -prev_cost)
        endif
    end function volpft_srch_refine

    function volpft_srch_costfun( fun_self, vec, D ) result( cost )
        class(*), intent(inout) :: fun_self
        integer,  intent(in)    :: D
        real,     intent(in)    :: vec(D)
        real                    :: cost, rmat(3,3)
        rmat = euler2m(vec(:3))
        if( doshift )then
            cost = -vpftcc%corr(rmat, shvec_glob)
        else
            cost = -vpftcc%corr(rmat)
        endif
    end function volpft_srch_costfun

    subroutine volpft_srch_kill
        integer :: ithr
        if( allocated(opt_eul) )then
            do ithr=1,size(opt_eul)
                call opt_eul(ithr)%ospec%kill
                if( associated(opt_eul(ithr)%opt_obj) )then
                    call opt_eul(ithr)%opt_obj%kill
                    nullify(opt_eul(ithr)%opt_obj)
                endif
            end do
            deallocate(opt_eul)
        endif
        call vpftcc%kill
    end subroutine volpft_srch_kill

end module simple_volpft_srch
