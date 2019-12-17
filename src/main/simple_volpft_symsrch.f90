module simple_volpft_symsrch
include 'simple_lib.f08'
!$ use omp_lib
!$ use omp_lib_kinds
use simple_opt_spec,        only: opt_spec
use simple_volpft_corrcalc, only: volpft_corrcalc
use simple_optimizer,       only: optimizer
use simple_oris,            only: oris
use simple_ori,             only: ori, euler2m
use simple_sym,             only: sym
implicit none

public :: volpft_symsrch_init, volpft_srch4symaxis
private
#include "simple_local_flags.inc"

integer, parameter :: NPROJ = 6000
integer, parameter :: NBEST = 20
logical, parameter :: SCOREFUN_JACOB = .false.

type opt4openMP
    type(opt_spec)            :: ospec              !< optimizer specification object
    class(optimizer), pointer :: nlopt => null()    !< optimizer object
end type opt4openMP

type(opt4openMP), allocatable :: opt_symaxes(:)     !< parallel optimisation, symmetry axes
real,             allocatable :: sym_rmats(:,:,:)   !< symmetry operations rotation matrices
type(volpft_corrcalc)         :: vpftcc             !< corr calculator
type(ori)                     :: saxis_glob         !< best symaxis solution found so far
integer                       :: nrestarts     = 3  !< simplex restarts (randomized bounds)
integer                       :: nsym          = 0  !< # symmetry ops
integer                       :: nspace        = 0  !< # complex vectors in vpftcc
integer                       :: nspace_nonred = 0  !< # nonredundant complex vectors in vpftcc
integer                       :: kfromto(2)         !< Fourier index range

contains

    subroutine volpft_symsrch_init( vol, pgrp, hp, lp, nrestarts_in )
        use simple_projector,   only: projector
        use simple_opt_factory, only: opt_factory
        class(projector),  intent(in) :: vol
        character(len=*),  intent(in) :: pgrp
        real,              intent(in) :: hp, lp
        integer, optional, intent(in) :: nrestarts_in
        type(opt_factory) :: ofac
        type(ori)         :: o
        type(sym)         :: symobj
        integer :: isym, ithr
        real    :: lims(3,2)
        call volpft_symsrch_kill
        ! create the correlator
        call vpftcc%new(vol, hp, lp, KBALPHA)
        nspace        = vpftcc%get_nspace()
        nspace_nonred = vpftcc%get_nspace_nonred()
        kfromto       = vpftcc%get_kfromto()
        ! create the symmetry object
        call symobj%new(pgrp, icorelion=.true.)
        nsym = symobj%get_nsym()
        ! allocate
        allocate(sym_rmats(nsym,3,3), source=0.0)
        ! extract the rotation matrices for the symops
        do isym=1,nsym
            call symobj%get_symori(isym, o)
            sym_rmats(isym,:,:) = o%get_mat()
        end do
        ! set nrestarts
        nrestarts = 3
        if( present(nrestarts_in) ) nrestarts = nrestarts_in
        ! make optimizer specs
        lims      = 0.
        lims(1,2) = 359.99
        lims(2,2) = 90. ! only northern hemisphere
        lims(3,2) = 359.99
        ! make parallel optimiser struct
        allocate(opt_symaxes(nthr_glob))
        do ithr=1,nthr_glob
            ! optimiser spec
            call opt_symaxes(ithr)%ospec%specify("simplex",ndim=3,limits=lims,nrestarts=nrestarts,maxits=30)
            ! point to costfun
            call opt_symaxes(ithr)%ospec%set_costfun(volpft_symsrch_costfun)
            ! generate optimizer object with the factory
            call ofac%new(opt_symaxes(ithr)%ospec, opt_symaxes(ithr)%nlopt)
        end do
        ! create global symaxis
        call saxis_glob%new
        call o%kill
        call symobj%kill
    end subroutine volpft_symsrch_init

    subroutine volpft_srch4symaxis( symaxis_best )
        class(ori), intent(inout) :: symaxis_best
        real,    allocatable :: rmats(:,:,:), corrs(:)
        integer, allocatable :: order(:)
        class(*), pointer    :: fun_self => null()
        type(ori)  :: symaxis
        type(oris) :: espace, cand_axes
        integer    :: fromto(2), ntot, iproj, iproj_best
        integer    :: istop, ithr, iloc,iinpl_idx
        integer, parameter :: Ninpl_idx = 10
        integer    :: Ncorrs, idx, iii
        real       :: eul(3), corr_best, cost, lims(3,2), best_cost, best_x(3), dInpl_deg, inpl_corr
        integer, allocatable :: coarse_srch_best_inpl(:)
        ! For inplane angle
        dInpl_deg = 360./10./real(Ninpl_idx)
        ! set range
        fromto(1) = 1
        fromto(2) = NPROJ
        ! container for candidate symmetry axes
        call cand_axes%new(NPROJ)
        call cand_axes%set_all2single('corr', -1.0) ! for later ordering
        ! projection directions in discrete search
        call espace%new(NPROJ)
        call espace%spiral
        ! only consider the northern hemisphere
        do iproj = 1, NPROJ
            if (espace%e2get(iproj) < 90.) then
                fromto(2) = iproj
                exit
            end if
        end do
        ntot   = fromto(2) - fromto(1) + 1
        Ncorrs = ntot*Ninpl_idx
        allocate(coarse_srch_best_inpl(fromto(1):fromto(2)))
        coarse_srch_best_inpl = -1
        ! allocate
        allocate(rmats(Ncorrs,3,3),corrs(Ncorrs))
        ! fill-in the rotation matrices
        do iproj=fromto(1),fromto(2)
            do iinpl_idx=1,Ninpl_idx
                eul    = espace%get_euler(iproj)
                eul(3) = (iinpl_idx-1)*dInpl_deg
                idx    = (iproj-1)*Ninpl_idx+iinpl_idx
                rmats(iproj,:,:) = euler2m(eul)
            enddo
        end do
        ! grid search using the spiral geometry
        corrs = -1.
        !$omp parallel do schedule(static) default(shared) private(iproj,idx) collapse(2) proc_bind(close)
        do iproj=fromto(1),fromto(2)
            do iinpl_idx=1,Ninpl_idx
                idx     = (iproj-1)*Ninpl_idx+iinpl_idx
                corrs(iproj) = volpft_symsrch_scorefun(rmats(idx,:,:))
            enddo
        end do
        !$omp end parallel do
        ! identify the best candidates (serial code)
        corr_best  = -1.0
        iproj_best = 0
        do iproj=fromto(1),fromto(2)
            inpl_corr = -HUGE(inpl_corr)
            do iinpl_idx = 1, Ninpl_idx
                idx = (iproj-1)*Ninpl_idx+iinpl_idx
                if( corrs(idx)  > inpl_corr ) then
                    ! for each solid angle, keep just the best inplane angle
                    coarse_srch_best_inpl(iproj) = iinpl_idx
                    inpl_corr = corrs(idx)
                endif
            enddo
            call espace%get_ori(iproj, symaxis)
            call symaxis%set('corr', inpl_corr)
            call cand_axes%set_ori(iproj,symaxis)
        end do
        ! refine local optima
        ! order the local optima according to correlation
        order = cand_axes%order_corr()
        ! determine end of range
        istop = min(fromto(2) - fromto(1) + 1,NBEST)
        nrestarts=10
        !$omp parallel do schedule(static) default(shared) private(best_x,iinpl_idx,iii,best_cost,iloc,ithr,cost,lims) proc_bind(close)
        do iloc=1,istop
            if( iloc == 1 ) cycle
            ithr = omp_get_thread_num() + 1
            best_cost = 10.
            opt_symaxes(ithr)%ospec%x = cand_axes%get_euler(order(iloc))
            opt_symaxes(ithr)%ospec%x(3) = (coarse_srch_best_inpl(order(iloc)) - 1) * dInpl_deg
            lims(1,:) = [opt_symaxes(ithr)%ospec%x(1)-10.,opt_symaxes(ithr)%ospec%x(1)+10.]
            lims(2,:) = [opt_symaxes(ithr)%ospec%x(2)-10.,opt_symaxes(ithr)%ospec%x(2)+10.]
            lims(3,:) = [opt_symaxes(ithr)%ospec%x(3)-dInpl_deg/3.,opt_symaxes(ithr)%ospec%x(3)+dInpl_deg/3.]
            ! optimiser spec
            call opt_symaxes(ithr)%ospec%specify("simplex",ndim=3,limits=lims,nrestarts=nrestarts,ftol=1e-7,maxits=200)
            ! point to costfun
            call opt_symaxes(ithr)%ospec%set_costfun(volpft_symsrch_costfun)
            opt_symaxes(ithr)%ospec%niter = 0
            call opt_symaxes(ithr)%nlopt%minimize(opt_symaxes(ithr)%ospec, fun_self, cost)
            if( cost < best_cost )then
                best_cost = cost
                best_x =  opt_symaxes(ithr)%ospec%x
            end if
            !end do
            call cand_axes%set_euler(order(iloc), best_x)!opt_symaxes(ithr)%ospec%x)
            call cand_axes%set(order(iloc), 'corr', -best_cost)
        end do
        !$omp end parallel do
        ! order the local optima according to correlation
        order = cand_axes%order_corr()
        ! update global symaxis
        call cand_axes%get_ori(order(1), saxis_glob)
        ! return best
        call cand_axes%get_ori(order(1), symaxis_best)
        call symaxis%kill
        call espace%kill
        call cand_axes%kill
    end subroutine volpft_srch4symaxis

    function volpft_symsrch_costfun( fun_self, vec, D ) result( cost )
        class(*), intent(inout) :: fun_self
        integer,  intent(in)    :: D
        real,     intent(in)    :: vec(D)
        real :: cost, rmat(3,3)
        rmat = euler2m(vec(1:3))
        cost = -volpft_symsrch_scorefun(rmat)
    end function volpft_symsrch_costfun

    function volpft_symsrch_scorefun( rmat_symaxis ) result( cc )
        real, intent(in) :: rmat_symaxis(3,3)
        real    :: cc, rmat(3,3)
        complex :: sym_targets(nsym,kfromto(1):kfromto(2),nspace_nonred)
        complex :: sum_of_sym_targets(kfromto(1):kfromto(2),nspace_nonred)
        real    :: sqsum_targets(nsym), sqsum_sum
        integer :: isym, k
        sum_of_sym_targets = cmplx(0.,0.)
        do isym=1,nsym
            ! extracts Fourier component distribution @ symaxis @ symop isym
            ! ROTATION MATRICES DO NOT COMMUTE
            ! this is the correct order
            rmat = matmul(sym_rmats(isym,:,:), rmat_symaxis)
            call vpftcc%extract_target(rmat, sym_targets(isym,:,:), sqsum_targets(isym))
            if( SCOREFUN_JACOB )then
                do k = kfromto(1),kfromto(2)
                    sym_targets(isym,k,:) = sym_targets(isym,k,:) * real(k)
                end do
                sqsum_targets(isym) = sum(csq(sym_targets(isym,:,:)))
            end if
            sum_of_sym_targets = sum_of_sym_targets + sym_targets(isym,:,:)
        end do
        ! correlate with the average to score the symmetry axis
        sqsum_sum = sum(csq(sum_of_sym_targets))
        cc = 0.
        do isym=1,nsym
            cc = cc + sum(real(sum_of_sym_targets * conjg(sym_targets(isym,:,:))))&
                &/ sqrt(sqsum_sum * sqsum_targets(isym))
        end do
        cc = cc / real(nsym)
    end function volpft_symsrch_scorefun

    subroutine volpft_symsrch_kill
        integer :: ithr
        if( allocated(sym_rmats)   ) deallocate(sym_rmats)
        if( allocated(opt_symaxes) )then
            do ithr=1,size(opt_symaxes)
                call opt_symaxes(ithr)%ospec%kill
                if( associated(opt_symaxes(ithr)%nlopt) )then
                    call opt_symaxes(ithr)%nlopt%kill
                    nullify(opt_symaxes(ithr)%nlopt)
                endif
            end do
            deallocate(opt_symaxes)
        endif
        call vpftcc%kill
        call saxis_glob%kill
    end subroutine volpft_symsrch_kill

end module simple_volpft_symsrch
