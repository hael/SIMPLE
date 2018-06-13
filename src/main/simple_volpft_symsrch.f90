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

logical, parameter :: DEBUG   = .false.
integer, parameter :: NPROJ   = 200
integer, parameter :: NBEST   = 20
integer, parameter :: ANGSTEP = 10

type opt4openMP
    type(opt_spec)            :: ospec              !< optimizer specification object
    class(optimizer), pointer :: nlopt => null()    !< optimizer object
end type opt4openMP

type(opt4openMP), allocatable :: opt_symaxes(:)     !< parallel optimisation, symmetry axes
real,             allocatable :: sym_rmats(:,:,:)   !< symmetry operations rotation matrices
complex,          allocatable :: sym_targets(:,:,:) !< symmetry target distributions of Fourier components
type(volpft_corrcalc)         :: vpftcc             !< corr calculator
type(ori)                     :: saxis_glob         !< best symaxis solution found so far
type(oris)                    :: espace             !< projection directions in discrete search
type(oris)                    :: cand_axes          !< candidate symmetry axes
type(sym)                     :: symobj             !< symmetry object
integer                       :: nrestarts = 3      !< simplex restarts (randomized bounds)
integer                       :: nsym      = 0      !< # symmetry ops
integer                       :: nspace    = 0      !< # complex vectors in vpftcc
integer                       :: kfromto(2)         !< Fourier index range

contains

    subroutine volpft_symsrch_init( vol, pgrp, hp, lp, nrestarts_in )
        use simple_projector,    only: projector
         use simple_opt_factory, only: opt_factory
        class(projector),  intent(in) :: vol
        character(len=*),  intent(in) :: pgrp
        real,              intent(in) :: hp, lp
        integer, optional, intent(in) :: nrestarts_in
        type(opt_factory) :: ofac
        type(ori)         :: o
        integer :: isym, nspace, kfromto(2), ithr
        real    :: lims(3,2)
        call volpft_symsrch_kill
        ! create the correlator
        call vpftcc%new(vol, hp, lp, KBALPHA)
        nspace  = vpftcc%get_nspace()
        kfromto = vpftcc%get_kfromto()
        ! create the symmetry object
        call symobj%new(pgrp)
        nsym = symobj%get_nsym()
        ! allocate
        allocate(sym_rmats(nsym,3,3), source=0.0)
        allocate(sym_targets(nsym,kfromto(1):kfromto(2),nspace), source=cmplx(0.,0.))
        ! extract the rotation matrices for the symops
        do isym=1,nsym
            o = symobj%get_symori(isym)
            sym_rmats(isym,:,:) = o%get_mat()
        end do
        ! set nrestarts
        nrestarts = 3
        if( present(nrestarts_in) ) nrestarts = nrestarts_in
        ! make optimizer specs
        lims      = 0.
        lims(1,2) = 359.99
        lims(2,2) = 180.
        lims(3,2) = 359.99
        ! make parallel optimiser struct
        allocate(opt_symaxes(nthr_glob))
        do ithr=1,nthr_glob
            ! optimiser spec
            call opt_symaxes(ithr)%ospec%specify('simplex', 3, ftol=1e-4,&
            &gtol=1e-4, limits=lims, nrestarts=nrestarts, maxits=30)
            ! point to costfun
            call opt_symaxes(ithr)%ospec%set_costfun(volpft_symsrch_costfun)
            ! generate optimizer object with the factory
            call ofac%new(opt_symaxes(ithr)%ospec, opt_symaxes(ithr)%nlopt)
        end do
        ! create global symaxis
        call saxis_glob%new
        ! container for candidate symmetry axes
        call cand_axes%new(NPROJ)
        ! projection directions in discrete search
        call espace%new(NPROJ)
        call espace%spiral
    end subroutine volpft_symsrch_init

    subroutine volpft_srch4symaxis( symaxis_best, fromto )
        class(ori),        intent(out) :: symaxis_best
        integer, optional, intent(in)  :: fromto(2)
        real,    allocatable :: corrs(:,:)
        integer, allocatable :: order(:)
        class(*), pointer    :: fun_self => null()
        type(ori) :: symaxis
        integer   :: ffromto(2), ntot, inpl, iproj, iproj_best
        integer   :: inpl_best, istart, istop, ithr, iloc
        real      :: eul(3), corr_best, cost
        logical   :: distr_exec
        ! flag distributed/shmem exec & set range
        distr_exec = present(fromto)
        if( distr_exec )then
            ffromto = fromto
        else
            ffromto(1) = 1
            ffromto(2) = NPROJ
        endif
        if( ffromto(1) < 1 .or. ffromto(2) > NPROJ )then
           stop 'range out of bound; simple_volpft_srch :: volpft_symsrch_gridsrch'
        endif
        ntot = ffromto(2) - ffromto(1) + 1
        allocate(corrs(ffromto(1):ffromto(2),0:359))
        ! grid search using the spiral geometry & ANGSTEP degree in-plane resolution
        corrs  = -1.
        !$omp parallel do schedule(static) default(shared) private(iproj,eul,inpl)&
        !$omp proc_bind(close)
        do iproj=ffromto(1),ffromto(2)
            eul = espace%get_euler(iproj)
            do inpl=0,359,ANGSTEP
                eul(3) = real(inpl)
                corrs(iproj,inpl) = volpft_symsrch_scorefun( eul )
            end do
        end do
        !$omp end parallel do
        ! identify the best candidates (serial code)
        call symaxis%new
        do iproj=ffromto(1),ffromto(2)
            corr_best  = -1.0
            iproj_best = 0
            inpl_best  = 0
            do inpl=0,359,ANGSTEP
                if( corrs(iproj,inpl) > corr_best )then
                    corr_best  = corrs(iproj,inpl)
                    iproj_best = iproj
                    inpl_best  = inpl
                endif
            end do
            symaxis = espace%get_ori(iproj_best)
            call symaxis%e3set(real(inpl_best))
            call symaxis%set('corr', corr_best)
            call cand_axes%set_ori(iproj,symaxis)
        end do
        call symaxis%kill
        ! refine local optima
        ! order the local optima according to correlation
        order = cand_axes%order_corr()
        ! determine range
        istart = 1
        istop  = min(ffromto(2) - ffromto(1) + 1,NBEST)
        !$omp parallel do schedule(static) default(shared) private(iloc,ithr) proc_bind(close)
        do iloc=istart,istop
            ithr = omp_get_thread_num() + 1
            opt_symaxes(ithr)%ospec%x = cand_axes%get_euler(order(iloc))
            call opt_symaxes(ithr)%nlopt%minimize(opt_symaxes(ithr)%ospec, fun_self, cost)
            call cand_axes%set_euler(order(iloc), opt_symaxes(ithr)%ospec%x)
            call cand_axes%set(order(iloc), 'corr', -cost)
        end do
        !$omp end parallel do
        ! order the local optima according to correlation
        order = cand_axes%order_corr()
        ! update global symaxis
        saxis_glob = cand_axes%get_ori(order(1))
        ! return best
        symaxis_best = cand_axes%get_ori(order(1))
    end subroutine volpft_srch4symaxis

    function volpft_symsrch_costfun( fun_self, vec, D ) result( cost )
        class(*), intent(inout) :: fun_self
        integer,  intent(in)    :: D
        real,     intent(in)    :: vec(D)
        real :: cost
        cost = -volpft_symsrch_scorefun(vec(1:3))
    end function volpft_symsrch_costfun

    function volpft_symsrch_scorefun( eul ) result( cc )
        real, intent(in) :: eul(3)
        real    :: cc, rmat_symaxis(3,3), rmat(3,3), sqsum_sum, sqsum_target
        complex :: sum_targets(kfromto(1):kfromto(2),nspace)
        integer :: isym
        rmat_symaxis = euler2m(eul)
        sum_targets  = cmplx(0.,0.)
        do isym=1,nsym
            ! extracts Fourier component distribution @ symaxis @ symop isym
            rmat = matmul(rmat_symaxis,sym_rmats(isym,:,:))
            call vpftcc%extract_target(rmat)
            ! stores isym distribution and updates sum over symops
            call vpftcc%get_target(sym_targets(isym,:,:))
            sum_targets = sum_targets + sym_targets(isym,:,:)
        end do
        ! correlate sum over symmetry related distributions with
        ! the individual distributions to score the symmetry axis
        sqsum_sum = sum(csq(sum_targets))
        cc = 0.
        do isym=1,nsym
            sqsum_target = sum(csq(sym_targets(isym,:,:)))
            cc = cc + sum(real(sum_targets * conjg(sym_targets(isym,:,:)))) / sqrt(sqsum_sum * sqsum_target)
        end do
        cc = cc / real(nsym)
    end function volpft_symsrch_scorefun

    subroutine volpft_symsrch_kill
        integer :: ithr
        if( allocated(sym_rmats)   ) deallocate(sym_rmats)
        if( allocated(sym_targets) ) deallocate(sym_targets)
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
        call symobj%kill
        call cand_axes%kill
        call espace%kill
    end subroutine volpft_symsrch_kill

end module simple_volpft_symsrch
