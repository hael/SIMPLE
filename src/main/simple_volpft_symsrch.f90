module simple_volpft_symsrch
!$ use omp_lib
!$ use omp_lib_kinds
include 'simple_lib.f08'
use simple_opt_spec,        only: opt_spec
use simple_volpft_corrcalc, only: volpft_corrcalc
use simple_opt_simplex,     only: opt_simplex
use simple_oris,            only: oris
use simple_ori,             only: ori, euler2m
use simple_sym,             only: sym
implicit none

public :: volpft_symsrch_init, volpft_symsrch_gridsrch
private

integer, parameter :: NPROJ   = 200
integer, parameter :: NBEST   = 20
integer, parameter :: ANGSTEP = 10

type(volpft_corrcalc) :: vpftcc             !< corr calculator
type(opt_spec)        :: ospec              !< optimizer specification object
type(opt_simplex)     :: nlopt              !< optimizer object
type(ori)             :: saxis_glob         !< best symaxis solution found so far
type(sym)             :: symobj             !< symmetry object
real,    allocatable  :: sym_rmats(:,:,:)   !< symmetry operations rotation matrices
real,    allocatable  :: corrs(:,:)         !< symmetry axis correlations
complex, allocatable  :: sym_targets(:,:,:) !< symmetry target distributions of Fourier components
integer               :: nrestarts = 3      !< simplex restarts (randomized bounds)
integer               :: nsym      = 0      !< # symmetry ops
integer               :: nspace    = 0      !< # complex vectors in vpftcc
integer               :: kfromto(2)         !< Fourier index range

contains

    subroutine volpft_symsrch_init( vol, pgrp, hp, lp, nrestarts_in )
        use simple_projector, only: projector
        class(projector),  intent(in) :: vol
        character(len=*),  intent(in) :: pgrp
        real,              intent(in) :: hp, lp
        integer, optional, intent(in) :: nrestarts_in
        integer   :: isym, nspace, kfromto(2)
        type(ori) :: o
        real      :: lims(3,2)
        ! create the correlator
        call vpftcc%new(vol, vol, hp, lp, KBALPHA)
        nspace  = vpftcc%get_nspace()
        kfromto = vpftcc%get_kfromto()
        ! create the symmetry object
        call symobj%new(pgrp)
        nsym = symobj%get_nsym()
        ! allocate
        if( allocated(sym_rmats)   ) deallocate(sym_rmats)
        if( allocated(corrs)       ) deallocate(corrs)
        if( allocated(sym_targets) ) deallocate(sym_targets)
        allocate(sym_rmats(nsym,3,3), corrs(NPROJ,0:359), source=0.0)
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
        call ospec%specify('simplex', 3, ftol=1e-4,&
        &gtol=1e-4, limits=lims, nrestarts=nrestarts, maxits=30)
        call ospec%set_costfun(volpft_symsrch_costfun)
        call saxis_glob%new
    end subroutine volpft_symsrch_init

    subroutine volpft_symsrch_gridsrch( symaxis_best, fromto )
        class(ori),        intent(out) :: symaxis_best
        integer, optional, intent(in)  :: fromto(2)
        integer    :: ffromto(2), ntot, inpl, iproj, isym, loc(2)
        integer    :: iproj_best, inpl_best
        real       :: eul(3), corr_best
        logical    :: distr_exec
        type(oris) :: espace, cand_axes
        type(ori)  :: symaxis
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
        ! create
        call symaxis%new
        call cand_axes%new(NPROJ)
        call espace%new(NPROJ)
        call espace%spiral
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
    end subroutine volpft_symsrch_gridsrch

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

end module simple_volpft_symsrch
