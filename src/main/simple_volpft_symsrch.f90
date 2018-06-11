module simple_volpft_symsrch
use simple_opt_spec,        only: opt_spec
use simple_volpft_corrcalc, only: volpft_corrcalc
use simple_opt_simplex,     only: opt_simplex
use simple_oris,            only: oris
use simple_ori,             only: ori, euler2m
use simple_sym,             only: sym
use simple_defs             ! use all in there
implicit none

public :: volpft_symsrch_init, volpft_symsrch_gridsrch
private

integer, parameter :: NPROJ   = 200
integer, parameter :: NBEST   = 20
integer, parameter :: ANGSTEP = 10

type(volpft_corrcalc) :: vpftcc             !< corr calculator
type(opt_spec)        :: ospec              !< optimizer specification object
type(opt_simplex)     :: nlopt              !< optimizer object
type(ori)             :: saxis_glob         !< best syumaxis solution found so far
type(sym)             :: symobj             !< symmetry object
real, allocatable     :: sym_rmats(:,:,:)   !< symmetry operations rotation matrices
real, allocatable     :: corrs(:,:)         !< symmetry axis correlations
integer               :: nrestarts = 3      !< simplex restarts (randomized bounds)
integer               :: nsym      = 0      !< # symmetry ops
logical               :: serial    = .true. !< controls shared-mem parallellization of corr calc

contains

    subroutine volpft_symsrch_init( vol, pgrp, hp, lp, nrestarts_in )
        use simple_projector, only: projector
        class(projector),  intent(in) :: vol
        character(len=*),  intent(in) :: pgrp
        real,              intent(in) :: hp, lp
        integer, optional, intent(in) :: nrestarts_in
        integer   :: isym
        type(ori) :: o
        real      :: lims(3,2)
        ! create the symmetry object
        call symobj%new(pgrp)
        nsym = symobj%get_nsym()
        ! allocate
        if( allocated(sym_rmats) ) deallocate(sym_rmats)
        if( allocated(corrs)     ) deallocate(corrs)
        allocate(sym_rmats(nsym,3,3), corrs(NPROJ,0:359), source=0.0)
        ! extract the rotation matrices for the symops
        do isym=1,nsym
            o = symobj%get_symori(isym)
            sym_rmats(isym,:,:) = o%get_mat()
        end do
        ! create the correlator
        call vpftcc%new( vol, vol, hp, lp, KBALPHA )
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
        real       :: eul(3), rmat_symaxis(3,3), rmat(3,3), corr_best
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
        serial = .true. ! since the below loop is parallel
        corrs  = -1.
        call vpftcc%extract_ref
        !$omp parallel do schedule(static) default(shared) private(iproj,eul,inpl,rmat_symaxis,isym,rmat)&
        !$omp proc_bind(close)
        do iproj=ffromto(1),ffromto(2)
            eul = espace%get_euler(iproj)
            do inpl=0,359,ANGSTEP
                eul(3) = real(inpl)
                rmat_symaxis = euler2m(eul)
                corrs(iproj,inpl) = 0.
                do isym=1,nsym
                    rmat = matmul(rmat_symaxis,sym_rmats(isym,:,:))
                    call vpftcc%extract_target(rmat, serial)
                    corrs(iproj,inpl) = corrs(iproj,inpl) + vpftcc%corr()
                end do
                corrs(iproj,inpl) = corrs(iproj,inpl) / real(nsym)
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

    ! cost fun, assumes reference extracted (call vpftcc%extract_ref)
    function volpft_symsrch_costfun( fun_self, vec, D ) result( cost )
        class(*), intent(inout) :: fun_self
        integer,  intent(in)    :: D
        real,     intent(in)    :: vec(D)
        real    :: cost, rmat_symaxis(3,3), rmat(3,3), cc
        integer :: isym
        rmat_symaxis = euler2m(vec(1:3))
        cc = 0.
        do isym=1,nsym
            rmat = matmul(rmat_symaxis,sym_rmats(isym,:,:))
            call vpftcc%extract_target(rmat, .true.)
            cc = cc + vpftcc%corr()
        end do
        cc = cc / real(nsym)
    end function volpft_symsrch_costfun

end module simple_volpft_symsrch
