!@descr: flex_pca paired engine FINAL STAGE -- "merge, don't refit" (proposal 1 par.7)
!!
!! The paired fits ARE the production fit; the final stage is a delivery replay, never a fresh
!! EM. Under SIMPLE_COV_PAIRED_MERGE=1 the paired driver stashes each fit's LAST-iteration raw
!! M-step sufficient statistics (numerators Y + packed coupled per-voxel densities rho, pre-ridge
!! pre-solve) plus the entry-frame basis those statistics are expressed in; this submodule then
!!   1. frame-aligns fit B into fit A: R = the orthogonal polar factor of the entry-frame
!!      cross-Gram M_BA(p,q) = <uB_p, uA_q> (align_basis_to_reference + the
!!      Procrustes precedent). polar(M) IS the signed permutation from matching composed with
!!      the in-span rotation: for any signed permutation P, P * polar(P^T M) = polar(M).
!!   2. rotates fit B's statistics into A's frame -- numerators linearly (Y' = Y.R over the
!!      cmat_exp grids), packed densities quadratically (rho' = R^T rho R in pair-index space;
!!      zB = M_BA zA, so E[zA zA'] = R^T E[zB zB'] R) -- sums the four quarter-sets (A-even,
!!      A-odd, B'-even, B'-odd) into a merged even/odd pair, harvests the per-shell sampling H
!!      from each half's pair diagonals and adds the ridge from the honest cross-fit FSC in the
!!      sampling-aware Gilles-Singer S.11 form with the SUMMED H (crossfsc_to_invtau2 +
!!      add_invtausq2rho_coupled -- the refine3D ml_reg add_invtausq2rho precedent: per-half
!!      invtau2 added to each half equals summed-H invtau2 added once to the sum), then runs ONE
!!      joint per-voxel coupled solve on the summed statistics and realizes the merged basis
!!      through the existing tail (gridcorr, band-limit at the inherited band, mask,
!!      orthonormalize). Merged eigenvolumes: flex_pca_merged_pc*.mrc; meta:
!!      flex_pca_probe_merged.txt. The S.13 fixed-point iteration of the conversion is NOT run
!!      (single-shot S.11).
!!   3. SIMPLE_COV_PAIRED_MERGE=2 is the volume-level pooled-Gram fallback (bag_basis_pool on
!!      the delivered bases) -- cruder (merges restored volumes, not statistics), kept as the
!!      floor the IgG bagging evidence backs.
!! The all-N embedding against the merged basis is the CALLER's stage (run_flex_pca diverts the
!! probe-only early-return to it); this submodule only delivers the merged model handles.
!! The merged realization applies the SAME mean-shaped/background deflation the per-fit tail
!! applies to every delivered basis (the joint solve's raw output re-enters the contrast/
!! background component exactly as the per-fit solves do -- measured 19-42% of basis energy;
!! skipping it read as a spurious span disagreement). Recorded deviation: no internal-FSC
!! Wiener post-filter exists here -- the cross-fit ridge IS the regularizer, exactly as in
!! refine3D's final ml_reg maps.
submodule (simple_flex_pca_em) simple_flex_pca_em_pairmerge
use simple_flex_pca_util, only: dilation_template
use simple_flex_reconstructor_latent_ops, only: solve_coupled_basis_exp, add_invtausq2rho_coupled, &
    &pair_index
use simple_flex_pca_crossfsc, only: crossfsc_to_invtau2, crossfsc_harvest_h, crossfsc_inband_mean
implicit none
#include "simple_local_flags.inc"

character(len=*), parameter :: MERGED_PC_FBODY  = 'flex_pca_merged_pc'
character(len=*), parameter :: MERGED_META      = 'flex_pca_probe_merged.txt'
character(len=*), parameter :: MERGED_EIG_FNAME = 'flex_pca_eigenvalues_merged.txt'
character(len=*), parameter :: PAIRED_MANIFEST  = 'flex_pca_paired.txt'

contains

    !> Snapshot one fit's LAST-iteration raw M-step sufficient statistics + entry frame.
    !! Called by the paired driver after the batch loop and reductions, BEFORE fit_iter_finish
    !! (which ridges rho, mutates the numerators in the coupled solve, and frees everything).
    !! Overwrites every iteration: convergence and the crossfsc stop are evaluated after the
    !! tails, so any iteration can turn out to be the last. prev_real at this point holds the
    !! PREVIOUS iteration's delivered basis == the frame the statistics' latents live in.
    module subroutine probe_fit_merge_stash( fit )
        type(probe_fit_t), intent(inout) :: fit
        integer :: q, es(3)
        if( .not. (allocated(fit%Yeven) .and. allocated(fit%rho_e)) ) &
            &THROW_HARD('probe_fit_merge_stash: no live M-step accumulators to stash')
        es = fit%es
        ! (re)size on rank or lattice change
        if( allocated(fit%mg_ye) )then
            if( size(fit%mg_ye,4) /= fit%ncomp .or. size(fit%mg_ye,1) /= es(1) )then
                deallocate(fit%mg_ye, fit%mg_yo, fit%mg_rhe, fit%mg_rho)
            endif
        endif
        if( .not. allocated(fit%mg_ye) )then
            allocate(fit%mg_ye (es(1),es(2),es(3),fit%ncomp), fit%mg_yo(es(1),es(2),es(3),fit%ncomp))
            allocate(fit%mg_rhe(fit%npairs,es(1),es(2),es(3)), fit%mg_rho(fit%npairs,es(1),es(2),es(3)))
        endif
        do q = 1, fit%ncomp
            fit%mg_ye(:,:,:,q) = fit%Yeven(q)%cmat_exp
            fit%mg_yo(:,:,:,q) = fit%Yodd(q)%cmat_exp
        end do
        fit%mg_rhe = fit%rho_e
        fit%mg_rho = fit%rho_o
        if( allocated(fit%mg_kpe) ) deallocate(fit%mg_kpe)
        if( allocated(fit%mg_kpo) ) deallocate(fit%mg_kpo)
        if( allocated(fit%mg_rpe) ) deallocate(fit%mg_rpe)
        if( allocated(fit%mg_rpo) ) deallocate(fit%mg_rpo)
        if( fit%l_pcg )then
            allocate(fit%mg_kpe, source=fit%kpk_e)
            allocate(fit%mg_kpo, source=fit%kpk_o)
            allocate(fit%mg_rpe, source=fit%rpk_e)
            allocate(fit%mg_rpo, source=fit%rpk_o)
        endif
        fit%mg_ncomp  = fit%ncomp
        fit%mg_npairs = fit%npairs
        ! entry-frame basis images (absent only at iteration 1; the driver requires >= 2 iterations)
        if( allocated(fit%mg_prev) )then
            do q = 1, size(fit%mg_prev)
                call fit%mg_prev(q)%kill
            end do
            deallocate(fit%mg_prev)
        endif
        if( allocated(fit%prev_real) )then
            if( size(fit%prev_real) /= fit%ncomp ) &
                &THROW_HARD('probe_fit_merge_stash: entry frame rank /= accumulator rank')
            allocate(fit%mg_prev(fit%ncomp))
            do q = 1, fit%ncomp
                call fit%mg_prev(q)%copy(fit%prev_real(q))
            end do
        endif
        fit%l_mg_stash = .true.
    end subroutine probe_fit_merge_stash

    !> THE MERGE (par.7 final-stage steps 1-2). Consumes the two fits' stashes and delivered
    !! state, delivers the merged model handles + the merged eigenvolumes/meta/manifest lines.
    !> Init-deflated matched cosines: project the SHARED deterministic init (the it000 stamps
    !! fit A wrote at initialisation) out of both entry-frame bases before measuring the match.
    !! Raw cross-half agreement at short budgets partly measures init memory -- both fits start
    !! IDENTICALLY -- which inflates every match above any prune threshold (measured: the rank
    !! cut is inert at n_probe_iters <= 3 while cutting 20->13 at 6). The deflated cosine is the
    !! honest reproducibility signal at ANY budget. ok=.false. when no stamps are on disk.
    !!
    !! With thr/sub_cos/pa_cos present it ALSO measures the reproducible SUBSPACE: the deflated
    !! bases are orthonormalised, the principal-angle cosines of their spans are pa_cos (sorted),
    !! the reproducible dimension is the count of pa_cos >= thr, and sub_cos(q) is the cosine of
    !! fit A's deflated component q with that subspace. The per-component matched cosine is
    !! rotation-sensitive (a near-flat spectrum lets components rotate between halves: measured
    !! 0.72 vs 0.48 for the same data under two summation orders while the principal angles agreed
    !! to 3 decimals), the subspace cosine is not, so the rank gate consumes sub_cos.
    subroutine init_deflated_matchcos( imgsA, ncA, imgsB, ncB, pstar, defl_cos, ok, thr, sub_cos, pa_cos )
        use simple_linalg, only: jacobi, eigsrt
        type(image), intent(in)  :: imgsA(:), imgsB(:)
        integer,     intent(in)  :: ncA, ncB, pstar(:)
        real(dp), allocatable, intent(out) :: defl_cos(:)
        logical,     intent(out) :: ok
        real(dp), optional, intent(in) :: thr
        real(dp), allocatable, optional, intent(out) :: sub_cos(:), pa_cos(:)
        type(image) :: ivol
        type(string) :: fn
        real(dp), allocatable :: VA(:,:), VB(:,:), VI(:,:)
        real(dp), allocatable :: QA(:,:), QB(:,:), G(:,:), C(:,:), ev(:), evec(:,:), u(:)
        real, pointer :: rmat(:,:,:) => null()
        real(dp) :: nrm, pj
        integer  :: nvox, q, k, ni, ldim(3), nrot, r, j
        ok = .false.
        ! collect the init stamps (fit A namespace; the init is shared by construction)
        ni = 0
        do q = 1, ncA + 8
            fn = 'flex_pca_pc'//'it000_'//int2str_pad(q,3)//MRC_EXT
            if( .not. file_exists(fn) )then
                call fn%kill
                exit
            endif
            ni = q
            call fn%kill
        end do
        if( ni < 1 ) return
        ldim = imgsA(1)%get_ldim()
        nvox = product(ldim)
        allocate(VA(nvox,ncA), VB(nvox,ncB), VI(nvox,ni))
        do q = 1, ncA
            call imgsA(q)%get_rmat_ptr(rmat)
            VA(:,q) = reshape(real(rmat(1:ldim(1),1:ldim(2),1:ldim(3)),dp), [nvox])
        end do
        do q = 1, ncB
            call imgsB(q)%get_rmat_ptr(rmat)
            VB(:,q) = reshape(real(rmat(1:ldim(1),1:ldim(2),1:ldim(3)),dp), [nvox])
        end do
        call ivol%new(ldim, imgsA(1)%get_smpd())
        do q = 1, ni
            fn = 'flex_pca_pc'//'it000_'//int2str_pad(q,3)//MRC_EXT
            call ivol%read(fn)
            call ivol%get_rmat_ptr(rmat)
            VI(:,q) = reshape(real(rmat(1:ldim(1),1:ldim(2),1:ldim(3)),dp), [nvox])
            call fn%kill
        end do
        call ivol%kill
        ! Gram-Schmidt the init, then deflate + renormalize both bases
        do q = 1, ni
            do k = 1, q - 1
                VI(:,q) = VI(:,q) - dot_product(VI(:,q), VI(:,k))*VI(:,k)
            end do
            nrm = sqrt(max(sum(VI(:,q)**2), 1.d-30))
            VI(:,q) = VI(:,q)/nrm
        end do
        do q = 1, ncA
            do k = 1, ni
                VA(:,q) = VA(:,q) - dot_product(VA(:,q), VI(:,k))*VI(:,k)
            end do
            nrm = sqrt(max(sum(VA(:,q)**2), 1.d-30))
            VA(:,q) = VA(:,q)/nrm
        end do
        do q = 1, ncB
            do k = 1, ni
                VB(:,q) = VB(:,q) - dot_product(VB(:,q), VI(:,k))*VI(:,k)
            end do
            nrm = sqrt(max(sum(VB(:,q)**2), 1.d-30))
            VB(:,q) = VB(:,q)/nrm
        end do
        allocate(defl_cos(ncA))
        do q = 1, ncA
            if( pstar(q) >= 1 .and. pstar(q) <= ncB )then
                defl_cos(q) = abs(dot_product(VA(:,q), VB(:,pstar(q))))
            else
                defl_cos(q) = 0.d0
            endif
        end do
        if( present(sub_cos) .and. present(pa_cos) .and. present(thr) )then
            ! orthonormal bases of the two deflated spans (Gram-Schmidt; deflation broke the
            ! orthogonality the fits delivered), cross-Gram G = QA^T QB, principal cosines =
            ! sqrt(eig(G G^T)), principal directions in A's span = QA evec
            allocate(QA(nvox,ncA), QB(nvox,ncB), source=0.d0)
            QA = VA
            do q = 1, ncA
                do k = 1, q - 1
                    QA(:,q) = QA(:,q) - dot_product(QA(:,q), QA(:,k))*QA(:,k)
                end do
                nrm = sqrt(max(sum(QA(:,q)**2), 1.d-30))
                QA(:,q) = QA(:,q)/nrm
            end do
            QB = VB
            do q = 1, ncB
                do k = 1, q - 1
                    QB(:,q) = QB(:,q) - dot_product(QB(:,q), QB(:,k))*QB(:,k)
                end do
                nrm = sqrt(max(sum(QB(:,q)**2), 1.d-30))
                QB(:,q) = QB(:,q)/nrm
            end do
            allocate(G(ncA,ncB), C(ncA,ncA), ev(ncA), evec(ncA,ncA), u(nvox), pa_cos(ncA), sub_cos(ncA))
            G = matmul(transpose(QA), QB)
            C = matmul(G, transpose(G))
            call jacobi(C, ncA, ncA, ev, evec, nrot)
            call eigsrt(ev, evec, ncA, ncA)
            do j = 1, ncA
                pa_cos(j) = sqrt(max(0.d0, min(1.d0, ev(j))))
            end do
            r = count(pa_cos >= thr)
            do q = 1, ncA
                sub_cos(q) = 0.d0
                do j = 1, r
                    u  = matmul(QA, evec(:,j))
                    pj = dot_product(VA(:,q), u)
                    sub_cos(q) = sub_cos(q) + pj*pj
                end do
                sub_cos(q) = sqrt(min(1.d0, sub_cos(q)))
            end do
            deallocate(QA, QB, G, C, ev, evec, u)
        endif
        deallocate(VA, VB, VI)
        ok = .true.
    end subroutine init_deflated_matchcos

    module subroutine probe_paired_merge( params, build, fits, merge_mode, m_basis_recs, &
        &m_eigvals, m_ncomp, m_sig2, m_matchcos )
        class(parameters),   intent(inout) :: params
        type(builder),       intent(inout) :: build
        type(probe_fit_t),   intent(inout) :: fits(2)
        integer,             intent(in)    :: merge_mode
        type(reconstructor), allocatable, intent(out) :: m_basis_recs(:)
        real(dp),            allocatable, intent(out) :: m_eigvals(:)
        integer,             intent(out)              :: m_ncomp
        real(dp),            intent(out)              :: m_sig2
        !> per merged component: |cos| of the A<-B match (the reproducibility signal the rank
        !! gate consumes; 1.0 padding when the mode-2 pooled path skips matching)
        real(dp), allocatable, optional, intent(out)  :: m_matchcos(:)
        type(reconstructor), allocatable :: Ymrg(:), utilde(:)
        type(image),         allocatable :: realvols(:), utilde_real(:), pooled(:)
        type(image)  :: mstep_gridcorr, imga, imgb
        type(string) :: fname
        real(dp), allocatable :: Mba(:,:), Rm(:,:), sv_ab(:), eigvals_m(:), eig_pooled(:)
        real,     allocatable :: rho_me(:,:,:,:), rho_mo(:,:,:,:), kpk_m(:,:)
        complex,  allocatable :: rpk_m(:,:)
        type(flex_pcg_outcome_t) :: pcg_out
        logical  :: l_pcg_m
        real,     allocatable :: h_e(:,:), h_o(:,:), h_sum(:,:), invtau2(:,:), fscq(:,:), corrs(:)
        integer,  allocatable :: cnt(:), pstar(:), psign(:)
        real,     pointer     :: rmatp(:,:,:)
        real(dp) :: rdev, gA, gB, paircos
        real     :: lp_m, fmean
        integer  :: ncm, ncb, npairs_m, es(3), lb(3), nyq, filtsz, klo, kfr(2), khi_m
        integer  :: q, p, d_new, ncA_del, ncB_del
        ! ---- gate-4 pair-algebra self-check: packed R^T rho R vs brute-force dense, K=2 and a
        ! rectangular 3->2 case, deterministic values; THROWs on mismatch ----
        call merge_pair_algebra_selfcheck
        m_ncomp = 0
        m_sig2  = (real(fits(1)%nptcls,dp)*fits(1)%sig2_eff + real(fits(2)%nptcls,dp)*fits(2)%sig2_eff) &
            &/ real(max(1, fits(1)%nptcls + fits(2)%nptcls),dp)
        ncA_del = fits(1)%ncomp
        ncB_del = fits(2)%ncomp
        filtsz  = max(1, fdim(params%box_crop) - 1)
        if( merge_mode >= 2 )then
            ! ---- pooled-Gram fallback (bag_basis_pool on the DELIVERED bases) ----
            write(logfhandle,'(A)') '>>> FLEX_PCA MERGE mode=2: volume-level pooled-Gram fallback &
                &(bag_basis_pool over both delivered bases; statistics not merged)'
            call flush(logfhandle)
            call bag_basis_pool(fits(1)%prev_real, ncA_del, fits(1)%eigvals, &
                &fits(2)%prev_real, ncB_del, fits(2)%eigvals, ncA_del, pooled, eig_pooled)
            d_new = size(pooled)
            ! same gauge as mode 1: the pooled-Gram eigenvectors rotate freely in a
            ! near-degenerate span; fix the frame to fit A's delivered components
            if( d_new == ncA_del ) call gauge_fix_to_frame(fits(1)%prev_real, ncA_del, pooled, d_new)
            do q = 1, d_new
                fname = MERGED_PC_FBODY//int2str_pad(q,3)//MRC_EXT
                call pooled(q)%write(fname, del_if_exists=.true.); call fname%kill
            end do
            call basis_recs_from_images(params, build, pooled, d_new, m_basis_recs)
            allocate(eigvals_m(d_new), source=eig_pooled(1:d_new))
            call merged_delivery(params, fits, 2, pooled, d_new, eigvals_m, m_sig2)
            call move_alloc(eigvals_m, m_eigvals)
            m_ncomp = d_new
            do q = 1, d_new
                call pooled(q)%kill
            end do
            deallocate(pooled, eig_pooled)
            return
        endif
        ! ================= mode 1: accumulator-level merge =================
        if( .not. (fits(1)%l_mg_stash .and. fits(2)%l_mg_stash) ) &
            &THROW_HARD('paired merge: no stashed M-step statistics (driver gate error)')
        if( .not. (allocated(fits(1)%mg_prev) .and. allocated(fits(2)%mg_prev)) ) &
            &THROW_HARD('paired merge: entry-frame bases missing; the merge needs n_probe_iters >= 2')
        ncm = fits(1)%mg_ncomp        ! merged rank = fit A's stash (entry) rank
        ncb = fits(2)%mg_ncomp
        npairs_m = (ncm*(ncm+1))/2
        if( size(fits(1)%mg_ye,1) /= size(fits(2)%mg_ye,1) .or. &
           &size(fits(1)%mg_ye,2) /= size(fits(2)%mg_ye,2) .or. &
           &size(fits(1)%mg_ye,3) /= size(fits(2)%mg_ye,3) ) &
            &THROW_HARD('paired merge: the two stashes live on different exp lattices')
        es(1) = size(fits(1)%mg_ye,1); es(2) = size(fits(1)%mg_ye,2); es(3) = size(fits(1)%mg_ye,3)
        if( ncb /= ncm ) write(logfhandle,'(A,I0,A,I0,A)') &
            &'>>> FLEX_PCA MERGE WARNING: per-fit stash ranks differ (A=',ncm,' B=',ncb, &
            &'); the frame map is semi-orthogonal and B contributes only its matched span'
        ! ---- step 1: the frame map R (entry frames; see the module header) ----
        call align_basis_to_reference(fits(2)%mg_prev, ncb, fits(1)%mg_prev, ncm, Mba, sv_ab)
        call polar_factor(Mba, ncb, ncm, Rm)
        rdev = 0.d0
        do q = 1, ncm
            do p = 1, ncm
                rdev = max(rdev, abs(sum(Rm(:,q)*Rm(:,p)) - merge(1.d0, 0.d0, q == p)))
            end do
        end do
        write(logfhandle,'(A,I0,A,I0,A,ES9.2)') '>>> FLEX_PCA MERGE frame map R: ',ncb,' x ',ncm, &
            &'  max |R^T R - I|=',rdev
        write(logfhandle,'(A)',advance='no') '>>> FLEX_PCA MERGE principal-angle cos(A,B) ='
        do q = 1, size(sv_ab)
            write(logfhandle,'(1X,F7.4)',advance='no') real(sv_ab(q))
        end do
        write(logfhandle,*)
        allocate(pstar(ncm), psign(ncm))
        if( present(m_matchcos) )then
            if( allocated(m_matchcos) ) deallocate(m_matchcos)
            allocate(m_matchcos(ncm), source=1.d0)
        endif
        do q = 1, ncm
            pstar(q) = maxloc(abs(Rm(:,q)), dim=1)
            psign(q) = merge(1, -1, Rm(pstar(q),q) >= 0.d0)
            paircos  = Mba(pstar(q),q)
            if( present(m_matchcos) ) m_matchcos(q) = abs(paircos)
            write(logfhandle,'(A,I3,A,I3,A,I2,A,F7.4,A,F7.4)') '>>> FLEX_PCA MERGE pair  A', q, &
                &' <- B', pstar(q), '  sign=', psign(q), '  |cos|=', abs(real(paircos)), &
                &'  |R|=', real(abs(Rm(pstar(q),q)))
            if( abs(paircos) < 0.5d0 ) write(logfhandle,'(A,I0,A)') &
                &'>>> FLEX_PCA MERGE WARNING: component ',q,' matched below |cos|=0.5; the rank &
                &gate should have retired this axis before the merge'
        end do
        call flush(logfhandle)
        ! ---- init-deflated match (rank-gate signal at ANY budget) ----
        if( present(m_matchcos) )then
            block
                real(dp), allocatable :: dcos(:), scos(:), pcos(:)
                real(dp) :: thr
                logical :: l_dok
                integer :: qd, rdim
                thr = 0.50d0    ! the reporting bar for a reproducible principal direction
                call init_deflated_matchcos(fits(1)%mg_prev, ncm, fits(2)%mg_prev, ncb, &
                    &pstar, dcos, l_dok, thr=thr, sub_cos=scos, pa_cos=pcos)
                if( l_dok )then
                    write(logfhandle,'(A)',advance='no') '>>> FLEX_PCA MERGE init-DEFLATED match |cos| ='
                    do qd = 1, ncm
                        write(logfhandle,'(1X,F6.3)',advance='no') real(dcos(qd))
                    end do
                    write(logfhandle,*)
                    write(logfhandle,'(A)',advance='no') '>>> FLEX_PCA MERGE init-DEFLATED principal-angle cos ='
                    do qd = 1, ncm
                        write(logfhandle,'(1X,F6.3)',advance='no') real(pcos(qd))
                    end do
                    write(logfhandle,*)
                    rdim = count(pcos >= thr)
                    write(logfhandle,'(A,I0,A,F4.2,A)',advance='no') '>>> FLEX_PCA MERGE reproducible subspace: dim=', &
                        &rdim, ' (principal cos >= ', real(thr), ');  fit-A component cos with it ='
                    do qd = 1, ncm
                        write(logfhandle,'(1X,F6.3)',advance='no') real(scos(qd))
                    end do
                    write(logfhandle,*)
                    call flush(logfhandle)
                    ! the rank gate and the axis weights consume the SUBSPACE cosine
                    m_matchcos(1:ncm) = scos(1:ncm)
                    deallocate(dcos, scos, pcos)
                else
                    write(logfhandle,'(A)') '>>> FLEX_PCA MERGE init-deflated match unavailable &
                        &(no it000 stamps); rank gate falls back to the RAW match -- inflated at &
                        &short budgets'
                    call flush(logfhandle)
                endif
            end block
        endif
        ! ---- step 2a: merged even/odd pair -- numerators (Y' = Y.R) and densities
        ! (rho' = R^T rho R), the four quarter-sets summed ----
        allocate(Ymrg(ncm))
        do q = 1, ncm
            call init_basis_reconstructor(params, build, Ymrg(q))
            call Ymrg(q)%reset; call Ymrg(q)%reset_exp
            ! A-even + A-odd
            Ymrg(q)%cmat_exp = fits(1)%mg_ye(:,:,:,q) + fits(1)%mg_yo(:,:,:,q)
            ! + rotated B-even + B-odd
            do p = 1, ncb
                if( abs(Rm(p,q)) < 1.d-8 ) cycle
                Ymrg(q)%cmat_exp = Ymrg(q)%cmat_exp &
                    &+ real(Rm(p,q))*(fits(2)%mg_ye(:,:,:,p) + fits(2)%mg_yo(:,:,:,p))
            end do
        end do
        allocate(rho_me(npairs_m,es(1),es(2),es(3)), rho_mo(npairs_m,es(1),es(2),es(3)), source=0.)
        rho_me = fits(1)%mg_rhe        ! A-even (npairs_m == fit A's stash npairs by construction)
        rho_mo = fits(1)%mg_rho        ! A-odd
        call rotate_rho4_packed_add(Rm, ncb, ncm, fits(2)%mg_rhe, rho_me)
        call rotate_rho4_packed_add(Rm, ncb, ncm, fits(2)%mg_rho, rho_mo)
        ! rec_backend=pcg: the merged pair kernels, the same congruence on the packed sums, both halves
        l_pcg_m = fits(1)%l_pcg .and. fits(2)%l_pcg
        if( l_pcg_m )then
            if( .not. (allocated(fits(1)%mg_kpe) .and. allocated(fits(2)%mg_kpe)) ) &
                &THROW_HARD('probe_paired_merge: PCG fits without stashed pair kernels')
            call fits(1)%pcg%new(params%box_crop, params%smpd_crop, ncm)
            call flex_pcg_install_window(fits(1)%pcg, params)
            call fits(1)%pcg%alloc_packed(kpk_m)
            kpk_m = fits(1)%mg_kpe + fits(1)%mg_kpo
            call rotate_rho_packed_add(Rm, ncb, ncm, fits(2)%mg_kpe, kpk_m)
            call rotate_rho_packed_add(Rm, ncb, ncm, fits(2)%mg_kpo, kpk_m)
            deallocate(fits(1)%mg_kpe, fits(1)%mg_kpo, fits(2)%mg_kpe, fits(2)%mg_kpo)
            ! the merged right-hand sides: the same congruence as the numerators, Y' = Y.R
            call fits(1)%pcg%alloc_rhs_packed(rpk_m)
            rpk_m = fits(1)%mg_rpe + fits(1)%mg_rpo
            call rotate_rhs_packed_add(Rm, ncb, ncm, fits(2)%mg_rpe, rpk_m)
            call rotate_rhs_packed_add(Rm, ncb, ncm, fits(2)%mg_rpo, rpk_m)
            deallocate(fits(1)%mg_rpe, fits(1)%mg_rpo, fits(2)%mg_rpe, fits(2)%mg_rpo)
        endif
        ! stashes are consumed; free before the solve allocations peak
        deallocate(fits(1)%mg_ye, fits(1)%mg_yo, fits(1)%mg_rhe, fits(1)%mg_rho)
        deallocate(fits(2)%mg_ye, fits(2)%mg_yo, fits(2)%mg_rhe, fits(2)%mg_rho)
        fits(1)%l_mg_stash = .false.; fits(2)%l_mg_stash = .false.
        ! step 2b: per-shell sampling H from the merged halves' pair diagonals; the ridge
        ! consumes the summed H (h_e + h_o), equivalent to per-half addition
        lb  = lbound(Ymrg(1)%cmat_exp)
        nyq = Ymrg(1)%get_lfny(1)
        allocate(h_e(filtsz,ncm), h_o(filtsz,ncm), h_sum(filtsz,ncm), cnt(filtsz))
        call crossfsc_harvest_h(rho_me, npairs_m, ncm, lb, nyq, filtsz, h_e, cnt)
        call crossfsc_harvest_h(rho_mo, npairs_m, ncm, lb, nyq, filtsz, h_o, cnt)
        h_sum = h_e + h_o
        ! step 2c: cross-fit FSC per merged axis (fit A entry component q vs sign * fit B entry
        ! component pstar(q)) -> sampling-aware ridge
        kfr = covariance_kfromto(params)
        klo = max(6, kfr(1))
        allocate(fscq(filtsz,ncm), invtau2(ncm,filtsz), corrs(filtsz))
        do q = 1, ncm
            call imga%copy(fits(1)%mg_prev(q))
            call imgb%copy(fits(2)%mg_prev(pstar(q)))
            if( psign(q) < 0 ) call imgb%mul(-1.0)
            call imga%fft; call imgb%fft
            call imga%fsc(imgb, corrs)
            fscq(:,q) = corrs
            call imga%kill; call imgb%kill
            call crossfsc_to_invtau2(fscq(:,q), h_sum(:,q), params%tau, klo, invtau2(q,:))
        end do
        ! ---- step 2d: sum the merged halves and run ONE joint per-voxel coupled solve on the
        ! total statistics (numerators already hold e+o; densities summed here). The ridge is
        ! added ONCE with the summed H == per-half ridge added to each half (see header). ----
        rho_me = rho_me + rho_mo
        deallocate(rho_mo)
        call add_invtausq2rho_coupled(Ymrg, rho_me, ncm, invtau2)
        lp_m = max(fits(1)%lp_it, fits(2)%lp_it)   ! inherited band: min(khi_A, khi_B)
        khi_m = filtsz
        if( lp_m > 2.0*params%smpd_crop + TINY ) &
            &khi_m = max(2, min(filtsz, int(fits(1)%dstep_ann/lp_m)))
        write(logfhandle,'(A,F7.2,A,I0,A,F8.4)') '>>> FLEX_PCA MERGE joint coupled solve: band lp=', &
            &lp_m,' A (k_hi=',khi_m,')  tau fudge=',params%tau
        do q = 1, ncm
            fmean = crossfsc_inband_mean(fscq(:,q), khi_m)
            write(logfhandle,'(A,I3,A,F8.4,A,ES10.3)') '>>> FLEX_PCA MERGE crossFSC comp ',q, &
                &'  in-band mean=',fmean,'  mean invtau2 in band=', &
                &sum(invtau2(q,2:khi_m))/real(max(1,khi_m-1))
        end do
        call flush(logfhandle)
        if( l_pcg_m )then
            call fits(1)%pcg%finalize(kpk_m)
            call fits(1)%pcg%set_ridge(invtau2)
            call fits(1)%pcg%solve(Ymrg, rho_me, rpk_m, params%maxits_pcg, params%rtol, pcg_out, 'FLEX_PCA PCG MERGE')
            call fits(1)%pcg%clear_ridge
            write(logfhandle,'(A,I0,A,ES10.3,A,ES10.3,A,ES10.3,A,A,A,F8.1)') '>>> FLEX_PCA PCG MERGE joint solve: iters=', &
                &pcg_out%iteration_count, '  init=', pcg_out%initial_rel_residual, '  resid=', pcg_out%final_rel_residual, &
                &'  update=', pcg_out%final_rel_update, '  stop=', trim(pcg_out%stop_reason), '  seconds=', pcg_out%seconds
            call flush(logfhandle)
            deallocate(kpk_m, rpk_m)
        else
            call solve_coupled_basis_exp(Ymrg, rho_me, ncm)
        endif
        deallocate(rho_me)
        ! realize through the existing tail: compress, deapodize, band-limit, mask
        allocate(realvols(ncm))
        mstep_gridcorr = prep3D_inv_kbenvelope4mul([params%box_crop,params%box_crop,params%box_crop], &
            &params%smpd_crop)
        do q = 1, ncm
            Ymrg(q)%rho_exp = 1.0
            call Ymrg(q)%compress_exp; call Ymrg(q)%ifft
            call realvols(q)%new([params%box_crop,params%box_crop,params%box_crop], params%smpd_crop)
            call Ymrg(q)%get_rmat_ptr(rmatp); call realvols(q)%set_rmat(rmatp, .false.)
            call realvols(q)%mul(mstep_gridcorr)
            if( lp_m > 2.0*params%smpd_crop + TINY )then
                call realvols(q)%fft; call realvols(q)%bp(0., lp_m); call realvols(q)%ifft
            endif
            call flex_window_apply(realvols(q), params)
            call Ymrg(q)%dealloc_rho; call Ymrg(q)%kill
        end do
        call mstep_gridcorr%kill
        deallocate(Ymrg)
        ! mean-shaped/background deflation as the per-fit tail applies it: the joint solve's raw
        ! output carries the same re-entering contrast/background component
        call merge_deflate_mean_shaped(params, fits(1), realvols, ncm)
        call orthonormalize_representatives(params, build, realvols, ncm, utilde, utilde_real, d_new)
        ! gauge fix: the joint solve + re-orthonormalisation return an arbitrary in-span frame, so
        ! the merged set is rotated onto fit A's entry frame (where the statistics and the
        ! index-aligned Gamma live); span, model and embedding are unchanged
        if( d_new == ncm )then
            call gauge_fix_to_frame(fits(1)%mg_prev, ncm, utilde_real, d_new)
            ! utilde was realized from the pre-gauge images; rebuild in the fixed gauge
            do q = 1, size(utilde)
                call utilde(q)%dealloc_rho; call utilde(q)%kill
            end do
            deallocate(utilde)
            call basis_recs_from_images(params, build, utilde_real, d_new, utilde)
        else
            write(logfhandle,'(A,I0,A,I0,A)') '>>> FLEX_PCA MERGE gauge fix skipped: the &
                &orthonormalisation dropped rank (',d_new,' of ',ncm,'); frame left as realized'
            call flush(logfhandle)
        endif
        ! merged Gamma: index-aligned average of fit A's Gamma and fit B's rotated Gamma
        ! (R^2 weights), the existing tail's idiom
        allocate(eigvals_m(d_new))
        do q = 1, d_new
            gA = fits(1)%eigvals(min(q, ncA_del))
            gB = 0.d0
            do p = 1, ncb
                gB = gB + Rm(p,min(q,ncm))**2 * fits(2)%eigvals(min(p, ncB_del))
            end do
            eigvals_m(q) = max(0.5d0*(gA + gB), DTINY)
        end do
        do q = 1, d_new
            fname = MERGED_PC_FBODY//int2str_pad(q,3)//MRC_EXT
            call utilde_real(q)%write(fname, del_if_exists=.true.); call fname%kill
        end do
        call merged_delivery(params, fits, 1, utilde_real, d_new, eigvals_m, m_sig2, &
            &Rm=Rm, Mba=Mba, sv_ab=sv_ab, pstar=pstar, psign=psign, fscq=fscq, khi_m=khi_m)
        ! hand the merged model to the caller (the all-N embedding stage)
        call move_alloc(utilde, m_basis_recs)
        call move_alloc(eigvals_m, m_eigvals)
        m_ncomp = d_new
        do q = 1, size(realvols)
            call realvols(q)%kill
        end do
        do q = 1, size(utilde_real)
            call utilde_real(q)%kill
        end do
        deallocate(realvols, utilde_real, Mba, Rm, sv_ab, pstar, psign)
        deallocate(h_e, h_o, h_sum, cnt, fscq, invtau2, corrs)
    end subroutine probe_paired_merge

    !> Shared delivery for both merge modes: merged probe meta + eigenvalue table, the
    !! merged-vs-fitA / merged-vs-fitB matched |cos| tables (log + manifest), and the merge
    !! provenance appended to the paired manifest.
    subroutine merged_delivery( params, fits, mode, mimgs, nm, eigvals_m, sig2_m, &
        &Rm, Mba, sv_ab, pstar, psign, fscq, khi_m )
        class(parameters), intent(in)    :: params
        type(probe_fit_t), intent(inout) :: fits(2)
        integer,           intent(in)    :: mode, nm
        type(image),       intent(inout) :: mimgs(:)
        real(dp),          intent(in)    :: eigvals_m(:)
        real(dp),          intent(in)    :: sig2_m
        real(dp), optional, intent(in)   :: Rm(:,:), Mba(:,:), sv_ab(:)
        integer,  optional, intent(in)   :: pstar(:), psign(:)
        real,     optional, intent(in)   :: fscq(:,:)
        integer,  optional, intent(in)   :: khi_m
        integer,  allocatable :: ma(:), mb(:)
        real(dp), allocatable :: mcos_a(:), mcos_b(:)
        integer  :: q, p, u
        ! merged meta + eigenvalue table
        call save_probe_state(nm, eigvals_m, sig2_m, fname=MERGED_META)
        call del_file(MERGED_EIG_FNAME)
        open(newunit=u, file=MERGED_EIG_FNAME, status='replace', action='write')
        write(u,'(A)') '# component eigenvalue'
        do q = 1, nm
            write(u,'(I6,1X,ES20.10)') q, eigvals_m(q)
        end do
        close(u)
        ! matched |cos| of the merged basis against each fit's DELIVERED basis
        call matched_abs_cos(mimgs, nm, fits(1)%prev_real, fits(1)%ncomp, ma, mcos_a)
        call matched_abs_cos(mimgs, nm, fits(2)%prev_real, fits(2)%ncomp, mb, mcos_b)
        do q = 1, nm
            write(logfhandle,'(A,I3,A,I3,A,F7.4,A,I3,A,F7.4)') '>>> FLEX_PCA MERGED comp ',q, &
                &'  vs fitA ',ma(q),' |cos|=',real(mcos_a(q)),'   vs fitB ',mb(q),' |cos|=', &
                &real(mcos_b(q))
        end do
        write(logfhandle,'(A,I0,A,F7.4,A,F7.4)') '>>> FLEX_PCA MERGED delivered: ncomp=',nm, &
            &'  min matched |cos| vs fitA=',real(minval(mcos_a)),'  vs fitB=',real(minval(mcos_b))
        call flush(logfhandle)
        ! merge provenance -> the paired manifest (append; the paired record was written first)
        open(newunit=u, file=PAIRED_MANIFEST, status='unknown', position='append', action='write')
        write(u,'(A,I0,A,I0,A,ES16.8,A,F7.2)') '# merge: mode=',mode,'  ncomp=',nm,'  sig2=', &
            &sig2_m,'  lp=',max(fits(1)%lp_it, fits(2)%lp_it)
        if( present(Rm) )then
            write(u,'(A)') '# merge frame map R (rows: fit B entry components; cols: fit A entry components)'
            do p = 1, size(Rm,1)
                write(u,'(A)',advance='no') '# merge R'
                do q = 1, size(Rm,2)
                    write(u,'(1X,F9.5)',advance='no') real(Rm(p,q))
                end do
                write(u,*)
            end do
            write(u,'(A)',advance='no') '# merge principal-angle cos(A,B):'
            do q = 1, size(sv_ab)
                write(u,'(1X,F7.4)',advance='no') real(sv_ab(q))
            end do
            write(u,*)
            write(u,'(A)') '# merge pair table: Aq  Bp  sign  |cos|  crossFSC_inband_mean'
            do q = 1, size(pstar)
                write(u,'(A,I4,1X,I4,1X,I3,1X,F7.4,1X,F8.4)') '# merge pair ', q, pstar(q), &
                    &psign(q), abs(real(Mba(pstar(q),q))), &
                    &crossfsc_inband_mean(fscq(:,q), khi_m)
            end do
        else
            write(u,'(A)') '# merge mode=2: pooled-Gram fallback (bag_basis_pool over the delivered bases)'
        endif
        write(u,'(A)') '# merge matched: comp  fitA_comp  |cos|A  fitB_comp  |cos|B'
        do q = 1, nm
            write(u,'(A,I4,1X,I4,1X,F7.4,1X,I4,1X,F7.4)') '# merge matched ', q, ma(q), &
                &real(mcos_a(q)), mb(q), real(mcos_b(q))
        end do
        close(u)
        deallocate(ma, mb, mcos_a, mcos_b)
    end subroutine merged_delivery

    !> Greedy matched |cos| between two real-space volume stacks (the xfsc_paired_record
    !! matching, factored): my(i) is the y-component matched to x-component i, without
    !! replacement while y-components last, then best-|cos| with replacement.
    subroutine matched_abs_cos( ximgs, nx, yimgs, ny, my, mcos )
        integer,     intent(in)    :: nx, ny
        type(image), intent(inout) :: ximgs(:), yimgs(:)
        integer,  allocatable, intent(out) :: my(:)
        real(dp), allocatable, intent(out) :: mcos(:)
        real,    pointer     :: prx(:,:,:), pry(:,:,:)
        real(dp), allocatable :: cosm(:,:)
        logical,  allocatable :: used_x(:), used_y(:)
        real(dp) :: nrm_x, nrm_y, c, best_c
        integer  :: q, p, k, best_q, best_p
        allocate(cosm(nx,ny), used_x(nx), used_y(ny), my(nx), mcos(nx))
        used_x = .false.; used_y = .false.
        my = 0; mcos = 0.d0
        do p = 1, ny
            call yimgs(p)%get_rmat_ptr(pry)
            nrm_y = sqrt(max(sum(real(pry,dp)**2), DTINY))
            do q = 1, nx
                call ximgs(q)%get_rmat_ptr(prx)
                nrm_x = sqrt(max(sum(real(prx,dp)**2), DTINY))
                cosm(q,p) = sum(real(prx,dp)*real(pry,dp)) / (nrm_x*nrm_y)
            end do
        end do
        do k = 1, min(nx, ny)
            best_c = -1.d0; best_q = 0; best_p = 0
            do q = 1, nx
                if( used_x(q) ) cycle
                do p = 1, ny
                    if( used_y(p) ) cycle
                    c = abs(cosm(q,p))
                    if( c > best_c )then
                        best_c = c; best_q = q; best_p = p
                    endif
                end do
            end do
            used_x(best_q) = .true.; used_y(best_p) = .true.
            my(best_q)   = best_p
            mcos(best_q) = best_c
        end do
        ! surplus x-components (nx > ny): best |cos| with replacement, honestly reported
        do q = 1, nx
            if( my(q) > 0 ) cycle
            best_c = -1.d0; best_p = 1
            do p = 1, ny
                c = abs(cosm(q,p))
                if( c > best_c )then
                    best_c = c; best_p = p
                endif
            end do
            my(q) = best_p; mcos(q) = best_c
        end do
        deallocate(cosm, used_x, used_y)
    end subroutine matched_abs_cos

    !> Gauge-fix an orthonormal volume set onto a reference frame of the same rank: the
    !! orthogonal Procrustes rotation T = polar(X^T Aref) applied in-span (new_q = sum_k
    !! T(k,q) X_k), so new_q approximates ref_q as closely as an orthogonal in-span map
    !! allows. Logs the principal-angle cosines of the two spans (the faithfulness of the
    !! merge itself, gauge-independent) before rotating.
    subroutine gauge_fix_to_frame( ref_imgs, nref, imgs, n )
        integer,     intent(in)    :: nref, n
        type(image), intent(inout) :: ref_imgs(:), imgs(:)
        type(image), allocatable :: rot(:)
        real(dp), allocatable :: Mra(:,:), sv(:), Tg(:,:), Mt(:,:)
        integer :: q, k
        if( n /= nref ) THROW_HARD('gauge_fix_to_frame: rank mismatch')
        ! Mra(i,j) = <ref_i, img_j>; svals = span principal-angle cosines (log: faithfulness)
        call align_basis_to_reference(ref_imgs, nref, imgs, n, Mra, sv)
        write(logfhandle,'(A)',advance='no') '>>> FLEX_PCA MERGE span principal-angle cos &
            &(merged span vs gauge frame):'
        do q = 1, size(sv)
            write(logfhandle,'(1X,F7.4)',advance='no') real(sv(q))
        end do
        write(logfhandle,*)
        call flush(logfhandle)
        ! T = polar(X^T Aref) = polar(Mra^T)
        allocate(Mt(n,nref))
        Mt = transpose(Mra)
        call polar_factor(Mt, n, nref, Tg)
        allocate(rot(n))
        do q = 1, n
            call rot(q)%copy(imgs(1))
            call rot(q)%zero_and_unflag_ft
            do k = 1, n
                call rot(q)%add(imgs(k), real(Tg(k,q)))
            end do
        end do
        do q = 1, n
            call imgs(q)%copy(rot(q))
            call rot(q)%kill
        end do
        deallocate(rot, Mra, sv, Tg, Mt)
        write(logfhandle,'(A)') '>>> FLEX_PCA MERGE gauge fixed: merged components rotated &
            &in-span onto the reference frame (span and model unchanged)'
        call flush(logfhandle)
    end subroutine gauge_fix_to_frame

    !> Orthogonal polar factor R = M (M^T M)^{-1/2} of the (ncb x ncm) cross-Gram -- the
    !! Procrustes frame-rotation precedent (jacobi on M^T M, inverse-sqrt with a 1e-12 floor).
    !! Semi-orthogonal (columns orthonormal) when ncb >= ncm.
    subroutine polar_factor( Mba, ncb, ncm, Rm )
        integer,  intent(in)  :: ncb, ncm
        real(dp), intent(in)  :: Mba(ncb,ncm)
        real(dp), allocatable, intent(out) :: Rm(:,:)
        real(dp) :: MtM(ncm,ncm), Vo(ncm,ncm), Wo(ncm,ncm), evo(ncm)
        integer  :: q, nrot
        MtM = matmul(transpose(Mba), Mba)
        call jacobi(MtM, ncm, ncm, evo, Vo, nrot)
        do q = 1, ncm
            evo(q) = 1.d0/sqrt(max(evo(q), 1.d-12))
        end do
        do q = 1, ncm
            Wo(:,q) = Vo(:,q)*evo(q)
        end do
        allocate(Rm(ncb,ncm))
        Rm = matmul(Mba, matmul(Wo, transpose(Vo)))
    end subroutine polar_factor

    !> One voxel of the packed-pair quadratic transform: dout = pack(R^T unpack(din) R).
    !! din is fit B's packed coupled density row-set at one voxel (ncb*(ncb+1)/2, pair_index
    !! order), dout the (ncm*(ncm+1)/2) packed result in fit A's frame. dp accumulation,
    !! sp storage -- matching the accumulators themselves.
    pure subroutine rotate_pair_voxel( Rm, ncb, ncm, din, dout )
        integer,  intent(in)  :: ncb, ncm
        real(dp), intent(in)  :: Rm(ncb,ncm)
        real,     intent(in)  :: din(:)
        real,     intent(out) :: dout(:)
        real(dp) :: Dd(ncb,ncb), T(ncb,ncm), s
        integer  :: p, s2, q, r
        do s2 = 1, ncb
            do p = 1, s2
                Dd(p,s2) = real(din(pair_index(p,s2)),dp)
                Dd(s2,p) = Dd(p,s2)
            end do
        end do
        ! T = D R
        do q = 1, ncm
            do p = 1, ncb
                s = 0.d0
                do s2 = 1, ncb
                    s = s + Dd(p,s2)*Rm(s2,q)
                end do
                T(p,q) = s
            end do
        end do
        ! dout(pair(q,r)) = (R^T T)(q,r), upper triangle q <= r
        do r = 1, ncm
            do q = 1, r
                s = 0.d0
                do p = 1, ncb
                    s = s + Rm(p,q)*T(p,r)
                end do
                dout(pair_index(q,r)) = real(s)
            end do
        end do
    end subroutine rotate_pair_voxel

    !> Rotate fit B's packed coupled density into fit A's frame and add it to the merged density:
    !! rho_acc += pack(R^T unpack(rho_b) R), per voxel.
    !> packed right-hand sides of fit B rotated into the merged frame and added: acc(q') += sum_p R(p,q') b(p)
    subroutine rotate_rhs_packed_add( Rm, ncb, ncm, rpk_b, rpk_acc )
        real(dp), intent(in)    :: Rm(:,:)
        integer,  intent(in)    :: ncb, ncm
        complex,  intent(in)    :: rpk_b(:,:)
        complex,  intent(inout) :: rpk_acc(:,:)
        integer :: t, p, q
        if( size(rpk_b,1) /= ncb .or. size(rpk_acc,1) /= ncm ) &
            &THROW_HARD('rotate_rhs_packed_add: leading extents do not match the ranks')
        if( size(rpk_b,2) /= size(rpk_acc,2) ) THROW_HARD('rotate_rhs_packed_add: band lists differ')
        do t = 1, size(rpk_b,2)
            do q = 1, ncm
                do p = 1, ncb
                    if( abs(Rm(p,q)) < 1.d-8 ) cycle
                    rpk_acc(q,t) = rpk_acc(q,t) + real(Rm(p,q)) * rpk_b(p,t)
                end do
            end do
        end do
    end subroutine rotate_rhs_packed_add

    !> the coupled rho of fit B (pair-leading on the crop lattice) rotated into the merged frame and added
    subroutine rotate_rho4_packed_add( Rm, ncb, ncm, rho_b, rho_acc )
        real(dp), intent(in)    :: Rm(:,:)
        integer,  intent(in)    :: ncb, ncm
        real,     intent(in)    :: rho_b(:,:,:,:)
        real,     intent(inout) :: rho_acc(:,:,:,:)
        real, allocatable :: dout(:)
        integer :: i1, i2, i3, n1, n2, n3, npb, npm
        npb = (ncb*(ncb+1))/2
        npm = (ncm*(ncm+1))/2
        if( size(rho_b,1) /= npb .or. size(rho_acc,1) /= npm ) &
            &THROW_HARD('rotate_rho4_packed_add: packed leading extents do not match the ranks')
        n1 = size(rho_b,2); n2 = size(rho_b,3); n3 = size(rho_b,4)
        allocate(dout(npm))
        do i3 = 1, n3
            do i2 = 1, n2
                do i1 = 1, n1
                    call rotate_pair_voxel(Rm, ncb, ncm, rho_b(:,i1,i2,i3), dout)
                    rho_acc(:,i1,i2,i3) = rho_acc(:,i1,i2,i3) + dout
                end do
            end do
        end do
        deallocate(dout)
    end subroutine rotate_rho4_packed_add

    !> the packed pair kernels of fit B rotated into the merged frame and added, slot by slot
    subroutine rotate_rho_packed_add( Rm, ncb, ncm, rho_b, rho_acc )
        real(dp), intent(in)    :: Rm(:,:)
        integer,  intent(in)    :: ncb, ncm
        real,     intent(in)    :: rho_b(:,:)
        real,     intent(inout) :: rho_acc(:,:)
        real, allocatable :: dout(:)
        integer :: t, npb, npm
        npb = (ncb*(ncb+1))/2
        npm = (ncm*(ncm+1))/2
        if( size(rho_b,1) /= npb .or. size(rho_acc,1) /= npm ) &
            &THROW_HARD('rotate_rho_packed_add: packed leading extents do not match the ranks')
        if( size(rho_b,2) /= size(rho_acc,2) ) THROW_HARD('rotate_rho_packed_add: band lists differ')
        allocate(dout(npm))
        do t = 1, size(rho_b,2)
            call rotate_pair_voxel(Rm, ncb, ncm, rho_b(:,t), dout)
            rho_acc(:,t) = rho_acc(:,t) + dout
        end do
        deallocate(dout)
    end subroutine rotate_rho_packed_add

    !> Gate-4 self-check of the packed-pair rotation: R^T rho R through rotate_pair_voxel
    !! (the exact code path the merge uses) against a brute-force dense rotation evaluated
    !! independently via matmul, on deterministic pseudo-random values. Square K=2 case and a
    !! rectangular ncb=3 -> ncm=2 case. Runs at every merge; THROWs on mismatch.
    subroutine merge_pair_algebra_selfcheck
        real(dp) :: R2(2,2), R32(3,2)
        real     :: din2(3), dout2(3), din3(6), dout32(3)
        real(dp) :: D2(2,2), D3(3,3), B2(2,2), B32(2,2)
        real(dp) :: dmax
        integer  :: p, q, k
        ! K=2 square case
        k = 0
        do q = 1, 2
            do p = 1, q
                k = k + 1
                din2(pair_index(p,q)) = real(sin(3.7d0*real(k,dp) + 0.31d0))
            end do
        end do
        do q = 1, 2
            do p = 1, 2
                R2(p,q) = cos(1.3d0*real(p,dp) + 2.1d0*real(q,dp))
                D2(p,q) = real(din2(pair_index(min(p,q),max(p,q))),dp)
            end do
        end do
        call rotate_pair_voxel(R2, 2, 2, din2, dout2)
        B2 = matmul(transpose(R2), matmul(D2, R2))
        dmax = 0.d0
        do q = 1, 2
            do p = 1, q
                dmax = max(dmax, abs(real(dout2(pair_index(p,q)),dp) - B2(p,q)))
            end do
        end do
        ! rectangular 3 -> 2 case (rank-mismatched fits)
        k = 0
        do q = 1, 3
            do p = 1, q
                k = k + 1
                din3(pair_index(p,q)) = real(sin(2.9d0*real(k,dp) - 0.77d0))
            end do
        end do
        do q = 1, 3
            do p = 1, 3
                D3(p,q) = real(din3(pair_index(min(p,q),max(p,q))),dp)
            end do
        end do
        do q = 1, 2
            do p = 1, 3
                R32(p,q) = sin(0.9d0*real(p,dp) + 1.7d0*real(q,dp))
            end do
        end do
        call rotate_pair_voxel(R32, 3, 2, din3, dout32)
        B32 = matmul(transpose(R32), matmul(D3, R32))
        do q = 1, 2
            do p = 1, q
                dmax = max(dmax, abs(real(dout32(pair_index(p,q)),dp) - B32(p,q)))
            end do
        end do
        if( dmax > 1.d-6 ) THROW_HARD('paired-merge pair-algebra self-check FAILED: packed R^T &
            &rho R disagrees with the dense brute force')
        write(logfhandle,'(A,ES9.2)') '>>> FLEX_PCA MERGE pair-algebra self-check PASSED &
            &(K=2 + 3->2 rectangular): max |packed - dense| = ', dmax
        call flush(logfhandle)
    end subroutine merge_pair_algebra_selfcheck

    !> Mean-shaped (contrast) + background (+ optional pose-derivative) deflation of the merged
    !! solve output, mirroring the per-fit tail's per-iteration deflation statement for
    !! statement (simple_flex_pca_em_iter::fit_iter_finish, SIMPLE_COV_EM_DEFLATE block):
    !! consensus shells split in resolution, optional flat-in-mask background template
    !! (SIMPLE_COV_DEFLATE_BG), modified Gram-Schmidt with the relative-norm floor, then
    !! projection out of every merged component. Config comes from the same sources the tail
    !! uses: fit%l_deflate_mean / fit%vdfl / fit%dstep_ann (per-fit stage config; identical
    !! across fits) and the two env gates.
    subroutine merge_deflate_mean_shaped( params, fit, realvols, nvols )
        class(parameters), intent(inout) :: params
        type(probe_fit_t), intent(in)    :: fit
        type(image),       intent(inout) :: realvols(:)
        integer,           intent(in)    :: nvols
        type(image), allocatable :: dfl_basis(:)
        type(image) :: mvol_dfl
        real, pointer :: rm_dfl(:,:,:), rv_dfl(:,:,:)
        real(dp) :: mm_dfl, mv_dfl, rem_dfl, tot_dfl, mnorm_dfl
        real     :: res_lo, res_hi
        integer  :: ndfl, ndfl_sh, idfl, jdfl, nkeep_dfl, kfr_dfl(2), q
        logical  :: l_dfl_bg, l_dfl_dil
        if( .not. fit%l_deflate_mean ) return
        ndfl = max(1, fit%vdfl)
        l_dfl_bg = .true.
        if( l_dfl_bg ) ndfl = ndfl + 1
        ndfl_sh = ndfl
        if( l_dfl_bg )   ndfl_sh = ndfl_sh - 1
        ! the DILATION of the consensus, (x-c).grad rho, is the breathing mode (magnification and
        ! defocus scatter): reproducible, peripheral, and the first axis a focused basis takes
        ! (10028 continuum, 2026-09-09). Deflated by default; SIMPLE_COV_DEFLATE_DILATION=0 opts out.
        l_dfl_dil = .true.
        if( l_dfl_dil ) ndfl = ndfl + 1
        allocate(dfl_basis(ndfl))
        call mvol_dfl%read_and_crop(params%vols(1), params%smpd, params%box_crop, params%smpd_crop)
        kfr_dfl = covariance_kfromto(params)
        if( l_dfl_bg )then
            ! after the shells: uniform density under the soft spherical mask
            call dfl_basis(ndfl_sh+1)%copy(mvol_dfl)
            call dfl_basis(ndfl_sh+1)%get_rmat_ptr(rm_dfl)
            rm_dfl = 0.
            rm_dfl(1:params%box_crop, 1:params%box_crop, 1:params%box_crop) = 1.
            call flex_window_apply(dfl_basis(ndfl_sh+1), params)
        endif
        if( l_dfl_dil )then
            call dfl_basis(ndfl)%copy(mvol_dfl)
            call dilation_template(dfl_basis(ndfl), params%box_crop)
            call flex_window_apply(dfl_basis(ndfl), params)
        endif
        do idfl = 1, ndfl_sh
            call dfl_basis(idfl)%copy(mvol_dfl)
            if( ndfl_sh > 1 )then
                if( idfl == 1 )then
                    res_lo = 0.
                else
                    res_lo = fit%dstep_ann / &
                        &max(1., real(kfr_dfl(1)) + real(idfl-1)*real(max(1,kfr_dfl(2)-kfr_dfl(1)))/real(ndfl_sh))
                endif
                res_hi = fit%dstep_ann / &
                    &max(1., real(kfr_dfl(1)) + real(idfl)  *real(max(1,kfr_dfl(2)-kfr_dfl(1)))/real(ndfl_sh))
                call dfl_basis(idfl)%fft
                call dfl_basis(idfl)%bp(res_lo, res_hi, width=1.0)
                call dfl_basis(idfl)%ifft
                call dfl_basis(idfl)%get_rmat_ptr(rm_dfl)
                rm_dfl(params%box_crop+1:,:,:) = 0.
                call flex_window_apply(dfl_basis(idfl), params)
            endif
        end do
        call mvol_dfl%get_rmat_ptr(rv_dfl)
        mnorm_dfl = sqrt(sum(real(rv_dfl,dp)**2))
        ! modified Gram-Schmidt; a shell the band edges or the mask emptied drops out
        nkeep_dfl = 0
        do idfl = 1, ndfl
            call dfl_basis(idfl)%get_rmat_ptr(rm_dfl)
            do jdfl = 1, nkeep_dfl
                call dfl_basis(jdfl)%get_rmat_ptr(rv_dfl)
                mv_dfl = sum(real(rm_dfl,dp)*real(rv_dfl,dp))
                rm_dfl = rm_dfl - real(mv_dfl)*rv_dfl
            end do
            mm_dfl = sqrt(sum(real(rm_dfl,dp)*real(rm_dfl,dp)))
            if( mm_dfl <= 1.d-3*mnorm_dfl ) cycle
            rm_dfl    = rm_dfl / real(mm_dfl)
            nkeep_dfl = nkeep_dfl + 1
            if( nkeep_dfl /= idfl ) call dfl_basis(nkeep_dfl)%copy(dfl_basis(idfl))
        end do
        if( nkeep_dfl > 0 )then
            rem_dfl = 0.d0; tot_dfl = 0.d0
            do q = 1, nvols
                call realvols(q)%get_rmat_ptr(rv_dfl)
                tot_dfl = tot_dfl + sum(real(rv_dfl,dp)**2)
                do idfl = 1, nkeep_dfl
                    call dfl_basis(idfl)%get_rmat_ptr(rm_dfl)
                    mv_dfl  = sum(real(rm_dfl,dp)*real(rv_dfl,dp))   ! unit-norm already
                    rem_dfl = rem_dfl + mv_dfl*mv_dfl
                    rv_dfl  = rv_dfl - real(mv_dfl)*rm_dfl
                end do
            end do
            write(logfhandle,'(A,I0,A,F6.2,A)') '>>> FLEX_PCA MERGE mean-shaped deflation (rank ', &
                &nkeep_dfl,') removed ',100.d0*rem_dfl/max(tot_dfl,DTINY),' % of merged basis energy'
            call flush(logfhandle)
        endif
        do idfl = 1, ndfl
            call dfl_basis(idfl)%kill
        end do
        deallocate(dfl_basis)
        call mvol_dfl%kill
    end subroutine merge_deflate_mean_shaped

    !> Soft inner exclusion (focused heterogeneity), matching the per-fit tail's annulus.

end submodule simple_flex_pca_em_pairmerge
