!@descr: flex_pca EM: per-particle latent embedding with contrast fitting
submodule (simple_flex_pca_em) simple_flex_pca_em_embed
use simple_flex_pca_distr,  only: flex_pca_nparts
use simple_flex_pca_parts,  only: flex_pca_part_fname, write_embed_stats_part,&
    &reduce_embed_zhalf_parts, read_embed_stats_part
use simple_matcher_3Drec,   only: init_rec, cleanup_rec_buffers
use simple_matcher_ptcl_io, only: discrete_read_imgbatch, prepimgbatch
use simple_flex_reconstructor_latent_ops, only: project_fplanes_mean_basis
use simple_flex_projected_latent_model, only: prep_imgs4projected_model
implicit none
#include "simple_local_flags.inc"

contains

    !> Contrast-aware MAP embedding (supplement S.E, eqs S.14-S.15).
    module subroutine embed_latents_with_contrast( params, build, mean_rec, basis_recs, ncomp, eigvals, sig2_eff, &
        &pinds, nptcls, z, contrast, precision, resid_energy, resid_mean_energy, rho_out, stats_only, &
        &from_parts )
        class(parameters),   intent(inout) :: params
        type(builder),       intent(inout) :: build
        type(reconstructor), intent(inout) :: mean_rec
        type(reconstructor), intent(inout) :: basis_recs(ncomp)
        integer,             intent(in)    :: ncomp, pinds(:), nptcls
        real(dp),            intent(in)    :: eigvals(ncomp)
        real(dp),            intent(in)    :: sig2_eff
        real(dp),            intent(out)   :: z(nptcls,ncomp), contrast(nptcls)
        real(dp),            intent(out)   :: precision(ncomp,ncomp,nptcls)
        real(dp),            intent(out)   :: resid_energy(nptcls), resid_mean_energy(nptcls)
        ! per-component split-half reliability, exported so the STATE stage can order the latent
        ! path by how well each component is measured rather than by how much variance it carries.
        ! A high-variance, low-rho nuisance component otherwise dominates target placement.
        real(dp), optional,  intent(out)   :: rho_out(ncomp)
        !> worker: run the image pass over THIS part's particles, ship the sufficient statistics, stop
        logical,  optional,  intent(in)    :: stats_only
        !> master: skip the image pass entirely, gather the parts, run the coupled phase
        logical,  optional,  intent(in)    :: from_parts
        real(dp), parameter :: A_LO = 0.1d0, A_HI = 5.0d0
        type(fplane_type), allocatable :: fpls(:)
        type(fplane_type), allocatable :: basis_fpls(:,:), mean_fpl(:), data_fpl(:)
        type(ori), allocatable :: orientations(:)
        real(dp), allocatable :: prior(:), rho(:), Gcache(:,:,:), bcache(:,:), ccache(:,:)
        real(dp), allocatable :: zhalf(:,:,:), Ghf(:,:,:,:), bhf(:,:,:), chf(:,:,:)
        real(dp), allocatable :: Gpart(:,:,:), bpart(:,:), cpart(:,:)
        integer,  allocatable :: prows(:)
        integer :: ipart, pn_part
        real(dp), parameter   :: RHO_FLOOR = 1.d-3
        real(dp) :: rho_max, rrel
        logical :: l_relprior, l_stats_only, l_from_parts, l_polar_embed, l_devprep, l_cache_stats
        integer :: ihf
        integer :: batchlims(2), batchsz, ibatch, i, q, r, ithr, nthr, ia, row
        integer, allocatable :: nzeroG_thr(:), nzeroR_thr(:), nzeroZ_thr(:)   ! dead-basis counters
        real(dp) :: a, a_best, a_keep, a_num, a_den, e_yy, e_mm, best_res, res, aa, sig2
        real(dp), allocatable :: Ainv_th(:,:,:)     ! posterior covariance, for the tr(G A^-1) term
        integer  :: icm
        integer(timer_int_kind) :: t_phase
        real(dp), allocatable :: Gth(:,:,:), Ath(:,:,:), zth(:,:), zbest(:,:), cth(:,:), bth(:,:), myth(:)
        real(dp), allocatable :: Gtilth(:,:,:)   ! per-thread noise-whitened projected Gram
        real(dp), allocatable :: gwork(:,:,:), gvec(:,:,:), gev(:,:), gspec_thr(:,:), gspec(:)
        integer,  allocatable :: gcnt_thr(:)
        integer :: nrot_t, gcnt, blk_rp
        real(dp) :: gsum
        nthr = nthr_glob
        sig2 = max(sig2_eff, DTINY)      ! whitened-noise variance for the MAP shrinkage
        allocate(nzeroG_thr(nthr), nzeroR_thr(nthr), nzeroZ_thr(nthr), source=0)   ! dead-basis counters
        ! DEFAULT FLIPPED 2026-08-19: the PLAIN prior wins on AMI at n=3 replications
        ! (+0.024/+0.022/+0.012 across three operating points incl. production nparts=4),
        ! so the reliability prior is now OPT-IN via SIMPLE_COV_RELPRIOR=1. cov_env_int
        ! cannot express "off", hence the explicit >0 read (the recurring env-int trap).
        blk_rp = 0
        call cov_env_int('SIMPLE_COV_RELPRIOR', blk_rp)
        l_relprior = blk_rp > 0
        l_stats_only = .false.
        l_from_parts = .false.
        if( present(stats_only) ) l_stats_only = stats_only
        if( present(from_parts) ) l_from_parts = from_parts
        ! RELPRIOR=0 no longer forces the stage in-process. The distributed flow ships G/b/c
        ! and the split-half solves so the master can re-solve under the reliability-scaled
        ! prior; with the prior PLAIN the same re-solve applies 1/Gamma_q instead, so the
        ! caches are kept and only the rho computation is skipped. Measured motivation: two
        ! nparts=1-matched pairs read +0.037/+0.024 and +0.025/+0.022 (ARI/AMI) for the plain
        ! prior on the EM arm -- the rho^2 rescaling over-shrinks the reproducible directions
        ! 5-8 whose rho sits at 0.46-0.55.
        l_cache_stats = l_relprior .or. l_stats_only .or. l_from_parts
        allocate(prior(ncomp))
        if( l_cache_stats )then
            ! the reducing master never holds Gcache at nptcls: it reads one part's blocks at a time
            ! in the re-solve below
            if( l_from_parts )then
                allocate(Gcache(0,0,0), bcache(0,0), ccache(0,0))
            else
                allocate(Gcache(ncomp,ncomp,nptcls), bcache(ncomp,nptcls), ccache(ncomp,nptcls), source=0.d0)
            endif
            allocate(zhalf(nptcls,ncomp,2), source=0.d0)
            allocate(Ghf(ncomp,ncomp,2,nthr), bhf(ncomp,2,nthr), chf(ncomp,2,nthr), source=0.d0)
        endif
        do q = 1, ncomp
            prior(q) = 1.d0 / max(eigvals(q), DTINY)
        end do
        allocate(Gth(ncomp,ncomp,nthr), Ath(ncomp,ncomp,nthr), zth(ncomp,nthr), zbest(ncomp,nthr), &
            &cth(ncomp,nthr), bth(ncomp,nthr), myth(nthr), Gtilth(ncomp,ncomp,nthr))
        allocate(gwork(ncomp,ncomp,nthr), gvec(ncomp,ncomp,nthr), gev(ncomp,nthr))
        allocate(Ainv_th(ncomp,ncomp,nthr), source=0.d0)
        allocate(gspec_thr(ncomp,nthr), source=0.d0)
        allocate(gcnt_thr(nthr), source=0)
        allocate(basis_fpls(ncomp,nthr), mean_fpl(nthr), data_fpl(nthr), orientations(MAXIMGBATCHSZ))
        call mean_rec%expand_exp
        do q = 1, ncomp
            call basis_recs(q)%expand_exp
        end do
        call init_rec(params, build, MAXIMGBATCHSZ, fpls)
        call prepimgbatch(params, build, MAXIMGBATCHSZ)
        z = 0.d0; contrast = 1.d0; precision = 0.d0; resid_energy = 0.d0; resid_mean_energy = 0.d0
        t_phase = tic()
        ! master reducing parts: the workers have already paid the image pass below, so the master
        ! only needs their sufficient statistics to run the coupled phase (rho over every particle,
        ! then the re-solve)
        if( l_from_parts )then
            ! pass 1 only: zhalf and the per-particle scalars are all rho needs, and rho has to exist
            ! before any particle can be solved. The Gram blocks stay on disk until the re-solve
            ! reads them one part at a time.
            call reduce_embed_zhalf_parts(params, flex_pca_nparts(), pinds, contrast, resid_energy, &
                &resid_mean_energy, zhalf, nptcls, ncomp)
            goto 200
        endif
        write(logfhandle,'(A)') '>>> FLEX_PCA CONTRAST-AWARE EMBEDDING'
        call flush(logfhandle)
        ! POLAR PATH: the shared-direction former replaces the per-particle batch loop below. It
        ! fills the same sufficient statistics, so the reliability prior, the re-solve, the precision
        ! matrices and the worker part-write are all untouched. It needs Gcache, which only exists on
        ! the reliability-prior path -- the same constraint the distributed path already carries.
        l_polar_embed = cov_polar_embed_enabled()
        if( l_polar_embed )then
            if( .not. l_relprior ) THROW_HARD('polar embedding requires the reliability prior; &
                &unset SIMPLE_COV_RELPRIOR')
            call embed_accumulate_polar(params, build, mean_rec, basis_recs, ncomp, eigvals, sig2, &
                &pinds, nptcls, Gcache, bcache, ccache, zhalf, contrast, resid_energy, &
                &resid_mean_energy, prior, nzeroG_thr(1))
            nzeroG_thr(1) = 0
        endif
        if( .not. l_polar_embed )then
        call cov_dev_prep_start(params, build, l_devprep)
        do ibatch = 1, nptcls, MAXIMGBATCHSZ
            batchlims = [ibatch, min(nptcls, ibatch + MAXIMGBATCHSZ - 1)]
            batchsz   = batchlims(2) - batchlims(1) + 1
            call discrete_read_imgbatch(params, build, nptcls, pinds, batchlims)
            call prep_imgs4projected_model(params, build, batchsz, build%imgbatch(:batchsz), &
                &pinds(batchlims(1):batchlims(2)), fpls(:batchsz), mskrad=cov_image_mask_radius(params))
            do i = 1, batchsz
                call build%spproj_field%get_ori(pinds(batchlims(1)+i-1), orientations(i))
            end do
            !$omp parallel do default(shared) schedule(dynamic) proc_bind(close) &
            !$omp& private(i,ithr,q,r,ia,a,a_best,a_keep,a_num,a_den,icm,aa,e_yy,e_mm,best_res,res,row)
            do i = 1, batchsz
                if( orientations(i)%isstatezero() ) cycle
                ithr = omp_get_thread_num() + 1
                row  = batchlims(1) + i - 1
                call project_fplanes_mean_basis(mean_rec, basis_recs, orientations(i), fpls(i), &
                    &mean_fpl(ithr), basis_fpls(:,ithr), apply_ctf_amp=.true.)
                ! data plane = whitened observation (fpls(i)); mean_fpl = T mu ; basis = T U
                e_yy = real(cov_herm_inner(fpls(i), fpls(i)), dp)
                e_mm = real(cov_herm_inner(mean_fpl(ithr), mean_fpl(ithr)), dp)
                myth(ithr) = real(cov_herm_inner(mean_fpl(ithr), fpls(i)), dp)
                do q = 1, ncomp
                    bth(q,ithr) = real(cov_herm_inner(basis_fpls(q,ithr), fpls(i)), dp)      ! (TU)* y
                    cth(q,ithr) = real(cov_herm_inner(basis_fpls(q,ithr), mean_fpl(ithr)), dp) ! (TU)* T mu
                    do r = q, ncomp
                        Gth(q,r,ithr) = real(cov_herm_inner(basis_fpls(q,ithr), basis_fpls(r,ithr)), dp)
                        Gth(r,q,ithr) = Gth(q,r,ithr)
                    end do
                end do
                ! split-half sufficient statistics, for the reliability-weighted prior below
                if( l_relprior )then
                    do ihf = 1, 2
                        do q = 1, ncomp
                            bhf(q,ihf,ithr) = real(cov_herm_inner(basis_fpls(q,ithr), fpls(i), ihf), dp)
                            chf(q,ihf,ithr) = real(cov_herm_inner(basis_fpls(q,ithr), mean_fpl(ithr), ihf), dp)
                            do r = q, ncomp
                                Ghf(q,r,ihf,ithr) = real(cov_herm_inner(basis_fpls(q,ithr), basis_fpls(r,ithr), ihf), dp)
                                Ghf(r,q,ihf,ithr) = Ghf(q,r,ihf,ithr)
                            end do
                        end do
                    end do
                endif
                resid_mean_energy(row) = e_yy - 2.d0*myth(ithr) + e_mm                       ! contrast=1 mean residual
                ! Contrast (S.E): for each a on the grid solve the fixed-a MAP
                ! (a^2 G/sig2 + Gamma^-1) z = (a b - a^2 c)/sig2 and keep the a with the lowest residual.
                ! With COV_EMBED_CONTRAST_GRID off this is a single pass at a = 1.
                best_res = huge(1.d0)
                a_best   = 1.d0
                a_keep   = 1.d0     ! a NaN residual would leave this unset and hand garbage downstream
                    aa = a_best*a_best
                    Ath(:,:,ithr) = (aa/sig2)*Gth(:,:,ithr)
                    do q = 1, ncomp
                        Ath(q,q,ithr) = Ath(q,q,ithr) + prior(q)
                        zth(q,ithr)   = (a_best*bth(q,ithr) - aa*cth(q,ithr))/sig2
                    end do
                        call spd_solve_dp(Ath(:,:,ithr), zth(:,ithr), ncomp)
                    a = a_best
                    aa  = a*a
                    res = e_yy + aa*e_mm - 2.d0*a*myth(ithr) + quad_form(Gth(:,:,ithr), zth(:,ithr), ncomp)*aa
                    do q = 1, ncomp
                        res = res + 2.d0*aa*zth(q,ithr)*cth(q,ithr) - 2.d0*a*zth(q,ithr)*bth(q,ithr)
                    end do
                    res = res/sig2
                    do q = 1, ncomp
                        res = res + prior(q)*zth(q,ithr)*zth(q,ithr)
                    end do
                    if( res < best_res )then
                        best_res      = res
                        a_keep        = a
                        zbest(:,ithr) = zth(:,ithr)
                    endif
                a_best = a_keep
                ! projected-Gram spectrum on a subsample, for the conditioning report below
                if( mod(row, GRAM_DIAG_STRIDE) == 0 )then
                    gwork(:,:,ithr) = Gth(:,:,ithr)
                    call jacobi(gwork(:,:,ithr), ncomp, ncomp, gev(:,ithr), gvec(:,:,ithr), nrot_t)
                    call eigsrt(gev(:,ithr), gvec(:,:,ithr), ncomp, ncomp)
                    do q = 1, ncomp
                        gspec_thr(q,ithr) = gspec_thr(q,ithr) + max(gev(q,ithr), 0.d0)
                    end do
                    gcnt_thr(ithr) = gcnt_thr(ithr) + 1
                endif
                ! count the particles whose projected basis or rhs came out numerically dead
                if( maxval(abs(Gth(:,:,ithr))) <= 0.d0 ) nzeroG_thr(ithr) = nzeroG_thr(ithr) + 1
                if( maxval(abs(bth(:,ithr) - cth(:,ithr))) <= 0.d0 ) nzeroR_thr(ithr) = nzeroR_thr(ithr) + 1
                if( maxval(abs(zbest(:,ithr))) <= 0.d0 ) nzeroZ_thr(ithr) = nzeroZ_thr(ithr) + 1
                contrast(row)     = a_best
                z(row,:)          = zbest(:,ithr)
                resid_energy(row) = best_res
                aa = contrast(row)*contrast(row)
                Gtilth(:,:,ithr) = (aa/sig2)*Gth(:,:,ithr)
                call map_sampling_precision(Gtilth(:,:,ithr), prior, ncomp, precision(:,:,row))
                if( l_cache_stats )then
                    ! cache the sufficient statistics so the master's re-solve can run in closed
                    ! form with no second pass over the images. Cached whenever stats FLOW
                    ! (reliability prior on, or distributed either side): with RELPRIOR=0 under
                    ! distribution the master re-solves with the PLAIN prior, so skipping the
                    ! cache here would ship all-zero blocks and silently zero every latent.
                    Gcache(:,:,row) = Gth(:,:,ithr)
                    bcache(:,row)   = bth(:,ithr)
                    ccache(:,row)   = cth(:,ithr)
                    ! and the two half-data solves at the chosen contrast
                    do ihf = 1, 2
                        Ath(:,:,ithr) = (aa/sig2)*Ghf(:,:,ihf,ithr)
                        do q = 1, ncomp
                            Ath(q,q,ithr) = Ath(q,q,ithr) + prior(q)
                            zth(q,ithr)   = (contrast(row)*bhf(q,ihf,ithr) - aa*chf(q,ihf,ithr))/sig2
                        end do
                        call spd_solve_dp(Ath(:,:,ithr), zth(:,ithr), ncomp)
                        zhalf(row,:,ihf) = zth(:,ithr)
                    end do
                endif
            end do
            !$omp end parallel do
            if( batchlims(2) == nptcls .or. mod(batchlims(2), 5*MAXIMGBATCHSZ) == 0 )then
                write(logfhandle,'(A,I0,A,I0)') '>>> FLEX_PCA CONTRAST EMBED PARTICLES: ',batchlims(2),' / ',nptcls
                call flush(logfhandle)
            endif
        end do
        call cov_dev_prep_stop(l_devprep)
        endif
        ! the reducing master lands here rather than after the diagnostics, so it still frees the
        ! per-thread Gram workspace; gcnt_thr is all zero when no batch loop ran, so the per-particle
        ! spectrum report below skips itself
200     continue
        allocate(gspec(ncomp), source=0.d0)
        gcnt = sum(gcnt_thr)
        if( gcnt > 0 )then
            do q = 1, ncomp
                gspec(q) = sum(gspec_thr(q,:)) / real(gcnt, dp)
            end do
            write(logfhandle,'(A,I0,A)') '>>> FLEX_PCA projected-Gram spectrum (mean over ', &
                &gcnt,' particles), largest first:'
            write(logfhandle,'(A,10(1X,ES9.2))') '>>>   ', (gspec(q), q=1,min(10,ncomp))
            write(logfhandle,'(A,ES11.4,A,ES11.4)') '>>>   Gram condition number lam1/lamN = ', &
                &gspec(1)/max(gspec(ncomp),DTINY), '   lam1/lam5 = ', gspec(1)/max(gspec(min(5,ncomp)),DTINY)
            ! the participation ratio is only a rank if the spectrum is normalised first
            gsum = sum(gspec)
            if( gsum > 0.d0 )then
                gspec = gspec / gsum
                write(logfhandle,'(A,F7.3,A,I0)') '>>>   effective rank (participation ratio) = ', &
                    &1.d0 / max(sum(gspec**2), 1.d-300), '  of ', ncomp
            endif
            call flush(logfhandle)
        endif
        deallocate(gwork, gvec, gev, gspec_thr, gcnt_thr, gspec)
        deallocate(Ainv_th)
        ! The distribution of the fitted contrast is the first falsification test for the whole
        ! contrast story, and it should be readable without a script: if SIMPLE's per-particle
        ! normalisation has already absorbed the amplitude then this spread is ~0 and freeing a_i
        ! cannot help anything. The count at each bound says whether the bracket is binding, which
        ! would mean the bracket is doing the fitting rather than the data.
        write(logfhandle,'(A,I0,A,I0,A,I0,A,I0)') '>>> FLEX_PCA D2 zeroG=',sum(nzeroG_thr), &
            &' zeroRHS=',sum(nzeroR_thr),' zeroZ=',sum(nzeroZ_thr),' of nptcls=',nptcls
        write(logfhandle,'(A,F8.1)') '>>> FLEX_PCA CONTRAST EMBED SECONDS: ', toc(t_phase)
        call flush(logfhandle)
        ! worker: the image pass is done for this part's particles and everything below is coupled
        ! across all of them, so ship the sufficient statistics and leave the rest to the master
        if( l_stats_only )then
            call write_embed_stats_part(flex_pca_part_fname('embedstats', params%part, params%numlen), &
                &pinds, contrast, resid_energy, resid_mean_energy, Gcache, bcache, ccache, zhalf, &
                &nptcls, ncomp)
            if( present(rho_out) ) rho_out = 1.d0
            goto 900
        endif
        ! ---- RELIABILITY-WEIGHTED PRIOR ---- The plain prior precision sig2/Gamma_q hands the LARGEST
        ! eigenvalue the WEAKEST prior, so a high-variance but poorly measured component becomes
        ! near-unregularized least squares. Rescale each prior by the component's split-half reliability.
        if( l_cache_stats )then
            if( l_relprior )then
            allocate(rho(ncomp))
            ! optional export of the half-data solves for calibrating the per-particle error model:
            ! rho below collapses the pair to one correlation per component, whereas the variance of
            ! their difference measures the error directly. Off unless SIMPLE_COV_ZHALF is set.
            call write_zhalf_replicates(zhalf, prior, nptcls, ncomp)
            do q = 1, ncomp
                rho(q) = corr_dp(zhalf(:,q,1), zhalf(:,q,2), nptcls)
                rho(q) = max(0.d0, rho(q))
                rho(q) = 2.d0*rho(q) / (1.d0 + rho(q))            ! Spearman-Brown to full length
            end do
            ! Scale rho RELATIVE to the most reliable component, not absolutely: an absolute rho^2 shrinks
            ! the informative components as well and compresses the latent spread state placement needs.
            rho_max = maxval(rho)
            if( rho_max <= DTINY ) rho_max = 1.d0
            do q = 1, ncomp
                rrel = (rho(q)*rho(q)) / (rho_max*rho_max)
                prior(q) = 1.d0 / max(max(rrel, RHO_FLOOR) * eigvals(q), DTINY)
            end do
            ! all components, not the leading 10: rho drives state-target placement via comp_rho, so
            ! its ranking has to be checkable over the whole basis
            write(logfhandle,'(A)') '>>> FLEX_PCA split-half reliability per component (rho, corrected):'
            do q = 1, ncomp
                write(logfhandle,'(A,I3,A,F7.4,A,ES11.3,A,ES11.3)') '>>>   z',q,'  rho=',rho(q), &
                    &'  eigval=',eigvals(q),'  prior_precision=',prior(q)
            end do
            call flush(logfhandle)
            else
                write(logfhandle,'(A)') '>>> FLEX_PCA re-solving with the PLAIN prior &
                    &(RELPRIOR=0, distributed): rho not computed, no rescaling'
                call flush(logfhandle)
            endif
            ! re-solve every particle in closed form from the cached sufficient statistics. Two
            ! routes to the same arithmetic: in process the blocks are already in Gcache, while a
            ! reducing master streams them back one part at a time to keep its footprint flat.
            if( l_from_parts )then
                do ipart = 1, flex_pca_nparts()
                    call read_embed_stats_part(params, ipart, pinds, prows, Gpart, bpart, cpart, &
                        &pn_part, nptcls, ncomp)
                    !$omp parallel do default(shared) private(i,row,q,aa,ithr) schedule(static) proc_bind(close)
                    do i = 1, pn_part
                        ithr = omp_get_thread_num() + 1
                        row  = prows(i)
                        aa   = contrast(row)*contrast(row)
                        Ath(:,:,ithr) = (aa/sig2)*Gpart(:,:,i)
                        do q = 1, ncomp
                            Ath(q,q,ithr) = Ath(q,q,ithr) + prior(q)
                            zth(q,ithr)   = (contrast(row)*bpart(q,i) - aa*cpart(q,i))/sig2
                        end do
                        call spd_solve_dp(Ath(:,:,ithr), zth(:,ithr), ncomp)
                        z(row,:) = zth(:,ithr)
                        Gtilth(:,:,ithr) = (aa/sig2)*Gpart(:,:,i)
                        call map_sampling_precision(Gtilth(:,:,ithr), prior, ncomp, precision(:,:,row))
                    end do
                    !$omp end parallel do
                    deallocate(prows, Gpart, bpart, cpart)
                end do
            else
                !$omp parallel do default(shared) private(row,q,aa,ithr) schedule(static) proc_bind(close)
                do row = 1, nptcls
                    ithr = omp_get_thread_num() + 1
                    aa   = contrast(row)*contrast(row)
                    Ath(:,:,ithr) = (aa/sig2)*Gcache(:,:,row)
                    do q = 1, ncomp
                        Ath(q,q,ithr) = Ath(q,q,ithr) + prior(q)
                        zth(q,ithr)   = (contrast(row)*bcache(q,row) - aa*ccache(q,row))/sig2
                    end do
                    call spd_solve_dp(Ath(:,:,ithr), zth(:,ithr), ncomp)
                    z(row,:) = zth(:,ithr)
                    Gtilth(:,:,ithr) = (aa/sig2)*Gcache(:,:,row)
                    call map_sampling_precision(Gtilth(:,:,ithr), prior, ncomp, precision(:,:,row))
                end do
                !$omp end parallel do
            endif
            write(logfhandle,'(A)') '>>> FLEX_PCA latents re-solved from the cached statistics'
            call flush(logfhandle)
            if( present(rho_out) )then
                if( l_relprior )then
                    rho_out = rho
                else
                    rho_out = 1.d0
                endif
            endif
            if( allocated(rho) ) deallocate(rho)
            deallocate(Gcache, bcache, ccache, zhalf, Ghf, bhf, chf)
        else
            ! no split-half statistics available; treat every component as equally measured
            if( present(rho_out) ) rho_out = 1.d0
        endif
900     continue
        do i = 1, size(orientations)
            call orientations(i)%kill
        end do
        do ithr = 1, nthr
            call cleanup_plane(mean_fpl(ithr)); call cleanup_plane(data_fpl(ithr))
            do q = 1, ncomp
                call cleanup_plane(basis_fpls(q,ithr))
            end do
        end do
        call cleanup_rec_buffers(build, fpls)
        deallocate(prior, Gth, Ath, zth, zbest, cth, bth, myth, Gtilth, basis_fpls, mean_fpl, data_fpl, orientations)
        deallocate(nzeroG_thr, nzeroR_thr, nzeroZ_thr)
        if( allocated(Gcache) ) deallocate(Gcache, bcache, ccache, zhalf)
        if( allocated(Ghf)    ) deallocate(Ghf, bhf, chf)
        if( allocated(gwork)  ) deallocate(gwork, gvec, gev, gspec_thr, gcnt_thr)
    end subroutine embed_latents_with_contrast

end submodule simple_flex_pca_em_embed
