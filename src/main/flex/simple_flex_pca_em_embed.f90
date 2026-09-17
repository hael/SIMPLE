!@descr: flex_pca EM: per-particle latent embedding with contrast fitting
submodule (simple_flex_pca_em) simple_flex_pca_em_embed
use simple_matcher_3Drec,   only: init_rec, cleanup_rec_buffers
use simple_matcher_ptcl_io, only: discrete_read_imgbatch, prepimgbatch
use simple_flex_reconstructor_latent_ops, only: project_fplanes_mean_basis
use simple_flex_pca_planes, only: planes_batch_load
use simple_flex_pca_plane_cache, only: plane_cache_in_use
use simple_flex_pca_util, only: cov_env_flag_on
implicit none
#include "simple_local_flags.inc"

contains

    !> Contrast-aware MAP embedding (supplement S.E, eqs S.14-S.15).
    module subroutine embed_latents_with_contrast( params, build, mean_rec, basis_recs, ncomp, eigvals, sig2_eff, &
        &pinds, nptcls, z, contrast, precision, resid_energy, resid_mean_energy, rho_out, stats_only, &
        &from_parts, rounds, zhalf_out)
        class(flex_pca_rounds), intent(inout) :: rounds
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
        !> the even/odd Fourier-half solutions of every particle (nptcls,ncomp,2), for the noise
        !! calibration of the latent deconvolution; requesting them keeps the sufficient statistics
        real(dp), optional,  intent(out)   :: zhalf_out(:,:,:)
        real(dp), parameter :: A_LO = 0.1d0, A_HI = 5.0d0
        type(fplane_type), allocatable :: fpls(:)
        type(fplane_type), allocatable :: basis_fpls(:,:), mean_fpl(:), data_fpl(:)
        type(ori), allocatable :: orientations(:)
        real(dp), allocatable :: prior(:), rho(:), Gcache(:,:,:), bcache(:,:), ccache(:,:)
        real(dp), allocatable :: zhalf(:,:,:), Ghf(:,:,:,:), bhf(:,:,:), chf(:,:,:)
        real(dp), allocatable :: myhf(:,:), emmhf(:,:)   ! per-half <Tmu,y>, <Tmu,Tmu> for the half-solve contrast
        real(dp) :: ah, aah
        logical  :: l_halfcontrast
        real(dp), allocatable :: Gpart(:,:,:), bpart(:,:), cpart(:,:)
        integer,  allocatable :: prows(:)
        integer :: ipart, pn_part
        real(dp), parameter   :: RHO_FLOOR = 1.d-3
        real(dp) :: rho_max, rrel
        logical :: l_relprior, l_stats_only, l_from_parts, l_devprep, l_cache_stats, l_pcache
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
        integer :: nrot_t, gcnt
        real(dp) :: gsum
        nthr = omp_get_max_threads()
        sig2 = max(sig2_eff, DTINY)      ! whitened-noise variance for the MAP shrinkage
        allocate(nzeroG_thr(nthr), nzeroR_thr(nthr), nzeroZ_thr(nthr), source=0)   ! dead-basis counters
        ! The PLAIN 1/Gamma prior is the prior: it won on AMI at n=3 replications
        ! (+0.024/+0.022/+0.012 across three operating points incl. production nparts=4), and the
        ! paired path's axis reliability is already in eigvals as the merge weight 2c/(1+c).
        l_relprior = .false.
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
        l_cache_stats = l_relprior .or. l_stats_only .or. l_from_parts .or. present(zhalf_out)
        ! Per-half fitted contrast in the split-half solves (SIMPLE_COV_HALF_CONTRAST=0 opts out).
        ! The delivered z keeps a=1. WHY: the basis is deflated against the mean, so the FULL-plane
        ! <TU,Tmu> nearly cancels while its two half-plane parts do not (they are +-D_i, pose
        ! dependent, and Tmu dwarfs TUz). With a fixed a=1 the residual carries (a_i-1)*Tmu, which
        ! therefore enters the two half solves with OPPOSITE signs and drives the split-half
        ! correlation negative (measured -0.13..-0.24 on 10028, 2026-09-07). Fitting a per half
        ! removes that term from the reliability estimate without touching the delivered latents.
        l_halfcontrast = .true.
        if( l_relprior )then
            if( l_halfcontrast )then
                write(logfhandle,'(A)') '>>> FLEX_PCA split-half solves: per-half fitted contrast'
            else
                write(logfhandle,'(A)') '>>> FLEX_PCA split-half solves: contrast fixed to the full-plane value'
            endif
            call flush(logfhandle)
        endif
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
            allocate(myhf(2,nthr), emmhf(2,nthr), source=0.d0)
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
        ! one read path for every pass of the run: the downscaled cache when it is in use
        l_pcache = plane_cache_in_use(params, build)
        call init_rec(params, build, MAXIMGBATCHSZ, fpls, cropped=l_pcache)
        if( l_pcache )then
            call prepimgbatch(params, build, MAXIMGBATCHSZ, box=params%box_crop, smpd=params%smpd_crop)
        else
            call prepimgbatch(params, build, MAXIMGBATCHSZ)
        endif
        z = 0.d0; contrast = 1.d0; precision = 0.d0; resid_energy = 0.d0; resid_mean_energy = 0.d0
        t_phase = tic()
        ! master reducing parts: the workers have already paid the image pass below, so the master
        ! only needs their sufficient statistics to run the coupled phase (rho over every particle,
        ! then the re-solve)
        if( l_from_parts )then
            ! pass 1 only: zhalf and the per-particle scalars are all rho needs, and rho has to exist
            ! before any particle can be solved. The Gram blocks stay on disk until the re-solve
            ! reads them one part at a time.
            call reduce_embed_zhalf_parts(params, rounds%nparts(), pinds, contrast, resid_energy, &
                &resid_mean_energy, zhalf, nptcls, ncomp)
            goto 200
        endif
        write(logfhandle,'(A)') '>>> FLEX_PCA CONTRAST-AWARE EMBEDDING'
        call flush(logfhandle)
        call cov_dev_prep_start(params, build, l_devprep)
        do ibatch = 1, nptcls, MAXIMGBATCHSZ
            batchlims = [ibatch, min(nptcls, ibatch + MAXIMGBATCHSZ - 1)]
            batchsz   = batchlims(2) - batchlims(1) + 1
            call planes_batch_load(params, build, nptcls, pinds, batchlims, fpls, &
                &cov_image_mask_radius(params), l_pcache)
            do i = 1, batchsz
                call build%spproj_field%get_ori(pinds(batchlims(1)+i-1), orientations(i))
            end do
            !$omp parallel do default(shared) schedule(dynamic) proc_bind(close) &
            !$omp& private(i,ithr,q,r,ia,a,a_best,a_keep,a_num,a_den,icm,aa,e_yy,e_mm,best_res,res,row,ah,aah)
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
                ! split-half sufficient statistics, for the reliability-weighted prior below and for
                ! the deconvolution's noise calibration (a worker always forms them: its part may be
                ! reduced by a master that needs either)
                if( l_cache_stats )then
                    do ihf = 1, 2
                        myhf(ihf,ithr)  = real(cov_herm_inner(mean_fpl(ithr), fpls(i), ihf), dp)
                        emmhf(ihf,ithr) = real(cov_herm_inner(mean_fpl(ithr), mean_fpl(ithr), ihf), dp)
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
                    ! and the two half-data solves, each at its OWN fitted contrast (see l_halfcontrast)
                    do ihf = 1, 2
                        if( l_halfcontrast )then
                            ah = myhf(ihf,ithr) / max(emmhf(ihf,ithr), DTINY)
                            ah = max(0.1d0, min(5.d0, ah))
                        else
                            ah = contrast(row)
                        endif
                        aah = ah*ah
                        Ath(:,:,ithr) = (aah/sig2)*Ghf(:,:,ihf,ithr)
                        do q = 1, ncomp
                            Ath(q,q,ithr) = Ath(q,q,ithr) + prior(q)
                            zth(q,ithr)   = (ah*bhf(q,ihf,ithr) - aah*chf(q,ihf,ithr))/sig2
                        end do
                        ! the half solves: their disagreement calibrates the DATA noise
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
                do ipart = 1, rounds%nparts()
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
            if( present(zhalf_out) ) zhalf_out(1:nptcls,1:ncomp,1:2) = zhalf
            deallocate(Gcache, bcache, ccache, zhalf, Ghf, bhf, chf, myhf, emmhf)
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
        if( allocated(myhf) ) deallocate(myhf, emmhf)
        if( allocated(Ghf)    ) deallocate(Ghf, bhf, chf)
        if( allocated(gwork)  ) deallocate(gwork, gvec, gev, gspec_thr, gcnt_thr)
    end subroutine embed_latents_with_contrast

    !> One part's embedding sufficient statistics.
    !!
    !! The embedding is not a clean partition, so a part cannot ship finished latents.
    !! The reliability prior comes from
    !! rho(q) = corr(zhalf(:,q,1), zhalf(:,q,2)) over every particle, and each particle's final z is
    !! re-solved against it; per-part rho would solve the parts against different priors.
    !!
    !! So a part ships what it can compute independently -- the per-particle sufficient statistics
    !! from the image pass plus its own rows of the split-half latents -- and the master does the
    !! coupled arithmetic: reduce zhalf, form rho and the prior once, then re-solve. The re-solve
    !! touches no images, so the stage that actually costs is the part that distributes.
    subroutine write_embed_stats_part( fname, pinds, contrast, resid_energy, resid_mean_energy, &
        &Gcache, bcache, ccache, zhalf, nptcls, ncomp )
        class(string), intent(in) :: fname
        integer,       intent(in) :: pinds(:), nptcls, ncomp
        real(dp),      intent(in) :: contrast(:), resid_energy(:), resid_mean_energy(:)
        real(dp),      intent(in) :: Gcache(:,:,:), bcache(:,:), ccache(:,:), zhalf(:,:,:)
        type(string) :: tmp_fname
        integer :: funit, io_stat, header(4)
        header = [FLEX_PCA_PART_MAGIC, EMBED_STATS_VERSION, nptcls, ncomp]
        tmp_fname = fname//'.tmp'
        call fopen(funit, file=tmp_fname, access='STREAM', action='WRITE', status='REPLACE', iostat=io_stat)
        call fileiochk('write_embed_stats_part; open', io_stat)
        write(funit, iostat=io_stat) header
        call fileiochk('write_embed_stats_part; header', io_stat)
        ! zhalf before Gcache, deliberately: the master forms rho from the split-half latents before
        ! it can solve anything, and consumes the much larger Gram blocks one part at a time. Small
        ! arrays first lets pass 1 stop reading at zhalf.
        write(funit, iostat=io_stat) pinds(1:nptcls)
        write(funit, iostat=io_stat) contrast(1:nptcls), resid_energy(1:nptcls), resid_mean_energy(1:nptcls)
        write(funit, iostat=io_stat) zhalf(1:nptcls,1:ncomp,1:2)
        write(funit, iostat=io_stat) Gcache(1:ncomp,1:ncomp,1:nptcls)
        write(funit, iostat=io_stat) bcache(1:ncomp,1:nptcls), ccache(1:ncomp,1:nptcls)
        call fileiochk('write_embed_stats_part; payload', io_stat)
        call fclose(funit)
        call simple_rename(tmp_fname, fname)
        call tmp_fname%kill
    end subroutine write_embed_stats_part

    !> Pass 1: gather what the master needs before it can solve anything -- the split-half latents
    !! that form rho, plus the per-particle scalars. Reads no Gram blocks and deletes nothing; the
    !! files are consumed by read_embed_stats_part below.
    !!
    !! Splitting the reduce in two keeps the master's footprint flat in dataset size: holding every
    !! part's Gcache at once costs ncomp^2 doubles per particle, one part at a time costs that
    !! divided by nparts.
    subroutine reduce_embed_zhalf_parts( params, nparts, gpinds, contrast, resid_energy, &
        &resid_mean_energy, zhalf, nptcls, ncomp )
        class(parameters), intent(in)    :: params
        integer,           intent(in)    :: nparts, gpinds(:), nptcls, ncomp
        real(dp),          intent(inout) :: contrast(:), resid_energy(:), resid_mean_energy(:)
        real(dp),          intent(inout) :: zhalf(:,:,:)
        integer,  allocatable :: ppinds(:)
        real(dp), allocatable :: pc(:), pre(:), prme(:), pzh(:,:,:)
        type(string) :: fname
        integer :: ipart, funit, io_stat, header(4), pn, i, hit, nfilled
        integer(timer_int_kind) :: t_red
        t_red   = tic()
        nfilled = 0
        do ipart = 1, nparts
            fname = flex_pca_part_fname('embedstats', ipart, params%numlen)
            if( .not. file_exists(fname) ) THROW_HARD('missing embed-stats part: '//fname%to_char())
            call fopen(funit, file=fname, access='STREAM', action='READ', status='OLD', iostat=io_stat)
            call fileiochk('reduce_embed_zhalf_parts; open '//fname%to_char(), io_stat)
            read(funit, iostat=io_stat) header
            call fileiochk('reduce_embed_zhalf_parts; header', io_stat)
            if( header(1) /= FLEX_PCA_PART_MAGIC ) THROW_HARD('bad embed-stats part magic')
            if( header(2) /= EMBED_STATS_VERSION ) THROW_HARD('bad embed-stats part version')
            if( header(4) /= ncomp               ) THROW_HARD('embed-stats part ncomp mismatch')
            pn = header(3)
            allocate(ppinds(pn), pc(pn), pre(pn), prme(pn), pzh(pn,ncomp,2))
            read(funit, iostat=io_stat) ppinds
            read(funit, iostat=io_stat) pc, pre, prme
            read(funit, iostat=io_stat) pzh
            call fileiochk('reduce_embed_zhalf_parts; payload', io_stat)
            call fclose(funit)
            ! match on pinds rather than assuming a contiguous layout, so a part boundary that does
            ! not line up cannot silently misplace rows
            do i = 1, pn
                hit = binsrch_int(gpinds, nptcls, ppinds(i))
                if( hit < 1 ) THROW_HARD('embed-stats part carries a particle not in the global set')
                contrast(hit)          = pc(i)
                resid_energy(hit)      = pre(i)
                resid_mean_energy(hit) = prme(i)
                zhalf(hit,:,:)         = pzh(i,:,:)
                nfilled = nfilled + 1
            end do
            deallocate(ppinds, pc, pre, prme, pzh)
            call fname%kill
        end do
        if( nfilled /= nptcls ) THROW_HARD('embed-stats parts did not cover every particle')
        write(logfhandle,'(A,I0,A,I0,A,F8.1)') '>>> FLEX_PCA reduced embed-stats (zhalf) parts=',nparts, &
            &'  particles=',nfilled,'  seconds=',toc(t_red)
        call flush(logfhandle)
    end subroutine reduce_embed_zhalf_parts

    !> Pass 2: one part's Gram blocks and its global row indices, so the caller can re-solve just
    !! those particles and free the buffer before reading the next. Deletes the part file.
    subroutine read_embed_stats_part( params, ipart, gpinds, rows, Gpart, bpart, cpart, pn, nptcls, ncomp )
        class(parameters),     intent(in)  :: params
        integer,               intent(in)  :: ipart, gpinds(:), nptcls, ncomp
        integer,  allocatable, intent(out) :: rows(:)
        real(dp), allocatable, intent(out) :: Gpart(:,:,:), bpart(:,:), cpart(:,:)
        integer,               intent(out) :: pn
        integer,  allocatable :: ppinds(:)
        real(dp), allocatable :: skip3(:), pzh(:,:,:)
        type(string) :: fname
        integer :: funit, io_stat, header(4), i, hit
        fname = flex_pca_part_fname('embedstats', ipart, params%numlen)
        if( .not. file_exists(fname) ) THROW_HARD('missing embed-stats part: '//fname%to_char())
        call fopen(funit, file=fname, access='STREAM', action='READ', status='OLD', iostat=io_stat)
        call fileiochk('read_embed_stats_part; open '//fname%to_char(), io_stat)
        read(funit, iostat=io_stat) header
        if( header(1) /= FLEX_PCA_PART_MAGIC ) THROW_HARD('bad embed-stats part magic')
        if( header(2) /= EMBED_STATS_VERSION ) THROW_HARD('bad embed-stats part version')
        if( header(4) /= ncomp               ) THROW_HARD('embed-stats part ncomp mismatch')
        pn = header(3)
        allocate(ppinds(pn), skip3(3*pn), pzh(pn,ncomp,2))
        allocate(Gpart(ncomp,ncomp,pn), bpart(ncomp,pn), cpart(ncomp,pn), rows(pn))
        read(funit, iostat=io_stat) ppinds
        read(funit, iostat=io_stat) skip3          ! contrast, resid_energy, resid_mean_energy: pass 1
        read(funit, iostat=io_stat) pzh            ! zhalf: pass 1
        read(funit, iostat=io_stat) Gpart
        read(funit, iostat=io_stat) bpart, cpart
        call fileiochk('read_embed_stats_part; payload', io_stat)
        call fclose(funit)
        do i = 1, pn
            hit = binsrch_int(gpinds, nptcls, ppinds(i))
            if( hit < 1 ) THROW_HARD('embed-stats part carries a particle not in the global set')
            rows(i) = hit
        end do
        deallocate(ppinds, skip3, pzh)
        call del_file(fname)
        call fname%kill
    end subroutine read_embed_stats_part

    !> Index of key in an ascending array, 0 if absent. pinds arrive in project order, so the scatter
    !! below can binary search rather than scan linearly, which is O(nptcls) per row.
    pure integer function binsrch_int( arr, n, key ) result( pos )
        integer, intent(in) :: n, arr(n), key
        integer :: lo, hi, mid
        pos = 0
        lo  = 1
        hi  = n
        do while( lo <= hi )
            mid = (lo + hi)/2
            if( arr(mid) == key )then
                pos = mid
                return
            else if( arr(mid) < key )then
                lo = mid + 1
            else
                hi = mid - 1
            endif
        end do
    end function binsrch_int

end submodule simple_flex_pca_em_embed
