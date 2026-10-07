!@descr: flex_pca: the MAP embedding of every particle with its per-particle contrast and statistics
module simple_flex_pca_embed
!$ use omp_lib, only: omp_get_max_threads, omp_get_thread_num
use simple_core_module_api, only: del_file, dp, dtiny, eigsrt, fclose, file_exists, fileiochk, fopen, fplane_type, &
    &jacobi, logfhandle, maximgbatchsz, ori, simple_exception, simple_rename, string, tic, timer_int_kind, toc
use simple_flex_pca_records,              only: flex_fit_model, flex_selection, flex_latent
use simple_srch_sort_loc,                 only: locate_2
use simple_builder,                       only: builder
use simple_parameters,                    only: parameters
use simple_linalg,                        only: jacobi, eigsrt
use simple_flex_reconstructor_latent_ops, only: project_fplanes_mean_basis, planes_batch_load
use simple_ori,                           only: ori
use simple_flex_pca_rounds,               only: flex_pca_rounds
use simple_flex_pca_stages,               only: flex_stage_request, PCA_STAGE_EMBED
use simple_flex_pca_artifacts,            only: FLEX_PCA_PART_MAGIC
use simple_flex_pca_posterior,            only: quad_form, spd_solve_dp, map_sampling_precision
use simple_flex_pca_basis,                only: cov_herm_inner, cov_image_mask_radius
use simple_flex_pca_util,                 only: corr_dp
use simple_flex_pca_fit_types,            only: cleanup_plane
use simple_matcher_3Drec,                 only: init_rec, cleanup_rec_buffers
use simple_matcher_ptcl_io,               only: prepimgbatch
use simple_flex_pca_planes,               only: flex_plane_store
implicit none
private
#include "simple_local_flags.inc"

public :: embed_latents_with_contrast

logical, parameter :: COV_EMBED_CONTRAST_GRID = .false.
integer, parameter :: GRAM_DIAG_STRIDE        = 200   ! subsample for the projected-Gram spectrum
integer, parameter :: EMBED_STATS_VERSION     = 1

contains

    !> Contrast-aware MAP embedding (supplement S.E, eqs S.14-S.15).
    subroutine embed_latents_with_contrast( params, build, plane_store, model, sel, latent, rounds, &
        &stats_only, from_parts, l_zhalf )
        type(flex_fit_model),   intent(inout) :: model
        type(flex_selection),   intent(in)    :: sel
        type(flex_latent),      intent(inout) :: latent  !< z, contrast, precision, residual energies, comp_rho (and zhalf when l_zhalf) written
        !> keep the even/odd half solutions in latent%zhalf (the latent deconvolution calibrates its noise on them)
        logical,                intent(in)    :: l_zhalf
        class(flex_pca_rounds),  intent(inout) :: rounds
        class(parameters),       intent(inout) :: params
        type(builder),           intent(inout) :: build
        class(flex_plane_store), intent(inout) :: plane_store
        !> worker: run the image pass over THIS part's particles, ship the sufficient statistics, stop
        logical, optional,      intent(in)    :: stats_only
        !> master: skip the image pass entirely, gather the parts, run the coupled phase
        logical, optional,      intent(in)    :: from_parts
        type(fplane_type), allocatable :: fpls(:)
        type(fplane_type), allocatable :: basis_fpls(:,:), mean_fpl(:), data_fpl(:)
        type(ori), allocatable :: orientations(:)
        real(dp), allocatable :: prior(:), rho(:), Gcache(:,:,:), bcache(:,:), ccache(:,:)
        real(dp), allocatable :: zhalf(:,:,:), Ghf(:,:,:,:), bhf(:,:,:), chf(:,:,:)
        real(dp), allocatable :: myhf(:,:), emmhf(:,:)   ! per-half <Tmu,y>, <Tmu,Tmu> for the half-solve contrast
        real(dp) :: ah, aah
        real(dp), allocatable :: Gpart(:,:,:), bpart(:,:), cpart(:,:)
        integer,  allocatable :: prows(:)
        type(string) :: part_fname
        integer :: ipart, pn_part
        real(dp), parameter :: RHO_FLOOR = 1.d-3
        real(dp) :: rho_max, rrel
        logical :: l_relprior, l_stats_only, l_from_parts, l_cache_stats, l_pcache
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
        if( allocated(latent%z) ) deallocate(latent%z)
        if( allocated(latent%contrast) ) deallocate(latent%contrast)
        if( allocated(latent%precision) ) deallocate(latent%precision)
        if( allocated(latent%resid_energy) ) deallocate(latent%resid_energy)
        if( allocated(latent%resid_mean_energy) ) deallocate(latent%resid_mean_energy)
        if( allocated(latent%comp_rho) ) deallocate(latent%comp_rho)
        allocate(latent%z(sel%nptcls,model%ncomp), latent%contrast(sel%nptcls), latent%precision(model%ncomp,model%ncomp,sel%nptcls))
        allocate(latent%resid_energy(sel%nptcls), latent%resid_mean_energy(sel%nptcls), latent%comp_rho(model%ncomp))
        nthr = omp_get_max_threads()
        sig2 = max(model%sig2_eff, DTINY)      ! whitened-noise variance for the MAP shrinkage
        allocate(nzeroG_thr(nthr), nzeroR_thr(nthr), nzeroZ_thr(nthr), source=0)   ! dead-basis counters
        ! The PLAIN 1/Gamma prior is the prior: it won on AMI at n=3 replications
        ! (+0.024/+0.022/+0.012 across three operating points incl. production nparts=4), and the
        ! paired path's axis reliability is already in eigvals as the merge weight 2c/(1+c).
        l_relprior = .false.
        l_stats_only = .false.
        l_from_parts = .false.
        if( present(stats_only) ) l_stats_only = stats_only
        if( present(from_parts) ) l_from_parts = from_parts
        ! the caches flow whenever stats do; under the plain prior the re-solve applies 1/Gamma_q, no rho
        l_cache_stats = l_relprior .or. l_stats_only .or. l_from_parts .or. l_zhalf
        if( l_relprior )then
            write(logfhandle,'(A)') '>>> FLEX_PCA split-half solves: per-half fitted contrast'
            call flush(logfhandle)
        endif
        allocate(prior(model%ncomp))
        if( l_cache_stats )then
            ! the reducing master never holds Gcache at nptcls: it reads one part's blocks at a time
            ! in the re-solve below
            if( l_from_parts )then
                allocate(Gcache(0,0,0), bcache(0,0), ccache(0,0))
            else
                allocate(Gcache(model%ncomp,model%ncomp,sel%nptcls), bcache(model%ncomp,sel%nptcls), ccache(model%ncomp,sel%nptcls), source=0.d0)
            endif
            allocate(zhalf(sel%nptcls,model%ncomp,2), source=0.d0)
            allocate(Ghf(model%ncomp,model%ncomp,2,nthr), bhf(model%ncomp,2,nthr), chf(model%ncomp,2,nthr), source=0.d0)
            allocate(myhf(2,nthr), emmhf(2,nthr), source=0.d0)
        endif
        do q = 1, model%ncomp
            prior(q) = 1.d0 / max(model%eigvals(q), DTINY)
        end do
        allocate(Gth(model%ncomp,model%ncomp,nthr), Ath(model%ncomp,model%ncomp,nthr), zth(model%ncomp,nthr), zbest(model%ncomp,nthr), &
            &cth(model%ncomp,nthr), bth(model%ncomp,nthr), myth(nthr), Gtilth(model%ncomp,model%ncomp,nthr))
        allocate(gwork(model%ncomp,model%ncomp,nthr), gvec(model%ncomp,model%ncomp,nthr), gev(model%ncomp,nthr))
        allocate(Ainv_th(model%ncomp,model%ncomp,nthr), source=0.d0)
        allocate(gspec_thr(model%ncomp,nthr), source=0.d0)
        allocate(gcnt_thr(nthr), source=0)
        allocate(basis_fpls(model%ncomp,nthr), mean_fpl(nthr), data_fpl(nthr), orientations(MAXIMGBATCHSZ))
        call model%mean_rec%expand_exp
        do q = 1, model%ncomp
            call model%basis_recs(q)%expand_exp
        end do
        ! one read path for every pass of the run: the downscaled cache when it is in use
        l_pcache = plane_store%cache_in_use()
        call init_rec(params, build, MAXIMGBATCHSZ, fpls, cropped=l_pcache)
        if( l_pcache )then
            call prepimgbatch(params, build, MAXIMGBATCHSZ, box=params%box_crop, smpd=params%smpd_crop)
        else
            call prepimgbatch(params, build, MAXIMGBATCHSZ)
        endif
        latent%z = 0.d0; latent%contrast = 1.d0; latent%precision = 0.d0; latent%resid_energy = 0.d0; latent%resid_mean_energy = 0.d0
        t_phase = tic()
        ! master reducing parts: the workers have already paid the image pass below, so the master
        ! only needs their sufficient statistics to run the coupled phase (rho over every particle,
        ! then the re-solve)
        if( l_from_parts )then
            ! pass 1 only: zhalf and the per-particle scalars are all rho needs, and rho has to exist
            ! before any particle can be solved. The Gram blocks stay on disk until the re-solve
            ! reads them one part at a time.
            call reduce_embed_zhalf_parts(params, rounds, sel%pinds, latent%contrast, latent%resid_energy, &
                &latent%resid_mean_energy, zhalf, sel%nptcls, model%ncomp)
            goto 200
        endif
        write(logfhandle,'(A)') '>>> FLEX_PCA CONTRAST-AWARE EMBEDDING'
        call flush(logfhandle)
        do ibatch = 1, sel%nptcls, MAXIMGBATCHSZ
            batchlims = [ibatch, min(sel%nptcls, ibatch + MAXIMGBATCHSZ - 1)]
            batchsz   = batchlims(2) - batchlims(1) + 1
            call planes_batch_load(plane_store, params, build, sel%nptcls, sel%pinds, batchlims, fpls, &
                &cov_image_mask_radius(params), l_pcache)
            do i = 1, batchsz
                call build%spproj_field%get_ori(sel%pinds(batchlims(1)+i-1), orientations(i))
            end do
            !$omp parallel do default(shared) schedule(dynamic) proc_bind(close) &
            !$omp& private(i,ithr,q,r,ia,a,a_best,a_keep,a_num,a_den,icm,aa,e_yy,e_mm,best_res,res,row,ah,aah)
            do i = 1, batchsz
                if( orientations(i)%isstatezero() ) cycle
                ithr = omp_get_thread_num() + 1
                row  = batchlims(1) + i - 1
                call project_fplanes_mean_basis(model%mean_rec, model%basis_recs, orientations(i), fpls(i), &
                    &mean_fpl(ithr), basis_fpls(:,ithr), apply_ctf_amp=.true.)
                ! data plane = whitened observation (fpls(i)); mean_fpl = T mu ; basis = T U
                e_yy = real(cov_herm_inner(fpls(i), fpls(i)), dp)
                e_mm = real(cov_herm_inner(mean_fpl(ithr), mean_fpl(ithr)), dp)
                myth(ithr) = real(cov_herm_inner(mean_fpl(ithr), fpls(i)), dp)
                do q = 1, model%ncomp
                    bth(q,ithr) = real(cov_herm_inner(basis_fpls(q,ithr), fpls(i)), dp)      ! (TU)* y
                    cth(q,ithr) = real(cov_herm_inner(basis_fpls(q,ithr), mean_fpl(ithr)), dp) ! (TU)* T mu
                    do r = q, model%ncomp
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
                        do q = 1, model%ncomp
                            bhf(q,ihf,ithr) = real(cov_herm_inner(basis_fpls(q,ithr), fpls(i), ihf), dp)
                            chf(q,ihf,ithr) = real(cov_herm_inner(basis_fpls(q,ithr), mean_fpl(ithr), ihf), dp)
                            do r = q, model%ncomp
                                Ghf(q,r,ihf,ithr) = real(cov_herm_inner(basis_fpls(q,ithr), basis_fpls(r,ithr), ihf), dp)
                                Ghf(r,q,ihf,ithr) = Ghf(q,r,ihf,ithr)
                            end do
                        end do
                    end do
                endif
                latent%resid_mean_energy(row) = e_yy - 2.d0*myth(ithr) + e_mm                       ! contrast=1 mean residual
                ! Contrast (S.E): for each a on the grid solve the fixed-a MAP
                ! (a^2 G/sig2 + Gamma^-1) z = (a b - a^2 c)/sig2 and keep the a with the lowest residual.
                ! With COV_EMBED_CONTRAST_GRID off this is a single pass at a = 1.
                best_res = huge(1.d0)
                a_best   = 1.d0
                a_keep   = 1.d0     ! a NaN residual would leave this unset and hand garbage downstream
                aa = a_best*a_best
                Ath(:,:,ithr) = (aa/sig2)*Gth(:,:,ithr)
                do q = 1, model%ncomp
                    Ath(q,q,ithr) = Ath(q,q,ithr) + prior(q)
                    zth(q,ithr)   = (a_best*bth(q,ithr) - aa*cth(q,ithr))/sig2
                end do
                call spd_solve_dp(Ath(:,:,ithr), zth(:,ithr), model%ncomp)
                a = a_best
                aa  = a*a
                res = e_yy + aa*e_mm - 2.d0*a*myth(ithr) + quad_form(Gth(:,:,ithr), zth(:,ithr), model%ncomp)*aa
                do q = 1, model%ncomp
                    res = res + 2.d0*aa*zth(q,ithr)*cth(q,ithr) - 2.d0*a*zth(q,ithr)*bth(q,ithr)
                end do
                res = res/sig2
                do q = 1, model%ncomp
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
                    call jacobi(gwork(:,:,ithr), model%ncomp, model%ncomp, gev(:,ithr), gvec(:,:,ithr), nrot_t)
                    call eigsrt(gev(:,ithr), gvec(:,:,ithr), model%ncomp, model%ncomp)
                    do q = 1, model%ncomp
                        gspec_thr(q,ithr) = gspec_thr(q,ithr) + max(gev(q,ithr), 0.d0)
                    end do
                    gcnt_thr(ithr) = gcnt_thr(ithr) + 1
                endif
                ! count the particles whose projected basis or rhs came out numerically dead
                if( maxval(abs(Gth(:,:,ithr))) <= 0.d0 ) nzeroG_thr(ithr) = nzeroG_thr(ithr) + 1
                if( maxval(abs(bth(:,ithr) - cth(:,ithr))) <= 0.d0 ) nzeroR_thr(ithr) = nzeroR_thr(ithr) + 1
                if( maxval(abs(zbest(:,ithr))) <= 0.d0 ) nzeroZ_thr(ithr) = nzeroZ_thr(ithr) + 1
                latent%contrast(row)     = a_best
                latent%z(row,:)          = zbest(:,ithr)
                latent%resid_energy(row) = best_res
                aa = latent%contrast(row)*latent%contrast(row)
                Gtilth(:,:,ithr) = (aa/sig2)*Gth(:,:,ithr)
                call map_sampling_precision(Gtilth(:,:,ithr), prior, model%ncomp, latent%precision(:,:,row))
                if( l_cache_stats )then
                    ! cache the sufficient statistics so the master's re-solve can run in closed
                    ! form with no second pass over the images. Cached whenever stats FLOW
                    ! (reliability prior on, or distributed either side): with RELPRIOR=0 under
                    ! distribution the master re-solves with the PLAIN prior, so skipping the
                    ! cache here would ship all-zero blocks and silently zero every latent.
                    Gcache(:,:,row) = Gth(:,:,ithr)
                    bcache(:,row)   = bth(:,ithr)
                    ccache(:,row)   = cth(:,ithr)
                    ! and the two half-data solves, each at its OWN fitted contrast (the delivered z keeps a=1):
                    ! at a=1 the residual (a_i-1)*Tmu enters the halves with opposite signs (basis deflated vs the mean)
                    do ihf = 1, 2
                        ah = myhf(ihf,ithr) / max(emmhf(ihf,ithr), DTINY)
                        ah = max(0.1d0, min(5.d0, ah))
                        aah = ah*ah
                        Ath(:,:,ithr) = (aah/sig2)*Ghf(:,:,ihf,ithr)
                        do q = 1, model%ncomp
                            Ath(q,q,ithr) = Ath(q,q,ithr) + prior(q)
                            zth(q,ithr)   = (ah*bhf(q,ihf,ithr) - aah*chf(q,ihf,ithr))/sig2
                        end do
                        ! the half solves: their disagreement calibrates the DATA noise
                        call spd_solve_dp(Ath(:,:,ithr), zth(:,ithr), model%ncomp)
                        zhalf(row,:,ihf) = zth(:,ithr)
                    end do
                endif
            end do
            !$omp end parallel do
            if( batchlims(2) == sel%nptcls .or. mod(batchlims(2), 5*MAXIMGBATCHSZ) == 0 )then
                write(logfhandle,'(A,I0,A,I0)') '>>> FLEX_PCA CONTRAST EMBED PARTICLES: ',batchlims(2),' / ',sel%nptcls
                call flush(logfhandle)
            endif
        end do
        ! the reducing master lands here rather than after the diagnostics, so it still frees the
        ! per-thread Gram workspace; gcnt_thr is all zero when no batch loop ran, so the per-particle
        ! spectrum report below skips itself
200     allocate(gspec(model%ncomp), source=0.d0)
        gcnt = sum(gcnt_thr)
        if( gcnt > 0 )then
            do q = 1, model%ncomp
                gspec(q) = sum(gspec_thr(q,:)) / real(gcnt, dp)
            end do
            write(logfhandle,'(A,I0,A)') '>>> FLEX_PCA projected-Gram spectrum (mean over ', &
                &gcnt,' particles), largest first:'
            write(logfhandle,'(A,10(1X,ES9.2))') '>>>   ', (gspec(q), q=1,min(10,model%ncomp))
            write(logfhandle,'(A,ES11.4,A,ES11.4)') '>>>   Gram condition number lam1/lamN = ', &
                &gspec(1)/max(gspec(model%ncomp),DTINY), '   lam1/lam5 = ', gspec(1)/max(gspec(min(5,model%ncomp)),DTINY)
            ! the participation ratio is only a rank if the spectrum is normalised first
            gsum = sum(gspec)
            if( gsum > 0.d0 )then
                gspec = gspec / gsum
                write(logfhandle,'(A,F7.3,A,I0)') '>>>   effective rank (participation ratio) = ', &
                    &1.d0 / max(sum(gspec**2), 1.d-300), '  of ', model%ncomp
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
            &' zeroRHS=',sum(nzeroR_thr),' zeroZ=',sum(nzeroZ_thr),' of nptcls=',sel%nptcls
        write(logfhandle,'(A,F8.1)') '>>> FLEX_PCA CONTRAST EMBED SECONDS: ', toc(t_phase)
        call flush(logfhandle)
        ! worker: the image pass is done for this part's particles and everything below is coupled
        ! across all of them, so ship the sufficient statistics and leave the rest to the master
        if( l_stats_only )then
            part_fname = rounds%part_fname('embedstats', params%part, params%numlen)
            call write_embed_stats_part(part_fname, &
                &sel%pinds, latent%contrast, latent%resid_energy, latent%resid_mean_energy, Gcache, bcache, ccache, zhalf, &
                &sel%nptcls, model%ncomp)
            call part_fname%kill
            latent%comp_rho = 1.d0
            goto 900
        endif
        ! ---- RELIABILITY-WEIGHTED PRIOR ---- The plain prior precision sig2/Gamma_q hands the LARGEST
        ! eigenvalue the WEAKEST prior, so a high-variance but poorly measured component becomes
        ! near-unregularized least squares. Rescale each prior by the component's split-half reliability.
        if( l_cache_stats )then
            if( l_relprior )then
                allocate(rho(model%ncomp))
                do q = 1, model%ncomp
                    rho(q) = corr_dp(zhalf(:,q,1), zhalf(:,q,2), sel%nptcls)
                    rho(q) = max(0.d0, rho(q))
                    rho(q) = 2.d0*rho(q) / (1.d0 + rho(q))            ! Spearman-Brown to full length
                end do
                ! Scale rho RELATIVE to the most reliable component, not absolutely: an absolute rho^2 shrinks
                ! the informative components as well and compresses the latent spread state placement needs.
                rho_max = maxval(rho)
                if( rho_max <= DTINY ) rho_max = 1.d0
                do q = 1, model%ncomp
                    rrel = (rho(q)*rho(q)) / (rho_max*rho_max)
                    prior(q) = 1.d0 / max(max(rrel, RHO_FLOOR) * model%eigvals(q), DTINY)
                end do
                ! all components, not the leading 10: rho drives state-target placement via comp_rho, so
                ! its ranking has to be checkable over the whole basis
                write(logfhandle,'(A)') '>>> FLEX_PCA split-half reliability per component (rho, corrected):'
                do q = 1, model%ncomp
                    write(logfhandle,'(A,I3,A,F7.4,A,ES11.3,A,ES11.3)') '>>>   z',q,'  rho=',rho(q), &
                        &'  eigval=',model%eigvals(q),'  prior_precision=',prior(q)
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
                    call read_embed_stats_part(params, rounds, ipart, sel%pinds, prows, Gpart, bpart, cpart, &
                        &pn_part, sel%nptcls, model%ncomp)
                    !$omp parallel do default(shared) private(i,row,q,aa,ithr) schedule(static) proc_bind(close)
                    do i = 1, pn_part
                        ithr = omp_get_thread_num() + 1
                        row  = prows(i)
                        aa   = latent%contrast(row)*latent%contrast(row)
                        Ath(:,:,ithr) = (aa/sig2)*Gpart(:,:,i)
                        do q = 1, model%ncomp
                            Ath(q,q,ithr) = Ath(q,q,ithr) + prior(q)
                            zth(q,ithr)   = (latent%contrast(row)*bpart(q,i) - aa*cpart(q,i))/sig2
                        end do
                        call spd_solve_dp(Ath(:,:,ithr), zth(:,ithr), model%ncomp)
                        latent%z(row,:) = zth(:,ithr)
                        Gtilth(:,:,ithr) = (aa/sig2)*Gpart(:,:,i)
                        call map_sampling_precision(Gtilth(:,:,ithr), prior, model%ncomp, latent%precision(:,:,row))
                    end do
                    !$omp end parallel do
                    deallocate(prows, Gpart, bpart, cpart)
                end do
            else
                !$omp parallel do default(shared) private(row,q,aa,ithr) schedule(static) proc_bind(close)
                do row = 1, sel%nptcls
                    ithr = omp_get_thread_num() + 1
                    aa   = latent%contrast(row)*latent%contrast(row)
                    Ath(:,:,ithr) = (aa/sig2)*Gcache(:,:,row)
                    do q = 1, model%ncomp
                        Ath(q,q,ithr) = Ath(q,q,ithr) + prior(q)
                        zth(q,ithr)   = (latent%contrast(row)*bcache(q,row) - aa*ccache(q,row))/sig2
                    end do
                    call spd_solve_dp(Ath(:,:,ithr), zth(:,ithr), model%ncomp)
                    latent%z(row,:) = zth(:,ithr)
                    Gtilth(:,:,ithr) = (aa/sig2)*Gcache(:,:,row)
                    call map_sampling_precision(Gtilth(:,:,ithr), prior, model%ncomp, latent%precision(:,:,row))
                end do
                !$omp end parallel do
            endif
            write(logfhandle,'(A)') '>>> FLEX_PCA latents re-solved from the cached statistics'
            call flush(logfhandle)
            if( l_relprior )then
                latent%comp_rho = rho
            else
                latent%comp_rho = 1.d0
            endif
            if( allocated(rho) ) deallocate(rho)
            if( l_zhalf )then
                if( allocated(latent%zhalf) ) deallocate(latent%zhalf)
                allocate(latent%zhalf(sel%nptcls,model%ncomp,2)); latent%zhalf = zhalf
            endif
            deallocate(Gcache, bcache, ccache, zhalf, Ghf, bhf, chf, myhf, emmhf)
        else
            ! no split-half statistics available; treat every component as equally measured
            latent%comp_rho = 1.d0
        endif
900     do i = 1, size(orientations)
            call orientations(i)%kill
        end do
        do ithr = 1, nthr
            call cleanup_plane(mean_fpl(ithr)); call cleanup_plane(data_fpl(ithr))
            do q = 1, model%ncomp
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

    !> One part's embedding sufficient statistics (per-particle scalars, split-half latents, G/b/c
    !! blocks). The master reduces zhalf over every particle, forms the prior once, then re-solves
    !! each part's rows from these blocks without touching an image.
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
    subroutine reduce_embed_zhalf_parts( params, rounds, gpinds, contrast, resid_energy, &
        &resid_mean_energy, zhalf, nptcls, ncomp )
        class(parameters), intent(in)    :: params
        class(flex_pca_rounds), intent(in) :: rounds
        integer,           intent(in)    :: gpinds(:), nptcls, ncomp
        real(dp),          intent(inout) :: contrast(:), resid_energy(:), resid_mean_energy(:)
        real(dp),          intent(inout) :: zhalf(:,:,:)
        integer,  allocatable :: ppinds(:)
        real(dp), allocatable :: pc(:), pre(:), prme(:), pzh(:,:,:)
        type(string) :: fname
        integer :: ipart, funit, io_stat, header(4), pn, i, hit, nfilled
        integer(timer_int_kind) :: t_red
        t_red   = tic()
        nfilled = 0
        do ipart = 1, rounds%nparts()
            fname = rounds%part_fname('embedstats', ipart, params%numlen)
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
                hit = locate_2(gpinds, nptcls, ppinds(i))
                if( hit > 0 )then
                    if( gpinds(hit) /= ppinds(i) ) hit = 0
                endif
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
        write(logfhandle,'(A,I0,A,I0,A,F8.1)') '>>> FLEX_PCA reduced embed-stats (zhalf) parts=',rounds%nparts(), &
            &'  particles=',nfilled,'  seconds=',toc(t_red)
        call flush(logfhandle)
    end subroutine reduce_embed_zhalf_parts

    !> Pass 2: one part's Gram blocks and its global row indices, so the caller can re-solve just
    !! those particles and free the buffer before reading the next. Deletes the part file.
    subroutine read_embed_stats_part( params, rounds, ipart, gpinds, rows, Gpart, bpart, cpart, pn, nptcls, ncomp )
        class(parameters),     intent(in)  :: params
        class(flex_pca_rounds), intent(in) :: rounds
        integer,               intent(in)  :: ipart, gpinds(:), nptcls, ncomp
        integer,  allocatable, intent(out) :: rows(:)
        real(dp), allocatable, intent(out) :: Gpart(:,:,:), bpart(:,:), cpart(:,:)
        integer,               intent(out) :: pn
        integer,  allocatable :: ppinds(:)
        real(dp), allocatable :: skip3(:), pzh(:,:,:)
        type(string) :: fname
        integer :: funit, io_stat, header(4), i, hit
        fname = rounds%part_fname('embedstats', ipart, params%numlen)
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
            hit = locate_2(gpinds, nptcls, ppinds(i))
            if( hit > 0 )then
                if( gpinds(hit) /= ppinds(i) ) hit = 0
            endif
            if( hit < 1 ) THROW_HARD('embed-stats part carries a particle not in the global set')
            rows(i) = hit
        end do
        deallocate(ppinds, skip3, pzh)
        call del_file(fname)
        call fname%kill
    end subroutine read_embed_stats_part

end module simple_flex_pca_embed
