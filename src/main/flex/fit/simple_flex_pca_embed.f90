!@descr: flex_pca: the MAP embedding of every particle with its per-particle contrast and statistics, and basis composition across runs
module simple_flex_pca_embed
use simple_core_module_api
use simple_flex_pca_records, only: flex_fit_model, flex_selection, flex_latent
use simple_srch_sort_loc, only: locate_2
use simple_builder, only: builder
use simple_image, only: image
use simple_parameters, only: parameters
use simple_reconstructor, only: reconstructor
use simple_linalg, only: jacobi, eigsrt
use simple_flex_reconstructor_latent_ops, only: project_fplanes_mean_basis
use simple_ori, only: ori
use simple_flex_pca_rounds, only: flex_pca_rounds
use simple_flex_pca_stages, only: flex_stage_request, PCA_STAGE_EMBED
use simple_flex_pca_artifacts, only: flex_pca_part_fname, FLEX_PCA_PART_MAGIC
use simple_flex_pca_run_types, only: flex_run_settings
use simple_flex_pca_posterior, only: quad_form, spd_solve_dp, map_sampling_precision
use simple_flex_pca_basis, only: cov_herm_inner, save_probe_state, load_probe_state, cov_image_mask_radius,&
    &basis_recs_from_images
use simple_flex_pca_util, only: corr_dp
use simple_flex_pca_fit_types, only: cleanup_plane
use simple_matcher_3Drec, only: init_rec, cleanup_rec_buffers
use simple_matcher_ptcl_io, only: prepimgbatch
use simple_flex_pca_planes, only: planes_batch_load
use simple_flex_pca_plane_cache, only: plane_cache_in_use
use simple_imghead, only: find_ldim_nptcls
use simple_flex_pca_deconv, only: calibrate_noise_scale, noise_and_projection, signal_subspace
implicit none
private
#include "simple_local_flags.inc"

public :: embed_latents_with_contrast, compose_basis_from_runs, compose_cut_reembed

logical,  parameter :: COV_EMBED_CONTRAST_GRID = .false.
integer,  parameter :: GRAM_DIAG_STRIDE = 200   ! subsample for the projected-Gram spectrum
integer, parameter :: EMBED_STATS_VERSION = 1

!! SIMPLE_COV_COMPOSE=<dir>[,<dir>...]: finished runs; per run polished > merged > plain namespace.
!! Columns are Fourier-padded to box_crop (zero beyond their band), Gram-Schmidt'ed coarse box first;
!! residuals below COMPOSE_R2_FLOOR drop, prior variances follow the rescaling. compose_cut_reembed
!! then cuts the union to its signal subspace (SIMPLE_COV_COMPOSE_CUT=0 skips it).

!> a finished run's delivered basis: namespace, rank, grid, prior variances
type compose_src_t
    character(len=LONGSTRLEN) :: dir  = ''
    character(len=32)         :: pfx  = ''
    character(len=64)         :: meta = ''
    integer  :: box   = 0
    integer  :: ncomp = 0
    real(dp) :: sig2  = 0.d0
    real(dp), allocatable :: eigvals(:)
end type compose_src_t

!> a composed column: a fine column whose residual after the coarse projection is below this
!! fraction of its norm is already spanned and is dropped
real(dp), parameter :: COMPOSE_R2_FLOOR = 0.05d0

contains

    !> Contrast-aware MAP embedding (supplement S.E, eqs S.14-S.15).
    subroutine embed_latents_with_contrast( params, build, model, sel, latent, rounds, stats_only, from_parts, l_zhalf )
        type(flex_fit_model), intent(inout) :: model
        type(flex_selection), intent(in)    :: sel
        type(flex_latent),    intent(inout) :: latent  !< z, contrast, precision, residual energies, comp_rho (and zhalf when l_zhalf) written
        !> keep the even/odd half solutions in latent%zhalf (the latent deconvolution calibrates its noise on them)
        logical,              intent(in)    :: l_zhalf
        class(flex_pca_rounds), intent(inout) :: rounds
        class(parameters),   intent(inout) :: params
        type(builder),       intent(inout) :: build
        !> worker: run the image pass over THIS part's particles, ship the sufficient statistics, stop
        logical,  optional,  intent(in)    :: stats_only
        !> master: skip the image pass entirely, gather the parts, run the coupled phase
        logical,  optional,  intent(in)    :: from_parts
        type(fplane_type), allocatable :: fpls(:)
        type(fplane_type), allocatable :: basis_fpls(:,:), mean_fpl(:), data_fpl(:)
        type(ori), allocatable :: orientations(:)
        real(dp), allocatable :: prior(:), rho(:), Gcache(:,:,:), bcache(:,:), ccache(:,:)
        real(dp), allocatable :: zhalf(:,:,:), Ghf(:,:,:,:), bhf(:,:,:), chf(:,:,:)
        real(dp), allocatable :: myhf(:,:), emmhf(:,:)   ! per-half <Tmu,y>, <Tmu,Tmu> for the half-solve contrast
        real(dp) :: ah, aah
        real(dp), allocatable :: Gpart(:,:,:), bpart(:,:), cpart(:,:)
        integer,  allocatable :: prows(:)
        integer :: ipart, pn_part
        real(dp), parameter   :: RHO_FLOOR = 1.d-3
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
        l_pcache = plane_cache_in_use(params, build)
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
            call reduce_embed_zhalf_parts(params, rounds%nparts(), sel%pinds, latent%contrast, latent%resid_energy, &
                &latent%resid_mean_energy, zhalf, sel%nptcls, model%ncomp)
            goto 200
        endif
        write(logfhandle,'(A)') '>>> FLEX_PCA CONTRAST-AWARE EMBEDDING'
        call flush(logfhandle)
        do ibatch = 1, sel%nptcls, MAXIMGBATCHSZ
            batchlims = [ibatch, min(sel%nptcls, ibatch + MAXIMGBATCHSZ - 1)]
            batchsz   = batchlims(2) - batchlims(1) + 1
            call planes_batch_load(params, build, sel%nptcls, sel%pinds, batchlims, fpls, &
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
200     continue
        allocate(gspec(model%ncomp), source=0.d0)
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
            call write_embed_stats_part(flex_pca_part_fname('embedstats', params%part, params%numlen), &
                &sel%pinds, latent%contrast, latent%resid_energy, latent%resid_mean_energy, Gcache, bcache, ccache, zhalf, &
                &sel%nptcls, model%ncomp)
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
                    call read_embed_stats_part(params, ipart, sel%pinds, prows, Gpart, bpart, cpart, &
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
900     continue
        do i = 1, size(orientations)
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
    subroutine compose_basis_from_runs( params, cfg, build, model )
        type(flex_fit_model), intent(inout) :: model   !< the composed basis, prior variances, rank and noise level
        class(parameters),   intent(inout) :: params
        type(flex_run_settings), intent(in)    :: cfg
        type(builder),       intent(inout) :: build
        type(compose_src_t), allocatable :: srcs(:)
        type(image),         allocatable :: imgs(:)
        type(image)  :: raw, fine
        type(string) :: fname
        real(dp), allocatable :: VV(:,:), qcol(:), ev(:), s2(:), r2(:)
        integer,  allocatable :: c_src(:), c_comp(:), c_box(:)
        real,     pointer     :: rmat(:,:,:) => null()
        real(dp) :: pj
        real     :: smpd_s
        integer  :: nsrc, isrc, ntot, bc, bs, nvox, j, k, q, kept
        integer  :: u
        if( allocated(model%basis_recs) ) deallocate(model%basis_recs)
        if( allocated(model%eigvals) ) deallocate(model%eigvals)
        bc   = params%box_crop
        nvox = bc*bc*bc
        if( .not. cfg%l_compose ) THROW_HARD('compose_basis_from_runs: SIMPLE_COV_COMPOSE is not set')
        call parse_sources(cfg%compose_dirs, srcs, nsrc)
        do isrc = 1, nsrc
            call resolve_source(params, srcs(isrc))
        end do
        call sort_sources_by_box(srcs, nsrc)
        ntot = sum(srcs(1:nsrc)%ncomp)
        write(logfhandle,'(A,I0,A,I0,A,I0)') '>>> FLEX_PCA COMPOSE: ', nsrc, ' finished run(s), ', ntot, &
            &' delivered components, composing box=', bc
        do isrc = 1, nsrc
            write(logfhandle,'(A,I0,A,I0,A,I0,A,ES10.3,A,A,A,A)') '>>>   run ', isrc, ': box ', &
                &srcs(isrc)%box, '  ncomp ', srcs(isrc)%ncomp, '  sig2 ', srcs(isrc)%sig2, '  namespace ', &
                &trim(srcs(isrc)%pfx), '  ', trim(srcs(isrc)%dir)
        end do
        call flush(logfhandle)
        allocate(VV(nvox,ntot), ev(ntot), s2(ntot), r2(ntot), c_src(ntot), c_comp(ntot), c_box(ntot))
        ev = 0.d0; s2 = 0.d0; r2 = 0.d0; c_src = 0; c_comp = 0; c_box = 0
        ! ---- load every delivered column on its own grid, weight, pad to the composing grid ----
        j = 0
        do isrc = 1, nsrc
            bs       = srcs(isrc)%box
            smpd_s   = params%smpd_crop*real(bc)/real(bs)
            call raw%new([bs,bs,bs], smpd_s)
            do q = 1, srcs(isrc)%ncomp
                j = j + 1
                c_src(j) = isrc; c_comp(j) = q; c_box(j) = bs
                ev(j)    = srcs(isrc)%eigvals(q)
                fname = trim(srcs(isrc)%dir)//trim(srcs(isrc)%pfx)//int2str_pad(q,3)//MRC_EXT
                call raw%read(fname)
                call fname%kill
                ! Fourier-pad to the composing grid: coefficients copied verbatim, zero beyond the
                ! source band (no antialiasing window -- the band edge is the source's own low-pass)
                if( bs == bc )then
                    call fine%copy(raw)
                else
                    call raw%fft
                    call fine%new([bc,bc,bc], params%smpd_crop)
                    call raw%pad(fine, antialiasing=.false.)
                    call fine%ifft
                    call raw%ifft
                endif
                call fine%set_smpd(params%smpd_crop)
                call fine%get_rmat_ptr(rmat)
                VV(:,j) = reshape(real(rmat(1:bc,1:bc,1:bc),dp), [nvox])
                call fine%kill
            end do
            call raw%kill
        end do
        call flush(logfhandle)
        ! ---- coarse-first Gram-Schmidt on the composing grid; the prior variance follows the
        ! column's rescaling (z_j u_j = (z_j s_j) u_j/s_j) and its retained fraction ----
        allocate(qcol(nvox))
        kept = 0
        do j = 1, ntot
            s2(j) = sum(VV(:,j)**2)
            if( s2(j) <= 1.d-30 )then
                write(logfhandle,'(A,I0,A,I0,A)') '>>>   run ', c_src(j), ' comp ', c_comp(j), ': empty column; dropped'
                r2(j) = 0.d0
                cycle
            endif
            qcol = VV(:,j)/sqrt(s2(j))
            do k = 1, kept
                pj   = dot_product(qcol, VV(:,k))
                qcol = qcol - pj*VV(:,k)
            end do
            r2(j) = sum(qcol**2)
            if( r2(j) < COMPOSE_R2_FLOOR )then
                write(logfhandle,'(A,I0,A,I0,A,F6.3,A)') '>>>   run ', c_src(j), ' comp ', c_comp(j), &
                    &': residual fraction ', real(r2(j)), ' already spanned by the coarser columns; dropped'
                cycle
            endif
            kept = kept + 1
            VV(:,kept)   = qcol/sqrt(r2(j))
            ev(kept)     = ev(j)*s2(j)*r2(j)
            c_src(kept)  = c_src(j); c_comp(kept) = c_comp(j); c_box(kept) = c_box(j)
            s2(kept) = s2(j); r2(kept) = r2(j)
        end do
        deallocate(qcol)
        if( kept < 2 ) THROW_HARD('compose_basis_from_runs: fewer than 2 independent columns survived')
        model%ncomp    = kept
        model%sig2_eff = srcs(maxloc(srcs(1:nsrc)%box, dim=1))%sig2   ! the finest run's noise level: its band is the composing band
        allocate(model%eigvals(model%ncomp))
        model%eigvals = ev(1:model%ncomp)
        write(logfhandle,'(A,I0,A,I0,A,ES10.3)') '>>> FLEX_PCA COMPOSE: ', model%ncomp, ' of ', ntot, &
            &' columns kept (orthonormal, coarse first); sig2_eff=', model%sig2_eff
        write(logfhandle,'(A)') '>>>   column  run  comp  box   s2 (fine-grid norm2)   residual   prior var'
        do j = 1, model%ncomp
            write(logfhandle,'(A,I5,1X,I4,1X,I5,1X,I5,3X,ES12.4,7X,F8.4,3X,ES12.4)') '>>>   ', j, &
                &c_src(j), c_comp(j), c_box(j), s2(j), r2(j), model%eigvals(j)
        end do
        call flush(logfhandle)
        ! ---- realise: real-space images -> column reconstructors; the composed basis goes to disk
        ! under the plain namespace (the distributed embed workers load flex_pca_pc*.mrc) ----
        allocate(imgs(model%ncomp))
        do j = 1, model%ncomp
            call imgs(j)%new([bc,bc,bc], params%smpd_crop)
            call imgs(j)%get_rmat_ptr(rmat)
            rmat(1:bc,1:bc,1:bc) = real(reshape(VV(:,j), [bc,bc,bc]))
            fname = 'flex_pca_pc'//int2str_pad(j,3)//MRC_EXT
            call imgs(j)%write(fname, del_if_exists=.true.)
            call fname%kill
        end do
        deallocate(VV)
        call basis_recs_from_images(params, build, imgs, model%ncomp, model%basis_recs)
        do j = 1, model%ncomp
            call imgs(j)%kill
        end do
        deallocate(imgs)
        call save_probe_state(model%ncomp, model%eigvals, model%sig2_eff)
        ! provenance table
        open(newunit=u, file='flex_pca_compose.txt', status='replace', action='write')
        write(u,'(A)') '# composed basis: column  run  source_comp  source_box  norm2_fine  residual_frac  prior_var'
        do j = 1, model%ncomp
            write(u,'(I5,1X,I4,1X,I5,1X,I5,1X,ES14.6,1X,F9.5,1X,ES14.6)') j, c_src(j), c_comp(j), &
                &c_box(j), s2(j), r2(j), model%eigvals(j)
        end do
        do isrc = 1, nsrc
            write(u,'(A,I0,A,I0,A,A,A,A)') '# run ', isrc, ' box ', srcs(isrc)%box, ' namespace ', &
                &trim(srcs(isrc)%pfx), ' dir ', trim(srcs(isrc)%dir)
        end do
        close(u)
        do isrc = 1, nsrc
            if( allocated(srcs(isrc)%eigvals) ) deallocate(srcs(isrc)%eigvals)
        end do
        deallocate(srcs, ev, s2, r2, c_src, c_comp, c_box)
    end subroutine compose_basis_from_runs

    !> comma-separated run directories -> sources (trailing slash normalised)
    subroutine parse_sources( str, srcs, nsrc )
        character(len=*), intent(in) :: str
        type(compose_src_t), allocatable, intent(out) :: srcs(:)
        integer, intent(out) :: nsrc
        character(len=LONGSTRLEN) :: tok
        integer :: i, i0, n, l
        n = 1
        do i = 1, len(str)
            if( str(i:i) == ',' ) n = n + 1
        end do
        allocate(srcs(n))
        nsrc = 0
        i0   = 1
        do i = 1, len(str) + 1
            if( i > len(str) )then
                tok = adjustl(str(i0:len(str)))
            else if( str(i:i) == ',' )then
                tok = adjustl(str(i0:i-1))
            else
                cycle
            endif
            i0 = i + 1
            l  = len_trim(tok)
            if( l < 1 ) cycle
            if( tok(l:l) /= '/' ) tok = tok(1:l)//'/'
            nsrc = nsrc + 1
            srcs(nsrc)%dir = tok
        end do
        if( nsrc < 1 ) THROW_HARD('SIMPLE_COV_COMPOSE names no run directory')
    end subroutine parse_sources

    !> pick the namespace the run embedded with (polished > merged > plain), read its meta,
    !! count its eigenvolumes, read the grid
    subroutine resolve_source( params, src )
        class(parameters),   intent(in)    :: params
        type(compose_src_t), intent(inout) :: src
        type(string) :: fn
        integer :: ldim(3), nvol, q, nfound
        if( file_exists(trim(src%dir)//'flex_pca_polished_pc001.mrc') .and. &
           &file_exists(trim(src%dir)//'flex_pca_probe_polished.txt') )then
            src%pfx = 'flex_pca_polished_pc'; src%meta = 'flex_pca_probe_polished.txt'
        else if( file_exists(trim(src%dir)//'flex_pca_merged_pc001.mrc') .and. &
               &file_exists(trim(src%dir)//'flex_pca_probe_merged.txt') )then
            src%pfx = 'flex_pca_merged_pc'; src%meta = 'flex_pca_probe_merged.txt'
        else if( file_exists(trim(src%dir)//'flex_pca_pc001.mrc') .and. &
               &file_exists(trim(src%dir)//'flex_pca_probe.txt') )then
            src%pfx = 'flex_pca_pc'; src%meta = 'flex_pca_probe.txt'
        else
            THROW_HARD('SIMPLE_COV_COMPOSE: no delivered basis + probe meta in '//trim(src%dir))
        endif
        call load_probe_state(src%ncomp, src%eigvals, src%sig2, fname=trim(src%dir)//trim(src%meta))
        nfound = 0
        do q = 1, src%ncomp
            if( .not. file_exists(trim(src%dir)//trim(src%pfx)//int2str_pad(q,3)//MRC_EXT) ) exit
            nfound = q
        end do
        if( nfound < 1 ) THROW_HARD('SIMPLE_COV_COMPOSE: no eigenvolumes under '//trim(src%dir)//trim(src%pfx))
        if( nfound < src%ncomp )then
            write(logfhandle,'(A,I0,A,I0,A,A)') '>>> FLEX_PCA COMPOSE: meta rank ', src%ncomp, &
                &' but only ', nfound, ' eigenvolumes on disk under ', trim(src%dir)//trim(src%pfx)
            src%ncomp = nfound
        endif
        fn = trim(src%dir)//trim(src%pfx)//'001'//MRC_EXT
        call find_ldim_nptcls(fn, ldim, nvol)
        call fn%kill
        src%box = ldim(1)
        if( src%box > params%box_crop ) THROW_HARD('SIMPLE_COV_COMPOSE: a source box exceeds box_crop; compose at the finest box')
    end subroutine resolve_source

    !> stable insertion sort, coarse box first
    subroutine sort_sources_by_box( srcs, nsrc )
        type(compose_src_t), intent(inout) :: srcs(:)
        integer,             intent(in)    :: nsrc
        type(compose_src_t) :: tmp
        integer :: i, j
        do i = 2, nsrc
            tmp = srcs(i)
            j   = i - 1
            do while( j >= 1 )
                if( srcs(j)%box <= tmp%box ) exit
                srcs(j+1) = srcs(j)
                j = j - 1
            end do
            srcs(j+1) = tmp
        end do
    end subroutine sort_sources_by_box

    !> Cut the composed union to its signal subspace BEFORE the deconvolution, by re-embedding:
    !! one all-N embedding with the union, the noise calibration and the noise-whitened
    !! population-variance criterion of the stage cut (SIMPLE_COV_CUT_SNR), then the basis is
    !! rotated onto the kept directions (U W+), re-orthonormalised and written back under the plain
    !! namespace so the run's embedding proper (and its distributed workers) see k columns. The
    !! deconvolution's cost grows with the square of the latent dimension, and an uncut union's
    !! per-particle precision matrices are mostly noise dimensions.
    subroutine compose_cut_reembed( params, cfg, build, model, sel, rounds )
        type(flex_fit_model), intent(inout) :: model
        type(flex_selection), intent(in)    :: sel
        type(flex_latent) :: lat
        class(flex_pca_rounds), intent(inout) :: rounds
        class(parameters),   intent(inout) :: params
        type(flex_run_settings), intent(in)    :: cfg
        type(builder),       intent(inout) :: build
        type(image), allocatable :: imgs(:)
        type(image)  :: vol
        type(string) :: fname
        real(dp), allocatable :: prior(:), a_comp(:)
        real(dp), allocatable :: R(:,:,:), Nz(:,:,:), mu(:,:), Sig(:,:,:), W(:,:), WWt(:,:), WWti(:,:), Wp(:,:)
        real(dp), allocatable :: VV(:,:), U2(:,:), zc(:,:), evc(:), qcol(:), s2(:), r2(:), nzc(:,:)
        real(dp) :: pik(1), a, snr_thr, pj
        real,     pointer :: rmat(:,:,:) => null()
        integer  :: d, k, i, j, q, n, bc, nvox, errflg, kept
        d    = model%ncomp
        n    = sel%nptcls
        bc   = params%box_crop
        nvox = bc*bc*bc
        write(logfhandle,'(A,I0,A)') '>>> FLEX_PCA COMPOSE CUT: embedding the ', d, &
            &'-column union once to find its signal subspace'
        call flush(logfhandle)
        allocate(lat%z(n,d), lat%contrast(n), lat%precision(d,d,n), lat%resid_energy(n), lat%resid_mean_energy(n), lat%zhalf(n,d,2), prior(d), a_comp(d))
        lat%zhalf = 0.d0
        do q = 1, d
            prior(q) = 1.d0/max(model%eigvals(q), DTINY)
        end do
        if( rounds%distributed() )then
            call save_probe_state(d, model%eigvals, model%sig2_eff)
            call rounds%run_stage(params, flex_stage_request(stage=PCA_STAGE_EMBED, label='embedding'))
            call embed_latents_with_contrast(params, build, model, sel, lat, rounds, from_parts=.true., l_zhalf=.true.)
        else
            call embed_latents_with_contrast(params, build, model, sel, lat, rounds, l_zhalf=.true.)
        endif
        deallocate(lat%contrast, lat%resid_energy, lat%resid_mean_energy)
        call calibrate_noise_scale(lat%zhalf, lat%precision, prior, n, d, a, a_comp)
        deallocate(lat%zhalf, a_comp)
        allocate(R(d,d,n), Nz(d,d,n))
        call noise_and_projection(lat%precision, prior, n, d, a, R, Nz)
        deallocate(R, lat%precision)
        ! method-of-moments population covariance: observed scatter minus the mean measurement noise
        allocate(mu(d,1), Sig(d,d,1))
        mu(:,1) = sum(lat%z, dim=1)/real(n,dp)
        Sig = 0.d0
        do i = 1, n
            do q = 1, d
                Sig(:,q,1) = Sig(:,q,1) + (lat%z(i,:) - mu(:,1))*(lat%z(i,q) - mu(q,1))
            end do
        end do
        Sig(:,:,1) = Sig(:,:,1)/real(max(n-1,1),dp)
        do i = 1, n
            Sig(:,:,1) = Sig(:,:,1) - Nz(:,:,i)/real(n,dp)
        end do
        Sig(:,:,1) = 0.5d0*(Sig(:,:,1) + transpose(Sig(:,:,1)))
        pik(1)  = 1.d0
        snr_thr = cfg%cut_snr
        call signal_subspace(mu, Sig, pik, 1, Nz, n, d, snr_thr, W, k)
        write(logfhandle,'(A,F7.3,A,I0,A,I0)') '>>> FLEX_PCA COMPOSE CUT: noise scale a=', a, &
            &'  signal subspace ', k, ' of ', d
        call flush(logfhandle)
        if( k >= d )then
            write(logfhandle,'(A)') '>>> FLEX_PCA COMPOSE CUT: every direction carries population variance; no cut'
            deallocate(lat%z, prior, Nz, mu, Sig, W)
            return
        endif
        ! cut coordinates and their population variance (observed minus noise along each direction)
        allocate(zc(n,k), evc(k), nzc(k,k))
        zc  = matmul(lat%z, transpose(W))
        nzc = 0.d0
        do i = 1, n
            nzc = nzc + matmul(W, matmul(Nz(:,:,i), transpose(W)))/real(n,dp)
        end do
        do q = 1, k
            evc(q) = sum((zc(:,q) - sum(zc(:,q))/real(n,dp))**2)/real(max(n-1,1),dp)
            evc(q) = max(evc(q) - nzc(q,q), 0.1d0*evc(q), DTINY)
        end do
        deallocate(lat%z, Nz, zc, nzc, mu, Sig, prior)
        ! U W+ : the composed columns combined onto the kept directions (z = W+ zc on the signal subspace)
        allocate(WWt(k,k), WWti(k,k), Wp(d,k))
        WWt = matmul(W, transpose(W))
        call matinv(WWt, WWti, k, errflg)
        if( errflg /= 0 ) THROW_HARD('compose_cut_reembed: singular W W^T')
        Wp = matmul(transpose(W), WWti)
        allocate(VV(nvox,d), U2(nvox,k))
        call vol%new([bc,bc,bc], params%smpd_crop)
        do j = 1, d
            fname = 'flex_pca_pc'//int2str_pad(j,3)//MRC_EXT
            call vol%read(fname)
            call fname%kill
            call vol%get_rmat_ptr(rmat)
            VV(:,j) = reshape(real(rmat(1:bc,1:bc,1:bc),dp), [nvox])
        end do
        U2 = matmul(VV, Wp)
        deallocate(VV, W, WWt, WWti, Wp)
        ! re-orthonormalise (the kept directions are orthogonal in the whitened frame, not in real space);
        ! the prior variance follows the rescaling as in the composition
        allocate(qcol(nvox), s2(k), r2(k))
        kept = 0
        do j = 1, k
            s2(j) = sum(U2(:,j)**2)
            if( s2(j) <= 1.d-30 ) cycle
            qcol = U2(:,j)/sqrt(s2(j))
            do q = 1, kept
                pj   = dot_product(qcol, U2(:,q))
                qcol = qcol - pj*U2(:,q)
            end do
            r2(j) = sum(qcol**2)
            if( r2(j) < 1.d-6 ) cycle
            kept = kept + 1
            U2(:,kept) = qcol/sqrt(r2(j))
            evc(kept)  = evc(j)*s2(j)*r2(j)
        end do
        deallocate(qcol)
        if( kept < 2 ) THROW_HARD('compose_cut_reembed: fewer than 2 columns survived the cut')
        ! realise: images -> reconstructors, plain namespace rewritten (k files; the surplus removed)
        allocate(imgs(kept))
        do j = 1, kept
            call imgs(j)%new([bc,bc,bc], params%smpd_crop)
            call imgs(j)%get_rmat_ptr(rmat)
            rmat(1:bc,1:bc,1:bc) = real(reshape(U2(:,j), [bc,bc,bc]))
            fname = 'flex_pca_pc'//int2str_pad(j,3)//MRC_EXT
            call imgs(j)%write(fname, del_if_exists=.true.)
            call fname%kill
        end do
        do j = kept+1, d
            fname = 'flex_pca_pc'//int2str_pad(j,3)//MRC_EXT
            call del_file(fname)
            call fname%kill
        end do
        deallocate(U2)
        do q = 1, size(model%basis_recs)
            call model%basis_recs(q)%dealloc_rho; call model%basis_recs(q)%kill
        end do
        deallocate(model%basis_recs, model%eigvals)
        call basis_recs_from_images(params, build, imgs, kept, model%basis_recs)
        do j = 1, kept
            call imgs(j)%kill
        end do
        deallocate(imgs)
        allocate(model%eigvals(kept))
        model%eigvals = evc(1:kept)
        model%ncomp   = kept
        call save_probe_state(model%ncomp, model%eigvals, model%sig2_eff)
        write(logfhandle,'(A,I0,A)') '>>> FLEX_PCA COMPOSE CUT: ', model%ncomp, ' columns re-embedded; prior variances:'
        write(logfhandle,'(A,12(1X,ES9.2))') '>>>   ', (model%eigvals(q), q=1,min(12,model%ncomp))
        call flush(logfhandle)
        deallocate(evc, s2, r2)
        call vol%kill
    end subroutine compose_cut_reembed

end module simple_flex_pca_embed
