!@descr: flex_pca EM: basis-fit driver, noise-prior calibration, probe-state I/O and band selection
submodule (simple_flex_pca_em) simple_flex_pca_em_fit
use simple_flex_pca_distr,  only: flex_pca_is_master, flex_pca_is_worker, PCA_STAGE_PROBE,&
    &PCA_STAGE_EMBED
use simple_matcher_3Drec,   only: init_rec, cleanup_rec_buffers
use simple_matcher_ptcl_io, only: discrete_read_imgbatch, prepimgbatch
use simple_flex_reconstructor_latent_ops, only: project_fplanes_mean_basis
use simple_flex_projected_latent_model, only: prep_imgs4projected_model
implicit none
#include "simple_local_flags.inc"

contains

    !> Full column-covariance eigenbasis pipeline.
    module subroutine build_covariance_eigenbasis( params, build, mean_rec, pinds, nptcls, &
        &col_sep, neigs_req, basis_recs, eigvals, ncomp_out, sig2_out, basis_imgs, fprefix )
        class(parameters),   intent(inout) :: params
        type(builder),       intent(inout) :: build
        type(reconstructor), intent(inout) :: mean_rec
        integer,             intent(in)    :: pinds(:), nptcls, col_sep, neigs_req
        type(reconstructor), allocatable, intent(out) :: basis_recs(:)
        real(dp),            allocatable, intent(out) :: eigvals(:)
        integer,             intent(out)   :: ncomp_out
        real(dp),            intent(out)   :: sig2_out
        !> optional clean real-space eigenvolumes + output-name prefix, used by the
        !! held-out (cross-halfset) embedding to align two independently fitted bases
        type(image), allocatable, optional, intent(out) :: basis_imgs(:)
        character(len=*),        optional, intent(in)  :: fprefix
        integer,   allocatable :: col_hkl(:,:), col_lookup(:,:,:)
        complex,   allocatable :: Bcol_e(:,:,:,:), Bcol_o(:,:,:,:), colvol(:,:,:,:)
        real,      allocatable :: Hcol_e(:,:,:,:), Hcol_o(:,:,:,:), col_fsc(:)
        type(image), allocatable :: realvols(:), utilde_real(:)
        type(reconstructor) :: work
        type(reconstructor), allocatable :: utilde(:)
        !> capped, halfset-balanced particle subset the column-subspace initialiser runs on. Equals
        !! pinds when the cap is off, which is the default.
        integer, allocatable :: bpinds(:)
        integer :: nbp, nparts_sub
        !> reproducibility-weighted column pruning (SIMPLE_COV_COLFSC_W / _MIN)
        integer :: vcolw, vcolmin, s_keep, s_src, ngood_col
        real    :: colfsc_thresh, colfsc_w
        !> probe-worker handoff, read back from the master's flex_pca_probe.txt
        real(dp),            allocatable :: eig_probe(:)
        real(dp),            allocatable :: zw(:,:), contrastw(:), precw(:,:,:), rew(:), rmew(:)
        real(dp) :: sig2_probe
        integer  :: ncomp_probe
        real(dp), allocatable :: vred(:,:)
        integer :: ncol, nreal, s, lb(3), ub(3), nyq_rec, d_tilde, q, directsvd
        !> optional solve-specific subset (SIMPLE_COV_SOLVE_MAX/_STRIDE); defaults to bpinds
        integer, allocatable :: spinds_solve(:)
        integer :: nsolve
        integer(timer_int_kind) :: t_blk
        logical :: l_cols_ok
        real(dp), allocatable :: svals(:)
        ! one work reconstructor defines the expanded lattice / Nyquist / grid correction
        call init_basis_reconstructor(params, build, work)
        lb      = lbound(work%cmat_exp)
        ub      = ubound(work%cmat_exp)
        nyq_rec = work%get_lfny(1)
        ! NOTE: do not restrict the column sample sweep to the lp band. Everything above the lp shell is
        ! hard-zeroed by realize_hermitian_volume, yet it still moves the eigenvalues by up to 1.2e-4
        ! relative for only ~1.09x -- some path out of band escapes that argument. Find it first.
        ! column selection. A WORKER in the COLS round must not re-run the SNR-greedy selection: its
        ! own estimate_snr_volume call only produced this part's contribution, so the selection would
        ! see a fraction of the variance. It reads the master's choice instead. A worker in the SNR
        ! round runs the selection path only to reach estimate_snr_volume, and exits below.
        ! PROBE WORKER. Like the SOLVE worker, nothing upstream is needed: the master refreshed
        ! flex_pca_pc*.mrc and flex_pca_probe.txt before scheduling this round, so the worker rebuilds
        ! the current basis from disk and contributes one EM half-pass over its own particle range.
        ! niters=1 -- the iteration loop lives on the master, one qsys round per iteration, because
        ! the basis the E-step projects against changes every iteration.
        if( flex_pca_is_worker() .and. params%stage == PCA_STAGE_PROBE )then
            call load_probe_state(ncomp_probe, eig_probe, sig2_probe)
            call load_probe_basis(params, build, ncomp_probe, basis_recs)
            call probe_subspace_iteration(params, build, mean_rec, basis_recs, eig_probe, sig2_probe, &
                &pinds, nptcls, ncomp_probe, 1)
            do s = 1, size(basis_recs)
                call basis_recs(s)%dealloc_rho; call basis_recs(s)%kill
            end do
            deallocate(basis_recs)
            if( allocated(eig_probe) ) deallocate(eig_probe)
            if( .not. allocated(eigvals) ) allocate(eigvals(0))
            allocate(basis_recs(0))
            ncomp_out = 0
            sig2_out  = 0._dp
            call work%dealloc_rho; call work%kill
            return
        endif
        ! embed worker: same handoff as the probe worker (basis on disk as flex_pca_pc*.mrc,
        ! dimension/prior variances/noise level in flex_pca_probe.txt) but only one round, since the
        ! basis is final by now and does not change under the workers
        if( flex_pca_is_worker() .and. params%stage == PCA_STAGE_EMBED )then
            call load_probe_state(ncomp_probe, eig_probe, sig2_probe)
            call load_probe_basis(params, build, ncomp_probe, basis_recs)
            allocate(zw(nptcls,ncomp_probe), contrastw(nptcls), precw(ncomp_probe,ncomp_probe,nptcls))
            allocate(rew(nptcls), rmew(nptcls))
            call embed_latents_with_contrast(params, build, mean_rec, basis_recs, ncomp_probe, &
                &eig_probe, sig2_probe, pinds, nptcls, zw, contrastw, precw, rew, rmew, &
                &stats_only=.true.)
            deallocate(zw, contrastw, precw, rew, rmew)
            do s = 1, size(basis_recs)
                call basis_recs(s)%dealloc_rho; call basis_recs(s)%kill
            end do
            deallocate(basis_recs)
            if( allocated(eig_probe) ) deallocate(eig_probe)
            if( .not. allocated(eigvals) ) allocate(eigvals(0))
            allocate(basis_recs(0))
            ncomp_out = 0
            sig2_out  = 0._dp
            call work%dealloc_rho; call work%kill
            return
        endif
        ! ---- EM PATH (the only estimator) ----
        ! Hand probe_subspace_iteration a data-free basis. The moment/covariance estimator that used
        ! to live here -- SNR volume, column selection and accumulation, reduced solve -- has been
        ! removed; probe_subspace_iteration is itself a PPCA EM and supersedes it. Only reachable on
        ! the master: a worker in a PROBE or EMBED round has already returned above with the basis it
        ! read off disk, so nothing here runs twice.
        ! basis_imgs is the PRE-refinement basis for the held-out/bagged arms. Under EM that basis is
        ! a data-free placeholder and returning it would silently compare initialisers instead of
        ! fits, so refuse rather than mislead. Two-fit reproducibility is measured by running two
        ! jobs and comparing their written eigenvolumes.
        if( present(basis_imgs) ) THROW_HARD('the held-out/bagged basis arms need a pre-refinement &
            &basis, which the EM initialiser cannot provide; run two jobs and compare flex_pca_pc*.mrc')
        call init_basis_datafree(params, build, mean_rec, pinds, nptcls, col_sep, neigs_req, &
            &basis_recs, eigvals, ncomp_out, sig2_out)
        call work%dealloc_rho; call work%kill
    end subroutine build_covariance_eigenbasis

    !> DATA-FREE EM INITIALISER (SIMPLE_COV_EM=1).
    !!
    !! Replaces the whole moment estimator -- SNR volume, column selection, column accumulation,
    !! merge, reduced solve -- as the thing that hands probe_subspace_iteration its starting basis.
    !! probe_subspace_iteration is ALREADY a PPCA EM (posterior precision A = (a^2/sig2)G + Gamma^-1,
    !! E[zz'] = zz' + A^-1, coupled per-voxel M-step, Gamma from the posterior second moment). What it
    !! has never been given is a start that is not the covariance eigenbasis, which is why running it
    !! today is evidence for neither the EM formulation nor against it: if the columns are noise, the
    !! probe starts inside a noise subspace and its Gamma^(0) are noise variances.
    !!
    !! The subspace here is the lowest-frequency Fourier modes the band admits, realised as cos/sin
    !! Friedel pairs through exactly the same band-limit + mask + deapodisation the column path uses,
    !! then orthonormalised. NO data moment is formed anywhere: select_frequencies_lowfreq is
    !! pure geometry (greedy smallest-|xi|, col_sep-separated) and the "columns" fed to it are unit
    !! impulses, not estimates. Deterministic, so two arms are comparable run to run.
    module subroutine init_basis_datafree( params, build, mean_rec, pinds, nptcls, col_sep, neigs_req, &
        &basis_recs, eigvals, ncomp_out, sig2_out )
        class(parameters),   intent(inout) :: params
        type(builder),       intent(inout) :: build
        type(reconstructor), intent(inout) :: mean_rec
        integer,             intent(in)    :: pinds(:), nptcls, col_sep, neigs_req
        type(reconstructor), allocatable, intent(out) :: basis_recs(:)
        real(dp),            allocatable, intent(out) :: eigvals(:)
        integer,             intent(out)   :: ncomp_out
        real(dp),            intent(out)   :: sig2_out
        type(reconstructor) :: work
        type(reconstructor), allocatable :: utilde(:)
        type(image),         allocatable :: realvols(:), utilde_real(:)
        integer,             allocatable :: col_hkl(:,:)
        complex,             allocatable :: colvol(:,:,:,:)
        real(dp),            allocatable :: svals(:)
        integer  :: lb(3), ub(3), ncols_req, ncol, nreal, d_tilde, s, q, h, k, l
        integer, allocatable :: cpinds(:)
        integer  :: ncal
        real(dp) :: gam0
        type(string) :: fname
        call init_basis_reconstructor(params, build, work)
        lb = lbound(work%cmat_exp)
        ub = ubound(work%cmat_exp)
        ! each impulse yields a cos AND a sin volume, so half as many impulses as components are
        ! needed; a small margin covers the pairs orthonormalisation drops at the energy floor
        ncols_req = max(1, (neigs_req + 1)/2 + 2)
        call select_frequencies_lowfreq(params, ncols_req, col_sep, col_hkl, ncol)
        allocate(colvol(lb(1):ub(1),lb(2):ub(2),lb(3):ub(3),ncol), source=cmplx(0.,0.))
        do s = 1, ncol
            h = col_hkl(1,s); k = col_hkl(2,s); l = col_hkl(3,s)
            if( h < lb(1) .or. h > ub(1) .or. k < lb(2) .or. k > ub(2) &
                &.or. l < lb(3) .or. l > ub(3) ) cycle
            colvol(h,k,l,s) = cmplx(1.,0.)
        end do
        call basis_to_real_representatives(params, work, colvol, ncol, lb, ub, realvols, nreal)
        deallocate(colvol, col_hkl)
        if( nreal < 1 ) THROW_HARD('flex_pca EM initialiser produced no basis representatives')
        call orthonormalize_representatives(params, build, realvols, nreal, utilde, utilde_real, &
            &d_tilde, svals, nptcls_basis=nptcls)
        do s = 1, nreal
            call realvols(s)%kill
        end do
        deallocate(realvols)
        ncomp_out = max(1, min(neigs_req, d_tilde))
        call basis_recs_from_images(params, build, utilde_real(1:ncomp_out), ncomp_out, basis_recs)
        ! The eigenvolume MRCs are the master->worker handoff for every distributed probe round; on
        ! the moment path form_eigenbasis_from_reduced writes them, and nothing else does, so the EM
        ! path has to write its own initial basis or the first PROBE round finds no flex_pca_pc001.mrc.
        if( flex_pca_is_master() .or. .not. flex_pca_is_worker() )then
            do q = 1, ncomp_out
                fname = string('flex_pca_pc')//int2str_pad(q,3)//MRC_EXT
                call utilde_real(q)%write(fname, del_if_exists=.true.); call fname%kill
            end do
        endif
        write(logfhandle,'(A,I0,A,I0,A,I0)') '>>> FLEX_PCA EM INIT (data-free): impulses=',ncol, &
            &'  representatives=',nreal,'  rank=',ncomp_out
        call flush(logfhandle)
        ! the ONLY data pass: the noise convention constant and the initial prior variance.
        ! Two scalars do not need a million particles, and at box_crop=96 on the full dataset an
        ! uncapped pass costs more than several EM iterations, so it shares the initialiser's
        ! budget SIMPLE_COV_CALIB_MAX (default COV_CALIB_MAX_PTCLS; SIMPLE_COV_BASIS_MAX used to
        ! be shared here but its default of 0 left this pass uncapped). Always master here -- a worker returned
        ! long before this point -- so nparts is 1 and the budget is not divided twice.
        call cov_stage_subsample(build, pinds, nptcls, 1, COV_CALIB_MAX_PTCLS, &
            &'SIMPLE_COV_CALIB_STRIDE', 'SIMPLE_COV_CALIB_MAX', 'EM CALIBRATION', cpinds, ncal)
        call em_calibrate_noise_prior(params, build, mean_rec, basis_recs, ncomp_out, cpinds, ncal, &
            &sig2_out, gam0)
        deallocate(cpinds)
        allocate(eigvals(ncomp_out), source=gam0)
        do s = 1, d_tilde
            call utilde(s)%dealloc_rho; call utilde(s)%kill
            call utilde_real(s)%kill
        end do
        deallocate(utilde, utilde_real)
        if( allocated(svals) ) deallocate(svals)
        call work%dealloc_rho; call work%kill
    end subroutine init_basis_datafree

    !> One pass over the particles to fix the two scalars the EM cannot invent: the whitened-noise
    !! convention constant sig2 (measured on the high-frequency shells, where conformational signal is
    !! negligible) and the initial prior variance Gamma^(0).
    !!
    !! Gamma^(0) is deliberately set from the mean-deflated residual with only the noise floor removed,
    !! i.e. an OVER-estimate of the conformational signal, which makes the iteration-1 prior WEAK. The
    !! asymmetry is on purpose: a prior that starts too tight shrinks z, which shrinks the Gamma update,
    !! which tightens the prior -- the self-reinforcing collapse that is the single most plausible cause
    !! of the previous EM's failure. Starting loose is corrected by the very first Gamma M-step; starting
    !! tight is not corrected by anything.
    module subroutine em_calibrate_noise_prior( params, build, mean_rec, basis_recs, ncomp, pinds, nptcls, &
        &sig2_out, gam0 )
        class(parameters),   intent(inout) :: params
        type(builder),       intent(inout) :: build
        type(reconstructor), intent(inout) :: mean_rec
        type(reconstructor), intent(inout) :: basis_recs(:)
        integer,             intent(in)    :: ncomp, pinds(:), nptcls
        real(dp),            intent(out)   :: sig2_out, gam0
        type(fplane_type), allocatable :: fpls(:), basis_fpls(:,:), mean_fpl(:)
        type(ori),         allocatable :: orientations(:)
        real(dp), allocatable :: res_thr(:), trg_thr(:), hfp_thr(:), hfc_thr(:), aa_thr(:)
        integer,  allocatable :: nval_thr(:)
        integer  :: nthr, ithr, i, q, ibatch, batchlims(2), batchsz, nyq_rec, nval
        real(dp) :: a, aa, e_mm, e_yy, myv, res, trg, pw, cnt, wcnt
        real(dp) :: res_sum, trg_sum, aa_sum, hfpw, hfcnt
        nthr    = nthr_glob
        nyq_rec = mean_rec%get_lfny(1)
        call mean_rec%expand_exp
        do q = 1, ncomp
            call basis_recs(q)%expand_exp
        end do
        call init_rec(params, build, MAXIMGBATCHSZ, fpls)
        call prepimgbatch(params, build, MAXIMGBATCHSZ)
        allocate(mean_fpl(nthr), basis_fpls(ncomp,nthr), orientations(MAXIMGBATCHSZ))
        allocate(res_thr(nthr), trg_thr(nthr), hfp_thr(nthr), hfc_thr(nthr), aa_thr(nthr), source=0.d0)
        allocate(nval_thr(nthr), source=0)
        wcnt = 0.d0
        do ibatch = 1, nptcls, MAXIMGBATCHSZ
            batchlims = [ibatch, min(nptcls, ibatch + MAXIMGBATCHSZ - 1)]
            batchsz   = batchlims(2) - batchlims(1) + 1
            call discrete_read_imgbatch(params, build, nptcls, pinds, batchlims)
            call prep_imgs4projected_model(params, build, batchsz, build%imgbatch(:batchsz), &
                &pinds(batchlims(1):batchlims(2)), fpls(:batchsz), mskrad=cov_image_mask_radius(params))
            do i = 1, batchsz
                call build%spproj_field%get_ori(pinds(batchlims(1)+i-1), orientations(i))
            end do
            if( wcnt <= 0.d0 ) call cov_herm_selfpower(fpls(1), pw, wcnt)
            !$omp parallel do default(shared) schedule(dynamic) proc_bind(close) &
            !$omp& private(i,ithr,q,a,aa,e_mm,e_yy,myv,res,trg,pw,cnt)
            do i = 1, batchsz
                if( orientations(i)%isstatezero() ) cycle
                ithr = omp_get_thread_num() + 1
                call project_fplanes_mean_basis(mean_rec, basis_recs, orientations(i), fpls(i), &
                    &mean_fpl(ithr), basis_fpls(:,ithr), apply_ctf_amp=.true.)
                e_mm = real(cov_herm_inner(mean_fpl(ithr), mean_fpl(ithr)), dp)
                myv  = real(cov_herm_inner(mean_fpl(ithr), fpls(i)), dp)
                e_yy = real(cov_herm_inner(fpls(i), fpls(i)), dp)
                a    = max(0.1d0, min(5.0d0, myv / max(e_mm, DTINY)))
                aa   = a*a
                res  = max(e_yy - 2.d0*a*myv + aa*e_mm, 0.d0)
                trg  = 0.d0
                do q = 1, ncomp
                    trg = trg + real(cov_herm_inner(basis_fpls(q,ithr), basis_fpls(q,ithr)), dp)
                end do
                call plane_hf_power(fpls(i), nyq_rec, 0.7, pw, cnt)
                res_thr(ithr)  = res_thr(ithr)  + res
                trg_thr(ithr)  = trg_thr(ithr)  + trg
                aa_thr(ithr)   = aa_thr(ithr)   + aa
                hfp_thr(ithr)  = hfp_thr(ithr)  + pw
                hfc_thr(ithr)  = hfc_thr(ithr)  + cnt
                nval_thr(ithr) = nval_thr(ithr) + 1
            end do
            !$omp end parallel do
        end do
        res_sum = sum(res_thr); trg_sum = sum(trg_thr); aa_sum = sum(aa_thr)
        hfpw    = sum(hfp_thr); hfcnt   = sum(hfc_thr); nval = sum(nval_thr)
        if( nval < 1 ) THROW_HARD('flex_pca EM calibration saw no valid particles')
        sig2_out = max(hfpw / max(hfcnt, 1.d0), DTINY)
        ! signal power = mean-deflated residual with the noise floor removed; floor it at a tenth of
        ! the residual so a mis-measured sig2 cannot drive the prior to zero
        res_sum = max(res_sum/real(nval,dp) - sig2_out*wcnt, 0.1d0*res_sum/real(nval,dp))
        trg_sum = max(trg_sum/real(nval,dp), DTINY)
        aa_sum  = max(aa_sum /real(nval,dp), DTINY)
        gam0    = max(res_sum / (aa_sum * trg_sum), DTINY)
        write(logfhandle,'(A,ES12.4,A,ES12.4,A,I0)') '>>> FLEX_PCA EM CALIBRATION: sig2=',sig2_out, &
            &'  gamma0=',gam0,'  particles=',nval
        write(logfhandle,'(A,ES12.4,A,ES12.4,A,F6.3)') '>>>   mean deflated signal=',res_sum, &
            &'  mean tr(G)=',trg_sum,'  mean a^2=',real(aa_sum)
        call flush(logfhandle)
        do ithr = 1, nthr
            call cleanup_plane(mean_fpl(ithr))
            do q = 1, ncomp
                call cleanup_plane(basis_fpls(q,ithr))
            end do
        end do
        do i = 1, size(orientations)
            call orientations(i)%kill
        end do
        call cleanup_rec_buffers(build, fpls)
        deallocate(mean_fpl, basis_fpls, orientations)
        deallocate(res_thr, trg_thr, hfp_thr, hfc_thr, aa_thr, nval_thr)
    end subroutine em_calibrate_noise_prior





    !>  Master -> probe-worker handoff: basis dimension, prior variances, whitened-noise level.
    module subroutine save_probe_state( ncomp, eigvals, sig2_eff )
        integer,  intent(in) :: ncomp
        real(dp), intent(in) :: eigvals(:), sig2_eff
        integer :: funit, io_stat, q
        call fopen(funit, file=string(COV_PROBE_META), action='WRITE', status='REPLACE', iostat=io_stat)
        call fileiochk('save_probe_state', io_stat)
        write(funit,*) ncomp
        write(funit,*) sig2_eff
        do q = 1, ncomp
            write(funit,*) eigvals(q)
        end do
        call fclose(funit)
    end subroutine save_probe_state

    module subroutine load_probe_state( ncomp, eigvals, sig2_eff )
        integer,               intent(out) :: ncomp
        real(dp), allocatable, intent(out) :: eigvals(:)
        real(dp),              intent(out) :: sig2_eff
        integer :: funit, io_stat, q
        if( .not. file_exists(string(COV_PROBE_META)) ) &
            &THROW_HARD('flex_pca probe worker found no '//COV_PROBE_META//' from the master')
        call fopen(funit, file=string(COV_PROBE_META), action='READ', status='OLD', iostat=io_stat)
        call fileiochk('load_probe_state', io_stat)
        read(funit,*) ncomp
        read(funit,*) sig2_eff
        if( ncomp < 1 ) THROW_HARD('invalid cached probe basis dimension')
        allocate(eigvals(ncomp))
        do q = 1, ncomp
            read(funit,*) eigvals(q)
        end do
        call fclose(funit)
    end subroutine load_probe_state

    !>  Rebuild the projection-ready basis a probe worker needs from the master's flex_pca_pc*.mrc.
    !!  Same idiom as load_utilde_stack and probe_external_basis: set_rmat then fft then expand_exp,
    !!  never add(), which would leave the reconstructor flagged Fourier and propagate an
    !!  untransformed grid.
    module subroutine load_probe_basis( params, build, ncomp, basis_recs )
        class(parameters),   intent(inout) :: params
        type(builder),       intent(inout) :: build
        integer,             intent(in)    :: ncomp
        type(reconstructor), allocatable, intent(out) :: basis_recs(:)
        type(image)  :: vol
        type(string) :: fname
        integer      :: q
        allocate(basis_recs(ncomp))
        call vol%new([params%box_crop,params%box_crop,params%box_crop], params%smpd_crop)
        do q = 1, ncomp
            fname = string('flex_pca_pc')//int2str_pad(q,3)//MRC_EXT
            if( .not. file_exists(fname) ) &
                &THROW_HARD('flex_pca probe worker found no '//fname%to_char()//' from the master')
            call vol%read(fname)
            call init_basis_reconstructor(params, build, basis_recs(q))
            call basis_recs(q)%set_rmat(vol%get_rmat(), .false.)
            call basis_recs(q)%fft
            call basis_recs(q)%expand_exp
            call fname%kill
        end do
        call vol%kill
    end subroutine load_probe_basis



    !> Deterministic low-frequency column selection: repeatedly take the lowest-|xi| candidate in the
    !! canonical Hermitian half that is at least col_sep away from every already-chosen column.
    module subroutine select_frequencies_lowfreq( params, ncols_req, col_sep, col_hkl, ncol )
        class(parameters),    intent(in)  :: params
        integer,              intent(in)  :: ncols_req, col_sep
        integer, allocatable, intent(out) :: col_hkl(:,:)
        integer,              intent(out) :: ncol
        integer, allocatable :: cand(:,:)
        integer :: kfromto(2), kmax, kmax_sq, h, k, l, sep, r_sq, ncand, i, target
        kfromto = covariance_kfromto(params)
        kmax    = max(2, kfromto(2))
        kmax_sq = kmax*(kmax+1)
        sep     = max(1, col_sep)
        target  = max(1, ncols_req)
        allocate(cand(3, (2*kmax+1)**3))
        ncand = 0
        do h = 0, kmax
            do k = -kmax, kmax
                do l = -kmax, kmax
                    r_sq = h*h + k*k + l*l
                    if( r_sq == 0 .or. r_sq > kmax_sq ) cycle
                    if( h == 0 )then
                        if( k < 0 ) cycle
                        if( k == 0 .and. l < 0 ) cycle
                    endif
                    ncand = ncand + 1
                    cand(:,ncand) = [h,k,l]
                end do
            end do
        end do
        allocate(col_hkl(3, target))
        ncol = 0
        do
            call pick_next_lowfreq(cand, ncand, col_hkl, ncol, sep, i)
            if( i == 0 ) exit
            ncol = ncol + 1
            col_hkl(:,ncol) = cand(:,i)
            cand(:,i) = huge(1)
            if( ncol >= target ) exit
        end do
        if( ncol < 1 ) THROW_HARD('flex_pca laid out no impulse columns for the EM initialiser; increase lp or neigs')
        deallocate(cand)
    end subroutine select_frequencies_lowfreq

    module subroutine pick_next_lowfreq( cand, ncand, chosen, nchosen, sep, best )
        integer, intent(in)  :: cand(:,:), ncand, chosen(:,:), nchosen, sep
        integer, intent(out) :: best
        integer :: i, j, r_sq, best_r, d(3)
        logical :: ok
        best   = 0
        best_r = huge(1)
        do i = 1, ncand
            if( cand(1,i) == huge(1) ) cycle
            r_sq = sum(cand(:,i)**2)
            if( r_sq >= best_r ) cycle
            ok = .true.
            do j = 1, nchosen
                d = cand(:,i) - chosen(:,j)
                if( sum(d**2) < sep*sep )then
                    ok = .false.; exit
                endif
            end do
            if( ok )then
                best   = i
                best_r = r_sq
            endif
        end do
    end subroutine pick_next_lowfreq

    module function covariance_kfromto( params ) result( kfromto )
        class(parameters), intent(in) :: params
        integer :: kfromto(2), kto_full
        real    :: dstep_crop
        kto_full   = max(1, fdim(params%box_crop) - 1)
        kfromto(1) = 1
        kfromto(2) = kto_full
        if( params%lp > 2.0 * params%smpd_crop + TINY )then
            dstep_crop = real(max(1, params%box_crop - 1)) * params%smpd_crop
            kfromto(2) = max(1, min(kto_full, int(dstep_crop / params%lp)))
        endif
    end function covariance_kfromto


    !> EXTERNAL-BASIS PROBE: embed the particles in a basis read from disk (the run's own eigenvolumes,
    !! optionally with extra probe volumes appended) and write the coefficients, so the embedding stage can
    !! be exercised against a known basis without re-fitting the covariance.
    module subroutine probe_external_basis( params, build, mean_rec, pinds, nptcls, eigdir, neigs, eigvals, &
        &sig2_eff, probe_prefix, nprobe )
        class(parameters),   intent(inout) :: params
        type(builder),       intent(inout) :: build
        type(reconstructor), intent(inout) :: mean_rec
        integer,             intent(in)    :: pinds(:), nptcls, neigs, nprobe
        character(len=*),    intent(in)    :: eigdir, probe_prefix
        real(dp),            intent(in)    :: eigvals(:), sig2_eff
        type(reconstructor), allocatable :: basis_recs(:)
        type(image)  :: vol
        type(string) :: fname
        real(dp), allocatable :: ev(:), z(:,:), contrast(:), precision(:,:,:)
        real(dp), allocatable :: resid_energy(:), resid_mean_energy(:), sorted(:)
        real     :: dummy
        integer  :: ncomb, k, u, i, q
        real(dp) :: evmed
        ! nprobe = 0 is legitimate: an appended probe volume can dominate the projected Gram spectrum, in
        ! which case the reported conditioning describes the probe rather than the basis.
        if( nprobe < 0 ) return
        ncomb = neigs + nprobe
        write(logfhandle,'(A,I0,A,I0,A)') '>>> FLEX_PCA EXTERNAL-BASIS PROBE: ', neigs, &
            &' eigenvolumes + ', nprobe, ' probe volumes in one joint solve'
        if( ncomb < 2 ) THROW_HARD('probe_external_basis: need at least 2 basis volumes')
        call flush(logfhandle)
        allocate(basis_recs(ncomb), ev(ncomb))
        call vol%new([params%box_crop,params%box_crop,params%box_crop], params%smpd_crop)
        do k = 1, ncomb
            if( k <= neigs )then
                fname = trim(eigdir)//'flex_pca_pc'//int2str_pad(k,3)//MRC_EXT
            else
                fname = trim(probe_prefix)//int2str_pad(k-neigs,3)//MRC_EXT
            endif
            if( .not. file_exists(fname%to_char()) )then
                write(logfhandle,'(A,A)') '>>> FLEX_PCA probe basis: missing ', fname%to_char()
                call fname%kill; call vol%kill
                do i = 1, k-1
                    call basis_recs(i)%dealloc_rho; call basis_recs(i)%kill
                end do
                deallocate(basis_recs, ev)
                return
            endif
            call vol%read(fname)
            call fname%kill
            if( params%msk_crop > TINY ) call vol%mask3D_soft(params%msk_crop, backgr=0.)
            ! the set_rmat + fft + expand_exp idiom -- NEVER add(), which silently yields a
            ! zero projected basis when the reconstructor is left flagged as Fourier
            call init_basis_reconstructor(params, build, basis_recs(k))
            call basis_recs(k)%set_rmat(vol%get_rmat(), .false.)
            call basis_recs(k)%fft
            call basis_recs(k)%expand_exp
        end do
        call vol%kill
        allocate(sorted(neigs), source=eigvals(1:neigs))
        evmed = sorted(max(1,neigs/2))
        do k = 1, ncomb
            if( k <= neigs )then
                ev(k) = max(eigvals(k), DTINY)
            else
                ev(k) = max(evmed, DTINY)
            endif
        end do
        write(logfhandle,'(A,ES12.4)') '>>> FLEX_PCA probe prior variance (median eigenvalue): ', evmed
        allocate(z(nptcls,ncomb), contrast(nptcls), precision(ncomb,ncomb,nptcls), &
            &resid_energy(nptcls), resid_mean_energy(nptcls))
        call embed_latents_with_contrast(params, build, mean_rec, basis_recs, ncomb, ev, sig2_eff, &
            &pinds, nptcls, z, contrast, precision, resid_energy, resid_mean_energy)
        call del_file('flex_pca_probe_coordinates.txt')
        open(newunit=u, file='flex_pca_probe_coordinates.txt', status='replace', action='write')
        write(u,'(A)',advance='no') '# particle'
        do q = 1, neigs
            write(u,'(A,I0)',advance='no') ' pc', q
        end do
        do q = 1, nprobe
            write(u,'(A,I0)',advance='no') ' probe', q
        end do
        write(u,*)
        do i = 1, nptcls
            write(u,'(I10)',advance='no') pinds(i)
            do q = 1, ncomb
                write(u,'(1X,ES16.8)',advance='no') z(i,q)
            end do
            write(u,*)
        end do
        close(u)
        write(logfhandle,'(A)') '>>> FLEX_PCA probe coefficients -> flex_pca_probe_coordinates.txt'
        call flush(logfhandle)
        do k = 1, ncomb
            call basis_recs(k)%dealloc_rho; call basis_recs(k)%kill
        end do
        deallocate(basis_recs, ev, z, contrast, precision, resid_energy, resid_mean_energy, sorted)
        dummy = 0.
    end subroutine probe_external_basis

    module subroutine init_basis_reconstructor( params, build, rec )
        class(parameters),   intent(inout) :: params
        type(builder),       intent(inout) :: build
        type(reconstructor), intent(inout) :: rec
        call rec%new([params%box_crop,params%box_crop,params%box_crop],params%smpd_crop)
        call rec%alloc_rho(params,build%spproj,expand=.true.)
        call rec%reset
        call rec%reset_exp
    end subroutine init_basis_reconstructor

    !> Particle-image mask radius in pixels at params%box, or 0 to disable.
    module function cov_image_mask_radius( params ) result( r )
        class(parameters), intent(in) :: params
        real :: r
        integer :: vmi
        ! Runtime override (SIMPLE_COV_MASK_IMAGES=1). The compile-time default is OFF, which is
        ! safe ONLY when the solvent is pure noise -- true for synthetic data, FALSE for real data,
        ! where the region outside the envelope carries ice-thickness gradients, neighbouring
        ! particles and carbon edges. Those are low-frequency AND reproducible between halfsets, so
        ! they enter the covariance as apparent signal and can dominate the leading eigenvectors.
        r = 0.
        vmi = 0
        call cov_env_int('SIMPLE_COV_MASK_IMAGES', vmi)
        if( .not. (COV_MASK_IMAGES .or. vmi > 0) ) return
        if( params%msk_crop <= 0. .or. params%box_crop <= 0 ) return
        r = COV_MASK_MARGIN * params%msk_crop * real(params%box) / real(params%box_crop)
        r = min(r, 0.5*real(params%box) - COSMSKHALFWIDTH - 1.)
    end function cov_image_mask_radius

    !> SIMPLE_COV_DUMP_ACC=1 dumps the raw stage accumulators (SNR var/dens, column B/H) and the
    !! column selection to the run dir, for byte-level A/B validation of threading restructures.
    !! Reuses the distributed part-file writers, so the dumps carry the versioned headers.
    logical module function cov_dump_acc_on() result( on )
        character(len=32) :: envval
        integer :: stat, ln, ival
        on = .false.
        call get_environment_variable('SIMPLE_COV_DUMP_ACC', envval, ln, stat)
        if( stat /= 0 .or. ln < 1 ) return
        read(envval(:ln), *, iostat=stat) ival
        if( stat == 0 ) on = ival /= 0
    end function cov_dump_acc_on

end submodule simple_flex_pca_em_fit
