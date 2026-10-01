!@descr: flex_pca: the mean and its scale, the data-free basis initialiser and its calibration, probe-state and basis I/O, basis pooling, deflation and cross-half angles
module simple_flex_pca_basis
use simple_core_module_api
use simple_flex_pca_records, only: flex_fit_model, flex_selection
use simple_builder, only: builder
use simple_image, only: image
use simple_parameters, only: parameters
use simple_reconstructor, only: reconstructor
use simple_gridding, only: prep3D_inv_kbenvelope4mul
use simple_linalg, only: jacobi, eigsrt
use simple_math, only: ceil_div, floor_div
use simple_flex_reconstructor_latent_ops, only: project_fplane_mean, project_fplanes_mean_basis,&
    &prep_imgs4projected_model
use simple_flex_pca_pcg, only: flex_window_apply_rec
use simple_ori, only: ori
use simple_flex_pca_rounds, only: flex_pca_rounds
use simple_flex_pca_artifacts, only: FLEX_PCA_PART_MAGIC
use simple_flex_pca_run_types, only: flex_run_settings
use simple_flex_pca_mstep, only: init_basis_reconstructor
use simple_flex_pca_util, only: COV_ATHR_BUDGET, cov_signal_rank, cov_stage_subsample, cov_accum_bytes,&
    &cov_dim_budget
use simple_flex_pca_fit_types, only: cleanup_plane
use simple_matcher_3Drec, only: init_rec, cleanup_rec_buffers
use simple_matcher_ptcl_io, only: discrete_read_imgbatch, prepimgbatch
use simple_flex_pca_planes, only: planes_batch_load
use simple_flex_pca_plane_cache, only: plane_cache_in_use
implicit none
private
#include "simple_local_flags.inc"

public :: estimate_covariance_mean, apply_cached_mean_scale, estimate_mean_scale, init_mean_reconstructor, cov_herm_inner
public :: init_basis_datafree, save_probe_state, load_probe_state, load_probe_basis
public :: covariance_kfromto, cov_image_mask_radius, orthonormalize_representatives, align_basis_to_reference
public :: basis_recs_from_images, deflate_against_basis, cross_half_subspace_angles, COV_EIG_REL_FLOOR, COV_MAX_DTILDE
public :: COV_DEFAULT_DTILDE, COV_SAMPLES_PER_PARAM, COV_UNIT_CONTRAST, COV_MEAN_FROM_DATA, COV_MASK_IMAGES
public :: COV_MASK_MARGIN, COV_PROBE_META

real(dp), parameter :: COV_EIG_REL_FLOOR = 1.0d-6
! Rank cap on the orthonormalised representative subspace (orthonormalize_representatives).
integer,  parameter :: COV_MAX_DTILDE    = 320
! Default column-subspace dimension, applied as a min against the memory budget so the rank follows
! the data rather than free RAM.
integer,  parameter :: COV_DEFAULT_DTILDE = 128
real(dp), parameter :: COV_SAMPLES_PER_PARAM = 10.0d0
logical,  parameter :: COV_UNIT_CONTRAST  = .true.
logical, parameter :: COV_MEAN_FROM_DATA = .false.
logical, parameter :: COV_MASK_IMAGES = .false.
real, parameter :: COV_MASK_MARGIN = 1.4
character(len=*), parameter :: COV_PROBE_META   = 'flex_pca_probe.txt'

character(len=*), parameter :: MEAN_SCALE_FNAME = 'flex_pca_mean_scale.bin'


contains

    !> Single entry point for the covariance mean.
    subroutine estimate_covariance_mean( params, build, mean_rec, pinds, nptcls , rounds)
        class(flex_pca_rounds), intent(inout) :: rounds
        class(parameters),   intent(inout) :: params
        type(builder),       intent(inout) :: build
        type(reconstructor), intent(inout) :: mean_rec
        integer,             intent(in)    :: pinds(:), nptcls
        write(logfhandle,'(A)') '>>> FLEX_PCA SPLIT-HALF: hashed lattice split (alias-free)'
        call flush(logfhandle)
        if( COV_MEAN_FROM_DATA )then
            call estimate_mean_from_data(params, build, mean_rec, pinds, nptcls)
        else
            call init_mean_reconstructor(params, build, mean_rec)
            if( rounds%is_worker() )then
                call apply_cached_mean_scale(params, mean_rec)
            else
                call estimate_mean_scale(params, build, mean_rec, pinds, nptcls, rounds=rounds)
            endif
        endif
    end subroutine estimate_covariance_mean

    !> Kernel-regression consensus mean (eq. S.1) estimated from the particles themselves, as an
    !! alternative to reading the supplied consensus volume. Selected by COV_MEAN_FROM_DATA.
    subroutine estimate_mean_from_data( params, build, mean_rec, pinds, nptcls )
        class(parameters),   intent(inout) :: params
        type(builder),       intent(inout) :: build
        type(reconstructor), intent(inout) :: mean_rec
        integer,             intent(in)    :: pinds(:), nptcls
        type(fplane_type), allocatable :: fpls(:)
        type(fplane_type) :: num_fpl
        type(ori)    :: orientation
        type(image)  :: gridcorr_img
        type(string) :: fname
        integer :: batchlims(2), batchsz, ibatch, i, iptcl, used
        logical :: l_pcache
        integer(timer_int_kind) :: t_phase
        call init_basis_reconstructor(params, build, mean_rec)
        ! one read path for every pass of the run: the downscaled cache when it is in use
        l_pcache = plane_cache_in_use(params, build)
        call init_rec(params, build, MAXIMGBATCHSZ, fpls, cropped=l_pcache)
        if( l_pcache )then
            call prepimgbatch(params, build, MAXIMGBATCHSZ, box=params%box_crop, smpd=params%smpd_crop)
        else
            call prepimgbatch(params, build, MAXIMGBATCHSZ)
        endif
        used    = 0
        t_phase = tic()
        write(logfhandle,'(A)') '>>> FLEX_PCA MEAN ESTIMATION (eq. S.1 kernel regression)'
        call flush(logfhandle)
        do ibatch = 1, nptcls, MAXIMGBATCHSZ
            batchlims = [ibatch, min(nptcls, ibatch + MAXIMGBATCHSZ - 1)]
            batchsz   = batchlims(2) - batchlims(1) + 1
            call planes_batch_load(params, build, nptcls, pinds, batchlims, fpls, &
                &cov_image_mask_radius(params), l_pcache)
            do i = 1, batchsz
                iptcl = pinds(batchlims(1)+i-1)
                call build%spproj_field%get_ori(iptcl, orientation)
                if( orientation%isstatezero() ) cycle
                call form_reconstruction_plane(fpls(i), num_fpl)
                call mean_rec%insert_plane_oversamp(build%pgrpsyms, orientation, num_fpl)
                used = used + 1
            end do
        end do
        call orientation%kill
        call cleanup_plane(num_fpl)
        call cleanup_rec_buffers(build, fpls)
        if( used < 1 ) THROW_HARD('flex_pca mean estimation found no valid particles')
        ! canonical gridding finalization, identical to reconstruct3D_reference
        call mean_rec%compress_exp
        call mean_rec%sampl_dens_correct
        call mean_rec%ifft
        gridcorr_img = prep3D_inv_kbenvelope4mul([params%box_crop,params%box_crop,params%box_crop], params%smpd_crop)
        call mean_rec%mul(gridcorr_img)
        call gridcorr_img%kill
        fname = 'flex_pca_mean'//MRC_EXT
        call mean_rec%write(fname, del_if_exists=.true.)
        call fname%kill
        ! back to the projectable expanded-Fourier state
        call mean_rec%fft
        call mean_rec%expand_exp
        write(logfhandle,'(A,I0,A,F8.1)') '>>> FLEX_PCA mean estimated from ',used, &
            &' particles, seconds=',toc(t_phase)
        call flush(logfhandle)
    end subroutine estimate_mean_from_data

    !> Worker-side mean scaling: apply the radial scale the MASTER fitted, rather than re-fitting it
    !! from this part's particles (which would use a different stride and hence a different subset).
    subroutine apply_cached_mean_scale( params, mean_rec, cache_fname )
        class(parameters),   intent(inout) :: params
        type(reconstructor), intent(inout) :: mean_rec
        !> per-fit namespace (paired engine); default flex_pca_mean_scale.bin
        character(len=*), optional, intent(in) :: cache_fname
        real, allocatable :: filt(:)
        integer :: nyq
        logical :: ok
        nyq = max(1, fdim(params%box_crop) - 1)
        allocate(filt(nyq))
        call read_mean_scale(nyq, filt, ok, cache_fname=cache_fname)
        if( .not. ok ) THROW_HARD('flex_pca worker found no mean-scale cache from the master')
        call mean_rec%apply_filter(filt)
        call mean_rec%expand_exp
        deallocate(filt)
    end subroutine apply_cached_mean_scale

    !> Self-estimate the amplitude scale of the consensus mean map relative to the whitened data, which
    !! carry SIMPLE's non-unitary gridding convention. A smoothed, clamped per-shell scale is applied to
    !! the mean so that y - T*mu is a residual rather than a difference of two amplitude conventions.
    subroutine estimate_mean_scale( params, build, mean_rec, pinds, nptcls, cache_fname , rounds)
        class(flex_pca_rounds), intent(inout) :: rounds
        class(parameters),   intent(inout) :: params
        type(builder),       intent(inout) :: build
        type(reconstructor), intent(inout) :: mean_rec
        integer,             intent(in)    :: pinds(:), nptcls
        !> per-fit cache namespace (paired distributed master writes one per fit)
        character(len=*), optional, intent(in) :: cache_fname
        integer, parameter :: NSAMPLE = 4000
        type(fplane_type), allocatable :: fpls(:)
        type(fplane_type), allocatable :: mean_fpl_t(:)
        type(ori),         allocatable :: ori_t(:)
        integer  :: batchlims(2), batchsz, ibatch, i, iptcl, stride, used, nyq, sh, j, nsub
        integer  :: nthr_here, ithr, t
        integer, allocatable :: sub_pinds(:), used_t(:)
        real(dp) :: s_my, s_mm, s, sm
        real(dp), allocatable :: smy_sh(:), smm_sh(:), sprof(:)
        real(dp), allocatable :: s_my_t(:), s_mm_t(:), smy_sh_t(:,:), smm_sh_t(:,:)
        real,     allocatable :: filt(:)
        logical  :: l_pcache
        nyq = max(1, fdim(params%box_crop) - 1)
        allocate(smy_sh(0:nyq), smm_sh(0:nyq), source=0.d0)
        stride = max(1, nptcls / NSAMPLE)
        s_my = 0.d0; s_mm = 0.d0; used = 0
        ! THREADED OVER PARTICLES: per-thread partial sums, folded in fixed thread order below.
        ! Reproducible at a given nthr; changing nthr moves the fitted scale only at rounding level.
        !$ call omp_set_max_active_levels(1)
        nthr_here = max(1, omp_get_max_threads())
        allocate(mean_fpl_t(nthr_here), ori_t(nthr_here), used_t(nthr_here))
        allocate(s_my_t(nthr_here), s_mm_t(nthr_here), source=0.d0)
        allocate(smy_sh_t(0:nyq,nthr_here), smm_sh_t(0:nyq,nthr_here), source=0.d0)
        used_t = 0
        call mean_rec%expand_exp
        ! one read path for every pass of the run: the downscaled cache when it is in use
        l_pcache = plane_cache_in_use(params, build)
        call init_rec(params, build, MAXIMGBATCHSZ, fpls, cropped=l_pcache)
        if( l_pcache )then
            call prepimgbatch(params, build, MAXIMGBATCHSZ, box=params%box_crop, smpd=params%smpd_crop)
        else
            call prepimgbatch(params, build, MAXIMGBATCHSZ)
        endif
        ! Select the strided sample UP FRONT, not inside the batch loop -- otherwise every particle is
        ! read, normalised, padded, FFT'd and CTF-evaluated before ~(1 - 1/stride) of that is discarded.
        ! Same particles in the same order as the serial code; only the summation grouping is per-thread.
        nsub = 0
        do j = 1, nptcls, stride
            nsub = nsub + 1
        end do
        allocate(sub_pinds(nsub))
        nsub = 0
        do j = 1, nptcls, stride
            nsub = nsub + 1
            sub_pinds(nsub) = pinds(j)
        end do
        do ibatch = 1, nsub, MAXIMGBATCHSZ
            batchlims = [ibatch, min(nsub, ibatch + MAXIMGBATCHSZ - 1)]
            batchsz   = batchlims(2) - batchlims(1) + 1
            call planes_batch_load(params, build, nsub, sub_pinds, batchlims, fpls, &
                &cov_image_mask_radius(params), l_pcache)
            !$omp parallel do default(shared) private(i,iptcl,ithr) schedule(static) proc_bind(close)
            do i = 1, batchsz
                ithr  = omp_get_thread_num() + 1
                iptcl = sub_pinds(batchlims(1)+i-1)
                call build%spproj_field%get_ori(iptcl, ori_t(ithr))
                if( ori_t(ithr)%isstatezero() ) cycle
                call project_fplane_mean(mean_rec, ori_t(ithr), fpls(i), mean_fpl_t(ithr), &
                    &apply_ctf_amp=.true.)
                s_my_t(ithr) = s_my_t(ithr) + real(cov_herm_inner(mean_fpl_t(ithr), fpls(i)), dp)
                s_mm_t(ithr) = s_mm_t(ithr) + real(cov_herm_inner(mean_fpl_t(ithr), mean_fpl_t(ithr)), dp)
                call plane_shell_cross_accum(mean_fpl_t(ithr), fpls(i), nyq, smy_sh_t(:,ithr), smm_sh_t(:,ithr))
                used_t(ithr) = used_t(ithr) + 1
            end do
            !$omp end parallel do
        end do
        deallocate(sub_pinds)
        do t = 1, nthr_here
            s_my   = s_my + s_my_t(t)
            s_mm   = s_mm + s_mm_t(t)
            smy_sh = smy_sh + smy_sh_t(:,t)
            smm_sh = smm_sh + smm_sh_t(:,t)
            used   = used + used_t(t)
            call ori_t(t)%kill
            call cleanup_plane(mean_fpl_t(t))
        end do
        deallocate(mean_fpl_t, ori_t, used_t, s_my_t, s_mm_t, smy_sh_t, smm_sh_t)
        call cleanup_rec_buffers(build, fpls)
        if( s_mm > DTINY )then
            s = s_my / s_mm
        else
            s = 1.d0
        endif
        if( s <= 0.d0 ) s = 1.d0
        write(logfhandle,'(A,ES12.4,A,I0,A)') '>>> FLEX_PCA mean amplitude self-scale s=',s, &
            &' (from ',used,' particles)'
        ! Per-shell mean/data amplitude scale
        allocate(sprof(0:nyq), filt(nyq))
        write(logfhandle,'(A)') '>>> FLEX_PCA per-shell mean scale s(sh) and s(sh)/s_global (D5):'
        do sh = 0, nyq
            if( smm_sh(sh) > DTINY )then
                sprof(sh) = smy_sh(sh) / smm_sh(sh)
            else
                sprof(sh) = s
            endif
            if( sh <= min(nyq,20) ) write(logfhandle,'(A,I3,A,ES11.3,A,F7.3)') '>>>   sh=',sh, &
                &'  s=',sprof(sh),'  ratio=',sprof(sh)/s
        end do
        call flush(logfhandle)
        ! 3-point-smoothed, clamped radial scale applied to the mean (FT state), then re-expand
        do sh = 1, nyq
            if( sh == 1 )then
                sm = 0.5d0*sprof(1) + 0.5d0*sprof(min(2,nyq))
            else if( sh == nyq )then
                sm = 0.5d0*sprof(nyq) + 0.5d0*sprof(nyq-1)
            else
                sm = 0.25d0*sprof(sh-1) + 0.5d0*sprof(sh) + 0.25d0*sprof(sh+1)
            endif
            filt(sh) = real(min(2.d0*s, max(0.5d0*s, sm)))
        end do
        call mean_rec%apply_filter(filt)
        call mean_rec%expand_exp
        if( rounds%is_master() ) call write_mean_scale(nyq, filt, cache_fname=cache_fname)
        deallocate(smy_sh, smm_sh, sprof, filt)
    end subroutine estimate_mean_scale

    !>  Accumulate per-shell mean/data cross power Re<T mu, y> and mean auto power |T mu|^2 over the
    !!  native k<=0 half. The per-shell ratio s(sh)=sum my_sh/sum mm_sh is the ML mean amplitude scale
    !!  at each shell; the k=0 double-count cancels in the ratio so no weighting is needed.
    subroutine plane_shell_cross_accum( mean_fpl, fpl, nyq, my_sh, mm_sh )
        type(fplane_type), intent(in)    :: mean_fpl, fpl
        integer,           intent(in)    :: nyq
        real(dp),          intent(inout) :: my_sh(0:), mm_sh(0:)
        integer     :: pf, h, k, hmin, hmax, kmin, kmax, sh
        complex(dp) :: m, y
        pf   = OSMPL_PAD_FAC
        hmin = pf*ceil_div(lbound(fpl%cmplx_plane,1),pf); hmax = pf*floor_div(ubound(fpl%cmplx_plane,1),pf)
        kmin = pf*ceil_div(lbound(fpl%cmplx_plane,2),pf); kmax = min(0, pf*floor_div(ubound(fpl%cmplx_plane,2),pf))
        do k = kmin, kmax, pf
            do h = hmin, hmax, pf
                sh = nint(sqrt(real((h/pf)**2 + (k/pf)**2)))
                if( sh > nyq ) cycle
                m = cmplx(mean_fpl%cmplx_plane(h,k), kind=dp)
                y = cmplx(fpl%cmplx_plane(h,k),      kind=dp)
                my_sh(sh) = my_sh(sh) + real(conjg(m)*y, dp)
                mm_sh(sh) = mm_sh(sh) + real(conjg(m)*m, dp)
            end do
        end do
    end subroutine plane_shell_cross_accum

    subroutine init_mean_reconstructor( params, build, mean_rec )
        class(parameters),  intent(inout) :: params
        type(builder),      intent(inout) :: build
        type(reconstructor),intent(inout) :: mean_rec
        type(image) :: meanvol
        ! alloc_rho() ends with reset(), which zeros the reconstructor's cmat (and, since rmat/cmat
        ! share the in-place FFT buffer, the real map too).
        call mean_rec%read_and_crop(params%vols(1),params%smpd,params%box_crop,params%smpd_crop)
        call mean_rec%alloc_rho(params,build%spproj,expand=.true.)
        call meanvol%read_and_crop(params%vols(1),params%smpd,params%box_crop,params%smpd_crop)
        call mean_rec%set_rmat(meanvol%get_rmat(), .false.)
        call meanvol%kill
        call mean_rec%fft
        call mean_rec%expand_exp
    end subroutine init_mean_reconstructor

    !> Reconstruction-mode plane from a whitened observation-model plane: numerator T*y and density |T|^2.
    subroutine form_reconstruction_plane( fpl, num )
        type(fplane_type), intent(in)    :: fpl
        type(fplane_type), intent(inout) :: num
        integer :: h1, h2, k1, k2
        h1 = lbound(fpl%cmplx_plane,1); h2 = ubound(fpl%cmplx_plane,1)
        k1 = lbound(fpl%cmplx_plane,2); k2 = ubound(fpl%cmplx_plane,2)
        if( .not. allocated(num%cmplx_plane) ) allocate(num%cmplx_plane(h1:h2,k1:k2))
        if( .not. allocated(num%ctfsq_plane) ) allocate(num%ctfsq_plane(h1:h2,k1:k2))
        num%cmplx_plane = conjg(fpl%transfer_plane) * fpl%cmplx_plane
        num%ctfsq_plane = fpl%ctfsq_plane
        num%frlims  = fpl%frlims
        num%nyq     = fpl%nyq
        num%shconst = fpl%shconst
    end subroutine form_reconstruction_plane

    !> Split half (1 or 2) of lattice point (ih,ik) of the unpadded plane, by an integer hash with no
    !! spatial structure, so half1 - half2 carries no coherent term. A structured split such as ih+ik
    !! parity acts as a half-box shift of y, a pose-dependent term entering the halves with opposite signs.
    pure integer function cov_half_parity( ih, ik ) result( par )
        integer, intent(in) :: ih, ik
        integer(kind=8) :: key
        key = int(ih,8)*73856093_8 + int(ik,8)*19349663_8 + 83492791_8
        key = iand(key, 2147483647_8)
        key = ieor(key, ishft(key, -15))
        key = iand(key*2654435761_8, 4294967295_8)
        key = ieor(key, ishft(key, -13))
        key = iand(key*97_8 + 13_8, 4294967295_8)
        par = int(iand(ishft(key, -9), 1_8)) + 1
    end function cov_half_parity

    !> Complex inner product over the native k<=0 half-plane (stored half-plane, k in [kmin,0]) inside the
    !! shared nyq disc; the optional half (1 or 2) restricts it to one cov_half_parity split half.
    function cov_herm_inner( lhs, rhs, half ) result( val )
        type(fplane_type), intent(in) :: lhs, rhs
        integer, optional, intent(in) :: half
        complex(dp) :: val, acc
        integer :: h, k, hmin, hmax, kmin, kmax, nyq_eff, pf, h_hi, hlf, par, nyq_disk, k_sq
        hlf = 0
        if( present(half) ) hlf = half
        acc = cmplx(0.d0,0.d0,dp)
        pf  = OSMPL_PAD_FAC
        nyq_eff = lhs%nyq
        if( rhs%nyq > 0 ) nyq_eff = min(nyq_eff, rhs%nyq)
        if( nyq_eff <= 0 ) nyq_eff = ubound(lhs%cmplx_plane,1)
        hmin = max(pf*ceil_div(lbound(lhs%cmplx_plane,1),pf), pf*ceil_div(-nyq_eff,pf))
        hmax = min(pf*floor_div(ubound(lhs%cmplx_plane,1),pf), pf*floor_div(nyq_eff,pf))
        kmin = max(pf*ceil_div(lbound(lhs%cmplx_plane,2),pf), pf*ceil_div(-nyq_eff,pf))
        kmax = min(0, pf*floor_div(nyq_eff,pf))
        ! Integer form of the shell test below. For integer x >= 0 and integer n >= 0,
        ! nint(sqrt(x)) > n  <=>  sqrt(x) >= n+0.5  <=>  x >= n^2+n+0.25  <=>  x > n*(n+1),
        ! so the disc gate selects exactly the same samples without a square root and a round per
        ! element. The embedding Gram alone calls this routine ncomp*(ncomp+1)/2 times per particle.
        nyq_disk = nyq_eff * (nyq_eff + 1)
        do k = kmin, kmax, pf
            ! the k=0 line is its own Friedel mate, so only h<=0 there, or it is counted twice
            h_hi = hmax
            if( k == 0 ) h_hi = 0
            k_sq = k*k
            if( k_sq > nyq_disk ) cycle
            do h = hmin, h_hi, pf
                if( h*h + k_sq > nyq_disk ) cycle
                if( hlf /= 0 )then
                    par = cov_half_parity(h/pf, k/pf)
                    if( par /= hlf ) cycle
                endif
                acc = acc + conjg(cmplx(lhs%cmplx_plane(h,k),kind=dp)) * cmplx(rhs%cmplx_plane(h,k),kind=dp)
            end do
        end do
        val = acc
    end function cov_herm_inner

    !> Closed-form per-particle ML contrast a = <T mu, y> / <T mu, T mu>, clamped to a sane range.
    real function particle_contrast( mean_fpl, fpl )
        type(fplane_type), intent(in) :: mean_fpl, fpl
        real(dp) :: emm, emy
        if( COV_UNIT_CONTRAST )then
            particle_contrast = 1.0
            return
        endif
        emm = real(cov_herm_inner(mean_fpl, mean_fpl), dp)
        emy = real(cov_herm_inner(mean_fpl, fpl), dp)
        particle_contrast = real(max(0.1d0, min(5.0d0, emy / max(emm, DTINY))))
    end function particle_contrast

    !> Whitened self-power and sample count of a plane in cov_herm_inner's index convention;
    !! em_calibrate_noise_prior uses the count to subtract the noise floor sig2*cnt.
    subroutine cov_herm_selfpower( fpl, pw, cnt )
        type(fplane_type), intent(in)  :: fpl
        real(dp),          intent(out) :: pw, cnt
        integer     :: h, k, hmin, hmax, kmin, kmax, nyq_eff, pf, h_hi, nyq_disk, k_sq
        complex(dp) :: c
        pf  = OSMPL_PAD_FAC
        pw  = 0.d0; cnt = 0.d0
        nyq_eff = fpl%nyq
        if( nyq_eff <= 0 ) nyq_eff = ubound(fpl%cmplx_plane,1)
        hmin = max(pf*ceil_div(lbound(fpl%cmplx_plane,1),pf), pf*ceil_div(-nyq_eff,pf))
        hmax = min(pf*floor_div(ubound(fpl%cmplx_plane,1),pf), pf*floor_div(nyq_eff,pf))
        kmin = max(pf*ceil_div(lbound(fpl%cmplx_plane,2),pf), pf*ceil_div(-nyq_eff,pf))
        kmax = min(0, pf*floor_div(nyq_eff,pf))
        nyq_disk = nyq_eff * (nyq_eff + 1)   ! see cov_herm_inner for why this replaces nint(sqrt(.))
        do k = kmin, kmax, pf
            h_hi = hmax
            if( k == 0 ) h_hi = 0
            k_sq = k*k
            if( k_sq > nyq_disk ) cycle
            do h = hmin, h_hi, pf
                if( h*h + k_sq > nyq_disk ) cycle
                c  = cmplx(fpl%cmplx_plane(h,k), kind=dp)
                pw = pw + real(c*conjg(c), dp)
                cnt= cnt + 1.d0
            end do
        end do
    end subroutine cov_herm_selfpower

    !> Signal-free whitened-noise variance per coefficient, from the high-frequency shells of a residual
    !! plane (where conformational signal is negligible).
    subroutine plane_hf_power( fpl, nyq, frac, pw, cnt )
        type(fplane_type), intent(in)  :: fpl
        integer,           intent(in)  :: nyq
        real,              intent(in)  :: frac
        real(dp),          intent(out) :: pw, cnt
        integer  :: pf, h, k, hmin, hmax, kmin, kmax, sh_lo, sh
        complex(dp) :: c
        pf   = OSMPL_PAD_FAC
        sh_lo= nint(frac*real(nyq))
        pw   = 0.d0; cnt = 0.d0
        hmin = pf*ceil_div(lbound(fpl%cmplx_plane,1),pf); hmax = pf*floor_div(ubound(fpl%cmplx_plane,1),pf)
        kmin = pf*ceil_div(lbound(fpl%cmplx_plane,2),pf); kmax = min(0, pf*floor_div(nyq,pf))
        do k = kmin, kmax, pf
            do h = hmin, hmax, pf
                sh = nint(sqrt(real((h/pf)**2 + (k/pf)**2)))
                if( sh < sh_lo .or. sh > nyq ) cycle
                c  = cmplx(fpl%cmplx_plane(h,k), kind=dp)
                pw = pw + real(c*conjg(c), dp)
                cnt= cnt + 1.d0
            end do
        end do
    end subroutine plane_hf_power

    !> The mean's radial amplitude scale. A worker MUST NOT re-fit this: estimate_mean_scale derives
    !! stride = max(1, nptcls/NSAMPLE) from the particle count it is given, so a worker holding a
    !! fraction of the particles would sample a different subset and fit a different scale. Shipping
    !! the nyq-length filter instead of the scaled volume keeps the handoff exact and tiny -- the
    !! worker rebuilds the mean deterministically from vol1 and applies the same array.
    subroutine write_mean_scale( nyq, filt, cache_fname )
        integer, intent(in) :: nyq
        real,    intent(in) :: filt(nyq)
        !> per-fit namespace (paired engine); default MEAN_SCALE_FNAME
        character(len=*), optional, intent(in) :: cache_fname
        type(string) :: fname, tmp_fname
        integer :: funit, io_stat
        fname     = string(MEAN_SCALE_FNAME)
        if( present(cache_fname) ) fname = trim(cache_fname)
        tmp_fname = fname//'.tmp'
        call fopen(funit, file=tmp_fname, access='STREAM', action='WRITE', status='REPLACE', iostat=io_stat)
        call fileiochk('write_mean_scale; open', io_stat)
        write(funit, iostat=io_stat) FLEX_PCA_PART_MAGIC, nyq
        call fileiochk('write_mean_scale; header', io_stat)
        write(funit, iostat=io_stat) filt
        call fileiochk('write_mean_scale; payload', io_stat)
        call fclose(funit)
        call simple_rename(tmp_fname, fname)
        call fname%kill; call tmp_fname%kill
    end subroutine write_mean_scale

    subroutine read_mean_scale( nyq, filt, ok, cache_fname )
        integer,           intent(in)  :: nyq
        real,              intent(out) :: filt(nyq)
        logical,           intent(out) :: ok
        character(len=*), optional, intent(in) :: cache_fname
        type(string) :: fname
        integer :: funit, io_stat, magic, nyq_in
        ok    = .false.
        fname = string(MEAN_SCALE_FNAME)
        if( present(cache_fname) ) fname = trim(cache_fname)
        if( .not. file_exists(fname) )then
            call fname%kill
            return
        endif
        call fopen(funit, file=fname, access='STREAM', action='READ', status='OLD', iostat=io_stat)
        call fileiochk('read_mean_scale; open', io_stat)
        read(funit, iostat=io_stat) magic, nyq_in
        call fileiochk('read_mean_scale; header', io_stat)
        if( magic /= FLEX_PCA_PART_MAGIC ) THROW_HARD('bad mean-scale magic')
        if( nyq_in /= nyq ) THROW_HARD('mean-scale band mismatch; master and worker disagree on box_crop')
        read(funit, iostat=io_stat) filt
        call fileiochk('read_mean_scale; payload', io_stat)
        call fclose(funit)
        ok = .true.
        call fname%kill
    end subroutine read_mean_scale

    !> Data-free EM start: lowest-|k| band lattice points (col_sep apart) realised as masked cos/sin
    !! pairs and orthonormalised; deterministic. One capped data pass calibrates sig2 and Gamma^0;
    !! the master also writes the initial eigenvolumes and the it000 copies the paired merge deflates against.
    subroutine init_basis_datafree( params, cfg, build, model, sel, col_sep, neigs_req, fprefix , rounds)
        type(flex_fit_model), intent(inout) :: model   !< mean in; basis, prior variances, rank, noise level out
        type(flex_selection), intent(in)    :: sel
        class(flex_pca_rounds), intent(inout) :: rounds
        class(parameters),   intent(inout) :: params
        type(flex_run_settings), intent(in)    :: cfg
        type(builder),       intent(inout) :: build
        integer,             intent(in) :: col_sep, neigs_req
        character(len=*), optional, intent(in) :: fprefix
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
        type(string) :: fname, pfx
        if( allocated(model%basis_recs) )then
            do q = 1, size(model%basis_recs)
                call model%basis_recs(q)%dealloc_rho; call model%basis_recs(q)%kill
            end do
            deallocate(model%basis_recs)
        endif
        if( allocated(model%eigvals) ) deallocate(model%eigvals)
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
            &d_tilde, svals, nptcls_basis=sel%nptcls)
        do s = 1, nreal
            call realvols(s)%kill
        end do
        deallocate(realvols)
        model%ncomp = max(1, min(neigs_req, d_tilde))
        call basis_recs_from_images(params, build, utilde_real(1:model%ncomp), model%ncomp, model%basis_recs)
        ! The eigenvolume MRCs are the master->worker handoff for every distributed probe round, so the
        ! initial basis is written here or the first PROBE round finds no flex_pca_pc001.mrc.
        pfx = 'flex_pca_pc'
        if( present(fprefix) ) pfx = trim(fprefix)
        if( .not. rounds%is_worker() )then
            do q = 1, model%ncomp
                fname = pfx//int2str_pad(q,3)//MRC_EXT
                call utilde_real(q)%write(fname, del_if_exists=.true.); call fname%kill
                ! stamp the INIT basis as it000: the merge's init-deflated matched cosines project
                ! the shared deterministic init out of both half bases before comparing
                fname = pfx//'it000_'//int2str_pad(q,3)//MRC_EXT
                call utilde_real(q)%write(fname, del_if_exists=.true.); call fname%kill
            end do
        endif
        call pfx%kill
        write(logfhandle,'(A,I0,A,I0,A,I0)') '>>> FLEX_PCA EM INIT (data-free): impulses=',ncol, &
            &'  representatives=',nreal,'  rank=',model%ncomp
        call flush(logfhandle)
        ! the only data pass (sig2 and Gamma^0), capped at COV_CALIB_MAX_PTCLS (SIMPLE_COV_CALIB_MAX);
        ! master only, so nparts=1 and the budget is not divided
        call cov_stage_subsample(build, sel%pinds, sel%nptcls, 1, cfg%calib_max, 'EM CALIBRATION', cpinds, ncal)
        call em_calibrate_noise_prior(params, build, model%mean_rec, model%basis_recs, model%ncomp, cpinds, ncal, &
            &model%sig2_eff, gam0)
        deallocate(cpinds)
        allocate(model%eigvals(model%ncomp), source=gam0)
        do s = 1, d_tilde
            call utilde(s)%dealloc_rho; call utilde(s)%kill
            call utilde_real(s)%kill
        end do
        deallocate(utilde, utilde_real)
        if( allocated(svals) ) deallocate(svals)
        call work%dealloc_rho; call work%kill
    end subroutine init_basis_datafree

    !> One data pass for the two scalars the EM cannot invent: sig2 from the shells above 0.7*Nyquist and
    !! Gamma^(0) from the mean-deflated residual minus that noise floor. Gamma^(0) over-estimates on purpose:
    !! a loose prior is corrected by the first M-step, a tight one shrinks z and so tightens itself.
    subroutine em_calibrate_noise_prior( params, build, mean_rec, basis_recs, ncomp, pinds, nptcls, &
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
        nthr    = omp_get_max_threads()
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

    !> Master -> probe-worker handoff of MODEL state only (dimension, noise level, prior variances).
    !! Round control (iteration, budget, fits, stage) travels in job_descr under registered keys.
    subroutine save_probe_state( ncomp, eigvals, sig2_eff, fname )
        integer,  intent(in) :: ncomp
        real(dp), intent(in) :: eigvals(:), sig2_eff
        character(len=*), optional, intent(in) :: fname
        type(string) :: fn
        integer :: funit, io_stat, q
        fn = COV_PROBE_META
        if( present(fname) ) fn = trim(fname)
        call fopen(funit, file=fn, action='WRITE', status='REPLACE', iostat=io_stat)
        call fileiochk('save_probe_state', io_stat)
        write(funit,*) ncomp
        write(funit,*) sig2_eff
        do q = 1, ncomp
            write(funit,*) eigvals(q)
        end do
        call fclose(funit)
        call fn%kill
    end subroutine save_probe_state

    subroutine load_probe_state( ncomp, eigvals, sig2_eff, fname )
        integer,               intent(out) :: ncomp
        real(dp), allocatable, intent(out) :: eigvals(:)
        real(dp),              intent(out) :: sig2_eff
        character(len=*), optional, intent(in) :: fname
        type(string) :: fn
        integer :: funit, io_stat, q
        fn = COV_PROBE_META
        if( present(fname) ) fn = trim(fname)
        if( .not. file_exists(fn) ) &
            &THROW_HARD('flex_pca probe worker found no '//fn%to_char()//' from the master')
        call fopen(funit, file=fn, action='READ', status='OLD', iostat=io_stat)
        call fileiochk('load_probe_state', io_stat)
        read(funit,*) ncomp
        read(funit,*) sig2_eff
        if( ncomp < 1 ) THROW_HARD('invalid cached probe basis dimension')
        allocate(eigvals(ncomp))
        do q = 1, ncomp
            read(funit,*) eigvals(q)
        end do
        call fclose(funit)
        call fn%kill
    end subroutine load_probe_state

    !>  Rebuild the projection-ready basis a probe worker needs from the master's flex_pca_pc*.mrc:
    !!  set_rmat then fft then expand_exp, never add(), which would leave the reconstructor flagged
    !!  Fourier and propagate an untransformed grid.
    subroutine load_probe_basis( params, build, ncomp, basis_recs, fprefix )
        class(parameters),   intent(inout) :: params
        type(builder),       intent(inout) :: build
        integer,             intent(in)    :: ncomp
        type(reconstructor), allocatable, intent(out) :: basis_recs(:)
        character(len=*), optional, intent(in) :: fprefix
        type(image)  :: vol
        type(string) :: fname, pfx
        integer      :: q
        pfx = 'flex_pca_pc'
        if( present(fprefix) ) pfx = trim(fprefix)
        allocate(basis_recs(ncomp))
        call vol%new([params%box_crop,params%box_crop,params%box_crop], params%smpd_crop)
        do q = 1, ncomp
            fname = pfx//int2str_pad(q,3)//MRC_EXT
            if( .not. file_exists(fname) ) &
                &THROW_HARD('flex_pca probe worker found no '//fname%to_char()//' from the master')
            call vol%read(fname)
            call init_basis_reconstructor(params, build, basis_recs(q))
            call basis_recs(q)%set_rmat(vol%get_rmat(), .false.)
            call basis_recs(q)%fft
            call basis_recs(q)%expand_exp
            call fname%kill
        end do
        call pfx%kill
        call vol%kill
    end subroutine load_probe_basis

    !> Deterministic low-frequency column selection: repeatedly take the lowest-|xi| candidate in the
    !! canonical Hermitian half that is at least col_sep away from every already-chosen column.
    subroutine select_frequencies_lowfreq( params, ncols_req, col_sep, col_hkl, ncol )
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

    subroutine pick_next_lowfreq( cand, ncand, chosen, nchosen, sep, best )
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

    function covariance_kfromto( params ) result( kfromto )
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

    !> Particle-image mask radius in pixels at params%box, or 0 to disable.
    function cov_image_mask_radius( params ) result( r )
        class(parameters), intent(in) :: params
        real :: r
        ! The compile-time default (COV_MASK_IMAGES) is OFF, which is
        ! safe ONLY when the solvent is pure noise -- true for synthetic data, FALSE for real data,
        ! where the region outside the envelope carries ice-thickness gradients, neighbouring
        ! particles and carbon edges. Those are low-frequency AND reproducible between halfsets, so
        ! they enter the covariance as apparent signal and can dominate the leading eigenvectors.
        r = 0.
        if( .not. COV_MASK_IMAGES ) return
        if( params%msk_crop <= 0. .or. params%box_crop <= 0 ) return
        r = COV_MASK_MARGIN * params%msk_crop * real(params%box) / real(params%box_crop)
        r = min(r, 0.5*real(params%box) - COSMSKHALFWIDTH - 1.)
    end function cov_image_mask_radius

    !> Convert each merged complex column C_q into its two real spatial representatives Re(ifft
    !! C_q)=Sigma*cos_q and Im(ifft C_q)=Sigma*sin_q.
    subroutine basis_to_real_representatives( params, work, colvol, ncol, lb, ub, realvols, nreal )
        class(parameters),   intent(inout) :: params
        type(reconstructor), intent(inout) :: work
        complex,             intent(in)    :: colvol(:,:,:,:)
        integer,             intent(in)    :: ncol, lb(3), ub(3)
        type(image), allocatable, intent(out) :: realvols(:)
        integer,                  intent(out) :: nreal
        type(image)  :: gridcorr_img
        complex, allocatable :: vr(:,:,:), vi(:,:,:)
        integer :: s, i1, i2, i3, n1, n2, n3, hn, kn, ln, h, k, l
        real    :: energy
        n1 = ub(1)-lb(1)+1; n2 = ub(2)-lb(2)+1; n3 = ub(3)-lb(3)+1
        allocate(vr(lb(1):ub(1),lb(2):ub(2),lb(3):ub(3)))
        allocate(vi(lb(1):ub(1),lb(2):ub(2),lb(3):ub(3)))
        gridcorr_img = prep3D_inv_kbenvelope4mul([params%box_crop,params%box_crop,params%box_crop], params%smpd_crop)
        allocate(realvols(2*ncol))
        nreal = 0
        do s = 1, ncol
            do i3 = 1, n3
                l = lb(3)+i3-1; ln = -l
                do i2 = 1, n2
                    k = lb(2)+i2-1; kn = -k
                    do i1 = 1, n1
                        h = lb(1)+i1-1; hn = -h
                        if( hn < lb(1) .or. hn > ub(1) .or. kn < lb(2) .or. kn > ub(2) &
                            &.or. ln < lb(3) .or. ln > ub(3) )then
                            vr(h,k,l) = 0.5*colvol(i1,i2,i3,s)
                            vi(h,k,l) = cmplx(0.,-0.5)*colvol(i1,i2,i3,s)
                        else
                            vr(h,k,l) = 0.5*(colvol(i1,i2,i3,s) + conjg(colvol(n1-i1+1,n2-i2+1,n3-i3+1,s)))
                            vi(h,k,l) = cmplx(0.,-0.5)*(colvol(i1,i2,i3,s) - conjg(colvol(n1-i1+1,n2-i2+1,n3-i3+1,s)))
                        endif
                    end do
                end do
            end do
            call realize_hermitian_volume(params, work, vr, gridcorr_img, energy)
            if( energy > 0. )then
                nreal = nreal + 1
                call realvols(nreal)%copy(work)
            endif
            call realize_hermitian_volume(params, work, vi, gridcorr_img, energy)
            if( energy > 0. )then
                nreal = nreal + 1
                call realvols(nreal)%copy(work)
            endif
        end do
        call gridcorr_img%kill
        deallocate(vr, vi)
    end subroutine basis_to_real_representatives

    !>  Load a Hermitian expanded Fourier volume into the work reconstructor, fold to
    !!  compressed storage, inverse-FFT to a real volume, deapodize, low-pass and mask.
    subroutine realize_hermitian_volume( params, work, vherm, gridcorr_img, energy )
        class(parameters),   intent(in)    :: params
        type(reconstructor), intent(inout) :: work
        complex,             intent(in)    :: vherm(:,:,:)
        type(image),         intent(inout) :: gridcorr_img
        real,                intent(out)   :: energy
        real, pointer :: rmat(:,:,:)
        integer       :: ldim_work(3)
        call work%reset
        call work%reset_exp
        work%cmat_exp = vherm
        call work%compress_exp
        ! Band-limit the covariance column to the signal band FIRST, in Fourier space, before any real-space
        ! operation, so out-of-band shells never enter the masking or the Gram products downstream.
        if( params%lp > 2.0*params%smpd_crop + TINY ) call work%bp(0., params%lp)
        call work%ifft
        ! deapodize on the native lattice (same correction as production gridding)
        call work%mul(gridcorr_img)
        ! the envelope window when pcg_mskfile is set (the sphere otherwise): the data-free representatives
        ! must live on the same support as every later basis, or the calibration Gram and the whole
        ! latent scale differ from the pcg line by the sphere-to-envelope volume ratio
        call flex_window_apply_rec(work, params)
        if( work%is_ft() ) call work%ifft
        call work%get_rmat_ptr(rmat)
        ldim_work = work%get_ldim()
        energy = sum(rmat(1:ldim_work(1),1:ldim_work(2),1:ldim_work(3))**2)
    end subroutine realize_hermitian_volume

    !> Orthonormalize the real column representatives into the column subspace Utilde by Gram
    !! eigendecomposition, keeping every direction above a relative energy floor.
    subroutine orthonormalize_representatives( params, build, realvols, nreal, utilde, utilde_real, d_tilde, svals, &
        &nptcls_basis )
        class(parameters),   intent(inout) :: params
        type(builder),       intent(inout) :: build
        type(image),         intent(inout) :: realvols(:)
        integer,             intent(in)    :: nreal
        type(reconstructor), allocatable, intent(out) :: utilde(:)
        type(image),         allocatable, intent(out) :: utilde_real(:)
        integer,             intent(out)   :: d_tilde
        !> squared singular values of the representative set.
        real(dp), allocatable, optional, intent(out) :: svals(:)
        !> particles that actually reached the basis stages. Only used to REPORT the
        !! samples-per-parameter rank bound; it selects nothing.
        integer, optional, intent(in) :: nptcls_basis
        real(dp), allocatable :: gram(:,:), evec(:,:), eval(:)
        real, pointer :: rmat_i(:,:,:), rmat_j(:,:,:)
        integer :: i, q, nrot, keep, d_budget, d_cap, d_signal, d_samples, nbasis
        real(dp) :: lam_max, nrm
        character(len=9) :: accum_model
        if( nreal < 1 ) THROW_HARD('flex_pca produced no covariance column representatives')
        allocate(gram(nreal,nreal), evec(nreal,nreal), eval(nreal))
        do i = 1, nreal
            call realvols(i)%get_rmat_ptr(rmat_i)
            do q = i, nreal
                call realvols(q)%get_rmat_ptr(rmat_j)
                gram(i,q) = sum(real(rmat_i,dp)*real(rmat_j,dp))
                gram(q,i) = gram(i,q)
            end do
        end do
        call jacobi(gram, nreal, nreal, eval, evec, nrot)
        call eigsrt(eval, evec, nreal, nreal)              ! descending
        lam_max = max(eval(1), DTINY)
        keep = 0
        do q = 1, nreal
            if( eval(q) > COV_EIG_REL_FLOOR*lam_max ) keep = keep + 1
        end do
        ! memory cap from COV_ATHR_BUDGET under the packed 8*[d(d+1)/2]^2-byte model (cov_dim_budget); no
        ! current solve forms that array, and with the shipped constants COV_DEFAULT_DTILDE binds first
        d_budget = cov_dim_budget()
        ! data-driven rank, REPORT ONLY: the energy floor and the memory budget never ask how many
        ! directions are real, so log what the data would support and let the discrepancy show
        d_signal = cov_signal_rank(eval, nreal)
        nbasis   = 0
        if( present(nptcls_basis) ) nbasis = nptcls_basis
        d_samples = 0
        if( nbasis > 0 ) d_samples = &
            &max(1, int((-1.d0 + sqrt(1.d0 + 8.d0*real(nbasis,dp)/COV_SAMPLES_PER_PARAM))/2.d0))
        ! memory budget is a GUARD, not the chooser: COV_DEFAULT_DTILDE sets the rank unless the energy
        ! floor or the budget is lower.
        d_cap = min(d_budget, COV_DEFAULT_DTILDE)
        d_tilde  = max(1, min(keep, COV_MAX_DTILDE, d_cap))
        if( d_samples > 0 )then
            write(logfhandle,'(A,I0,A,I0,A,F4.1,A)') '>>> FLEX_PCA d_samples=',d_samples, &
                &'  (samples-per-parameter bound from N=',nbasis,' at R=',COV_SAMPLES_PER_PARAM, &
                &') -- REPORT ONLY'
        endif
        write(logfhandle,'(A,I0,A)') '>>> FLEX_PCA d_signal=',d_signal, &
            &' (spectrum noise-bulk estimate; report only)'
        accum_model = 'packed+CG'
        write(logfhandle,'(A,I0,A,I0,A,I0,A,I0,A,I0,A)') '>>> FLEX_PCA d_tilde=',d_tilde, &
            &'  (above energy floor=',keep,', memory cap=',d_budget,', rank cap=',COV_MAX_DTILDE, &
            &', default=',COV_DEFAULT_DTILDE,')'
        write(logfhandle,'(A,A,A,F8.3,A,F6.3,A)') '>>> FLEX_PCA d_tilde memory-cap model: ', &
            &trim(accum_model),', ',cov_accum_bytes(d_tilde)/1.d9, &
            &' GB at this d_tilde (budget ',COV_ATHR_BUDGET/1.d9,' GB)'
        if( d_tilde == d_budget .and. keep > d_budget )then
            write(logfhandle,'(A,I0,A)') '>>> FLEX_PCA NOTE: the column subspace is limited by the &
                &d_tilde memory budget, not by the data; ',keep,' directions cleared the energy floor.'
        else if( d_tilde == COV_DEFAULT_DTILDE .and. d_budget > COV_DEFAULT_DTILDE )then
            write(logfhandle,'(A,I0,A)') '>>> FLEX_PCA NOTE: d_tilde is the measured default; the &
                &memory budget would have allowed ',d_budget,'.'
        endif
        call flush(logfhandle)
        allocate(utilde(d_tilde), utilde_real(d_tilde))
        do q = 1, d_tilde
            nrm = sqrt(max(eval(q), DTINY))
            ! Build the unit-norm orthonormal basis vector as a CLEAN plain image via image arithmetic on
            ! the (verified band-limited) representatives, rather than raw rmat-pointer math on padded
            ! reconstructor buffers.
            call utilde_real(q)%copy(realvols(1))
            call utilde_real(q)%zero_and_unflag_ft
            do i = 1, nreal
                call utilde_real(q)%add(realvols(i), real(evec(i,q)/nrm))
            end do
            ! set_rmat(...,.false.) then fft then expand_exp -- NEVER add(), which leaves the reconstructor
            ! flagged as Fourier and silently propagates an untransformed grid
            call init_basis_reconstructor(params, build, utilde(q))
            call utilde(q)%set_rmat(utilde_real(q)%get_rmat(), .false.)
            call utilde(q)%fft
            call utilde(q)%expand_exp
        end do
        if( present(svals) )then
            allocate(svals(d_tilde))
            do q = 1, d_tilde
                svals(q) = max(eval(q), DTINY)
            end do
        endif
        deallocate(gram, evec, eval)
    end subroutine orthonormalize_representatives

    !>  M(i,j) = <U_ref_i, U_tgt_j> over unit-normed volumes (per-vector only, so both sets should be
    !!  orthonormal); z_ref = M z_tgt maps a target-basis latent into the reference frame. svals are
    !!  the singular values of M, i.e. the principal-angle cosines between the two subspaces.
    subroutine align_basis_to_reference( ref_imgs, nref_c, tgt_imgs, ntgt_c, M, svals )
        integer,     intent(in)    :: nref_c, ntgt_c
        type(image), intent(inout) :: ref_imgs(nref_c), tgt_imgs(ntgt_c)
        real(dp), allocatable, intent(out) :: M(:,:), svals(:)
        real, pointer :: rmat_i(:,:,:), rmat_j(:,:,:)
        real(dp), allocatable :: nrm_r(:), nrm_t(:), Mwork(:,:), V(:,:), ev(:)
        integer  :: i, j, nrot, nsv
        allocate(M(nref_c,ntgt_c), source=0.d0)
        allocate(nrm_r(nref_c), nrm_t(ntgt_c), source=0.d0)
        do i = 1, nref_c
            call ref_imgs(i)%get_rmat_ptr(rmat_i)
            nrm_r(i) = sqrt(max(sum(real(rmat_i,dp)*real(rmat_i,dp)), DTINY))
        end do
        do j = 1, ntgt_c
            call tgt_imgs(j)%get_rmat_ptr(rmat_j)
            nrm_t(j) = sqrt(max(sum(real(rmat_j,dp)*real(rmat_j,dp)), DTINY))
        end do
        do i = 1, nref_c
            call ref_imgs(i)%get_rmat_ptr(rmat_i)
            do j = 1, ntgt_c
                call tgt_imgs(j)%get_rmat_ptr(rmat_j)
                M(i,j) = sum(real(rmat_i,dp)*real(rmat_j,dp)) / (nrm_r(i)*nrm_t(j))
            end do
        end do
        ! principal-angle cosines = singular values of M, via the eigenvalues of M^T M
        nsv = min(nref_c, ntgt_c)
        allocate(Mwork(ntgt_c,ntgt_c), V(ntgt_c,ntgt_c), ev(ntgt_c), svals(nsv))
        Mwork = matmul(transpose(M), M)
        call jacobi(Mwork, ntgt_c, ntgt_c, ev, V, nrot)
        call eigsrt(ev, V, ntgt_c, ntgt_c)
        do i = 1, nsv
            svals(i) = sqrt(max(ev(i), 0.d0))
        end do
        deallocate(nrm_r, nrm_t, Mwork, V, ev)
    end subroutine align_basis_to_reference

    !> Turn a set of real-space basis volumes into embedding-ready column reconstructors.
    subroutine basis_recs_from_images( params, build, imgs, ncomp, basis_recs )
        class(parameters),   intent(inout) :: params
        type(builder),       intent(inout) :: build
        integer,             intent(in)    :: ncomp
        type(image),         intent(inout) :: imgs(ncomp)
        type(reconstructor), allocatable, intent(out) :: basis_recs(:)
        integer :: q
        allocate(basis_recs(ncomp))
        do q = 1, ncomp
            call init_basis_reconstructor(params, build, basis_recs(q))
            call basis_recs(q)%set_rmat(imgs(q)%get_rmat(), .false.)
            call basis_recs(q)%fft
            call basis_recs(q)%expand_exp
        end do
    end subroutine basis_recs_from_images

    !> Project a stack of volumes onto the orthogonal complement of an ORTHONORMAL basis stack.
    !! Used to isolate what an EM update ADDS to the current basis, so two half-updates can be
    !! compared without the basis they share dominating the comparison.
    subroutine deflate_against_basis( imgs, n, basis, nb )
        integer,     intent(in)    :: n, nb
        type(image), intent(inout) :: imgs(n), basis(nb)
        real, pointer :: rv(:,:,:), rb(:,:,:)
        real(dp) :: c, bb
        integer  :: i, k
        do k = 1, nb
            call basis(k)%get_rmat_ptr(rb)
            bb = sum(real(rb,dp)*real(rb,dp))
            if( bb <= DTINY ) cycle
            do i = 1, n
                call imgs(i)%get_rmat_ptr(rv)
                c  = sum(real(rv,dp)*real(rb,dp)) / bb
                rv = rv - real(c)*rb
            end do
        end do
    end subroutine deflate_against_basis

    !> Principal angles between the subspaces spanned by two NON-orthonormal volume stacks.
    !!
    !! align_basis_to_reference only per-vector normalises, which is exact only when both inputs
    !! are already orthonormal -- it says so itself. The even/odd bases coming out of the coupled
    !! solve are not, so this computes the honest quantity: for spans E and O the principal-angle
    !! cosines are the singular values of (E'E)^-1/2 (E'O) (O'O)^-1/2. Everything is done in the
    !! n x n Gram algebra, so the volume work is three symmetric Gram products and nothing else.
    subroutine cross_half_subspace_angles( eimgs, oimgs, n, svals )
        integer,     intent(in)    :: n
        type(image), intent(inout) :: eimgs(n), oimgs(n)
        real(dp), allocatable, intent(out) :: svals(:)
        real, pointer :: ri(:,:,:), rj(:,:,:)
        real(dp), allocatable :: Gee(:,:), Goo(:,:), Geo(:,:), We(:,:), Wo(:,:)
        real(dp), allocatable :: M(:,:), MtM(:,:), V2(:,:), ev2(:)
        integer  :: i, j, nrot
        real(dp) :: lam
        allocate(Gee(n,n), Goo(n,n), Geo(n,n), source=0.d0)
        do i = 1, n
            call eimgs(i)%get_rmat_ptr(ri)
            do j = i, n
                call eimgs(j)%get_rmat_ptr(rj)
                Gee(i,j) = sum(real(ri,dp)*real(rj,dp)); Gee(j,i) = Gee(i,j)
            end do
            do j = 1, n
                call oimgs(j)%get_rmat_ptr(rj)
                Geo(i,j) = sum(real(ri,dp)*real(rj,dp))
            end do
        end do
        do i = 1, n
            call oimgs(i)%get_rmat_ptr(ri)
            do j = i, n
                call oimgs(j)%get_rmat_ptr(rj)
                Goo(i,j) = sum(real(ri,dp)*real(rj,dp)); Goo(j,i) = Goo(i,j)
            end do
        end do
        ! inverse square roots by symmetric eigendecomposition, with a relative floor so a
        ! collapsed half-basis direction cannot blow the whitening up
        allocate(We(n,n), Wo(n,n), source=0.d0)
        call inv_sqrt_sym(Gee, n, We)
        call inv_sqrt_sym(Goo, n, Wo)
        allocate(M(n,n), MtM(n,n), V2(n,n), ev2(n), svals(n))
        M   = matmul(We, matmul(Geo, Wo))
        MtM = matmul(transpose(M), M)
        call jacobi(MtM, n, n, ev2, V2, nrot)
        call eigsrt(ev2, V2, n, n)
        do i = 1, n
            svals(i) = sqrt(max(0.d0, min(1.d0, ev2(i))))
        end do
        deallocate(Gee, Goo, Geo, We, Wo, M, MtM, V2, ev2)
      contains
        subroutine inv_sqrt_sym( A, m, Ainvsq )
            integer,  intent(in)  :: m
            real(dp), intent(in)  :: A(m,m)
            real(dp), intent(out) :: Ainvsq(m,m)
            real(dp), allocatable :: Aw(:,:), Vv(:,:), ee(:)
            integer :: k, nr
            real(dp) :: mx
            allocate(Aw(m,m), Vv(m,m), ee(m))
            Aw = A
            call jacobi(Aw, m, m, ee, Vv, nr)
            mx = maxval(ee)
            Ainvsq = 0.d0
            do k = 1, m
                lam = ee(k)
                if( lam > 1.d-10*max(mx,DTINY) )then
                    Ainvsq = Ainvsq + matmul(reshape(Vv(:,k),[m,1]), &
                        &reshape(Vv(:,k),[1,m])) / sqrt(lam)
                endif
            end do
            deallocate(Aw, Vv, ee)
        end subroutine inv_sqrt_sym
    end subroutine cross_half_subspace_angles

end module simple_flex_pca_basis
