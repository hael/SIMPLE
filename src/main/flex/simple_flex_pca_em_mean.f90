!@descr: flex_pca EM: consensus mean estimation, mean scale and Fourier-plane accumulation
submodule (simple_flex_pca_em) simple_flex_pca_em_mean
use simple_flex_pca_distr,  only: flex_pca_is_master, flex_pca_is_worker
use simple_flex_pca_parts,  only: write_mean_scale, read_mean_scale
use simple_matcher_3Drec,   only: init_rec, cleanup_rec_buffers
use simple_matcher_ptcl_io, only: discrete_read_imgbatch, prepimgbatch
use simple_flex_reconstructor_latent_ops, only: project_fplane_mean
use simple_flex_projected_latent_model, only: prep_imgs4projected_model
implicit none
#include "simple_local_flags.inc"

contains

    !> Single entry point for the covariance mean.
    module subroutine estimate_covariance_mean( params, build, mean_rec, pinds, nptcls )
        class(parameters),   intent(inout) :: params
        type(builder),       intent(inout) :: build
        type(reconstructor), intent(inout) :: mean_rec
        integer,             intent(in)    :: pinds(:), nptcls
        integer :: vfc
        ! one-time runtime init of the fitted-contrast override (single-threaded here; every stage
        ! of every process passes through the mean before touching particles)
        vfc = 0
        call cov_env_int('SIMPLE_COV_CONTRAST', vfc)
        cov_fit_contrast_rt = vfc > 0
        if( cov_fit_contrast_rt )then
            write(logfhandle,'(A)') '>>> FLEX_PCA PER-PARTICLE CONTRAST ON (deviations against &
                &fitted a_i; kills the rank-one contrast mode -- RECOVAR A.4)'
            call flush(logfhandle)
        endif
        if( COV_MEAN_FROM_DATA )then
            call estimate_mean_from_data(params, build, mean_rec, pinds, nptcls)
        else
            call init_mean_reconstructor(params, build, mean_rec)
            if( flex_pca_is_worker() )then
                call apply_cached_mean_scale(params, mean_rec)
            else
                call estimate_mean_scale(params, build, mean_rec, pinds, nptcls)
            endif
        endif
    end subroutine estimate_covariance_mean

    !> Kernel-regression consensus mean (eq. S.1) estimated from the particles themselves, as an
    !! alternative to reading the supplied consensus volume. Selected by COV_MEAN_FROM_DATA.
    module subroutine estimate_mean_from_data( params, build, mean_rec, pinds, nptcls )
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
        logical :: l_devprep
        integer(timer_int_kind) :: t_phase
        call init_basis_reconstructor(params, build, mean_rec)
        call init_rec(params, build, MAXIMGBATCHSZ, fpls)
        call prepimgbatch(params, build, MAXIMGBATCHSZ)
        call cov_dev_prep_start(params, build, l_devprep)
        used    = 0
        t_phase = tic()
        write(logfhandle,'(A)') '>>> FLEX_PCA MEAN ESTIMATION (eq. S.1 kernel regression)'
        call flush(logfhandle)
        do ibatch = 1, nptcls, MAXIMGBATCHSZ
            batchlims = [ibatch, min(nptcls, ibatch + MAXIMGBATCHSZ - 1)]
            batchsz   = batchlims(2) - batchlims(1) + 1
            call discrete_read_imgbatch(params, build, nptcls, pinds, batchlims)
            call prep_imgs4projected_model(params, build, batchsz, build%imgbatch(:batchsz), &
                &pinds(batchlims(1):batchlims(2)), fpls(:batchsz), mskrad=cov_image_mask_radius(params))
            do i = 1, batchsz
                iptcl = pinds(batchlims(1)+i-1)
                call build%spproj_field%get_ori(iptcl, orientation)
                if( orientation%isstatezero() ) cycle
                call form_reconstruction_plane(fpls(i), num_fpl)
                call mean_rec%insert_plane_oversamp(build%pgrpsyms, orientation, num_fpl)
                used = used + 1
            end do
        end do
        call cov_dev_prep_stop(l_devprep)
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
    module subroutine apply_cached_mean_scale( params, mean_rec )
        class(parameters),   intent(inout) :: params
        type(reconstructor), intent(inout) :: mean_rec
        real, allocatable :: filt(:)
        integer :: nyq
        logical :: ok
        nyq = max(1, fdim(params%box_crop) - 1)
        allocate(filt(nyq))
        call read_mean_scale(nyq, filt, ok)
        if( .not. ok ) THROW_HARD('flex_pca worker found no flex_pca_mean_scale.bin from the master')
        call mean_rec%apply_filter(filt)
        call mean_rec%expand_exp
        deallocate(filt)
    end subroutine apply_cached_mean_scale

    !> Self-estimate the amplitude scale of the consensus mean map relative to the whitened data, which
    !! carry SIMPLE's non-unitary gridding convention. A smoothed, clamped per-shell scale is applied to
    !! the mean so that y - T*mu is a residual rather than a difference of two amplitude conventions.
    module subroutine estimate_mean_scale( params, build, mean_rec, pinds, nptcls )
        class(parameters),   intent(inout) :: params
        type(builder),       intent(inout) :: build
        type(reconstructor), intent(inout) :: mean_rec
        integer,             intent(in)    :: pinds(:), nptcls
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
        logical  :: l_devprep
        nyq = max(1, fdim(params%box_crop) - 1)
        allocate(smy_sh(0:nyq), smm_sh(0:nyq), source=0.d0)
        stride = max(1, nptcls / NSAMPLE)
        s_my = 0.d0; s_mm = 0.d0; used = 0
        ! THREADED OVER PARTICLES: per-thread partial sums, folded in fixed thread order below.
        ! Reproducible at a given nthr; changing nthr moves the fitted scale only at rounding level.
        !$ call omp_set_max_active_levels(1)
        nthr_here = max(1, params%nthr)
        allocate(mean_fpl_t(nthr_here), ori_t(nthr_here), used_t(nthr_here))
        allocate(s_my_t(nthr_here), s_mm_t(nthr_here), source=0.d0)
        allocate(smy_sh_t(0:nyq,nthr_here), smm_sh_t(0:nyq,nthr_here), source=0.d0)
        used_t = 0
        call mean_rec%expand_exp
        call init_rec(params, build, MAXIMGBATCHSZ, fpls)
        call prepimgbatch(params, build, MAXIMGBATCHSZ)
        call cov_dev_prep_start(params, build, l_devprep)
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
            call discrete_read_imgbatch(params, build, nsub, sub_pinds, batchlims)
            call prep_imgs4projected_model(params, build, batchsz, build%imgbatch(:batchsz), &
                &sub_pinds(batchlims(1):batchlims(2)), fpls(:batchsz), mskrad=cov_image_mask_radius(params))
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
        call cov_dev_prep_stop(l_devprep)
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
        if( flex_pca_is_master() ) call write_mean_scale(nyq, filt)
        deallocate(smy_sh, smm_sh, sprof, filt)
    end subroutine estimate_mean_scale
    !>  Accumulate per-shell mean/data cross power Re<T mu, y> and mean auto power |T mu|^2 over the
    !!  native k<=0 half. The per-shell ratio s(sh)=sum my_sh/sum mm_sh is the ML mean amplitude scale
    !!  at each shell; the k=0 double-count cancels in the ratio so no weighting is needed.
    module subroutine plane_shell_cross_accum( mean_fpl, fpl, nyq, my_sh, mm_sh )
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
    module subroutine init_mean_reconstructor( params, build, mean_rec )
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
    module subroutine form_reconstruction_plane( fpl, num )
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

    !> Complex Hermitian-half plane inner product over the native k<=0 half (the planes are stored
    !! half-plane, k in [kmin,0]). Properness is handled by the caller: reduced_covariance_solve
    !! accumulates Re(b_q)Re(b_r), not |b|^2, and pairs it with the 0.5*sig2 noise scaling, since
    !! E[Re(b_q)Re(b_r)]_noise = 0.5*sig2*G. The optional half selects a checkerboard sub-half.
    !> The (h,k) sample sequence cov_herm_inner visits, in exactly its order, for a plane pair that
    !! shares nyq and bounds with `ref` and is called without the `half` selector. Returns nsamp = -1
    !! if the sequence does not fit `slist`, in which case the caller must use cov_herm_inner directly.
    module subroutine cov_herm_sample_list( ref, slist, nsamp )
        type(fplane_type), intent(in)  :: ref
        integer,           intent(out) :: slist(:,:)
        integer,           intent(out) :: nsamp
        integer :: h, k, hmin, hmax, kmin, kmax, nyq_eff, pf, h_hi, nyq_disk, k_sq
        pf      = OSMPL_PAD_FAC
        nyq_eff = ref%nyq
        if( nyq_eff <= 0 ) nyq_eff = ubound(ref%cmplx_plane,1)
        hmin = max(pf*ceil_div(lbound(ref%cmplx_plane,1),pf), pf*ceil_div(-nyq_eff,pf))
        hmax = min(pf*floor_div(ubound(ref%cmplx_plane,1),pf), pf*floor_div(nyq_eff,pf))
        kmin = max(pf*ceil_div(lbound(ref%cmplx_plane,2),pf), pf*ceil_div(-nyq_eff,pf))
        kmax = min(0, pf*floor_div(nyq_eff,pf))
        nyq_disk = nyq_eff * (nyq_eff + 1)
        nsamp    = 0
        do k = kmin, kmax, pf
            h_hi = hmax
            if( k == 0 ) h_hi = 0
            k_sq = k*k
            if( k_sq > nyq_disk ) cycle
            do h = hmin, h_hi, pf
                if( h*h + k_sq > nyq_disk ) cycle
                nsamp = nsamp + 1
                if( nsamp > size(slist,2) )then
                    nsamp = -1
                    return
                endif
                slist(1,nsamp) = h
                slist(2,nsamp) = k
            end do
        end do
    end subroutine cov_herm_sample_list

    module function cov_herm_inner( lhs, rhs, half ) result( val )
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
        ! element. This routine is called d_tilde*(d_tilde+1)/2 times per particle in the reduced solve.
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
                    par = modulo((h/pf) + (k/pf), 2) + 1
                    if( par /= hlf ) cycle
                endif
                acc = acc + conjg(cmplx(lhs%cmplx_plane(h,k),kind=dp)) * cmplx(rhs%cmplx_plane(h,k),kind=dp)
            end do
        end do
        val = acc
    end function cov_herm_inner

    !> Closed-form per-particle ML contrast a = <T mu, y> / <T mu, T mu>, clamped to a sane range.
    real module function particle_contrast( mean_fpl, fpl )
        type(fplane_type), intent(in) :: mean_fpl, fpl
        real(dp) :: emm, emy
        if( COV_UNIT_CONTRAST .and. .not. cov_fit_contrast_rt )then
            particle_contrast = 1.0
            return
        endif
        emm = real(cov_herm_inner(mean_fpl, mean_fpl), dp)
        emy = real(cov_herm_inner(mean_fpl, fpl), dp)
        particle_contrast = real(max(0.1d0, min(5.0d0, emy / max(emm, DTINY))))
    end function particle_contrast

    !> Whitened self-power of a plane in cov_herm_inner's index convention. The reduced-solve debias
    !! assumes E[Re(b_q)Re(b_r)]_noise = 0.5*sig2*G with b and G from cov_herm_inner, so the sig2 fed to
    !! it has to be the noise variance measured in THAT convention; this is what the log compares against.
    module subroutine cov_herm_selfpower( fpl, pw, cnt )
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
    module subroutine plane_hf_power( fpl, nyq, frac, pw, cnt )
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

    !> Accumulate whitened plane power and voxel count per radial shell (native indices) into the passed
    !! shell-indexed arrays.
    module subroutine plane_shell_power_accum( fpl, nyq, pw_sh, cnt_sh )
        type(fplane_type), intent(in)    :: fpl
        integer,           intent(in)    :: nyq
        real(dp),          intent(inout) :: pw_sh(0:), cnt_sh(0:)
        integer  :: pf, h, k, hmin, hmax, kmin, kmax, sh
        complex(dp) :: c
        pf   = OSMPL_PAD_FAC
        hmin = pf*ceil_div(lbound(fpl%cmplx_plane,1),pf); hmax = pf*floor_div(ubound(fpl%cmplx_plane,1),pf)
        kmin = pf*ceil_div(lbound(fpl%cmplx_plane,2),pf); kmax = min(0, pf*floor_div(nyq,pf))
        do k = kmin, kmax, pf
            do h = hmin, hmax, pf
                sh = nint(sqrt(real((h/pf)**2 + (k/pf)**2)))
                if( sh > nyq ) cycle
                c = cmplx(fpl%cmplx_plane(h,k), kind=dp)
                pw_sh(sh)  = pw_sh(sh)  + real(c*conjg(c), dp)
                cnt_sh(sh) = cnt_sh(sh) + 1.d0
            end do
        end do
    end subroutine plane_shell_power_accum

    !>  Whitened residual plane r_i = data - transfer * P_i mu (mean_fpl already carries
    !!  the CTF-amplitude-weighted mean projection).
    module subroutine form_residual_plane( fpl, mean_fpl, resid, contrast )
        type(fplane_type), intent(in)    :: fpl, mean_fpl
        type(fplane_type), intent(inout) :: resid
        real, optional,    intent(in)    :: contrast
        integer :: h1, h2, k1, k2
        real    :: a
        a = 1.0; if( present(contrast) ) a = contrast
        h1 = lbound(fpl%cmplx_plane,1); h2 = ubound(fpl%cmplx_plane,1)
        k1 = lbound(fpl%cmplx_plane,2); k2 = ubound(fpl%cmplx_plane,2)
        if( .not. allocated(resid%cmplx_plane) ) allocate(resid%cmplx_plane(h1:h2,k1:k2))
        resid%cmplx_plane = fpl%cmplx_plane - a*mean_fpl%cmplx_plane
        resid%frlims = fpl%frlims
        resid%nyq    = fpl%nyq
    end subroutine form_residual_plane


    module subroutine cleanup_plane( fpl )
        type(fplane_type), intent(inout) :: fpl
        if( allocated(fpl%cmplx_plane)    ) deallocate(fpl%cmplx_plane)
        if( allocated(fpl%ctfsq_plane)    ) deallocate(fpl%ctfsq_plane)
        if( allocated(fpl%transfer_plane) ) deallocate(fpl%transfer_plane)
    end subroutine cleanup_plane

end submodule simple_flex_pca_em_mean
