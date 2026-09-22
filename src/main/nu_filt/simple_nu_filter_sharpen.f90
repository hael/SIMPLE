!@descr: NU-evidence nonuniform postprocessing, v2 (classical pipeline, local)
!
! v2 of the postprocess_nu gain stack (nu_evidence_local_sharpening.md S3b/S3c).
! v1 (confidence-as-Wiener band gains + amplitude-ratio restoration) measurably
! over-sharpened the high-resolution core on PfCRT: band support confidence is
! a calibrated SUPPORT PROBABILITY, not an SSNR -- it saturates to 1 wherever
! evidence exists, so the finest bands of a raw unfiltered input passed at full
! noise power and the ratio restoration then boosted local dips toward a
! noise-level target. v2 follows the classical postprocess recipe instead --
! shrinkage and sharpening coupled, in that order -- localized by the evidence:
!
!   1. sharpen: one classical Guinier B-factor, estimated from the merged map
!      between HPLIM_GUINIER and the finest evidenced local cutoff (only when
!      that cutoff is finer than NU_SHARP_BFAC_FINEST_A, exactly like the
!      standard postprocess gate);
!   2. filter (the local Wiener surrogate): per-voxel Butterworth low-pass at
!      the frozen state's evidenced local cutoff (selected_cutoff, from the
!      Potts-smoothed label optimization), composed exactly like the
!      production NU filter -- so sharpening never extends beyond the local
!      passband, which is the classical lesson v1 skipped; and, when the
!      evidence pair's FSC is supplied, its 2FSC/(1+FSC) weighting inside
!      each local passband, stretched so that its FSC=0.143 crossing sits at
!      the voxel's cutoff (2026-09-22): the Butterworth alone let the
!      low-SNR shoulder pass at full weight under the B-factor, the same
!      failure as the isotropic postprocess on exp_gate;
!   3. solvent: voxels the calibrated null claims (cutoff 0) and voxels
!      outside the spherical evidence support flatten to the half-map mean --
!      the evidence-derived envelope behavior validated in the v1 run.
!
! The global FSC-derived optlp is NOT applied as such: the global FSC
! averages over the map and would erase the evidenced local extension in the
! core. It is applied stretched to each local cutoff instead, so the core
! keeps its extension with the global SNR falloff's shape and the poorly
! resolved regions get it compressed.
!
! Discipline unchanged from v1: one frozen evidence identity; the shipped
! product is the single sharpened MERGED volume (classical-postprocess style)
! and is display/interpretation only -- the unregularized base half pair
! keeps sole resolution authority.
submodule (simple_nu_filter) simple_nu_filter_sharpen
implicit none
#include "simple_local_flags.inc"

contains

    !> The policy in one signature: everything is ESTIMATED on the evidence
    !! pair (vol_even/odd, the unregularized base pair: the evidence state it
    !! matches, the FSC whose 2FSC/(1+FSC) weighting is applied inside every
    !! local passband stretched so that its FSC=0.143 crossing sits at the
    !! voxel's evidenced cutoff, and the Guinier B-factor) and APPLIED to the
    !! apply pair (apply_even/odd, the solvent-prior'd pair when there is one;
    !! absent, the evidence pair is sharpened). 2026-09-22: the Butterworth
    !! alone let the low-SNR shoulder pass at full weight under the B-factor,
    !! the same failure as the isotropic postprocess on exp_gate.
    module subroutine nu_evidence_sharpen_vol( state, vol_even, vol_odd, fsc, vol_sharp, apply_even, apply_odd )
        type(nu_evidence_state), intent(in)    :: state
        type(image),             intent(in)    :: vol_even, vol_odd
        real,                    intent(in)    :: fsc(:)
        type(image),             intent(inout) :: vol_sharp
        type(image), optional,   intent(in)    :: apply_even, apply_odd
        type(nu_evidence_summary) :: summ
        type(image) :: vol_merged, vol_work, vol_filt, vol_supp
        real,    allocatable :: cutoffs(:), cutoff_map(:,:,:), distinct(:), filt(:), optlp(:), res(:)
        real,    allocatable :: rtmp(:,:,:), rout(:,:,:)
        logical, allocatable :: supp_lmask(:,:,:)
        real    :: finest, bfac, mean_half, c_here, fsc05, fsc0143, kscaled, frac
        integer :: ldim(3), i, j, k, imask, nyq, ndist, icut, cutoff_find, nnull, k0143, kk, klo
        logical :: l_apply
        if( .not.nu_evidence_state_is_valid(state) ) &
            &THROW_HARD('cannot sharpen from an invalid NU evidence state; nu_evidence_sharpen_vol')
        call get_nu_evidence_summary(state, summ)
        ldim = summ%ldim
        if( any(vol_even%get_ldim() /= ldim) .or. any(vol_odd%get_ldim() /= ldim) ) &
            &THROW_HARD('half-map dimensions do not match the frozen evidence state; nu_evidence_sharpen_vol')
        nyq = vol_even%get_filtsz()
        l_apply = present(apply_even) .and. present(apply_odd)
        if( present(apply_even) .neqv. present(apply_odd) ) &
            &THROW_HARD('the apply pair needs both halves; nu_evidence_sharpen_vol')
        if( l_apply )then
            if( any(apply_even%get_ldim() /= ldim) .or. any(apply_odd%get_ldim() /= ldim) ) &
                &THROW_HARD('apply pair dimensions do not match the frozen evidence state; nu_evidence_sharpen_vol')
        endif
        if( size(fsc) < nyq ) THROW_HARD('FSC has fewer shells than the volume; nu_evidence_sharpen_vol')
        allocate(optlp(nyq), source=0.)
        where( fsc(1:nyq) > 0. ) optlp = 2. * fsc(1:nyq) / (fsc(1:nyq) + 1.)
        where( fsc(1:nyq) < 0.05 ) optlp = 0.
        optlp = min(optlp, 0.99999)
        res = get_resarr(ldim(1), summ%smpd)
        call get_resolution(fsc(1:nyq), res, fsc05, fsc0143)
        k0143 = max(1, min(nyq, calc_fourier_index(fsc0143, ldim(1), summ%smpd)))
        ! expand the mask-packed evidenced local cutoffs to the grid (same
        ! spherical-support recreation and packing order as
        ! expand_nu_evidence_band_weights); 0 = null verdict or outside support
        call unpack_nu_evidence_state(state, selected_cutoff=cutoffs)
        call vol_supp%disc(ldim, summ%smpd, 0.5 * summ%mskdiam / summ%smpd, supp_lmask)
        call vol_supp%kill
        if( count(supp_lmask) /= summ%n_support ) &
            &THROW_HARD('NU evidence support geometry cannot be recreated from the frozen state; nu_evidence_sharpen_vol')
        allocate(cutoff_map(ldim(1),ldim(2),ldim(3)), source=0.0)
        imask = 0
        do k = 1, ldim(3)
            do j = 1, ldim(2)
                do i = 1, ldim(1)
                    if( .not.supp_lmask(i,j,k) ) cycle
                    imask = imask + 1
                    cutoff_map(i,j,k) = cutoffs(imask)
                enddo
            enddo
        enddo
        ! distinct evidenced cutoffs, coarse to fine
        allocate(distinct(0))
        do imask = 1, size(cutoffs)
            if( cutoffs(imask) <= TINY ) cycle
            if( size(distinct) > 0 )then
                if( any(abs(distinct - cutoffs(imask)) < 1.e-3) ) cycle
            endif
            distinct = [distinct, cutoffs(imask)]
        enddo
        ndist = size(distinct)
        if( ndist < 1 ) THROW_HARD('NU evidence state carries no evidenced cutoffs; nu_evidence_sharpen_vol')
        do i = 1, ndist - 1
            do j = i + 1, ndist
                if( distinct(j) > distinct(i) )then
                    c_here      = distinct(i)
                    distinct(i) = distinct(j)
                    distinct(j) = c_here
                endif
            enddo
        enddo
        finest = distinct(ndist)
        nnull  = count(cutoff_map <= TINY .and. supp_lmask)
        ! classical Guinier B-factor of the evidence (unregularized) pair
        ! average, anchored to the finest evidenced cutoff, gated exactly like
        ! the standard postprocess; the merged map that is sharpened is the
        ! apply pair's when there is one
        call vol_work%copy(vol_even)
        call vol_work%add(vol_odd)
        call vol_work%mul(0.5)
        if( l_apply )then
            call vol_merged%copy(apply_even)
            call vol_merged%add(apply_odd)
            call vol_merged%mul(0.5)
        else
            call vol_merged%copy(vol_work)
        endif
        bfac = 0.0
        if( finest < NU_SHARP_BFAC_FINEST_A )then
            bfac = vol_work%guinier_bfac(HPLIM_GUINIER, finest)
            write(logfhandle,'(A,1X,F8.2)') '>>> NU SHARPENING B-FACTOR (UNREGULARIZED PAIR) DETERMINED TO:', bfac
        else
            write(logfhandle,'(A,F6.2,A)') '>>> NU SHARPENING B-FACTOR SKIPPED (finest evidenced cutoff ', &
                &finest, ' A too coarse)'
        endif
        write(logfhandle,'(A)')       '>>> NU EVIDENCE SHARPENING (postprocess_nu v2): B-sharpen, FSC-weighted (2FSC/(1+FSC) '//&
            &'stretched to each local cutoff), then local low-pass'
        if( l_apply ) write(logfhandle,'(A)') '    estimated on the unregularized pair, applied to the solvent-prior pair'
        write(logfhandle,'(A,I0,A)')  '    evidenced local cutoffs: ', ndist, ' distinct'
        write(logfhandle,'(A,I0)')    '    null-claimed voxels flattened to the mean: ', nnull
        ! sharpen the MERGED map: every v2 operation is linear, so this is
        ! identical to averaging two identically-processed halves, and the
        ! shipped product is the single merged volume, classical-postprocess
        ! style
        allocate(filt(nyq))
        call vol_work%copy(vol_merged)
        rtmp      = vol_work%get_rmat()
        mean_half = real(sum(real(rtmp,dp)) / real(product(ldim),dp))
        if( abs(bfac) > TINY ) call vol_work%apply_bfac(bfac)
        call vol_work%fft()
        ! compose the output from the per-cutoff filtered versions,
        ! production NU-filter style; default is the flat solvent mean
        allocate(rout(ldim(1),ldim(2),ldim(3)), source=mean_half)
        do icut = 1, ndist
            cutoff_find = min(nyq, calc_fourier_index(distinct(icut), ldim(1), summ%smpd))
            call butterworth_filter(cutoff_find, filt)
            ! the global FSC weighting stretched so that its FSC=0.143
            ! crossing sits at this cutoff: shell k of the local passband
            ! takes the global weight at k * k0143 / cutoff_find
            do kk = 1, nyq
                kscaled = real(kk) * real(k0143) / real(max(1, cutoff_find))
                klo     = floor(kscaled)
                frac    = kscaled - real(klo)
                if( klo >= nyq )then
                    filt(kk) = 0.
                else if( klo < 1 )then
                    filt(kk) = filt(kk) * optlp(1)
                else
                    filt(kk) = filt(kk) * ((1. - frac) * optlp(klo) + frac * optlp(klo + 1))
                endif
            enddo
            call vol_filt%copy(vol_work)
            call vol_filt%apply_filter(filt)
            call vol_filt%ifft()
            rtmp = vol_filt%get_rmat()
            where( abs(cutoff_map - distinct(icut)) < 1.e-3 ) rout = rtmp
            write(logfhandle,'(A,F8.3,A,I12)') '    cutoff_A ', distinct(icut), ' voxels ', &
                &count(abs(cutoff_map - distinct(icut)) < 1.e-3)
        enddo
        call vol_sharp%new(ldim, summ%smpd)
        call vol_sharp%set_rmat(rout, .false.)
        deallocate(rout)
        call vol_merged%kill
        write(logfhandle,'(A)') '>>> NU EVIDENCE SHARPENING: display/interpretation product; the unregularized'
        write(logfhandle,'(A)') '>>> half pair keeps sole resolution authority (never feed _nu_sharp maps to FSC)'
        ! destruct
        call vol_work%kill
        call vol_filt%kill
        deallocate(cutoffs, cutoff_map, distinct, filt, optlp, res)
        if( allocated(rtmp)       ) deallocate(rtmp)
        if( allocated(supp_lmask) ) deallocate(supp_lmask)
    end subroutine nu_evidence_sharpen_vol

end submodule simple_nu_filter_sharpen
