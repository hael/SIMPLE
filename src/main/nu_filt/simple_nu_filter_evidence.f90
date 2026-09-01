!@descr: compact immutable evidence state for the direct NU-conditioned PCG replay
submodule (simple_nu_filter) simple_nu_filter_evidence
implicit none
#include "simple_local_flags.inc"

contains

    !> Build a phase-free compact state from the exact base-half unary bank.
    !! The explicit null candidate predicts zero from each opposite half.  Its
    !! systematic loss offset relative to the best member of the exact signal
    !! bank is robustly calibrated on the broad support.  This includes the
    !! best-of-bank advantage available under the null.  A signal candidate wins
    !! only when its cross-half predictive improvement exceeds that calibrated
    !! distribution.  The null participates in a separate ordered-label field;
    !! it never changes the established NU filtering label map.
    module subroutine build_nu_evidence_state( vol_even, vol_odd, state )
        class(image),            intent(in)  :: vol_even, vol_odd
        type(nu_evidence_state), intent(out) :: state
        type(image) :: vol_zero
        type(string) :: identity_seed_string, identity_hash
        integer(kind=NU_LABEL_KIND), allocatable :: evidence_map(:,:,:)
        real, allocatable :: null_full(:,:,:), smooth_tmp(:,:,:), null_cost(:), coords(:), signal_lps(:)
        real, allocatable :: gaps(:), probs(:), band_support_tmp(:,:)
        real, allocatable :: band_limits_active(:)
        real    :: nextb
        real(kind=8) :: fingerprint(6), cutoff_checksum, uncertainty_checksum, support_checksum
        real(kind=8) :: whitening_checksum
        real :: beta, temperature, best_e, second_e, e, prob_sum, entropy, null_lp
        real :: null_bias_median, null_bias_mad, null_bias_threshold
        integer :: n_signal, n_candidates, imask, icand, iband, i, j, k, label, n_uncertain, k25
        integer :: nb_active
        character(len=32) :: value_text
        character(len=XLONGSTRLEN) :: identity_seed

        if( .not.nu_evidence_requested ) &
            &THROW_HARD('NU evidence setup was not requested with an evidence source')
        if( trim(nu_evidence_source) /= NU_EVIDENCE_SOURCE_BASE .and. &
            &trim(nu_evidence_source) /= NU_EVIDENCE_SOURCE_PREV )then
            THROW_HARD('NU evidence source is not base_unfil or previous_shipped')
        endif
        if( nu_aux_replacement_label > 0 ) &
            &THROW_HARD('NU evidence bank contains an auxiliary replacement candidate')
        if( any(vol_even%get_ldim() /= ldim) .or. any(vol_odd%get_ldim() /= ldim) ) &
            &THROW_HARD('NU evidence half-map dimensions differ from setup')
        if( abs(vol_even%get_smpd() - smpd) > TINY .or. abs(vol_odd%get_smpd() - smpd) > TINY ) &
            &THROW_HARD('NU evidence half-map sampling differs from setup')
        if( .not.allocated(dmats_mask) .or. .not.allocated(candidate_coords) ) &
            &THROW_HARD('NU unary bank was released before evidence compaction')
        if( .not.allocated(nu_noise_profile_cached) ) &
            &THROW_HARD('NU whitening profile is unavailable for evidence compaction')
        if( n_nu_mask < 1 ) THROW_HARD('NU evidence support is empty')

        call calculate_nu_source_fingerprint(vol_even, vol_odd, fingerprint)
        if( any(abs(fingerprint - nu_evidence_source_fingerprint) > &
            &1.d-11 * max(1.d0, abs(nu_evidence_source_fingerprint))) )then
            THROW_HARD('NU evidence half maps differ from the evidence-source setup pair')
        endif

        n_signal = size(candidate_coords)
        if( n_signal < 1 .or. n_signal > NU_DMAT_CANDIDATE_CAP ) &
            &THROW_HARD('invalid NU evidence signal-candidate count')
        if( size(dmats_mask,1) /= n_nu_mask .or. size(dmats_mask,2) < n_signal ) &
            &THROW_HARD('NU unary bank shape is incompatible with evidence compaction')
        do icand = 2, n_signal
            if( candidate_coords(icand) <= candidate_coords(icand-1) ) &
                &THROW_HARD('NU evidence candidate ordering is not strictly coarse-to-fine')
        enddo

        n_candidates = n_signal + 1
        if( any(.not.ieee_is_finite(dmats_mask(:,:n_signal))) .or. any(dmats_mask(:,:n_signal) < 0.) .or. &
            &any(dmats_mask(:,:n_signal) >= NU_EVIDENCE_INVALID) ) &
            &THROW_HARD('NU evidence signal unaries must be finite and nonnegative')

        ! Candidate zero: no reproducible signal.  Smooth it at the same scale
        ! as the adjacent 20-A label so the null/coarsest comparison is not
        ! manufactured by unequal spatial averaging.  The raw zero predictor
        ! has an objective offset relative to smoothed predictors even for
        ! independent noise, and choosing the best of several signal candidates
        ! adds a multiple-comparison advantage.  Calibrate both effects from the
        ! distribution of C_zero-min(C_signal) over the deliberately generous
        ! support.  The subtracted offset is the NULL COMPONENT'S CENTER ONLY,
        ! estimated by the lower quartile of the gap distribution: the gaps are
        ! a solvent/signal mixture in which signal voxels sit strictly higher,
        ! so the lower quartile tracks the null component even on a
        ! majority-molecule support where the median is contaminated.  A
        ! center-plus-3MAD offset is a DETECTION threshold, not a likelihood
        ! offset -- on streptavidin (2026-08-27) it made the null unbeatable
        ! everywhere (null_fraction=1.0, global over-smoothing).  With the
        ! center-only offset the softmax stays a graded competition and the
        ! ordered-label spatial model does the consolidation; median and MAD
        ! are retained as recorded diagnostics.
        allocate(null_full(ldim(1),ldim(2),ldim(3)), source=0.)
        allocate(smooth_tmp(ldim(1),ldim(2),ldim(3)), source=0.)
        allocate(null_cost(n_nu_mask), source=0.)
        call vol_zero%new(ldim, smpd)
        call vol_zero%zero
        call vol_even%nu_objective(vol_zero, vol_odd, vol_zero, null_full, nu_lmask, &
            &nu_noise_profile_cached, nu_noise_rmax_cached)
        null_lp = nu_label_lowpass_limit(1)
        call smooth_nu_objective(null_full, smooth_tmp, null_lp)
        !$omp parallel do schedule(static) default(shared) private(imask,i,j,k) proc_bind(close)
        do imask = 1, n_nu_mask
            i = nu_mask_vox(1,imask)
            j = nu_mask_vox(2,imask)
            k = nu_mask_vox(3,imask)
            null_cost(imask) = null_full(i,j,k)
        enddo
        !$omp end parallel do
        call vol_zero%kill
        deallocate(null_full, smooth_tmp)
        call release_nu_smooth_norm
        if( any(.not.ieee_is_finite(null_cost)) .or. any(null_cost < 0.) .or. &
            &any(null_cost >= NU_EVIDENCE_INVALID) ) &
            &THROW_HARD('raw NU evidence null unary must be finite and nonnegative')

        allocate(gaps(n_nu_mask), source=0.)
        !$omp parallel do schedule(static) default(shared) private(imask) proc_bind(close)
        do imask = 1, n_nu_mask
            gaps(imask) = null_cost(imask) - minval(dmats_mask(imask,:n_signal))
        enddo
        !$omp end parallel do
        null_bias_median = median_nocopy(gaps)
        null_bias_mad = mad_gau(gaps, null_bias_median)
        ! lower-quartile center of the gap mixture; selec reorders but
        ! preserves the values, so it composes with the diagnostics above
        k25 = max(1, nint(0.25 * real(n_nu_mask)))
        null_bias_threshold = selec(k25, n_nu_mask, gaps)
        deallocate(gaps)
        if( .not.ieee_is_finite(null_bias_median) .or. .not.ieee_is_finite(null_bias_mad) .or. &
            &.not.ieee_is_finite(null_bias_threshold) .or. null_bias_mad < 0. ) &
            &THROW_HARD('NU evidence null-bias calibration is invalid')
        null_cost = null_cost - null_bias_threshold
        if( any(.not.ieee_is_finite(null_cost)) .or. any(abs(null_cost) >= NU_EVIDENCE_INVALID) ) &
            &THROW_HARD('calibrated NU evidence null score is invalid')

        allocate(coords(n_candidates), source=0.)
        coords(2:) = candidate_coords
        allocate(signal_lps(n_signal), source=0.)
        do icand = 1, n_signal
            signal_lps(icand) = nu_label_lowpass_limit(icand)
        enddo

        beta = estimate_evidence_beta(null_cost, dmats_mask(:,:n_signal))
        call regularize_evidence_labels(null_cost, dmats_mask(:,:n_signal), coords, beta, evidence_map)

        ! Calibrate confidence against the final spatial-model energy gap, not
        ! against raw Huber values.  This is a deterministic temperature for
        ! the exact bank, whitening profile, smoothing radii, and Potts model.
        allocate(gaps(n_nu_mask), source=0.)
        !$omp parallel do schedule(static) default(shared) &
        !$omp private(imask,i,j,k,icand,best_e,second_e,e) proc_bind(close)
        do imask = 1, n_nu_mask
            i = nu_mask_vox(1,imask)
            j = nu_mask_vox(2,imask)
            k = nu_mask_vox(3,imask)
            best_e = huge(best_e)
            second_e = huge(second_e)
            do icand = 1, n_candidates
                e = evidence_site_energy(imask, icand, null_cost, dmats_mask(:,:n_signal), coords, evidence_map, beta)
                if( e < best_e )then
                    second_e = best_e
                    best_e = e
                else if( e < second_e )then
                    second_e = e
                endif
            enddo
            gaps(imask) = max(0., second_e - best_e)
        enddo
        !$omp end parallel do
        temperature = median_nocopy(gaps)
        if( temperature <= TINY )then
            temperature = sum(gaps, mask=gaps > TINY) / real(max(1, count(gaps > TINY)))
        endif
        if( temperature <= TINY ) temperature = sqrt(epsilon(1.)) * &
            &max(1., (sum(null_cost) + sum(dmats_mask(:,:n_signal))) / real(n_nu_mask * n_candidates))
        if( .not.ieee_is_finite(temperature) .or. temperature <= TINY ) &
            &THROW_HARD('NU evidence calibration temperature is invalid')
        deallocate(gaps)

        ! Stage 6.6 band ladder from the ACTUAL bank: the static four bands,
        ! extended geometrically only while an accepted candidate is at least
        ! as fine as the next boundary. With the static bank (abinitio3D's
        ! discrete-ladder mode, nu_refine=no) this is exactly the pre-6.6
        ! four-band behavior; when the nu_refine shell walk has extended the
        ! bank (refine3D_auto, mirroring the gridding path), the ladder grows
        ! over the ACCEPTED candidates only, so no band can exist without a
        ! covering, challenger-validated probe.
        allocate(band_limits_active(NU_EVIDENCE_NBANDS), source=NU_EVIDENCE_BAND_LIMITS)
        do
            nb_active = size(band_limits_active)
            if( nb_active >= NU_EVIDENCE_MAX_NBANDS ) exit
            nextb = band_limits_active(nb_active) * NU_EVIDENCE_BAND_RATIO
            if( signal_lps(n_signal) > nextb + TINY ) exit
            band_limits_active = [band_limits_active, nextb]
        enddo
        nb_active = size(band_limits_active)
        allocate(state%selected_label(n_nu_mask), source=0_NU_LABEL_KIND)
        allocate(state%selected_cutoff(n_nu_mask), source=0.)
        allocate(state%uncertainty(n_nu_mask), source=0.)
        allocate(state%band_support(n_nu_mask,nb_active), source=0.)
        allocate(probs(n_candidates), source=0.)
        n_uncertain = 0
        do imask = 1, n_nu_mask
            i = nu_mask_vox(1,imask)
            j = nu_mask_vox(2,imask)
            k = nu_mask_vox(3,imask)
            if( nu_l_solvent_clamp )then
                ! solvent-constraint clamp: outside the envelope there is no
                ! evidenced detail by construction -- null assignment, zero
                ! band support (maximum lack-of-evidence in every band), zero
                ! uncertainty. Applied after the spherical-support statistics
                ! (calibration, whitening, spatial pass) are finalized.
                if( nu_solvent_lmask(i,j,k) )then
                    state%selected_label(imask)  = 0_NU_LABEL_KIND
                    state%selected_cutoff(imask) = 0.
                    state%uncertainty(imask)     = 0.
                    state%band_support(imask,:)  = 0.
                    cycle
                endif
            endif
            best_e = huge(best_e)
            do icand = 1, n_candidates
                probs(icand) = evidence_site_energy(imask, icand, null_cost, dmats_mask(:,:n_signal), &
                    &coords, evidence_map, beta)
                best_e = min(best_e, probs(icand))
            enddo
            do icand = 1, n_candidates
                e = (probs(icand) - best_e) / temperature
                if( e >= 80. )then
                    probs(icand) = 0.
                else
                    probs(icand) = exp(-e)
                endif
            enddo
            prob_sum = sum(probs)
            if( .not.ieee_is_finite(prob_sum) .or. prob_sum <= TINY ) &
                &THROW_HARD('NU evidence candidate probability normalization failed')
            probs = probs / prob_sum
            entropy = 0.
            do icand = 1, n_candidates
                if( probs(icand) > TINY ) entropy = entropy - probs(icand) * log(probs(icand))
            enddo
            entropy = entropy / log(real(n_candidates))
            state%uncertainty(imask) = min(1., max(0., entropy))
            if( state%uncertainty(imask) >= NU_EVIDENCE_UNCERTAIN_ENTROPY ) n_uncertain = n_uncertain + 1

            label = int(evidence_map(i,j,k)) - 1
            state%selected_label(imask) = int(label, kind=NU_LABEL_KIND)
            if( label > 0 ) state%selected_cutoff(imask) = signal_lps(label)
            do iband = 1, nb_active
                do icand = 1, n_signal
                    if( signal_lps(icand) <= band_limits_active(iband) + TINY ) &
                        &state%band_support(imask,iband) = state%band_support(imask,iband) + probs(icand+1)
                enddo
                state%band_support(imask,iband) = min(1., max(0., state%band_support(imask,iband)))
            enddo
        enddo
        deallocate(probs, evidence_map)

        state%summary%valid = .true.
        state%summary%ldim = ldim
        state%summary%smpd = smpd
        state%summary%mskdiam = nu_support_mskdiam
        state%summary%n_support = n_nu_mask
        state%summary%n_candidates = n_candidates
        ! Stage 6.6 evidence-gated retention: an appended (frontier-tracked)
        ! band is KEPT only if it actually earns support; otherwise the
        ! subdivision would replace the fine tail's partially evidenced
        ! weight with the full penalty and over-suppress genuine signal
        ! (measured on 1WCM, pcg_priors.md S6.6 run record). Pruning
        ! truncates finest-first; the static bands are never pruned, so
        ! pre-6.6 behavior is the guaranteed floor.
        do while( nb_active > NU_EVIDENCE_NBANDS )
            if( sum(state%band_support(:,nb_active)) / real(n_nu_mask) >= NU_EVIDENCE_MIN_BAND_SUPPORT ) exit
            nb_active = nb_active - 1
        enddo
        if( nb_active < size(state%band_support,2) )then
            if( nu_l_report ) write(logfhandle,'(A,I0,A)') &
                &'>>> NU ADAPTIVE EVIDENCE BANDS: pruned to ', nb_active, &
                &' band(s); appended band(s) earned no support'
            band_support_tmp = state%band_support(:,1:nb_active)
            call move_alloc(band_support_tmp, state%band_support)
        endif
        state%summary%n_bands = nb_active
        allocate(state%summary%band_limits(nb_active), source=band_limits_active(1:nb_active))
        allocate(state%summary%supported_fraction(nb_active), source=0.)
        state%summary%source = nu_evidence_source
        state%summary%null_fraction = real(count(state%selected_label == 0_NU_LABEL_KIND)) / real(n_nu_mask)
        state%summary%uncertain_fraction = real(n_uncertain) / real(n_nu_mask)
        state%summary%calibration_temperature = temperature
        state%summary%spatial_beta = beta
        state%summary%null_cost_mean = sum(null_cost) / real(n_nu_mask)
        state%summary%null_bias_median = null_bias_median
        state%summary%null_bias_mad = null_bias_mad
        state%summary%null_bias_threshold = null_bias_threshold
        do iband = 1, nb_active
            state%summary%supported_fraction(iband) = min(1., max(0., &
                &sum(state%band_support(:,iband)) / real(n_nu_mask)))
        enddo

        write(value_text,'(F10.4)') null_lp
        state%summary%provenance = 'algorithm='//NU_EVIDENCE_ALGORITHM//';source='//trim(nu_evidence_source)//&
            &';null=calibrated_zero_cross_half_prediction;null_reference=best_exact_signal_bank;'//&
            &'null_offset=lower_quartile_center;'//&
            &'null_smooth_A='//trim(adjustl(value_text))//&
            &';confidence=spatial_softmax_gap_temperature;'//&
            &'uncertainty=normalized_entropy;candidate_order=coarse_to_fine_validated;'//&
            &'smooth_awf=3;smooth_max_A=30;bands_A='
        do iband = 1, nb_active
            write(value_text,'(F10.4)') band_limits_active(iband)
            if( iband > 1 ) state%summary%provenance = trim(state%summary%provenance)//','
            state%summary%provenance = trim(state%summary%provenance)//trim(adjustl(value_text))
        enddo
        state%summary%provenance = trim(state%summary%provenance)//';candidates_A='
        do icand = 1, n_signal
            write(value_text,'(F10.4)') signal_lps(icand)
            if( icand > 1 ) state%summary%provenance = trim(state%summary%provenance)//','
            state%summary%provenance = trim(state%summary%provenance)//trim(adjustl(value_text))
        enddo
        write(value_text,'(ES14.6)') temperature
        state%summary%provenance = trim(state%summary%provenance)//';temperature='//trim(adjustl(value_text))
        write(value_text,'(ES14.6)') beta
        state%summary%provenance = trim(state%summary%provenance)//';beta='//trim(adjustl(value_text))
        ! the clamp is part of the evidence identity: a frozen state built
        ! with the solvent constraint must never be mistaken for one without
        if( nu_l_solvent_clamp )then
            state%summary%provenance = trim(state%summary%provenance)//';solvent_clamp=density_envelope'
        endif
        write(value_text,'(ES14.6)') null_bias_median
        state%summary%provenance = trim(state%summary%provenance)//';null_bias_median='//trim(adjustl(value_text))
        write(value_text,'(ES14.6)') null_bias_mad
        state%summary%provenance = trim(state%summary%provenance)//';null_bias_mad='//trim(adjustl(value_text))
        write(value_text,'(ES14.6)') null_bias_threshold
        state%summary%provenance = trim(state%summary%provenance)//';null_bias_threshold='//trim(adjustl(value_text))
        write(value_text,'(I0)') size(nu_noise_profile_cached)
        state%summary%provenance = trim(state%summary%provenance)//';whitening_shells='//trim(value_text)
        whitening_checksum = 0.d0
        do i = 1, size(nu_noise_profile_cached)
            whitening_checksum = whitening_checksum + real(i,8) * real(nu_noise_profile_cached(i),8)
        enddo
        write(value_text,'(ES22.14)') whitening_checksum
        state%summary%provenance = trim(state%summary%provenance)//';whitening_checksum='//trim(adjustl(value_text))
        write(value_text,'(ES14.6)') nu_noise_rmax_cached
        state%summary%provenance = trim(state%summary%provenance)//';whitening_rmax='//trim(adjustl(value_text))
        write(value_text,'(ES14.6)') nu_support_mskdiam
        state%summary%provenance = trim(state%summary%provenance)//';mskdiam_A='//trim(adjustl(value_text))

        cutoff_checksum = 0.d0
        uncertainty_checksum = 0.d0
        support_checksum = 0.d0
        do imask = 1, n_nu_mask
            cutoff_checksum = cutoff_checksum + real(1 + mod(imask,104729),8) * real(state%selected_cutoff(imask),8)
            uncertainty_checksum = uncertainty_checksum + &
                &real(1 + mod(imask,104723),8) * real(state%uncertainty(imask),8)
            do iband = 1, nb_active
                support_checksum = support_checksum + real(iband + mod(imask,104717),8) * &
                    &real(state%band_support(imask,iband),8)
            enddo
        enddo
        write(identity_seed,'(A,6(ES24.16,1X),3(ES24.16,1X))') trim(state%summary%provenance), &
            &nu_evidence_source_fingerprint, cutoff_checksum, uncertainty_checksum, support_checksum
        identity_seed_string = trim(identity_seed)
        identity_hash = identity_seed_string%to_fnv1a_hash64()
        state%summary%identity = identity_hash%to_char()

        call validate_compact_evidence_state(state)
        deallocate(null_cost, coords, signal_lps)
    end subroutine build_nu_evidence_state

    module subroutine calculate_nu_source_fingerprint( vol_even, vol_odd, fingerprint )
        class(image), intent(in) :: vol_even, vol_odd
        real(kind=8), intent(out) :: fingerprint(6)
        real(kind=8) :: even_val, odd_val, weight
        integer :: imask, i, j, k
        if( .not.allocated(nu_mask_vox) ) THROW_HARD('NU support is unavailable for source fingerprint')
        fingerprint = 0.d0
        do imask = 1, n_nu_mask
            i = nu_mask_vox(1,imask)
            j = nu_mask_vox(2,imask)
            k = nu_mask_vox(3,imask)
            even_val = real(vol_even%get([i,j,k]),8)
            odd_val  = real(vol_odd%get([i,j,k]),8)
            weight = real(1 + mod(imask,104729),8)
            fingerprint(1) = fingerprint(1) + even_val
            fingerprint(2) = fingerprint(2) + odd_val
            fingerprint(3) = fingerprint(3) + even_val * even_val
            fingerprint(4) = fingerprint(4) + odd_val * odd_val
            fingerprint(5) = fingerprint(5) + weight * even_val
            fingerprint(6) = fingerprint(6) + weight * odd_val
        enddo
    end subroutine calculate_nu_source_fingerprint

    module logical function nu_evidence_state_is_valid( state )
        type(nu_evidence_state), intent(in) :: state
        integer :: iband
        nu_evidence_state_is_valid = .false.
        if( .not.state%summary%valid ) return
        if( any(state%summary%ldim < 1) .or. state%summary%smpd <= TINY ) return
        if( state%summary%n_support < 1 ) return
        if( state%summary%mskdiam <= TINY ) return
        ! band count is adaptive (Stage 6.6): the static ladder is the floor,
        ! frontier-tracked extensions are bounded by the declared cap
        if( state%summary%n_candidates < 2 ) return
        if( state%summary%n_bands < NU_EVIDENCE_NBANDS .or. &
            &state%summary%n_bands > NU_EVIDENCE_MAX_NBANDS ) return
        if( .not.allocated(state%summary%band_limits) .or. &
            &.not.allocated(state%summary%supported_fraction) ) return
        if( size(state%summary%band_limits)        /= state%summary%n_bands ) return
        if( size(state%summary%supported_fraction) /= state%summary%n_bands ) return
        if( any(.not.ieee_is_finite(state%summary%band_limits)) .or. &
            &any(state%summary%band_limits <= 0.) ) return
        do iband = 2, state%summary%n_bands
            if( state%summary%band_limits(iband) >= state%summary%band_limits(iband-1) ) return
        enddo
        if( any(abs(state%summary%band_limits(1:NU_EVIDENCE_NBANDS) - NU_EVIDENCE_BAND_LIMITS) > 1.e-4) ) return
        if( trim(state%summary%source) /= NU_EVIDENCE_SOURCE_BASE .and. &
            &trim(state%summary%source) /= NU_EVIDENCE_SOURCE_PREV ) return
        if( len_trim(state%summary%identity) /= 16 ) return
        if( len_trim(state%summary%provenance) < 1 ) return
        if( .not.ieee_is_finite(state%summary%null_fraction) .or. &
            &state%summary%null_fraction < 0. .or. state%summary%null_fraction > 1. ) return
        if( .not.ieee_is_finite(state%summary%uncertain_fraction) .or. &
            &state%summary%uncertain_fraction < 0. .or. state%summary%uncertain_fraction > 1. ) return
        if( any(.not.ieee_is_finite(state%summary%supported_fraction)) .or. &
            &any(state%summary%supported_fraction < 0.) .or. any(state%summary%supported_fraction > 1.) ) return
        if( .not.ieee_is_finite(state%summary%null_bias_median) .or. &
            &.not.ieee_is_finite(state%summary%null_bias_mad) .or. &
            &.not.ieee_is_finite(state%summary%null_bias_threshold) .or. state%summary%null_bias_mad < 0. ) return
        if( .not.allocated(state%selected_label) .or. .not.allocated(state%selected_cutoff) .or. &
            &.not.allocated(state%uncertainty) .or. .not.allocated(state%band_support) ) return
        if( size(state%selected_label) /= state%summary%n_support ) return
        if( size(state%selected_cutoff) /= state%summary%n_support ) return
        if( size(state%uncertainty) /= state%summary%n_support ) return
        if( any(shape(state%band_support) /= [state%summary%n_support,state%summary%n_bands]) ) return
        if( any(state%selected_label < 0_NU_LABEL_KIND) .or. &
            &any(int(state%selected_label) >= state%summary%n_candidates) ) return
        if( any(.not.ieee_is_finite(state%selected_cutoff)) .or. any(state%selected_cutoff < 0.) ) return
        if( any(state%selected_cutoff <= TINY .and. state%selected_label /= 0_NU_LABEL_KIND) ) return
        if( any(state%selected_cutoff > TINY .and. state%selected_label == 0_NU_LABEL_KIND) ) return
        if( any(.not.ieee_is_finite(state%uncertainty)) .or. &
            &any(state%uncertainty < 0.) .or. any(state%uncertainty > 1.) ) return
        if( any(.not.ieee_is_finite(state%band_support)) .or. &
            &any(state%band_support < 0.) .or. any(state%band_support > 1.) ) return
        do iband = 2, state%summary%n_bands
            if( any(state%band_support(:,iband) > state%band_support(:,iband-1) + 1.e-6) ) return
        enddo
        nu_evidence_state_is_valid = .true.
    end function nu_evidence_state_is_valid

    module subroutine get_nu_evidence_summary( state, summary )
        type(nu_evidence_state),   intent(in)  :: state
        type(nu_evidence_summary), intent(out) :: summary
        if( .not.nu_evidence_state_is_valid(state) ) THROW_HARD('cannot summarize invalid NU evidence state')
        summary = state%summary
    end subroutine get_nu_evidence_summary

    module subroutine unpack_nu_evidence_state( state, selected_label, selected_cutoff, uncertainty, band_support )
        type(nu_evidence_state), intent(in) :: state
        integer, allocatable, optional, intent(out) :: selected_label(:)
        real,    allocatable, optional, intent(out) :: selected_cutoff(:), uncertainty(:), band_support(:,:)
        if( .not.nu_evidence_state_is_valid(state) ) THROW_HARD('cannot unpack invalid NU evidence state')
        if( present(selected_label) )then
            allocate(selected_label(size(state%selected_label)))
            selected_label = int(state%selected_label)
        endif
        if( present(selected_cutoff) )then
            allocate(selected_cutoff(size(state%selected_cutoff)), source=state%selected_cutoff)
        endif
        if( present(uncertainty) )then
            allocate(uncertainty(size(state%uncertainty)), source=state%uncertainty)
        endif
        if( present(band_support) )then
            allocate(band_support(size(state%band_support,1),size(state%band_support,2)), source=state%band_support)
        endif
    end subroutine unpack_nu_evidence_state

    !> Matching low-pass handoff with assignment support: the finest selected
    !! cutoff such that at least min_pct percent of the assigned (non-null)
    !! support voxels selected that cutoff or a finer one — the fine-end
    !! percentile of the assigned cutoffs. Robust against a single-voxel
    !! selection pinning the matching bandwidth, and against crop-Nyquist
    !! candidate aliasing, because duplicate finds pool naturally in the
    !! cumulative tail. min_pct=0 reproduces the raw finest-selected value.
    module real function nu_evidence_finest_supported_lp( state, min_pct )
        type(nu_evidence_state), intent(in) :: state
        real,                    intent(in) :: min_pct
        real, allocatable :: vals(:)
        integer :: n_assigned, imask, cnt, k
        nu_evidence_finest_supported_lp = 0.
        if( .not.nu_evidence_state_is_valid(state) ) return
        n_assigned = count(state%selected_label > 0_NU_LABEL_KIND)
        if( n_assigned == 0 ) return
        allocate(vals(n_assigned))
        cnt = 0
        do imask = 1, size(state%selected_label)
            if( state%selected_label(imask) > 0_NU_LABEL_KIND )then
                cnt = cnt + 1
                vals(cnt) = state%selected_cutoff(imask)
            endif
        enddo
        k = min(n_assigned, max(1, ceiling(real(n_assigned) * max(0., min_pct) / 100.)))
        nu_evidence_finest_supported_lp = selec(k, n_assigned, vals)
        deallocate(vals)
    end function nu_evidence_finest_supported_lp

    !> Replay-readiness contract (pcg_priors.md S6.2): a valid compact state
    !! is necessary but not sufficient to parameterize the replay precision.
    !! The spherical evidence support is deliberately generous and always
    !! contains BOTH substantial solvent and substantial molecule, so the null
    !! fraction is gated from both sides: below NU_EVIDENCE_MIN_NULL_FRAC is
    !! the zero-null failure mode (nothing reads as solvent, Q_NU regularizes
    !! nothing it should); above NU_EVIDENCE_MAX_NULL_FRAC is the
    !! saturated-null failure mode observed on streptavidin (everything reads
    !! as solvent, Q_NU degenerates into a global detail penalty). Hard error,
    !! no fallback: fallback policy is a workflow decision, never an inferred
    !! substitute prior.
    module subroutine assert_nu_evidence_replay_ready( state )
        type(nu_evidence_state), intent(in) :: state
        if( .not.nu_evidence_state_is_valid(state) ) &
            &THROW_HARD('NU evidence state is invalid; cannot parameterize the replay precision')
        if( state%summary%null_fraction < NU_EVIDENCE_MIN_NULL_FRAC .or. &
            &state%summary%null_fraction > NU_EVIDENCE_MAX_NULL_FRAC )then
            call print_nu_evidence_summary(state)
            write(logfhandle,'(A,F10.6,A,F10.6,A,F10.6,A)') '>>> NU REPLAY EVIDENCE REJECTED: null_fraction=', &
                &state%summary%null_fraction, ' outside [', NU_EVIDENCE_MIN_NULL_FRAC, ',', &
                &NU_EVIDENCE_MAX_NULL_FRAC, ']'
            THROW_HARD('NU evidence null population outside its readiness bounds (failed null calibration?)')
        endif
    end subroutine assert_nu_evidence_replay_ready

    !> Expand the frozen compact state into the per-band lack-of-evidence
    !! weight fields the PCG replay precision consumes: w_b = 1 - a_b inside
    !! the spherical evidence support, 1 outside it (no evidence there; the
    !! solver's own support grading confines the effect to the soft rim).
    !! The packed order is recreated from the frozen geometry alone -- the
    !! lexicographic traversal of the spherical support disc, exactly
    !! setup_nu_mask_voxels -- so this works after mutable NU state is
    !! released and never touches module-level filter state.
    module subroutine expand_nu_evidence_band_weights( state, band_w )
        type(nu_evidence_state), intent(in)  :: state
        real, allocatable,       intent(out) :: band_w(:,:,:,:)
        type(image) :: vol_supp
        logical, allocatable :: supp_lmask(:,:,:)
        integer :: i, j, k, iband, imask
        if( .not.nu_evidence_state_is_valid(state) ) &
            &THROW_HARD('cannot expand invalid NU evidence state')
        call vol_supp%disc(state%summary%ldim, state%summary%smpd, &
            &0.5 * state%summary%mskdiam / state%summary%smpd, supp_lmask)
        call vol_supp%kill
        if( count(supp_lmask) /= state%summary%n_support ) &
            &THROW_HARD('NU evidence support geometry cannot be recreated from the frozen state')
        allocate(band_w(state%summary%ldim(1),state%summary%ldim(2),state%summary%ldim(3),&
            &state%summary%n_bands), source=1.0)
        imask = 0
        do k = 1, state%summary%ldim(3)
            do j = 1, state%summary%ldim(2)
                do i = 1, state%summary%ldim(1)
                    if( .not.supp_lmask(i,j,k) ) cycle
                    imask = imask + 1
                    do iband = 1, state%summary%n_bands
                        band_w(i,j,k,iband) = min(1., max(0., 1. - state%band_support(imask,iband)))
                    enddo
                enddo
            enddo
        enddo
        if( imask /= state%summary%n_support ) &
            &THROW_HARD('NU evidence packed order mismatch during band-weight expansion')
        deallocate(supp_lmask)
    end subroutine expand_nu_evidence_band_weights

    module subroutine print_nu_evidence_summary( state )
        type(nu_evidence_state), intent(in) :: state
        integer :: iband
        if( .not.nu_evidence_state_is_valid(state) ) THROW_HARD('cannot print invalid NU evidence state')
        write(logfhandle,'(A)') '>>> NU REPLAY EVIDENCE'
        write(logfhandle,'(A,A)')    '    pcg_nu_evidence_source=', trim(state%summary%source)
        write(logfhandle,'(A,A)')    '    pcg_nu_evidence_identity=', trim(state%summary%identity)
        write(logfhandle,'(A,F10.6)') '    pcg_nu_null_fraction=', state%summary%null_fraction
        do iband = 1, state%summary%n_bands
            write(logfhandle,'(A,I2.2,A,F10.6)') '    pcg_nu_supported_fraction_band', iband, '=', &
                &state%summary%supported_fraction(iband)
        enddo
        if( .not. (NU_DEV_OUTPUT .and. nu_l_report) ) return
        write(logfhandle,'(A,I0)')   '    pcg_nu_candidate_count=', state%summary%n_candidates
        write(logfhandle,'(A,F10.6)') '    pcg_nu_uncertain_fraction=', state%summary%uncertain_fraction
        write(logfhandle,'(A,ES14.6)') '    pcg_nu_calibration_temperature=', &
            &state%summary%calibration_temperature
        write(logfhandle,'(A,ES14.6)') '    pcg_nu_null_bias_median=', state%summary%null_bias_median
        write(logfhandle,'(A,ES14.6)') '    pcg_nu_null_bias_mad=', state%summary%null_bias_mad
        write(logfhandle,'(A,ES14.6)') '    pcg_nu_null_bias_threshold=', state%summary%null_bias_threshold
        write(logfhandle,'(A,A)') '    pcg_nu_evidence_provenance=', trim(state%summary%provenance)
    end subroutine print_nu_evidence_summary

    subroutine validate_compact_evidence_state( state )
        type(nu_evidence_state), intent(in) :: state
        integer :: iband
        if( nu_evidence_state_is_valid(state) ) return
        write(logfhandle,'(A)') '>>> INVALID NU EVIDENCE STATE'
        write(logfhandle,'(A,L1,2X,A,A,2X,A,I0)') '    summary_valid=', state%summary%valid, &
            &'source=', trim(state%summary%source), 'identity_length=', len_trim(state%summary%identity)
        write(logfhandle,'(A,3(I0,1X),2X,A,I0,2X,A,I0)') '    ldim=', state%summary%ldim, &
            &'support=', state%summary%n_support, 'candidates=', state%summary%n_candidates
        if( allocated(state%selected_label) ) &
            &write(logfhandle,'(A,I0,A,I0)') '    label range=', minval(state%selected_label), ':', maxval(state%selected_label)
        if( allocated(state%selected_cutoff) ) &
            &write(logfhandle,'(A,ES12.4,A,ES12.4)') '    cutoff range=', minval(state%selected_cutoff), ':', &
                &maxval(state%selected_cutoff)
        if( allocated(state%uncertainty) ) &
            &write(logfhandle,'(A,ES12.4,A,ES12.4)') '    uncertainty range=', minval(state%uncertainty), ':', &
                &maxval(state%uncertainty)
        if( allocated(state%band_support) )then
            do iband = 1, size(state%band_support,2)
                write(logfhandle,'(A,I0,A,ES12.4,A,ES12.4)') '    band ', iband, ' range=', &
                    &minval(state%band_support(:,iband)), ':', maxval(state%band_support(:,iband))
            enddo
        endif
        THROW_HARD('constructed NU evidence state failed its immutable-state contract')
    end subroutine validate_compact_evidence_state

    real function estimate_evidence_beta( null_cost, signal_costs ) result( beta )
        real, intent(in) :: null_cost(:), signal_costs(:,:)
        real :: best_e, second_e, e
        integer :: imask, icand
        beta = 0.
        do imask = 1, size(signal_costs,1)
            best_e = huge(best_e)
            second_e = huge(second_e)
            do icand = 1, size(signal_costs,2) + 1
                e = evidence_unary(imask, icand, null_cost, signal_costs)
                if( e < best_e )then
                    second_e = best_e
                    best_e = e
                else if( e < second_e )then
                    second_e = e
                endif
            enddo
            beta = beta + max(0., second_e - best_e)
        enddo
        beta = NU_LABEL_SMOOTH_BETA_FRAC * beta / real(size(signal_costs,1))
    end function estimate_evidence_beta

    subroutine regularize_evidence_labels( null_cost, signal_costs, coords, beta, candmap )
        real, intent(in) :: null_cost(:), signal_costs(:,:), coords(:), beta
        integer(kind=NU_LABEL_KIND), allocatable, intent(out) :: candmap(:,:,:)
        integer :: imask, icand, i, j, k, iter, color, current, best, nchanged
        real :: e, best_e
        allocate(candmap(ldim(1),ldim(2),ldim(3)), source=1_NU_LABEL_KIND)
        !$omp parallel do schedule(static) default(shared) private(imask,icand,i,j,k,best,best_e) proc_bind(close)
        do imask = 1, n_nu_mask
            i = nu_mask_vox(1,imask)
            j = nu_mask_vox(2,imask)
            k = nu_mask_vox(3,imask)
            best = 1
            best_e = null_cost(imask)
            do icand = 2, size(signal_costs,2) + 1
                if( signal_costs(imask,icand-1) < best_e )then
                    best = icand
                    best_e = signal_costs(imask,icand-1)
                endif
            enddo
            candmap(i,j,k) = int(best,kind=NU_LABEL_KIND)
        enddo
        !$omp end parallel do
        if( beta <= TINY ) return
        do iter = 1, NU_LABEL_SMOOTH_MAXITS
            nchanged = 0
            do color = 0, NU_LABEL_SMOOTH_NCOLORS - 1
                !$omp parallel do schedule(static) default(shared) &
                !$omp private(imask,icand,i,j,k,current,best,e,best_e) reduction(+:nchanged) proc_bind(close)
                do imask = 1, n_nu_mask
                    i = nu_mask_vox(1,imask)
                    j = nu_mask_vox(2,imask)
                    k = nu_mask_vox(3,imask)
                    if( nu_label_smooth_color(i,j,k) /= color ) cycle
                    current = int(candmap(i,j,k))
                    best = current
                    best_e = evidence_site_energy(imask, current, null_cost, signal_costs, coords, candmap, beta)
                    do icand = 1, size(signal_costs,2) + 1
                        if( icand == current ) cycle
                        e = evidence_site_energy(imask, icand, null_cost, signal_costs, coords, candmap, beta)
                        if( nu_label_smooth_is_better(e, best_e) )then
                            best = icand
                            best_e = e
                        endif
                    enddo
                    if( best /= current )then
                        candmap(i,j,k) = int(best,kind=NU_LABEL_KIND)
                        nchanged = nchanged + 1
                    endif
                enddo
                !$omp end parallel do
            enddo
            if( nchanged == 0 ) exit
        enddo
    end subroutine regularize_evidence_labels

    real function evidence_site_energy( imask, icand, null_cost, signal_costs, coords, candmap, beta ) result( energy )
        integer, intent(in) :: imask, icand
        real, intent(in) :: null_cost(:), signal_costs(:,:), coords(:), beta
        integer(kind=NU_LABEL_KIND), intent(in) :: candmap(:,:,:)
        integer :: i, j, k, neigh(3,NU_LABEL_SMOOTH_NNEIGH), nsz, ineigh, ni, nj, nk, degree, other
        real :: pair_sum
        i = nu_mask_vox(1,imask)
        j = nu_mask_vox(2,imask)
        k = nu_mask_vox(3,imask)
        call neigh_8_3D(ldim, [i,j,k], neigh, nsz)
        pair_sum = 0.
        degree = 0
        do ineigh = 1, nsz
            ni = neigh(1,ineigh)
            nj = neigh(2,ineigh)
            nk = neigh(3,ineigh)
            if( .not.nu_lmask(ni,nj,nk) ) cycle
            degree = degree + 1
            other = int(candmap(ni,nj,nk))
            pair_sum = pair_sum + nu_label_smooth_coord_pair_cost(coords(icand),coords(other))
        enddo
        if( degree > 0 ) pair_sum = pair_sum / real(degree)
        energy = evidence_unary(imask, icand, null_cost, signal_costs) + beta * pair_sum
    end function evidence_site_energy

    real function evidence_unary( imask, icand, null_cost, signal_costs ) result( unary )
        integer, intent(in) :: imask, icand
        real, intent(in) :: null_cost(:), signal_costs(:,:)
        if( icand == 1 )then
            unary = null_cost(imask)
        else
            unary = signal_costs(imask,icand-1)
        endif
    end function evidence_unary

end submodule simple_nu_filter_evidence
