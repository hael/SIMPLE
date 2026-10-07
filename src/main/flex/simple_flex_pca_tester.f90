!@descr: unit tests for flex_pca numerical kernels, state inference, caches and codecs
! Posterior, hybrid former, coupled M-step, cross-FSC, probe codec, cache, derived settings,
! population placement, weights and latent deconvolution. Fixed seeds. The fast suite uses the
! compact numerical fixtures; the library suite deconvolves 20000 particles at realistic noise.
module simple_flex_pca_tester
use simple_core_module_api,               only: dp, DPI, CMPLX_ZERO, OSMPL_PAD_FAC, KBWINSZ, KBALPHA, &
    &fdim, fplane_type, string
use simple_defs_flex,                     only: FLEX_MAX_BW_GROW
use simple_syslib,                        only: del_file
use simple_parameters,                    only: parameters
use simple_sp_project,                    only: sp_project
use simple_reconstructor,                 only: reconstructor
use simple_image,                         only: image
use simple_ori,                           only: ori
use simple_sym,                           only: sym
use simple_kbinterpol,                    only: kbinterpol
use simple_math,                          only: ceil_div, floor_div
use simple_linalg,                        only: jacobi
use simple_gridding,                      only: prep3D_inv_kbenvelope4mul
use simple_rnd,                           only: seed_rnd_fixed
use simple_flex_pca_embedding_io,         only: write_embedding_cache, read_embedding_cache
use simple_flex_pca_records,              only: flex_selection, flex_fit_model, flex_latent, flex_state_set
use simple_flex_pca_state_service,        only: place_states_with_population_floor, &
    &auto_box_crop, auto_min_neff, auto_state_count, FLEX_AUTO_K_MIN, FLEX_AUTO_K_START
use simple_flex_pca_util,                 only: kernel_weights_at_bandwidth
use simple_flex_pca_weights,              only: build_covariance_state_weights
use simple_flex_pca_deconv,               only: calibrate_noise_scale, deconvolve_latent
use simple_flex_pca_fit_types,            only: flex_fit, flex_probe_part, xfsc_ctx_t, cleanup_plane
use simple_flex_pca_posterior,            only: probe_solve_plain
use simple_flex_pca_polar,                only: polar_grid_build, polar_project_recs
use simple_flex_pca_plane_cache,          only: flex_pca_plane_cache, plane_cache_contract_header, &
    &plane_cache_master_matches, plane_cache_worker_matches
use simple_flex_reconstructor_latent_ops, only: project_fplanes_mean_basis, &
    &insert_planes_oversamp_coupled_batch_scaled, solve_coupled_basis_exp, LATENT_WDIM
use simple_flex_pca_pcg,                  only: flex_pcg_t, flex_pcg_outcome_t
use simple_flex_pca_basis,                only: cross_half_subspace_angles
use simple_flex_pca_rounds,               only: flex_pca_rounds_shmem
use simple_flex_pca_crossfsc,             only: COV_XFSC_FNAME
use simple_flex_probe_fit,                only: flex_probe_fit, fit_estep_former_polar, write_probe_part, &
    &reduce_probe_parts, xfsc_paired_record, polar_ring_gram, polar_ring_selfpower
use simple_test_utils
implicit none
private
public :: run_all_flex_pca_tests, run_all_flex_pca_lib_tests

contains

    subroutine run_all_flex_pca_tests()
        write(*,'(A)') '**** running all flex_pca tests ****'
        call test_posterior_solve()
        call test_estep_former()
        call test_crossfit_fsc()
        call test_probe_codec()
        call test_embedding_cache_io()
        call test_plane_cache_contract()
        call test_auto_settings()
        call test_population_floor()
        call test_kernel_bandwidth()
        call test_state_weights()
        call test_deconvolution(4000, 0.5d0, 5.d0)
        call test_mstep_toy()
    end subroutine run_all_flex_pca_tests

    subroutine run_all_flex_pca_lib_tests()
        write(*,'(A)') '**** running all flex_pca library tests ****'
        call test_deconvolution(20000, 2.d0, 20.d0)
    end subroutine run_all_flex_pca_lib_tests

    subroutine test_estep_former()
        real(dp) :: full_rel, ring_rel, low_rel, linear_rel, cart_linear_rel, quad_rel, ring_fraction
        real(dp) :: full_tol, ring_tol, low_tol, linear_tol, ring_fraction_floor
        logical  :: full_pass, ring_pass, low_pass, cart_linear_pass, ring_power_pass
        write(*,'(A)') 'test_estep_former'
        call test_flex_estep_former(full_pass, ring_pass, low_pass, cart_linear_pass, ring_power_pass, &
            &full_rel, ring_rel, low_rel, linear_rel, cart_linear_rel, quad_rel, ring_fraction, &
            &full_tol, ring_tol, low_tol, linear_tol, ring_fraction_floor)
        call assert_true(full_pass, 'hybrid E-step statistics agree with the same-sample brute-force oracle')
        call assert_true(ring_power_pass, 'hybrid E-step fixture carries resolved mean and basis power on its rings')
        call assert_true(ring_pass, 'hybrid E-step ring statistics agree relative to the ring-only oracle')
        call assert_true(low_pass, 'hybrid exact-low-k statistics agree with the restricted direct oracle')
        call assert_true(cart_linear_pass, 'the noise-free Cartesian fixture preserves the sufficient-statistic identity')
        write(*,'(A,ES11.3,A,ES11.3)') '  same-sample relative error=',full_rel,' tolerance=',full_tol
        write(*,'(A,ES11.3,A,ES11.3)') '  minimum ring-power fraction=',ring_fraction, &
            &' floor=',ring_fraction_floor
        write(*,'(A,ES11.3,A,ES11.3)') '  ring-only relative error=',ring_rel,' tolerance=',ring_tol
        write(*,'(A,ES11.3,A,ES11.3)') '  exact-low-k relative error=',low_rel,' tolerance=',low_tol
        write(*,'(A,ES11.3,A,ES11.3)') '  Cartesian linearity relative error=',cart_linear_rel, &
            &' tolerance=',linear_tol
        write(*,'(A,ES11.3)') '  polar-linearity phase-0 observation=',linear_rel
        write(*,'(A,ES11.3)') '  polar-vs-Cartesian phase-0 observation=',quad_rel
    end subroutine test_estep_former

    subroutine test_crossfit_fsc()
        logical :: signed_pass, noise_pass
        real(dp) :: noise_mean, noise_floor
        write(*,'(A)') 'test_crossfit_fsc'
        call test_flex_crossfsc(signed_pass, noise_pass, noise_mean, noise_floor)
        call assert_true(signed_pass, 'signed-permuted half bases give identity cross-fit FSC curves')
        call assert_true(noise_pass, 'independent random half bases stay at the cross-fit FSC noise floor')
        write(*,'(A,ES11.3,A,ES11.3)') '  noise mean=',noise_mean,' derived floor=',noise_floor
    end subroutine test_crossfit_fsc

    subroutine test_probe_codec()
        logical :: pass
        write(*,'(A)') 'test_probe_codec'
        call test_probe_part_codec(pass)
        call assert_true(pass, 'probe-part codec round trip is bit exact for every payload block')
    end subroutine test_probe_codec

    !> The PPCA posterior against Gaussian elimination independent of the Cholesky implementation.
    subroutine test_posterior_solve()
        integer, parameter :: NC = 4
        class(flex_fit), allocatable :: fit
        real(dp) :: L(NC,NC), G(NC,NC), precision(NC,NC), precision_inv(NC,NC)
        real(dp) :: b(NC), zref(NC), contrast, ldA, quad, rel_z, rel_cov
        logical  :: lok
        integer  :: q
        write(*,'(A)') 'test_posterior_solve'
        allocate(fit)
        fit%model%ncomp = NC
        fit%model%sig2  = 1.d0
        allocate(fit%iter%Gth(NC,NC,1), fit%iter%bth(NC,1), fit%iter%cth(NC,1))
        allocate(fit%iter%zth(NC,1), fit%iter%Ainvth(NC,NC,1), fit%iter%prior(NC))
        L = 0.d0
        L(1,1) = 1.2d0
        L(2,1:2) = [0.3d0, 1.1d0]
        L(3,1:3) = [-0.2d0, 0.4d0, 0.9d0]
        L(4,1:4) = [0.1d0, -0.3d0, 0.2d0, 1.3d0]
        G = matmul(L, transpose(L))
        b = [0.7d0, -1.1d0, 0.4d0, 1.8d0]
        fit%iter%Gth(:,:,1) = G
        fit%iter%bth(:,1)   = b
        fit%iter%cth(:,1)   = 0.d0
        fit%iter%prior      = [0.5d0, 0.8d0, 1.3d0, 2.d0]
        contrast  = 1.4d0
        precision = contrast*contrast*G
        do q = 1, NC
            precision(q,q) = precision(q,q) + fit%iter%prior(q)
        end do
        call invert_gauss_jordan(precision, precision_inv)
        zref = matmul(precision_inv, contrast*b)
        call probe_solve_plain(fit, 1, contrast, ldA, lok, quad)
        rel_z   = maxval(abs(fit%iter%zth(:,1) - zref)) / max(1.d0, maxval(abs(zref)))
        rel_cov = maxval(abs(fit%iter%Ainvth(:,:,1) - precision_inv)) / &
            &max(1.d0, maxval(abs(precision_inv)))
        call assert_true(lok, 'posterior precision is positive definite')
        call assert_true(rel_z <= 1.d-10, 'posterior mean agrees with the closed-form Gaussian solve')
        call assert_true(rel_cov <= 1.d-10, 'posterior covariance agrees with the inverse precision')
        contrast = 0.d0
        call probe_solve_plain(fit, 1, contrast, ldA, lok, quad)
        call assert_true(maxval(abs(fit%iter%zth(:,1))) <= 1.d-14, &
            &'zero contrast returns the zero prior mean')
        call assert_true(maxval(abs(diagonal(fit%iter%Ainvth(:,:,1)) - 1.d0/fit%iter%prior)) <= 1.d-12, &
            &'zero contrast returns the prior covariance')
        call fit%kill
        deallocate(fit)
    end subroutine test_posterior_solve

    subroutine invert_gauss_jordan( A, Ainv )
        real(dp), intent(in)  :: A(:,:)
        real(dp), intent(out) :: Ainv(size(A,1),size(A,2))
        real(dp) :: aug(size(A,1),2*size(A,1)), pivot_row(2*size(A,1)), scale
        integer  :: i, j, p, n
        n = size(A,1)
        aug = 0.d0
        aug(:,1:n) = A
        do i = 1, n
            aug(i,n+i) = 1.d0
        end do
        do i = 1, n
            p = i - 1 + maxloc(abs(aug(i:n,i)), dim=1)
            if( p /= i )then
                pivot_row = aug(i,:)
                aug(i,:) = aug(p,:)
                aug(p,:) = pivot_row
            endif
            scale = aug(i,i)
            aug(i,:) = aug(i,:) / scale
            do j = 1, n
                if( j == i ) cycle
                aug(j,:) = aug(j,:) - aug(j,i)*aug(i,:)
            end do
        end do
        Ainv = aug(:,n+1:2*n)
    end subroutine invert_gauss_jordan

    pure function diagonal( A ) result( d )
        real(dp), intent(in) :: A(:,:)
        real(dp) :: d(min(size(A,1),size(A,2)))
        integer :: i
        do i = 1, size(d)
            d(i) = A(i,i)
        end do
    end function diagonal

    !> the resume path: every payload survives a write/read round trip bit for bit
    subroutine test_embedding_cache_io()
        integer,          parameter :: NP = 37, NC = 4
        character(len=*), parameter :: FN = 'test_flex_pca_cache.bin'
        integer  :: pinds(NP), i, q, r, ncomp_rd
        real(dp) :: z(NP,NC), eigvals(NC), contrast(NP), re(NP), rme(NP)
        real(dp) :: prec(NC,NC,NP), sig2, sig2_rd
        real(dp), allocatable :: z_rd(:,:), eig_rd(:), con_rd(:), re_rd(:), rme_rd(:), prec_rd(:,:,:)
        type(flex_selection) :: sel
        type(flex_fit_model) :: model, model_rd
        type(flex_latent)    :: latent, latent_rd
        write(*,'(A)') 'test_embedding_cache_io'
        do i = 1, NP
            pinds(i)    = 3*i + 1                   ! non-contiguous, as a real selection is
            contrast(i) = 0.5d0 + 0.01d0*real(i,dp)
            re(i)       = real(i,dp)
            rme(i)      = 2.d0*real(i,dp)
            do q = 1, NC
                z(i,q) = sin(real(i*q,dp))
                do r = 1, NC
                    prec(q,r,i) = 1.d0/real(q+r+i,dp)
                end do
            end do
        end do
        do q = 1, NC
            eigvals(q) = 10.d0/real(q,dp)
        end do
        sig2 = 0.137d0
        sel%pinds = pinds; sel%nptcls = NP
        model%ncomp = NC; model%eigvals = eigvals; model%sig2_eff = sig2
        latent%z = z; latent%contrast = contrast; latent%resid_energy = re; latent%resid_mean_energy = rme; latent%precision = prec
        call write_embedding_cache(FN, 0, 0., sel, model, latent)
        call read_embedding_cache(FN, 0, 0., sel, model_rd, latent_rd)
        ncomp_rd = model_rd%ncomp; sig2_rd = model_rd%sig2_eff
        call move_alloc(latent_rd%z, z_rd); call move_alloc(model_rd%eigvals, eig_rd); call move_alloc(latent_rd%contrast, con_rd)
        call move_alloc(latent_rd%resid_energy, re_rd); call move_alloc(latent_rd%resid_mean_energy, rme_rd); call move_alloc(latent_rd%precision, prec_rd)
        call assert_int(NC, ncomp_rd, 'cache round trip keeps the component count')
        call assert_true(all(z == z_rd),             'cache round trip: latents bit-exact')
        call assert_true(all(eigvals == eig_rd),     'cache round trip: eigenvalues bit-exact')
        call assert_true(all(contrast == con_rd),    'cache round trip: contrast bit-exact')
        call assert_true(all(re == re_rd),           'cache round trip: residual energy bit-exact')
        call assert_true(all(rme == rme_rd),         'cache round trip: mean residual energy bit-exact')
        call assert_true(all(prec == prec_rd),       'cache round trip: precision bit-exact')
        call assert_true(sig2 == sig2_rd,            'cache round trip: sig2_eff bit-exact')
        call del_file(FN)
    end subroutine test_embedding_cache_io

    !> The persistent plane cache is exact for a master selection and coverage-based for workers.
    subroutine test_plane_cache_contract()
        integer, parameter :: PBASE(4) = [2, 4, 6, 8], PCHANGED(4) = [2, 4, 7, 8]
        integer, parameter :: PWORKER_OK(2) = [2, 8], PWORKER_BAD(2) = [2, 9]
        integer(8), allocatable :: base(:), changed_selection(:), changed_project(:), changed_stamp(:)
        class(flex_pca_plane_cache), allocatable :: first, second
        type(string) :: first_name, second_name
        write(*,'(A)') 'test_plane_cache_contract'
        base = plane_cache_contract_header('/data/a.simple', 101_8, 12, 24, 12, 1.5, PBASE)
        changed_selection = plane_cache_contract_header('/data/a.simple', 101_8, 12, 24, 12, 1.5, PCHANGED)
        changed_project = plane_cache_contract_header('/data/b.simple', 101_8, 12, 24, 12, 1.5, PBASE)
        changed_stamp = plane_cache_contract_header('/data/a.simple', 102_8, 12, 24, 12, 1.5, PBASE)
        allocate(first, second)
        call assert_true(plane_cache_master_matches(base, base), &
            &'plane cache master accepts its exact project, geometry and selection')
        call assert_true(.not. plane_cache_master_matches(base, changed_selection), &
            &'plane cache master rebuilds after a changed selection')
        call assert_true(.not. plane_cache_master_matches(base, changed_project), &
            &'plane cache master rebuilds after a changed project path')
        call assert_true(.not. plane_cache_master_matches(base, changed_stamp), &
            &'plane cache master rebuilds after an in-place project replacement')
        call assert_true(plane_cache_worker_matches(base, base, PWORKER_OK), &
            &'plane cache worker adopts a partition covered by the completed cache')
        call assert_true(.not. plane_cache_worker_matches(base, base, PWORKER_BAD), &
            &'plane cache worker refuses a partition beyond the completed cache')
        first_name = 'first-cache.bin'; second_name = 'second-cache.bin'
        call first%new(first_name); call second%new(second_name)
        call first%kill
        call assert_true(first%pristine() .and. .not. second%pristine(), &
            &'plane cache kill resets one session without touching the next')
        call second%kill
        call second%new(second_name); call first%new(first_name)
        call second%kill
        call assert_true(second%pristine() .and. .not. first%pristine(), &
            &'plane cache sessions are pristine in the opposite lifecycle order')
        call first%kill
        call first_name%kill; call second_name%kill
        deallocate(first, second)
        deallocate(base, changed_selection, changed_project, changed_stamp)
    end subroutine test_plane_cache_contract

    !> derived settings reproduce the validated runs and move in the right direction
    subroutine test_auto_settings()
        write(*,'(A)') 'test_auto_settings'
        ! box_crop: IgG-RL and Ribosembly are box 128 at 3.0 A/px run at lp=15 with box_crop=64
        call assert_int(64, auto_box_crop(128, 3.0, 15.0), 'auto box_crop reproduces the validated 64')
        call assert_true(auto_box_crop(128, 3.0,  8.0) > 64,   'auto box_crop grows with finer lp')
        call assert_true(auto_box_crop(128, 3.0, 30.0) < 64,   'auto box_crop shrinks with coarser lp')
        call assert_true(auto_box_crop(128, 3.0, 15.0) <= 128, 'auto box_crop never exceeds the native box')
        ! min_neff: Ribosembly is occupancy-limited, IgG is SNR-limited
        call assert_true(abs(auto_min_neff(335240, 16, 0.d0) - 2095) <= 50, 'auto min_neff reproduces the Ribosembly scale')
        call assert_true(auto_min_neff(100000, 20, 0.0178d0) >= 56, 'auto min_neff meets the IgG SNR requirement')
        call assert_true(auto_min_neff(10000, 20, 1.d-3) > auto_min_neff(10000, 20, 1.d-1), &
            &'auto min_neff grows as conformational SNR falls')
        ! state count: over-provision, never below the floor nor above the validated level
        call assert_int(32, auto_state_count(335240, 2000), 'auto state count over-provisions to the validated level')
        call assert_int(FLEX_AUTO_K_MIN,   auto_state_count(1000, 2000),     'auto state count keeps the small-dataset floor')
        call assert_int(FLEX_AUTO_K_START, auto_state_count(100000000, 100), 'auto state count keeps the over-provision cap')
    end subroutine test_auto_settings

    !> two clusters plus a far outlier group: every particle gets a state, every state the floor
    subroutine test_population_floor()
        integer, parameter :: NA = 300, NB = 250, NO = 20, NP = NA+NB+NO, NC = 2, NST = 2
        real,    parameter :: FRAC = 0.2
        integer  :: i, q, state, nmin, nlab(NST)
        real(dp) :: z(NP,NC), eigvals(NC), prec(NC,NC,NP)
        real,     allocatable :: weights(:,:), targets(:,:), bandwidths(:), neff(:)
        integer,  allocatable :: labels(:)
        logical  :: indicators, neff_match, floor_met
        type(flex_latent)    :: latent
        type(flex_fit_model) :: model
        type(flex_state_set) :: st
        write(*,'(A)') 'test_population_floor'
        call set_fixed_seed(20260923)
        do i = 1, NA
            z(i,1) = -3.d0 + 0.01d0*real(mod(i,7),dp)
            z(i,2) =  0.02d0*real(mod(i,5),dp)
        end do
        do i = 1, NB
            z(NA+i,1) = 3.d0 + 0.01d0*real(mod(i,7),dp)
            z(NA+i,2) = 0.02d0*real(mod(i,5),dp)
        end do
        do i = 1, NO
            z(NA+NB+i,1) = 60.d0 + 0.05d0*real(mod(i,3),dp)
            z(NA+NB+i,2) = 60.d0 + 0.05d0*real(mod(i,4),dp)
        end do
        eigvals = 1.d0
        prec    = 0.d0
        do i = 1, NP
            do q = 1, NC
                prec(q,q,i) = 1.d0
            end do
        end do
        latent%z = z; latent%precision = prec; model%eigvals = eigvals; st%nstates = NST
        call place_states_with_population_floor(latent, model, NC, 0, 10, FRAC, st)
        call move_alloc(st%weights, weights); call move_alloc(st%targets, targets); call move_alloc(st%bandwidths, bandwidths)
        call move_alloc(st%neff, neff); call move_alloc(st%labels, labels)
        nmin = max(1, nint(FRAC*real(NP)))
        call assert_int(NP, size(labels), 'one label per particle')
        call assert_true(all(labels >= 1) .and. all(labels <= NST), 'every particle has a delivered state')
        call assert_true(size(weights,1) == NP .and. size(weights,2) == NST, 'weights are particles x states')
        indicators = .true.
        do i = 1, NP
            if( abs(sum(weights(i,:)) - 1.) > 1.e-6 ) indicators = .false.
        end do
        call assert_true(indicators, 'delivered weights are hard-label indicators')
        nlab = 0
        do i = 1, NP
            nlab(labels(i)) = nlab(labels(i)) + 1
        end do
        floor_met  = .true.
        neff_match = .true.
        do state = 1, NST
            if( nlab(state) < nmin ) floor_met = .false.
            if( nint(neff(state)) /= nlab(state) ) neff_match = .false.
        end do
        call assert_true(floor_met,  'every delivered state is at or above the population floor')
        call assert_true(neff_match, 'neff reports the delivered population')
    end subroutine test_population_floor

    !> compact support, unit peak, neff bounds and bounded widening towards min_neff
    subroutine test_kernel_bandwidth()
        integer, parameter :: NP = 400
        integer  :: i, nsupp
        real(dp) :: dist(NP), h_out, h_wide
        real     :: w(NP), neff, neff_wide
        write(*,'(A)') 'test_kernel_bandwidth'
        do i = 1, NP
            dist(i) = real(i,dp)                    ! squared distances 1..NP
        end do
        ! h^2 = 100: exactly 99 particles strictly inside the support (dist < 100)
        call kernel_weights_at_bandwidth(dist, NP, 10.d0, 1, w, h_out, neff)
        nsupp = count(w > 0.)
        call assert_int(99, nsupp, 'support count for h^2 = 100')
        call assert_true(abs(h_out - 10.d0) <= 1.d-12, 'bandwidth does not grow when the support suffices')
        call assert_true(all(w(100:) == 0.), 'the kernel is compactly supported')
        call assert_real(1., maxval(w), 1.e-6, 'the kernel peak is normalised to 1')
        call assert_true(neff <= real(nsupp) .and. neff >= 1., 'neff lies between 1 and the raw support count')
        call kernel_weights_at_bandwidth(dist, NP, 10.d0, 200, w, h_wide, neff_wide)
        call assert_true(h_wide > 10.d0, 'the bandwidth widens when the support falls short of min_neff')
        call assert_true(h_wide <= 10.d0*1.3d0**FLEX_MAX_BW_GROW + 1.d-9, 'the widening stays within its cap')
        call assert_true(count(w > 0.) > nsupp, 'widening increases the support')
        call assert_true(neff_wide > neff, 'neff increases with the bandwidth')
    end subroutine test_kernel_bandwidth

    !> covariance state weights on a clearly bimodal embedding
    subroutine test_state_weights()
        integer, parameter :: NPC = 150, NP = 2*NPC, NC = 2, NST = 2
        integer  :: i, q, state, nlab(NST)
        real(dp) :: z(NP,NC), eigvals(NC), prec(NC,NC,NP)
        real,     allocatable :: weights(:,:), targets(:,:), bandwidths(:), neff(:)
        integer,  allocatable :: labels(:)
        logical  :: in_hull
        type(flex_latent)    :: latent
        type(flex_state_set) :: st
        write(*,'(A)') 'test_state_weights'
        do i = 1, NPC
            z(i,     1) = -5.d0 + 0.01d0*real(mod(i,7),dp)
            z(i,     2) =  0.02d0*real(mod(i,5),dp)
            z(NPC+i, 1) =  5.d0 + 0.01d0*real(mod(i,7),dp)
            z(NPC+i, 2) =  0.02d0*real(mod(i,5),dp)
        end do
        eigvals = 1.d0
        prec    = 0.d0
        do i = 1, NP
            do q = 1, NC
                prec(q,q,i) = 1.d0                  ! identity posterior precision
            end do
        end do
        latent%z = z; latent%precision = prec; st%nstates = NST
        call build_covariance_state_weights(latent, NC, 0, 10, st)
        call move_alloc(st%weights, weights); call move_alloc(st%targets, targets); call move_alloc(st%bandwidths, bandwidths)
        call move_alloc(st%neff, neff); call move_alloc(st%labels, labels)
        call assert_true(size(weights,1) == NP .and. size(weights,2) == NST, 'weights are particles x states')
        call assert_int(NP, size(labels), 'one label per particle')
        call assert_true(all(labels >= 0) .and. all(labels <= NST), 'labels lie in 0..nstates')
        call assert_true(all(weights >= 0.), 'kernel weights are non-negative')
        in_hull = .true.
        do q = 1, NC
            do state = 1, NST
                if( real(targets(q,state),dp) < minval(z(:,q)) - 1.d-6 .or. &
                    &real(targets(q,state),dp) > maxval(z(:,q)) + 1.d-6 ) in_hull = .false.
            end do
        end do
        call assert_true(in_hull, 'state targets lie inside the occupied latent range')
        nlab = 0
        do i = 1, NP
            if( labels(i) >= 1 ) nlab(labels(i)) = nlab(labels(i)) + 1
        end do
        call assert_true(all(nlab >= 1), 'every state draws particles from a bimodal embedding')
        call assert_true(all(neff >= 1.), 'every state has a positive effective count')
        call assert_true(abs(real(targets(1,1),dp) - real(targets(1,2),dp)) >= 5.d0, &
            &'the state targets separate along the bimodal component')
    end subroutine test_state_weights

    !> a two-component truth observed with heteroscedastic noise through a weak prior:
    !! the noise scale calibrates to 1, the held-out rule picks K = 2 and the posterior
    !! means are closer to the truth than the raw latents
    subroutine test_deconvolution( n, smin, smax )
        integer,  intent(in) :: n            !< particles; the K ladder runs to min(4, n/2000)
        real(dp), intent(in) :: smin, smax   !< per-axis noise variance drawn uniformly in [smin,smax]
        integer, parameter :: D = 4
        real(dp), allocatable :: x(:,:), z(:,:), prec(:,:,:), zhalf(:,:,:), z0(:,:)
        real(dp) :: prior(D), a, a_comp(D), s, g(D), mse_z, mse_x, mu_true(D,2), u
        integer  :: i, q, k_out, kt
        character(len=24) :: tag
        write(tag,'(A,I0,A)') 'n=', n, ': '
        write(*,'(A,I0)') 'test_deconvolution n=', n
        allocate(x(n,D), z(n,D), prec(D,D,n), zhalf(n,D,2), z0(n,D))
        mu_true = 0.d0
        mu_true(1,1) = -1.5d0; mu_true(1,2) =  1.5d0
        mu_true(2,1) =  0.5d0; mu_true(2,2) = -0.5d0
        prior = 1.d0/4.d0                           ! prior variance 4 per axis (weak)
        call set_fixed_seed(20260924)
        do i = 1, n
            kt = merge(1, 2, mod(i,3) == 0)         ! weights 1/3, 2/3
            do q = 1, D
                x(i,q) = mu_true(q,kt) + 0.3d0*gauss()
            end do
            call random_number(u)
            s = smin + (smax - smin)*u              ! noise variance smin..smax per axis
            prec(:,:,i) = 0.d0
            do q = 1, D
                prec(q,q,i) = 1.d0/s + prior(q)     ! posterior precision = data + prior
                g(q)        = gauss()*sqrt(s)
            end do
            do q = 1, D
                z(i,q) = ((x(i,q) + g(q))/s)/prec(q,q,i)
                ! halves: each with half the data precision and independent noise
                zhalf(i,q,1) = ((x(i,q) + gauss()*sqrt(2.d0*s))/(2.d0*s))/(1.d0/(2.d0*s) + prior(q))
                zhalf(i,q,2) = ((x(i,q) + gauss()*sqrt(2.d0*s))/(2.d0*s))/(1.d0/(2.d0*s) + prior(q))
            end do
        end do
        z0 = z
        call calibrate_noise_scale(zhalf, prec, prior, n, D, a, a_comp)
        call assert_true(abs(a - 1.d0) <= 0.15d0, trim(tag)//' the noise scale calibrates to 1 (within 0.15)')
        call deconvolve_latent(z, prec, prior, n, D, a, 4, k_out)
        call assert_int(2, k_out, trim(tag)//' the held-out rule picks K = 2')
        mse_z = sum((z0 - x)**2)/real(n*D,dp)
        mse_x = sum((z  - x)**2)/real(n*D,dp)
        call assert_true(mse_x <= 0.6d0*mse_z, trim(tag)//' posterior means are closer to the truth than the raw latents')
        write(*,'(A,F8.4,A,F8.4,A,F6.3)') '  mse(z)=', mse_z, '  mse(xhat)=', mse_x, '  a=', a

    contains

        real(dp) function gauss()
            real(dp) :: u1, u2
            call random_number(u1); call random_number(u2)
            gauss = sqrt(-2.d0*log(max(u1, 1.d-12)))*cos(2.d0*DPI*u2)
        end function gauss

    end subroutine test_deconvolution

    !> The production hybrid former against two distinct references.  The strict oracle sums
    !! the same polar sample positions directly, without the ring Gram/GEMV factorisation.  The
    !! Cartesian comparison is a separately labelled phase-0 observation of quadrature error.
    !! Both references intentionally gather from expanded cmat_exp with the same normalised KB
    !! stencil as production: they guard ring algebra, indexing, half-plane and Friedel handling,
    !! but do not independently validate expansion or the interpolation model.
    subroutine test_flex_estep_former( full_pass, ring_pass, low_pass, cart_linear_pass, ring_power_pass, &
            &full_rel, ring_rel, low_rel, linear_rel, cart_linear_rel, quad_rel, ring_fraction, &
            &full_tol, ring_tol, low_tol, linear_tol, ring_fraction_floor )
        logical,  intent(out) :: full_pass, ring_pass, low_pass, cart_linear_pass, ring_power_pass
        real(dp), intent(out) :: full_rel, ring_rel, low_rel, linear_rel, cart_linear_rel, quad_rel, ring_fraction
        real(dp), intent(out) :: full_tol, ring_tol, low_tol, linear_tol, ring_fraction_floor
        integer,  parameter :: BOX = 24, NC = 2, BAND = 6, RHYB = 4, NTHR = 1
        real(dp), parameter :: COEFF(NC)               = [0.65d0, -0.40d0]
        real(dp), parameter :: MIN_RING_POWER_FRACTION = 0.10_dp
        type(flex_probe_fit), allocatable :: fit
        class(parameters), allocatable, target :: params
        class(sp_project), allocatable :: project
        type(fplane_type) :: fpl
        type(ori) :: o
        real, allocatable :: r0(:,:,:), rq(:,:,:,:)
        complex, allocatable :: bank(:,:)
        real(dp) :: G(NC,NC), b(NC), c(NC), e_mm, myv, a
        real(dp) :: Gref(NC,NC), bref(NC), cref(NC), eref, myref
        real(dp) :: Ghi(NC,NC), bhi(NC), chi(NC), ehi, myhi
        real(dp) :: Gsame(NC,NC), bsame(NC), csame(NC), esame, mysame
        real(dp) :: Gring(NC,NC), bring(NC), cring(NC), ering, myring
        real(dp) :: Glow(NC,NC), blow(NC), clow(NC), elow, mylow
        real(dp) :: Glref(NC,NC), blref(NC), clref(NC), elref, mylref
        real :: ctr, xyz(3), rot(3,3)
        integer :: i, j, k, q, r, ir, nlow, nops, lim, pf
        allocate(fit, params, project, r0(BOX,BOX,BOX), rq(BOX,BOX,BOX,NC))
        allocate(fit%estep)
        params%box_crop  = BOX
        params%smpd_crop = 1.0
        params%oritype   = 'cls3D'
        ctr = 0.5*real(BOX + 1)
        do k = 1, BOX
            xyz(3) = real(k) - ctr
            do j = 1, BOX
                xyz(2) = real(j) - ctr
                do i = 1, BOX
                    xyz(1) = real(i) - ctr
                    r0(i,j,k) = blob(xyz, [2.2,-1.6,0.7], 0.9) + &
                        &0.70*blob(xyz, [-2.4,1.2,-1.1], 1.1) + 0.45*blob(xyz, [0.1,2.5,1.4], 0.8)
                    rq(i,j,k,1) = blob(xyz, [1.4,2.0,-0.8], 0.8) - &
                        &0.80*blob(xyz, [-1.8,-1.3,1.0], 1.0)
                    rq(i,j,k,2) = blob(xyz, [-2.0,1.6,0.4], 0.9) + &
                        &0.65*blob(xyz, [2.3,-1.1,-1.2], 0.8) - 0.40*blob(xyz, [0.2,2.5,1.3], 1.1)
                end do
            end do
        end do
        fit%model%ncomp = NC
        allocate(fit%model%basis_recs(NC))
        call prepare_reconstructor(fit%model%mean_rec, r0)
        do q = 1, NC
            call prepare_reconstructor(fit%model%basis_recs(q), rq(:,:,:,q))
        end do
        call o%new(.false.)
        call o%set_euler([0.0, 0.0, 0.0])

        ! A no-CTF padded plane carrying exactly mean + U*COEFF.  Only the native-grid
        ! multiples of OSMPL_PAD_FAC are consumed by either path.
        pf  = OSMPL_PAD_FAC
        lim = pf*(BOX/2)
        allocate(fpl%cmplx_plane(-lim:lim,-lim:lim), source=CMPLX_ZERO)
        allocate(fpl%transfer_plane(-lim:lim,-lim:lim), source=cmplx(1.0,0.0))
        allocate(fpl%ctfsq_plane(-lim:lim,-lim:lim), source=1.0)
        fpl%frlims = 0
        fpl%frlims(1,:) = [-lim,lim]
        fpl%frlims(2,:) = [-lim,lim]
        fpl%nyq = pf*BAND
        allocate(fit%iter%mean_fpl(NTHR), fit%iter%basis_fpls(NC,NTHR))
        call project_fplanes_mean_basis(fit%model%mean_rec, fit%model%basis_recs, o, fpl, &
            &fit%iter%mean_fpl(1), fit%iter%basis_fpls(:,1), apply_ctf_amp=.true.)
        fpl%cmplx_plane = fit%iter%mean_fpl(1)%cmplx_plane
        do q = 1, NC
            fpl%cmplx_plane = fpl%cmplx_plane + real(COEFF(q))*fit%iter%basis_fpls(q,1)%cmplx_plane
        end do

        ! One exact bank direction, identity relative in-plane rotation.
        fit%estep%l_pol_hyb = .true.
        fit%estep%rhyb_es   = RHYB
        fit%estep%nyqb_es   = BAND
        fit%estep%nyqr_es   = BAND
        fit%estep%ph0_es    = lbound(fpl%cmplx_plane,1)
        fit%estep%pk0_es    = lbound(fpl%cmplx_plane,2)
        fit%estep%hlo_es    = ceil_div (lbound(fpl%cmplx_plane,1), pf)
        fit%estep%hhi_es    = floor_div(ubound(fpl%cmplx_plane,1), pf)
        fit%estep%klo_es    = ceil_div (lbound(fpl%cmplx_plane,2), pf)
        call polar_grid_build(fit%estep%pg_es, RHYB+1, BAND, &
            &fit%estep%hlo_es, fit%estep%hhi_es, fit%estep%klo_es, &
            &fit%estep%ph0_es, fit%estep%pk0_es, gate_lo=RHYB*(RHYB+1))
        fit%estep%nsamp_es  = fit%estep%pg_es%nsamp
        fit%estep%nsamp2_es = 2*fit%estep%nsamp_es
        fit%estep%nk_es     = fit%estep%pg_es%nk
        allocate(fit%estep%dir_es(1), fit%estep%cae(1), fit%estep%sae(1))
        fit%estep%dir_es = 1; fit%estep%cae = 1.0; fit%estep%sae = 0.0
        call build_exact_positions(fit%estep%hex_es, fit%estep%kex_es, fit%estep%npos_es)
        allocate(bank(fit%estep%nsamp_es,0:NC))
        allocate(fit%estep%UsallE(fit%estep%nsamp2_es,0:NC,1))
        allocate(fit%estep%CfE(NC*NC,fit%estep%nk_es,1))
        allocate(fit%estep%Cm0E(NC,fit%estep%nk_es,1), fit%estep%c00E(fit%estep%nk_es,1))
        allocate(fit%estep%CspE(0:NC,0:NC,NTHR))
        allocate(fit%estep%xws_es(fit%estep%nsamp2_es,NTHR))
        allocate(fit%estep%wr_es(fit%estep%nk_es,NTHR), fit%estep%wrd_es(fit%estep%nk_es,NTHR))
        allocate(fit%estep%Reb_es(0:NC,NTHR))
        rot = o%get_mat()
        call polar_project_recs(fit%model%mean_rec, fit%model%basis_recs, NC, rot, fit%estep%pg_es, bank)
        do q = 0, NC
            do r = 1, fit%estep%nsamp_es
                fit%estep%UsallE(2*r-1,q,1) = fit%estep%pg_es%sqwq(r)*real(bank(r,q))
                fit%estep%UsallE(2*r,  q,1) = fit%estep%pg_es%sqwq(r)*aimag(bank(r,q))
            end do
        end do
        do ir = 1, fit%estep%nk_es
            call polar_ring_gram(fit%estep%UsallE(1,0,1), fit%estep%nsamp2_es, NC, &
                &fit%estep%pg_es%rbeg(ir), fit%estep%pg_es%rend(ir)-fit%estep%pg_es%rbeg(ir)+1, &
                &fit%estep%CspE(0,0,1), fit%estep%CfE(1,ir,1), fit%estep%Cm0E(1,ir,1))
            fit%estep%c00E(ir,1) = polar_ring_selfpower(fit%estep%UsallE(1,0,1), &
                &fit%estep%nsamp2_es, fit%estep%pg_es%rbeg(ir), &
                &fit%estep%pg_es%rend(ir)-fit%estep%pg_es%rbeg(ir)+1)
        end do
        allocate(fit%iter%Gth(NC,NC,NTHR), fit%iter%bth(NC,NTHR), fit%iter%cth(NC,NTHR))
        allocate(fit%diag%sec_ring_thr(NTHR), fit%diag%sec_exact_thr(NTHR), &
            &fit%diag%sec_solve_thr(NTHR), source=0.d0)
        call fit_estep_former_polar(fit, fit%model%mean_rec, o, fpl, 1, 1, a, e_mm, myv)
        G = fit%iter%Gth(:,:,1); b = fit%iter%bth(:,1); c = fit%iter%cth(:,1)
        call direct_polar_stats(fit%estep%UsallE(:,:,1), fit%estep%xws_es(:,1), &
            &fit%estep%wrd_es(:,1), Ghi, bhi, chi, ehi, myhi)
        call direct_cartesian_stats(fit%model%mean_rec, fit%model%basis_recs, o, fpl, RHYB, &
            &Glref, blref, clref, elref, mylref)
        Gsame = Ghi + Glref; bsame = bhi + blref; csame = chi + clref
        esame = ehi + elref; mysame = myhi + mylref
        ring_fraction = abs(ehi)/max(abs(esame),tiny(1.0_dp))
        do q = 1, NC
            ring_fraction = min(ring_fraction, &
                &abs(Ghi(q,q))/max(abs(Gsame(q,q)),tiny(1.0_dp)))
        end do
        call direct_cartesian_stats(fit%model%mean_rec, fit%model%basis_recs, o, fpl, BAND, &
            &Gref, bref, cref, eref, myref)

        full_rel = stats_relerr(G, b, c, e_mm, myv, Gsame, bsame, csame, esame, mysame)
        Gring = G - Glref; bring = b - blref; cring = c - clref
        ering = e_mm - elref; myring = myv - mylref
        ring_rel = stats_relerr(Gring, bring, cring, ering, myring, Ghi, bhi, chi, ehi, myhi)
        Glow = G - Ghi; blow = b - bhi; clow = c - chi
        elow = e_mm - ehi; mylow = myv - myhi
        low_rel  = stats_relerr(Glow, blow, clow, elow, mylow, Glref, blref, clref, elref, mylref)
        quad_rel = stats_relerr(G, b, c, e_mm, myv, Gref, bref, cref, eref, myref)
        linear_rel      = linearity_relerr(G, b, c, e_mm, myv)
        cart_linear_rel = linearity_relerr(Gref, bref, cref, eref, myref)
        nlow  = fit%estep%npos_es
        nops  = nlow + 2*fit%estep%nsamp_es
        ! Dot products accumulate single-precision projected samples.  gamma_n <= n*eps to
        ! first order; the constants cover the two-stage ring sum/GEMV and the final low-k add.
        full_tol   = 32.d0*real(epsilon(1.0),dp)*real(max(1,nops),dp)
        ! The isolated pieces each include one subtraction from the full former result.
        ring_tol   = 2.d0*full_tol
        low_tol    = 2.d0*full_tol
        linear_tol = 64.d0*real(epsilon(1.0),dp)*real(max(1,nops)*(NC+1),dp)
        ! This is a fixture-design guard: every tested diagonal must carry material ring power.
        ring_fraction_floor = MIN_RING_POWER_FRACTION
        full_pass   = full_rel   <= full_tol
        ring_power_pass = ring_fraction >= ring_fraction_floor
        ring_pass   = ring_rel   <= ring_tol
        low_pass    = low_rel    <= low_tol
        cart_linear_pass = cart_linear_rel <= linear_tol

        deallocate(bank, r0, rq)
        if( allocated(fpl%cmplx_plane) )    deallocate(fpl%cmplx_plane)
        if( allocated(fpl%transfer_plane) ) deallocate(fpl%transfer_plane)
        if( allocated(fpl%ctfsq_plane) )    deallocate(fpl%ctfsq_plane)
        call o%kill
        call fit%kill
        deallocate(fit, project, params)

    contains

        pure real function blob( xyz, centre, sigma2 ) result( val )
            real, intent(in) :: xyz(3), centre(3)
            real, intent(in) :: sigma2
            val = exp(-sum((xyz-centre)**2)/(2.0*sigma2))
        end function blob

        subroutine prepare_reconstructor( rec, rmat )
            type(reconstructor), intent(inout) :: rec
            real,                intent(in)    :: rmat(BOX,BOX,BOX)
            call rec%new_accumulator(params, project, expand=.true., wthreads=.false.)
            call rec%set_rmat(rmat, .false.)
            call rec%fft
            call rec%expand_exp
        end subroutine prepare_reconstructor

        subroutine build_exact_positions( hs, ks, npos )
            integer, allocatable, intent(out) :: hs(:), ks(:)
            integer,              intent(out) :: npos
            integer :: h, kk, p
            npos = count_halfplane(RHYB)
            allocate(hs(npos), ks(npos))
            p = 0
            do kk = -RHYB, 0
                do h = -RHYB, merge(0,RHYB,kk==0)
                    if( h*h + kk*kk > RHYB*(RHYB+1) ) cycle
                    p = p + 1; hs(p) = h; ks(p) = kk
                end do
            end do
        end subroutine build_exact_positions

        integer function count_halfplane( khi ) result( n )
            integer, intent(in) :: khi
            integer :: h, kk
            n = 0
            do kk = -khi, 0
                do h = -khi, merge(0,khi,kk==0)
                    if( h*h + kk*kk <= khi*(khi+1) ) n = n + 1
                end do
            end do
        end function count_halfplane

        !> Brute-force sum over the exact polar rows consumed by the former.  Unlike
        !! polar_ring_gram + dgemv, this never forms ring Gram tables.
        subroutine direct_polar_stats( Us, xws, wr, Go, bo, co, eo, myo )
            real,     intent(in)  :: Us(:,0:), xws(:)
            real(dp), intent(in)  :: wr(:)
            real(dp), intent(out) :: Go(NC,NC), bo(NC), co(NC), eo, myo
            real(dp) :: u0r, u0i, uqr, uqi, urr, uri, yr, yi, wri
            integer :: ir_, j_, iq, jq
            Go = 0.d0; bo = 0.d0; co = 0.d0; eo = 0.d0; myo = 0.d0
            do ir_ = 1, fit%estep%nk_es
                wri = wr(ir_)
                do j_ = fit%estep%pg_es%rbeg(ir_), fit%estep%pg_es%rend(ir_)
                    u0r = real(Us(2*j_-1,0),dp); u0i = real(Us(2*j_,0),dp)
                    yr  = real(xws(2*j_-1),dp);  yi  = real(xws(2*j_),dp)
                    eo  = eo  + wri*(u0r*u0r + u0i*u0i)
                    myo = myo + u0r*yr + u0i*yi
                    do iq = 1, NC
                        uqr = real(Us(2*j_-1,iq),dp); uqi = real(Us(2*j_,iq),dp)
                        bo(iq) = bo(iq) + uqr*yr + uqi*yi
                        co(iq) = co(iq) + wri*(uqr*u0r + uqi*u0i)
                        do jq = iq, NC
                            urr = real(Us(2*j_-1,jq),dp); uri = real(Us(2*j_,jq),dp)
                            Go(iq,jq) = Go(iq,jq) + wri*(uqr*urr + uqi*uri)
                        end do
                    end do
                end do
            end do
            do jq = 1, NC
                do iq = jq+1, NC
                    Go(iq,jq) = Go(jq,iq)
                end do
            end do
        end subroutine direct_polar_stats

        subroutine direct_cartesian_stats( rec0, recs, ori_here, plane, khi, &
                &Go, bo, co, eo, myo )
            type(reconstructor), intent(in)    :: rec0, recs(NC)
            class(ori),          intent(inout) :: ori_here
            type(fplane_type),   intent(in)    :: plane
            integer,             intent(in)    :: khi
            real(dp),            intent(out)   :: Go(NC,NC), bo(NC), co(NC), eo, myo
            type(kbinterpol) :: kb
            real :: rmat(3,3), loc(3), wx(LATENT_WDIM), wy(LATENT_WDIM), wz(LATENT_WDIM)
            integer :: h, kk, iq, jq, win(2,3), lb(3), ub(3), hp, kp
            logical :: l_conjg
            complex :: u0, uq(NC), val, tf, yv
            complex(dp) :: u0d, yd
            Go = 0.d0; bo = 0.d0; co = 0.d0; eo = 0.d0; myo = 0.d0
            kb   = kbinterpol(KBWINSZ, KBALPHA)
            rmat = ori_here%get_mat()
            lb   = lbound(rec0%cmat_exp); ub = ubound(rec0%cmat_exp)
            do kk = -khi, 0
                do h = -khi, merge(0,khi,kk==0)
                    if( h*h + kk*kk > khi*(khi+1) ) cycle
                    loc = real(h)*rmat(1,:) + real(kk)*rmat(2,:)
                    l_conjg = loc(1) < 0.
                    if( l_conjg ) loc = -loc
                    call oracle_weights(kb, loc, win, wx, wy, wz)
                    if( any(win(1,:) < lb) .or. any(win(2,:) > ub) ) cycle
                    hp = OSMPL_PAD_FAC*h; kp = OSMPL_PAD_FAC*kk
                    tf = plane%transfer_plane(hp,kp)
                    yv = plane%cmplx_plane(hp,kp)
                    val = oracle_gather(rec0, win, wx, wy, wz)
                    if( l_conjg ) val = conjg(val)
                    u0 = tf*val
                    do iq = 1, NC
                        val = oracle_gather(recs(iq), win, wx, wy, wz)
                        if( l_conjg ) val = conjg(val)
                        uq(iq) = tf*val
                    end do
                    u0d = cmplx(u0,kind=dp); yd = cmplx(yv,kind=dp)
                    eo  = eo  + real(conjg(u0d)*u0d,dp)
                    myo = myo + real(conjg(u0d)*yd,dp)
                    do iq = 1, NC
                        bo(iq) = bo(iq) + real(conjg(cmplx(uq(iq),kind=dp))*yd,dp)
                        co(iq) = co(iq) + real(conjg(cmplx(uq(iq),kind=dp))*u0d,dp)
                        do jq = iq, NC
                            Go(iq,jq) = Go(iq,jq) + real(conjg(cmplx(uq(iq),kind=dp))* &
                                &cmplx(uq(jq),kind=dp),dp)
                        end do
                    end do
                end do
            end do
            do jq = 1, NC
                do iq = jq+1, NC
                    Go(iq,jq) = Go(jq,iq)
                end do
            end do
        end subroutine direct_cartesian_stats

        subroutine oracle_weights( kb, loc, win, wx, wy, wz )
            type(kbinterpol), intent(in)  :: kb
            real,             intent(in)  :: loc(3)
            integer,          intent(out) :: win(2,3)
            real,             intent(out) :: wx(LATENT_WDIM), wy(LATENT_WDIM), wz(LATENT_WDIM)
            real :: base(3), w3(3), s
            integer :: d, iw
            iw = ceiling(KBWINSZ - 0.5)
            win(1,:) = nint(loc) - iw
            win(2,:) = nint(loc) + iw
            base = real(win(1,:)) - loc
            do d = 1, LATENT_WDIM
                w3 = kb%apod(base + real(d-1))
                wx(d) = w3(1); wy(d) = w3(2); wz(d) = w3(3)
            end do
            s = sum(wx)
            if( abs(s) > epsilon(1.0) )then; wx = wx/s; else; wx = 1.0/real(LATENT_WDIM); endif
            s = sum(wy)
            if( abs(s) > epsilon(1.0) )then; wy = wy/s; else; wy = 1.0/real(LATENT_WDIM); endif
            s = sum(wz)
            if( abs(s) > epsilon(1.0) )then; wz = wz/s; else; wz = 1.0/real(LATENT_WDIM); endif
        end subroutine oracle_weights

        complex function oracle_gather( rec, win, wx, wy, wz ) result( val )
            type(reconstructor), intent(in) :: rec
            integer,             intent(in) :: win(2,3)
            real,                intent(in) :: wx(LATENT_WDIM), wy(LATENT_WDIM), wz(LATENT_WDIM)
            integer :: ix, iy, iz
            val = CMPLX_ZERO
            do iz = 1, LATENT_WDIM
                do iy = 1, LATENT_WDIM
                    do ix = 1, LATENT_WDIM
                        val = val + rec%cmat_exp(win(1,1)+ix-1,win(1,2)+iy-1,win(1,3)+iz-1) * &
                            &(wx(ix)*wy(iy)*wz(iz))
                    end do
                end do
            end do
        end function oracle_gather

        real(dp) function stats_relerr( Ga, ba, ca, ea, ma, Ge, be, ce, ee, me ) result( err )
            real(dp), intent(in) :: Ga(NC,NC), ba(NC), ca(NC), ea, ma
            real(dp), intent(in) :: Ge(NC,NC), be(NC), ce(NC), ee, me
            err = max(maxval(abs(Ga-Ge))/max(maxval(abs(Ge)),tiny(1.0_dp)), &
                &maxval(abs(ba-be))/max(maxval(abs(be)),tiny(1.0_dp)), &
                &maxval(abs(ca-ce))/max(maxval(abs(ce)),tiny(1.0_dp)), &
                &abs(ea-ee)/max(abs(ee),tiny(1.0_dp)), abs(ma-me)/max(abs(me),tiny(1.0_dp)))
        end function stats_relerr

        real(dp) function linearity_relerr( Gin, bin, cin, ein, mval ) result( err )
            real(dp), intent(in) :: Gin(NC,NC), bin(NC), cin(NC), ein, mval
            real(dp) :: bscale, mscale, lhs(NC), lin_m
            lhs   = bin - cin - matmul(Gin, COEFF)
            lin_m = mval - ein - dot_product(COEFF, cin)
            bscale = max(maxval(abs(bin)), maxval(abs(cin)), maxval(abs(matmul(Gin,COEFF))))
            mscale = max(abs(mval), abs(ein), abs(dot_product(COEFF,cin)))
            err = max(maxval(abs(lhs))/max(bscale,tiny(1.0_dp)), &
                &abs(lin_m)/max(mscale,tiny(1.0_dp)))
        end function linearity_relerr

    end subroutine test_flex_estep_former

    subroutine test_flex_crossfsc( signed_pass, noise_pass, noise_mean, noise_floor )
        logical,  intent(out) :: signed_pass, noise_pass
        real(dp), intent(out) :: noise_mean, noise_floor
        integer, parameter :: BOX = 16, NC = 2
        type(flex_probe_fit), allocatable :: fits(:)
        class(parameters), allocatable :: params
        type(xfsc_ctx_t) :: ctx
        real, allocatable :: ra(:,:,:,:), rb(:,:,:,:)
        real(dp) :: proj
        integer, allocatable :: shell_count(:)
        integer :: f, q, fsz, lims(3,2), h, k, l, sh
        allocate(params, fits(2), ra(BOX,BOX,BOX,NC), rb(BOX,BOX,BOX,NC))
        fsz = max(1, fdim(BOX) - 1)
        params%box_crop = BOX
        params%smpd_crop = 1.0
        do f = 1, 2
            fits(f)%model%ncomp = NC
            fits(f)%spec%khi_full = fsz
            allocate(fits(f)%history%prev_real(NC))
            allocate(fits(f)%diag%xf_fscq(fsz,NC), fits(f)%diag%xf_h_e(fsz,NC), &
                &fits(f)%diag%xf_h_o(fsz,NC), fits(f)%diag%xf_cnt(fsz), fits(f)%diag%xf_gam(NC))
            fits(f)%diag%xf_fscq = 1.0
            fits(f)%diag%xf_h_e  = 1.0
            fits(f)%diag%xf_h_o  = 1.0
            fits(f)%diag%xf_cnt  = 1
            fits(f)%diag%xf_gam  = 1.d0
            do q = 1, NC
                call fits(f)%history%prev_real(q)%new([BOX,BOX,BOX], 1.0)
            end do
        end do

        call seed_rnd_fixed(20261006)
        call random_number(ra)
        ra = ra - 0.5
        proj = sum(real(ra(:,:,:,1),dp)*real(ra(:,:,:,2),dp)) / &
            &sum(real(ra(:,:,:,1),dp)**2)
        ra(:,:,:,2) = ra(:,:,:,2) - real(proj)*ra(:,:,:,1)
        rb(:,:,:,1) = -ra(:,:,:,2)
        rb(:,:,:,2) =  ra(:,:,:,1)
        call install_bases(fits, ra, rb)
        call init_context(ctx, fsz)
        call xfsc_paired_record(ctx, params, fits, 1)
        signed_pass = ctx%xf%nrec == 1
        if( signed_pass )then
            signed_pass = all(ctx%xf%recs(1)%match_a == [1,2]) .and. &
                &all(ctx%xf%recs(1)%match_b == [2,1]) .and. &
                &all(ctx%xf%recs(1)%match_sign == [1,-1]) .and. &
                &all(ctx%xf%recs(1)%match_cos > 1.0 - 2.0e-6) .and. &
                &all(abs(ctx%xf%recs(1)%fsc_cross - 1.0) < 2.0e-5)
        endif
        call ctx%kill

        call seed_rnd_fixed(20261007)
        call random_number(ra)
        call random_number(rb)
        ra = ra - 0.5
        rb = rb - 0.5
        call install_bases(fits, ra, rb)
        call init_context(ctx, fsz)
        call xfsc_paired_record(ctx, params, fits, 1)
        noise_mean = sum(abs(real(ctx%xf%recs(1)%fsc_cross(2:fsz,:),dp))) / &
            &real(max(1,(fsz-1)*NC),dp)
        allocate(shell_count(fsz), source=0)
        lims = fits(1)%history%prev_real(1)%loop_lims(2)
        do k = lims(2,1), lims(2,2)
            do h = lims(1,1), lims(1,2)
                do l = lims(3,1), lims(3,2)
                    sh = nint(sqrt(real(h*h + k*k + l*l)))
                    if( sh >= 1 .and. sh <= fsz ) shell_count(sh) = shell_count(sh) + 1
                end do
            end do
        end do
        ! For independent Fourier fields, shell FSC fluctuates on the 1/sqrt(N_shell)
        ! scale.  Three times the mean scale is the declared finite-sample floor.
        noise_floor = 3.d0*sum(1.d0/sqrt(real(max(1,shell_count(2:fsz)),dp))) / &
            &real(max(1,fsz-1),dp)
        noise_pass = noise_mean < noise_floor

        call ctx%kill
        call del_file(COV_XFSC_FNAME)
        call del_file('flex_pca_crossfsc_series.txt')
        do f = 1, 2
            call fits(f)%kill
        end do
        deallocate(fits, params, ra, rb, shell_count)

    contains

        subroutine install_bases( fits_, a, b )
            type(flex_probe_fit), intent(inout) :: fits_(2)
            real,                 intent(in)    :: a(BOX,BOX,BOX,NC), b(BOX,BOX,BOX,NC)
            integer :: iq
            do iq = 1, NC
                call fits_(1)%history%prev_real(iq)%set_rmat(a(:,:,:,iq), .false.)
                call fits_(2)%history%prev_real(iq)%set_rmat(b(:,:,:,iq), .false.)
            end do
        end subroutine install_bases

        subroutine init_context( ctx_, filtsz )
            type(xfsc_ctx_t), intent(inout) :: ctx_
            integer,          intent(in)    :: filtsz
            call ctx_%new
            ctx_%l_loaded   = .false.
            ctx_%pairing_id = 1
            ctx_%reg_active = 1
            ctx_%filtsz     = filtsz
        end subroutine init_context

    end subroutine test_flex_crossfsc

    subroutine test_probe_part_codec( pass )
        logical, intent(out) :: pass
        type(flex_probe_part) :: src, dst
        type(flex_pca_rounds_shmem) :: rounds
        class(parameters), allocatable :: params
        type(string) :: fname
        integer :: i, j, k, q, p
        real :: xr, xi
        src%ncomp  = 2
        src%nval   = 17
        src%nll_sum = 1.25d0
        allocate(src%cmat_e(3,2,2,2), src%cmat_o(3,2,2,2))
        allocate(src%rho_ex(3,2,2,2), src%rho_ox(3,2,2,2))
        allocate(src%rho_e(3,3,2,2), src%rho_o(3,3,2,2), src%gam_sum(2))
        allocate(src%mix_sr(2), src%mix_sm(2,2), src%mix_smm(2,2,2), src%mix_sainv(2,2))
        allocate(src%z_sub(5,2)); src%nz = 5
        allocate(src%kpk_e(3,4), src%kpk_o(3,4), src%rpk_e(2,4), src%rpk_o(2,4))
        dst%ncomp = src%ncomp
        allocate(dst%cmat_e(3,2,2,2), dst%cmat_o(3,2,2,2), source=cmplx(0.,0.))
        allocate(dst%rho_ex(3,2,2,2), dst%rho_ox(3,2,2,2), source=0.)
        allocate(dst%rho_e(3,3,2,2), dst%rho_o(3,3,2,2), source=0.)
        allocate(dst%gam_sum(2), source=0.d0)
        allocate(dst%mix_sr(2), dst%mix_sm(2,2), dst%mix_smm(2,2,2), dst%mix_sainv(2,2), source=0.d0)
        allocate(dst%z_sub(5,2), source=0.d0)
        allocate(dst%kpk_e(3,4), dst%kpk_o(3,4), source=0.)
        allocate(dst%rpk_e(2,4), dst%rpk_o(2,4), source=cmplx(0.,0.))
        call seed_rnd_fixed(20261008)
        do q = 1, 2
            do k = 1, 2
                do j = 1, 2
                    do i = 1, 3
                        call random_number(xr); call random_number(xi)
                        src%cmat_e(i,j,k,q) = cmplx(xr,xi)
                        call random_number(xr); call random_number(xi)
                        src%cmat_o(i,j,k,q) = cmplx(xr,xi)
                        call random_number(src%rho_ex(i,j,k,q))
                        call random_number(src%rho_ox(i,j,k,q))
                    end do
                end do
            end do
        end do
        call random_number(src%rho_e); call random_number(src%rho_o)
        call random_number(src%gam_sum); call random_number(src%mix_sr)
        call random_number(src%mix_sm); call random_number(src%mix_smm)
        call random_number(src%mix_sainv); call random_number(src%z_sub)
        call random_number(src%kpk_e); call random_number(src%kpk_o)
        do p = 1, 4
            do q = 1, 2
                call random_number(xr); call random_number(xi)
                src%rpk_e(q,p) = cmplx(xr,xi)
                call random_number(xr); call random_number(xi)
                src%rpk_o(q,p) = cmplx(xr,xi)
            end do
        end do
        allocate(params)
        params%numlen = 4
        fname = rounds%part_fname('probe', 1, params%numlen)
        call write_probe_part(fname, src)
        call reduce_probe_parts(params, rounds, dst)
        ! The reducer owns normal cleanup; keep the fixture idempotent if that lifecycle
        ! changes or an implementation returns after a successful read without deleting.
        call del_file(fname)
        pass = src%ncomp == dst%ncomp .and. src%nval == dst%nval .and. src%nll_sum == dst%nll_sum
        pass = pass .and. src%nz == dst%nz
        pass = pass .and. all(src%cmat_e == dst%cmat_e) .and. all(src%cmat_o == dst%cmat_o)
        pass = pass .and. all(src%rho_ex == dst%rho_ex) .and. all(src%rho_ox == dst%rho_ox)
        pass = pass .and. all(src%rho_e == dst%rho_e) .and. all(src%rho_o == dst%rho_o)
        pass = pass .and. all(src%gam_sum == dst%gam_sum)
        pass = pass .and. all(src%mix_sr == dst%mix_sr) .and. all(src%mix_sm == dst%mix_sm)
        pass = pass .and. all(src%mix_smm == dst%mix_smm) .and. all(src%mix_sainv == dst%mix_sainv)
        pass = pass .and. all(src%z_sub == dst%z_sub)
        pass = pass .and. all(src%kpk_e == dst%kpk_e) .and. all(src%kpk_o == dst%kpk_o)
        pass = pass .and. all(src%rpk_e == dst%rpk_e) .and. all(src%rpk_o == dst%rpk_o)
        call fname%kill
        call src%kill
        call dst%kill
        deallocate(params)
    end subroutine test_probe_part_codec

    !> Coupled two-component M-step on known noiseless planes, with a dense PCG reference.
    subroutine test_mstep_toy()
        integer,  parameter :: BOX = 8, NC = 2, NPROJ = 40, BAND = 3
        integer,  parameter :: NPAIRS         = (NC*(NC+1))/2
        real(dp), parameter :: SUBSPACE_FLOOR = 0.90_dp
        class(parameters), allocatable, target :: params
        class(sp_project), allocatable :: project
        class(sym), allocatable :: c1
        class(flex_pcg_t), allocatable :: pcg
        type(reconstructor), allocatable :: zero_rec, truth_recs(:), grid_recs(:), pcg_recs(:)
        type(image), allocatable :: truth_imgs(:), grid_imgs(:), pcg_imgs(:), gridcorr, mask_img
        type(fplane_type), allocatable :: fpls(:), basis_fpls(:), mean_fpl
        type(ori), allocatable :: orientations(:)
        type(flex_pcg_outcome_t) :: outcome
        real, allocatable :: truth(:,:,:,:), mask(:,:,:), rho(:,:,:,:)
        real, allocatable :: kacc(:,:), kpk(:,:), rhs_zero(:,:,:,:), precond_zero(:,:,:,:), bvol(:,:,:,:)
        real, allocatable :: ridge_probe(:,:,:,:), h_unreg(:,:,:,:), h_reg(:,:,:,:)
        complex, allocatable :: racc(:,:), rpk(:,:)
        real(dp), allocatable :: z(:,:), dens(:,:,:), sv_grid(:), sv_pcg(:)
        real(dp), allocatable :: hdense(:,:), bdense(:), xdense(:), xpcg(:)
        integer, allocatable :: dof(:,:)
        real(dp) :: grid_cos, pcg_cos, dense_rel, dense_resid, resid_tol, condition
        real(dp) :: theta, proj, nrm, roundoff, ridge_lift(NC), ridge_rel, ridge_tol
        real(dp) :: u(3)
        real :: xyz(3), ctr, rr
        real, pointer :: rp(:,:,:)
        integer :: i, j, k, q, r, n, lim, ndof, mid
        logical, allocatable :: valid(:)
        logical :: dense_ok
        write(*,'(A)') 'test_mstep_toy'
        allocate(params, project, c1, pcg, zero_rec, mask_img)
        allocate(truth_recs(NC), grid_recs(NC), pcg_recs(NC))
        allocate(truth_imgs(NC), grid_imgs(NC), pcg_imgs(NC))
        allocate(fpls(NPROJ), basis_fpls(NC), mean_fpl, orientations(NPROJ))
        allocate(truth(BOX,BOX,BOX,NC), mask(BOX,BOX,BOX), z(NC,NPROJ), dens(NC,NC,NPROJ))
        allocate(valid(NPROJ), source=.true.)
        params%box        = BOX
        params%box_crop   = BOX
        params%box_croppd = OSMPL_PAD_FAC*BOX
        params%smpd_crop  = 1.0
        params%oritype    = 'cls3D'
        call c1%new('c1')

        ctr = 0.5*real(BOX+1)
        do k = 1, BOX
            xyz(3) = real(k)-ctr
            do j = 1, BOX
                xyz(2) = real(j)-ctr
                do i = 1, BOX
                    xyz(1) = real(i)-ctr
                    rr = sum(xyz*xyz)
                    mask(i,j,k) = merge(1.0, 0.0, rr <= 2.8)
                    truth(i,j,k,1) = mask(i,j,k) * (exp(-sum((xyz-[0.8,-0.3,0.2])**2)/0.9) - &
                        &0.55*exp(-sum((xyz-[-0.7,0.4,-0.3])**2)/1.1))
                    truth(i,j,k,2) = mask(i,j,k) * (exp(-sum((xyz-[-0.2,0.8,0.5])**2)/0.8) - &
                        &0.50*exp(-sum((xyz-[0.4,-0.7,-0.5])**2)/1.0))
                end do
            end do
        end do
        nrm = sqrt(sum(real(truth(:,:,:,1),dp)**2))
        truth(:,:,:,1) = truth(:,:,:,1)/real(nrm)
        proj = sum(real(truth(:,:,:,1),dp)*real(truth(:,:,:,2),dp))
        truth(:,:,:,2) = truth(:,:,:,2) - real(proj)*truth(:,:,:,1)
        nrm = sqrt(sum(real(truth(:,:,:,2),dp)**2))
        truth(:,:,:,2) = truth(:,:,:,2)/real(nrm)

        call new_empty_reconstructor(zero_rec)
        do q = 1, NC
            call new_empty_reconstructor(truth_recs(q))
            call truth_recs(q)%set_rmat(truth(:,:,:,q), .false.)
            call truth_recs(q)%fft
            call truth_recs(q)%expand_exp
            call new_empty_reconstructor(grid_recs(q))
            call new_empty_reconstructor(pcg_recs(q))
            call truth_imgs(q)%new([BOX,BOX,BOX], 1.0, wthreads=.false.)
            call truth_imgs(q)%set_rmat(truth(:,:,:,q), .false.)
        end do
        call mask_img%new([BOX,BOX,BOX], 1.0, wthreads=.false.)
        call mask_img%set_rmat(mask, .false.)

        do n = 1, NPROJ
            theta = 2.0_dp*DPI*real(n-1,dp)/real(NPROJ,dp)
            z(1,n) = cos(theta) + 0.25_dp*cos(3.0_dp*theta)
            z(2,n) = sin(theta) - 0.20_dp*sin(2.0_dp*theta)
            do r = 1, NC
                do q = 1, NC
                    dens(q,r,n) = z(q,n)*z(r,n)
                end do
            end do
        end do

        call seed_rnd_fixed(20261009)
        lim = OSMPL_PAD_FAC*BAND
        do n = 1, NPROJ
            call random_number(u)
            call orientations(n)%new(.false.)
            call orientations(n)%set_euler(real([360.0_dp*u(1), &
                &acos(max(-1.0_dp,min(1.0_dp,2.0_dp*u(2)-1.0_dp)))*180.0_dp/DPI, 360.0_dp*u(3)]))
            call init_reference_plane(fpls(n), lim)
            call project_fplanes_mean_basis(zero_rec, truth_recs, orientations(n), fpls(n), mean_fpl, &
                &basis_fpls, apply_ctf_amp=.true.)
            fpls(n)%cmplx_plane = CMPLX_ZERO
            do q = 1, NC
                fpls(n)%cmplx_plane = fpls(n)%cmplx_plane + real(z(q,n))*basis_fpls(q)%cmplx_plane
            end do
        end do

        allocate(rho(NPAIRS,size(grid_recs(1)%cmat_exp,1),size(grid_recs(1)%cmat_exp,2), &
            &size(grid_recs(1)%cmat_exp,3)), source=0.0)
        call insert_planes_oversamp_coupled_batch_scaled(grid_recs, rho, c1, orientations, fpls, &
            &z, dens, valid, NPROJ)
        do q = 1, NC
            pcg_recs(q)%cmat_exp = grid_recs(q)%cmat_exp
        end do

        call pcg%new(BOX, 1.0, NC)
        call pcg%set_band(BAND)
        call pcg%set_window_volume(mask_img)
        call pcg%alloc_accum(kacc)
        call pcg%alloc_rhs_accum(racc)
        call pcg%accumulate(kacc, c1, orientations, fpls, dens, valid, NPROJ)
        call pcg%accumulate_rhs(racc, c1, orientations, fpls, z, valid, NPROJ)
        call pcg%alloc_packed(kpk)
        call pcg%alloc_rhs_packed(rpk)
        call pcg%fold_accum(kacc, kpk)
        call pcg%fold_rhs(racc, rpk)
        call pcg%finalize(kpk)
        allocate(bvol(BOX,BOX,BOX,NC), rhs_zero(BOX,BOX,BOX,NC), &
            &precond_zero(BOX,BOX,BOX,NC), source=0.0)
        call pcg%finalize_rhs(rpk, bvol)
        allocate(ridge_probe(BOX,BOX,BOX,NC), h_unreg(BOX,BOX,BOX,NC), &
            &h_reg(BOX,BOX,BOX,NC), source=0.0)
        mid = BOX/2
        ridge_probe(mid,mid,mid,:) = 1.0
        call pcg%apply_operator(ridge_probe, h_unreg)
        ! A zero-vector preconditioner application finalizes the density-derived floors and lam.
        call pcg%apply_precond(rhs_zero, rho, lbound(pcg_recs(1)%cmat_exp), precond_zero)
        call pcg%apply_operator(ridge_probe, h_reg)
        do q = 1, NC
            ridge_lift(q) = real(h_reg(mid,mid,mid,q)-h_unreg(mid,mid,mid,q),dp)
        end do
        ridge_tol = 64.0_dp*real(epsilon(1.0),dp) * max(1.0_dp, &
            &(sqrt(sum(real(h_unreg,dp)**2))+sqrt(sum(real(h_reg,dp)**2))) / &
            &max(sqrt(sum(ridge_lift**2)),tiny(1.0_dp)))
        h_reg = h_reg - h_unreg
        do q = 1, NC
            h_reg(mid,mid,mid,q) = h_reg(mid,mid,mid,q) - real(ridge_lift(q))
        end do
        ridge_rel = sqrt(sum(real(h_reg,dp)**2)) / max(sqrt(sum(ridge_lift**2)),tiny(1.0_dp))
        call assemble_dense_system(pcg, bvol, mask, hdense, bdense, xdense, dof, condition, dense_ok)
        ndof = size(xdense)

        call pcg%solve(pcg_recs, rho, rpk, ndof, 0.0, outcome, 'FLEX PCG MSTEP TOY')
        call solve_coupled_basis_exp(grid_recs, rho, NC)
        allocate(gridcorr, source=prep3D_inv_kbenvelope4mul([BOX,BOX,BOX], 1.0))
        call finalize_reconstructors(grid_recs, gridcorr, mask, grid_imgs)
        call finalize_reconstructors(pcg_recs,  gridcorr, mask, pcg_imgs)
        call cross_half_subspace_angles(truth_imgs, grid_imgs, NC, sv_grid)
        call cross_half_subspace_angles(truth_imgs, pcg_imgs,  NC, sv_pcg)
        grid_cos = minval(sv_grid)
        pcg_cos  = minval(sv_pcg)

        allocate(xpcg(ndof))
        do n = 1, ndof
            q = dof(1,n)
            call pcg_imgs(q)%get_rmat_ptr(rp)
            xpcg(n) = real(rp(dof(2,n),dof(3,n),dof(4,n)),dp)
        end do
        dense_resid = sqrt(sum((matmul(hdense,xpcg)-bdense)**2)) / &
            &max(sqrt(sum(bdense**2)), tiny(1.0_dp))
        dense_rel = sqrt(sum((xpcg-xdense)**2)) / max(sqrt(sum(xdense**2)), tiny(1.0_dp))
        roundoff = 256.0_dp*real(epsilon(1.0),dp)*sqrt(real(ndof,dp))
        ! The residual is the guard.  A factor two covers order/compiler variation around
        ! the single-precision accumulation estimate without condition-number amplification.
        resid_tol = 2.0_dp*roundoff

        call assert_true(all(ridge_lift > 0.0_dp), &
            &'the M-step toy dense operator includes a positive density-derived Tikhonov lift')
        call assert_true(ridge_rel <= ridge_tol, &
            &'the M-step toy Tikhonov lift is a diagonal shift on the supported probe')
        call assert_true(dense_ok, 'the M-step toy dense PCG operator is positive definite')
        call assert_true(grid_cos >= SUBSPACE_FLOOR, &
            &'the coupled gridding M-step recovers both directions of the known blob span')
        call assert_true(pcg_cos >= SUBSPACE_FLOOR, &
            &'the coupled PCG M-step recovers both directions of the known blob span')
        call assert_true(dense_resid <= resid_tol, &
            &'the coupled PCG M-step satisfies the independently assembled dense system')
        write(*,'(A,2(ES11.3,1X),A,ES11.3,A,ES11.3)') '  Tikhonov lift=', ridge_lift, &
            &' off-diagonal relative error=', ridge_rel, ' tolerance=', ridge_tol
        write(*,'(A,ES11.3,A,ES11.3)') '  min subspace cosine: gridding=', grid_cos, ' PCG=', pcg_cos
        write(*,'(A,ES11.3,A,ES11.3)') '  PCG/dense relative-error observation=', dense_rel, &
            &' condition=', condition
        write(*,'(A,ES11.3,A,ES11.3,A,ES11.3)') '  dense residual=', dense_resid, &
            &' tolerance=', resid_tol, ' accumulation estimate=', roundoff

        do n = 1, NPROJ
            call cleanup_plane(fpls(n))
            call orientations(n)%kill
        end do
        call cleanup_plane(mean_fpl)
        do q = 1, NC
            call cleanup_plane(basis_fpls(q))
            call truth_recs(q)%dealloc_rho; call truth_recs(q)%kill
            call grid_recs(q)%dealloc_rho;  call grid_recs(q)%kill
            call pcg_recs(q)%dealloc_rho;   call pcg_recs(q)%kill
            call truth_imgs(q)%kill; call grid_imgs(q)%kill; call pcg_imgs(q)%kill
        end do
        call zero_rec%dealloc_rho; call zero_rec%kill
        call gridcorr%kill; call mask_img%kill
        call pcg%kill; call c1%kill; call project%kill
        deallocate(params, project, c1, pcg, zero_rec)

    contains

        subroutine new_empty_reconstructor( rec )
            type(reconstructor), intent(inout) :: rec
            call rec%new_accumulator(params, project, expand=.true., wthreads=.false.)
            call rec%reset
            call rec%reset_exp
        end subroutine new_empty_reconstructor

        subroutine init_reference_plane( plane, bound )
            type(fplane_type), intent(inout) :: plane
            integer,           intent(in)    :: bound
            allocate(plane%cmplx_plane(-bound:bound,-bound:bound), source=CMPLX_ZERO)
            allocate(plane%transfer_plane(-bound:bound,-bound:bound), source=cmplx(1.0,0.0))
            allocate(plane%ctfsq_plane(-bound:bound,-bound:bound), source=1.0)
            plane%frlims = 0
            plane%frlims(1,:) = [-bound,bound]
            plane%frlims(2,:) = [-bound,bound]
            plane%nyq = bound
        end subroutine init_reference_plane

        subroutine finalize_reconstructors( recs, correction, support, imgs )
            type(reconstructor), intent(inout) :: recs(NC)
            type(image),         intent(in)    :: correction
            real,                intent(in)    :: support(BOX,BOX,BOX)
            type(image),         intent(inout) :: imgs(NC)
            real, pointer :: rmat(:,:,:)
            integer :: iq
            do iq = 1, NC
                recs(iq)%rho_exp = 1.0
                call recs(iq)%compress_exp
                call recs(iq)%ifft
                call recs(iq)%get_rmat_ptr(rmat)
                call imgs(iq)%new([BOX,BOX,BOX], 1.0, wthreads=.false.)
                call imgs(iq)%set_rmat(rmat(1:BOX,1:BOX,1:BOX), .false.)
                call imgs(iq)%mul(correction)
                call imgs(iq)%get_rmat_ptr(rmat)
                rmat(1:BOX,1:BOX,1:BOX) = rmat(1:BOX,1:BOX,1:BOX)*support
            end do
        end subroutine finalize_reconstructors

        subroutine assemble_dense_system( op, rhs, support, amat, bvec, xref, indices, cond, ok )
            class(flex_pcg_t),     intent(inout) :: op
            real,                  intent(in)    :: rhs(BOX,BOX,BOX,NC), support(BOX,BOX,BOX)
            real(dp), allocatable, intent(out)   :: amat(:,:), bvec(:), xref(:)
            integer,  allocatable, intent(out)   :: indices(:,:)
            real(dp),              intent(out)   :: cond
            logical,               intent(out)   :: ok
            real, allocatable :: unit(:,:,:,:), hunit(:,:,:,:)
            real(dp), allocatable :: awork(:,:), eval(:), evec(:,:)
            real(dp) :: lmax, lmin
            integer :: col, row, iq, ii, jj, kk, id, nrot, nsys
            nsys = NC*count(support > 0.0)
            allocate(indices(4,nsys))
            id = 0
            do iq = 1, NC
                do kk = 1, BOX
                    do jj = 1, BOX
                        do ii = 1, BOX
                            if( support(ii,jj,kk) <= 0.0 ) cycle
                            id = id + 1
                            indices(:,id) = [iq,ii,jj,kk]
                        end do
                    end do
                end do
            end do
            allocate(amat(nsys,nsys), bvec(nsys), xref(nsys), source=0.0_dp)
            allocate(unit(BOX,BOX,BOX,NC), hunit(BOX,BOX,BOX,NC), source=0.0)
            do col = 1, nsys
                unit = 0.0
                unit(indices(2,col),indices(3,col),indices(4,col),indices(1,col)) = 1.0
                call op%apply_operator(unit, hunit)
                do row = 1, nsys
                    amat(row,col) = real(hunit(indices(2,row),indices(3,row),indices(4,row),indices(1,row)),dp)
                end do
            end do
            amat = 0.5_dp*(amat+transpose(amat))
            do row = 1, nsys
                bvec(row) = real(rhs(indices(2,row),indices(3,row),indices(4,row),indices(1,row)),dp)
            end do
            allocate(awork, source=amat)
            allocate(eval(nsys), evec(nsys,nsys))
            call jacobi(awork, nsys, nsys, eval, evec, nrot)
            lmax = maxval(eval)
            lmin = minval(eval)
            ok = lmin > 128.0_dp*epsilon(1.0_dp)*max(1.0_dp,lmax)
            cond = huge(1.0_dp)
            if( ok ) cond = lmax/lmin
            call dense_cholesky(amat, bvec, xref, ok)
            deallocate(unit, hunit, awork, eval, evec)
        end subroutine assemble_dense_system

        subroutine dense_cholesky( amat, bvec, xvec, ok )
            real(dp), intent(in)    :: amat(:,:), bvec(:)
            real(dp), intent(out)   :: xvec(:)
            logical,  intent(inout) :: ok
            real(dp), allocatable :: lower(:,:), y(:)
            real(dp) :: pivot, tol
            integer :: ii, jj, nsys
            if( .not. ok ) return
            nsys = size(bvec)
            allocate(lower(nsys,nsys), y(nsys), source=0.0_dp)
            tol = 128.0_dp*epsilon(1.0_dp)*max(1.0_dp,maxval(abs(amat)))
            do ii = 1, nsys
                do jj = 1, ii
                    pivot = amat(ii,jj) - dot_product(lower(ii,1:jj-1),lower(jj,1:jj-1))
                    if( ii == jj )then
                        if( pivot <= tol )then
                            ok = .false.
                            deallocate(lower, y)
                            return
                        endif
                        lower(ii,jj) = sqrt(pivot)
                    else
                        lower(ii,jj) = pivot/lower(jj,jj)
                    endif
                end do
            end do
            do ii = 1, nsys
                y(ii) = (bvec(ii)-dot_product(lower(ii,1:ii-1),y(1:ii-1)))/lower(ii,ii)
            end do
            do ii = nsys, 1, -1
                xvec(ii) = (y(ii)-dot_product(lower(ii+1:nsys,ii),xvec(ii+1:nsys)))/lower(ii,ii)
            end do
            deallocate(lower, y)
        end subroutine dense_cholesky

    end subroutine test_mstep_toy

end module simple_flex_pca_tester
