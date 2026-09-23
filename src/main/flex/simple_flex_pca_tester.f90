!@descr: unit tests for the flex_pca latent model, state weights and deconvolution (simple_flex_pca_model, _weights, _util, _deconv)
! The resume cache (bit-exact round trip of all seven payloads), the derived settings
! (box_crop, min_neff and state count against the validated IgG and Ribosembly scales),
! state placement with a population floor on two clusters plus outliers, kernel weights
! at a bandwidth with bounded widening, covariance state weights on a bimodal embedding,
! and latent deconvolution: noise-scale calibration, the held-out choice of K and
! posterior means closer to the truth than the raw latents. Random draws come from a
! fixed seed; the statistical checks have margins far outside their sampling noise.
! The fast suite deconvolves 4000 particles at noise variance 0.5..5 per axis (the K
! ladder stops at n/2000 = 2); the library suite repeats it on 20000 particles at the
! realistic 2..20, where the ladder runs to K = 4 and must still stop at 2.
module simple_flex_pca_tester
use simple_core_module_api,  only: dp, DPI
use simple_syslib,           only: del_file
use simple_flex_pca_model,   only: write_embedding_cache, read_embedding_cache, place_states_with_population_floor, &
    &auto_box_crop, auto_min_neff, auto_state_count, FLEX_AUTO_K_MIN, FLEX_AUTO_K_START
use simple_flex_pca_util,    only: kernel_weights_at_bandwidth, COV_MAX_BW_GROW
use simple_flex_pca_weights, only: build_covariance_state_weights
use simple_flex_pca_deconv,  only: calibrate_noise_scale, deconvolve_latent
use simple_test_utils
implicit none
private
public :: run_all_flex_pca_tests, run_all_flex_pca_lib_tests

contains

    subroutine run_all_flex_pca_tests()
        write(*,'(A)') '**** running all flex_pca tests ****'
        call test_embedding_cache_io()
        call test_auto_settings()
        call test_population_floor()
        call test_kernel_bandwidth()
        call test_state_weights()
        call test_deconvolution(4000, 0.5d0, 5.d0)
    end subroutine run_all_flex_pca_tests

    subroutine run_all_flex_pca_lib_tests()
        write(*,'(A)') '**** running all flex_pca library tests ****'
        call test_deconvolution(20000, 2.d0, 20.d0)
    end subroutine run_all_flex_pca_lib_tests

    !> the resume path: every payload survives a write/read round trip bit for bit
    subroutine test_embedding_cache_io()
        integer,  parameter :: NP = 37, NC = 4
        character(len=*), parameter :: FN = 'test_flex_pca_cache.bin'
        integer  :: pinds(NP), i, q, r, ncomp_rd
        real(dp) :: z(NP,NC), eigvals(NC), contrast(NP), re(NP), rme(NP)
        real(dp) :: prec(NC,NC,NP), sig2, sig2_rd
        real(dp), allocatable :: z_rd(:,:), eig_rd(:), con_rd(:), re_rd(:), rme_rd(:), prec_rd(:,:,:)
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
        call write_embedding_cache(FN, pinds, NP, NC, z, eigvals, contrast, re, rme, prec, sig2)
        call read_embedding_cache(FN, pinds, NP, ncomp_rd, z_rd, eig_rd, con_rd, re_rd, rme_rd, prec_rd, sig2_rd)
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
        integer,  parameter :: NA = 300, NB = 250, NO = 20, NP = NA+NB+NO, NC = 2, NST = 2
        real,     parameter :: FRAC = 0.2
        integer  :: i, q, state, nmin, nlab(NST)
        real(dp) :: z(NP,NC), eigvals(NC), prec(NC,NC,NP)
        real,     allocatable :: weights(:,:), targets(:,:), bandwidths(:), neff(:)
        integer,  allocatable :: labels(:)
        logical  :: indicators, neff_match, floor_met
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
        call place_states_with_population_floor(z, NP, NC, NC, NST, 0, 10, FRAC, eigvals, prec, &
            &weights, targets, bandwidths, neff, labels)
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
        integer,  parameter :: NP = 400
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
        call assert_true(h_wide <= 10.d0*1.3d0**COV_MAX_BW_GROW + 1.d-9, 'the widening stays within its cap')
        call assert_true(count(w > 0.) > nsupp, 'widening increases the support')
        call assert_true(neff_wide > neff, 'neff increases with the bandwidth')
    end subroutine test_kernel_bandwidth

    !> covariance state weights on a clearly bimodal embedding
    subroutine test_state_weights()
        integer,  parameter :: NPC = 150, NP = 2*NPC, NC = 2, NST = 2
        integer  :: i, q, state, nlab(NST)
        real(dp) :: z(NP,NC), eigvals(NC), prec(NC,NC,NP)
        real,     allocatable :: weights(:,:), targets(:,:), bandwidths(:), neff(:)
        integer,  allocatable :: labels(:)
        logical  :: in_hull
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
        call build_covariance_state_weights(z, NP, NC, NC, NST, 0, 10, eigvals, prec, &
            &weights, targets, bandwidths, neff, labels)
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
        integer,  parameter  :: D = 4
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
        call calibrate_noise_scale(z, zhalf, prec, prior, n, D, a, a_comp)
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

end module simple_flex_pca_tester
