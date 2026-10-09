!@descr: unit tests for simple_nano_species: intensity mixture and class choice, enclosed Gaussian fraction, threshold calibration
! Samples are drawn here at fixed seeds from the class means, spreads and sizes each test states, and the
! fits are judged against those generating parameters.
module simple_nano_species_tester
use simple_test_utils,   only: assert_true, assert_false, assert_int, assert_real, set_fixed_seed
use simple_defs,         only: dp, DPI
use simple_rnd,          only: gasdev
use simple_nano_species, only: fit_species_mixture, class_separation, enclosed_fraction, maxima_shape,&
    &calibrate_threshold, expected_false_count, fit_gauss_width, gauss_filter3D, local_maxima, robust_spread, MAX_NSPECIES
implicit none
private
public :: run_all_nano_species_tests

contains

    subroutine run_all_nano_species_tests()
        write(*,'(A)') '**** running all nano species tests ****'
        call test_enclosed_fraction()
        call test_one_class()
        call test_separation()
        call test_size_floor()
        call test_three_classes()
        call test_four_classes()
        call test_calibration()
        call test_gauss_width()
        call test_gauss_filter()
        call test_local_maxima()
        call test_robust_spread()
    end subroutine run_all_nano_species_tests

    ! x(i0+1:i0+m) drawn from N(mean, sdev**2)
    subroutine draw( x, i0, m, mean, sdev )
        real,    intent(inout) :: x(:)
        integer, intent(in)    :: i0, m
        real,    intent(in)    :: mean, sdev
        integer :: i
        do i = i0+1,i0+m
            x(i) = gasdev(mean, sdev)
        enddo
    end subroutine draw

    ! F(t) against Simpson integration of the chi density with three degrees of freedom
    subroutine test_enclosed_fraction()
        integer, parameter :: NSTEP = 2000
        real,    parameter :: TS(5) = [0.5, 1., 2., 3., 4.]
        real(dp) :: h, s, f
        integer  :: it, i
        logical  :: ok
        write(*,'(A)') 'test_enclosed_fraction'
        call assert_real(0., enclosed_fraction(0.), 0., 'enclosed_fraction(0) is 0')
        ok = .true.
        do it = 1,size(TS)
            h = real(TS(it), dp) / real(NSTEP, dp)
            s = 0._dp
            do i = 0,NSTEP
                f = sqrt(2._dp / DPI) * (i * h)**2 * exp(-0.5_dp * (i * h)**2)
                if( i == 0 .or. i == NSTEP )then
                    s = s + f
                elseif( mod(i,2) == 1 )then
                    s = s + 4._dp * f
                else
                    s = s + 2._dp * f
                endif
            enddo
            s = s * h / 3._dp
            if( abs(real(s) - enclosed_fraction(TS(it))) > 1.e-5 ) ok = .false.
        enddo
        call assert_true(ok, 'enclosed_fraction equals the integral of the 3D Gaussian radial density to 1e-5')
        call assert_real(0.198748, enclosed_fraction(1.), 1.e-5, 'enclosed_fraction(1) = erf(1/sqrt 2) - sqrt(2/pi) exp(-1/2)')
    end subroutine test_enclosed_fraction

    ! 400 intensities of one class: K = 1 with the closed-form BIC; nspecies forces K = 2
    subroutine test_one_class()
        integer, parameter :: N = 400
        real,    parameter :: MU0 = 1., SD0 = 0.05
        real, allocatable  :: post(:,:), mu(:), var(:)
        real,    allocatable :: bic(:)
        logical, allocatable :: adm(:)
        real     :: x(N)
        real(dp) :: m, v, lnl, bic1
        integer  :: labels(N), K
        write(*,'(A)') 'test_one_class'
        call set_fixed_seed(20261001)
        call draw(x, 0, N, MU0, SD0)
        call fit_species_mixture(x, SD0, 0, K, labels, post, mu, var, bic, adm)
        call assert_int(1, K, 'one class: K = 1')
        call assert_true(all(labels == 1), 'one class: every atom in class 1')
        call assert_true(all(abs(post(:,1) - 1.) < 1.e-6), 'one class: posteriors are 1')
        ! four standard errors of the mean
        call assert_real(MU0, mu(1), 4. * SD0 / sqrt(real(N)), 'one class: the mean is recovered')
        m    = sum(real(x, dp)) / N
        v    = max(sum((real(x, dp) - m)**2) / N, real(SD0, dp)**2)
        lnl  = -0.5_dp * N * log(2._dp * DPI * v) - 0.5_dp * sum((real(x, dp) - m)**2) / v
        bic1 = -2._dp * lnl + 2._dp * log(real(N, dp))
        call assert_real(real(bic1), bic(1), 1.e-4 * abs(real(bic1)), 'one class: BIC(1) is the closed form -2 lnL + 2 ln N')
        call assert_true(bic(1) < bic(2) .and. bic(1) < bic(3), 'one class: K = 1 has the lowest BIC')
        call assert_true(adm(1), 'K = 1 is always admissible')
        call fit_species_mixture(x, SD0, 2, K, labels, post, mu, var)
        call assert_int(2, K, 'nspecies = 2 fixes K')
        call assert_true(size(post,2) == 2 .and. mu(1) >= mu(2), 'nspecies = 2: two classes by decreasing mean')
        call assert_true(all(abs(sum(post, dim=2) - 1.) < 1.e-5), 'posteriors of an atom sum to 1')
    end subroutine test_one_class

    ! two classes of 200, spread SD0 each: at D = 2.5 the two-class fit is inadmissible and K = 1;
    ! at D = 4 K = 2 and the labels follow the generating classes up to the Bayes error (2.3%)
    subroutine test_separation()
        integer, parameter :: NC = 200, N = 2*NC
        real,    parameter :: SD0 = 0.05
        real, allocatable  :: post(:,:), mu(:), var(:)
        real,    allocatable :: bic(:)
        logical, allocatable :: adm(:)
        real    :: x(N)
        integer :: labels(N), truth(N), K
        write(*,'(A)') 'test_separation'
        truth(1:NC)   = 1
        truth(NC+1:N) = 2
        call set_fixed_seed(20261002)
        call draw(x,  0, NC, 1.,             SD0)
        call draw(x, NC, NC, 1. - 2.5 * SD0, SD0)
        call fit_species_mixture(x, SD0, 0, K, labels, post, mu, var, bic, adm)
        call assert_false(adm(2), 'D = 2.5: the two-class fit is inadmissible')
        call assert_int(1, K, 'D = 2.5: K = 1')
        call set_fixed_seed(20261003)
        call draw(x,  0, NC, 1.,           SD0)
        call draw(x, NC, NC, 1. - 4. * SD0, SD0)
        call fit_species_mixture(x, SD0, 0, K, labels, post, mu, var, bic, adm)
        call assert_true(adm(2), 'D = 4: the two-class fit is admissible')
        call assert_int(2, K, 'D = 4: K = 2')
        call assert_true(bic(2) < bic(1), 'D = 4: BIC(2) < BIC(1)')
        if( K == 2 )then
            call assert_true(count(labels == truth) >= nint(0.95 * N), 'D = 4: at least 95% of the labels are right')
            call assert_real(1.,           mu(1), 0.02, 'D = 4: the mean of class 1 is recovered')
            call assert_real(1. - 4. * SD0, mu(2), 0.02, 'D = 4: the mean of class 2 is recovered')
            call assert_true(class_separation(mu(1), var(1), mu(2), var(2)) >= 3., 'D = 4: the fitted separation passes the floor')
        endif
    end subroutine test_separation

    ! N = 400 sets the size floor at max(8, 0.02 N) = 8 atoms: a distinct class of 5 is not admitted,
    ! one of 8 is, and its atoms are labelled 2
    subroutine test_size_floor()
        integer, parameter :: N = 400
        real,    parameter :: SD0 = 0.05
        real, allocatable  :: post(:,:), mu(:), var(:)
        real,    allocatable :: bic(:)
        logical, allocatable :: adm(:)
        real    :: x(N)
        integer :: labels(N), K
        write(*,'(A)') 'test_size_floor'
        call set_fixed_seed(20261004)
        call draw(x,     0, N-5, 1.0, SD0)
        call draw(x, N-5,     5, 0.5, SD0)
        call fit_species_mixture(x, SD0, 0, K, labels, post, mu, var, bic, adm)
        call assert_false(adm(2), 'a class of 5 of 400 atoms is below the size floor')
        call assert_int(1, K, 'a class below the size floor gives K = 1')
        call set_fixed_seed(20261005)
        call draw(x,     0, N-8, 1.0, SD0)
        call draw(x, N-8,     8, 0.5, SD0)
        call fit_species_mixture(x, SD0, 0, K, labels, post, mu, var, bic, adm)
        call assert_true(adm(2), 'a class of 8 of 400 atoms meets the size floor')
        call assert_int(2, K, 'a class at the size floor gives K = 2')
        call assert_true(all(labels(1:N-8) == 1) .and. all(labels(N-7:N) == 2), 'a class at the size floor: every label right')
    end subroutine test_size_floor

    ! case 12 of the plan: 264, 132 and 132 atoms at intensities 1, 0.5 and 0.167, spread 0.04 (D >= 8)
    subroutine test_three_classes()
        integer, parameter :: N1 = 264, N2 = 132, N3 = 132, N = N1+N2+N3
        real,    parameter :: SD0 = 0.04, MUS(3) = [1., 0.5, 0.167]
        real, allocatable  :: post(:,:), mu(:), var(:)
        real,    allocatable :: bic(:)
        logical, allocatable :: adm(:)
        real    :: x(N)
        integer :: labels(N), truth(N), K
        write(*,'(A)') 'test_three_classes'
        call set_fixed_seed(20261006)
        ! drawn in the order 3, 1, 2 so that class numbers cannot follow the input order
        call draw(x,       0, N3, MUS(3), SD0)
        call draw(x,      N3, N1, MUS(1), SD0)
        call draw(x, N3 + N1, N2, MUS(2), SD0)
        truth(1:N3)          = 3
        truth(N3+1:N3+N1)    = 1
        truth(N3+N1+1:N)     = 2
        call fit_species_mixture(x, SD0, 0, K, labels, post, mu, var, bic, adm)
        call assert_int(3, K, 'three classes: K = 3')
        call assert_true(adm(3) .and. bic(3) < bic(2) .and. bic(2) < bic(1), 'three classes: admissible, BIC(3) < BIC(2) < BIC(1)')
        if( K == 3 )then
            call assert_true(all(labels == truth), 'three classes: every label right, classes by decreasing intensity')
            call assert_true(all(abs(mu - MUS) < 0.02), 'three classes: the class intensities are recovered')
        endif
    end subroutine test_three_classes

    ! four classes of 100 at intensities 1, 0.6, 0.35 and 0.2, spread 0.02: nspecies = 4 sizes the fit to four
    ! classes; the automatic search stops at its ceiling MAX_NSPECIES
    subroutine test_four_classes()
        integer, parameter :: NC = 100, N = 4*NC
        real,    parameter :: SD0 = 0.02, MUS(4) = [1., 0.6, 0.35, 0.2]
        integer, parameter :: ORDER(4) = [4, 2, 1, 3]
        real,    allocatable :: post(:,:), mu(:), var(:), bic(:)
        logical, allocatable :: adm(:)
        real    :: x(N)
        integer :: labels(N), truth(N), K, ic
        write(*,'(A)') 'test_four_classes'
        call set_fixed_seed(20261007)
        ! drawn in the order 4, 2, 1, 3 so that class numbers cannot follow the input order
        do ic = 1,4
            call draw(x, (ic-1) * NC, NC, MUS(ORDER(ic)), SD0)
            truth((ic-1)*NC+1:ic*NC) = ORDER(ic)
        enddo
        call fit_species_mixture(x, SD0, 4, K, labels, post, mu, var, bic, adm)
        call assert_int(4, K, 'four classes, nspecies = 4: K = 4')
        call assert_int(4, size(bic), 'four classes, nspecies = 4: a BIC per K up to 4')
        call assert_int(4, size(adm), 'four classes, nspecies = 4: an admissibility per K up to 4')
        call assert_true(size(post,1) == N .and. size(post,2) == 4, 'four classes, nspecies = 4: four posterior columns')
        if( K == 4 )then
            call assert_true(all(labels == truth), 'four classes: every label right, classes by decreasing intensity')
            call assert_true(all(abs(mu - MUS) < 0.02), 'four classes: the class intensities are recovered within 0.02')
            call assert_true(adm(4), 'four classes: the four-class fit is admissible')
        endif
        call fit_species_mixture(x, SD0, 0, K, labels, post, mu, var, bic, adm)
        call assert_int(MAX_NSPECIES, K, 'four classes, automatic: K is the ceiling of the search')
        call assert_int(MAX_NSPECIES, size(bic), 'four classes, automatic: a BIC per K up to the ceiling')
    end subroutine test_four_classes

    ! white-noise region counts of section 8.8 of the plan (154.5, 38.0, 7.2 maxima above 2.5, 3.0, 3.5;
    ! search-to-region volume ratio 3.628); the plan's independent script gave C = 2105 and k = 4.99
    subroutine test_calibration()
        real, parameter :: U(3) = [2.5, 3.0, 3.5], COUNTS(3) = [154.5, 38.0, 7.2]
        real, parameter :: RATIO = 3.628, TARGET = 0.2, C0 = 1000.
        real :: k, c_search, k2, c2
        write(*,'(A)') 'test_calibration'
        call calibrate_threshold(U, COUNTS, RATIO, TARGET, k, c_search)
        call assert_real(2105.3, c_search, 2., 'calibration: Poisson estimate of C from the section 8.8 counts')
        call assert_real(5.0, k, 0.1, 'calibration: k = 5.0 within 0.1 for the section 8.8 white-noise counts')
        call assert_real(TARGET, expected_false_count(c_search, k), 1.e-4, 'calibration: target false maxima expected at k')
        ! counts that follow the model exactly return its C
        call calibrate_threshold(U, C0 * maxima_shape(U), 1., TARGET, k2, c2)
        call assert_real(C0, c2, 1.e-3 * C0, 'calibration: counts from the model return its C')
        ! a larger search volume needs a higher threshold
        call calibrate_threshold(U, COUNTS, 2. * RATIO, TARGET, k2, c2)
        call assert_true(k2 > k, 'calibration: doubling the search volume raises k')
    end subroutine test_calibration

    ! samples of A exp(-4 pi**2 r**2 / B) on the 0.358 A grid within 1.38 A of a voxel centre
    subroutine gauss_samples( amp, bfac, y, r2 )
        real,              intent(in)  :: amp, bfac
        real, allocatable, intent(out) :: y(:), r2(:)
        real, parameter :: SMPD = 0.358, RAD = 1.38
        real    :: d2(15**3)
        integer :: i, j, k, n
        n = 0
        do k = -7,7
            do j = -7,7
                do i = -7,7
                    if( real(i*i + j*j + k*k) * SMPD**2 > RAD**2 ) cycle
                    n     = n + 1
                    d2(n) = real(i*i + j*j + k*k) * SMPD**2
                enddo
            enddo
        enddo
        r2 = d2(:n)
        y  = amp * exp(-4. * real(DPI)**2 * r2 / bfac)
    end subroutine gauss_samples

    ! noiseless samples return their amplitude and B; with no information in the data the prior or the class
    ! tie decides; the tied amplitude integrates to the class intensity
    subroutine test_gauss_width()
        real, allocatable :: y(:), r2(:)
        real :: amp, bfac, ic
        write(*,'(A)') 'test_gauss_width'
        call gauss_samples(2.0, 15.0, y, r2)
        call fit_gauss_width(y, r2, 1.e-6, log(13.9), 0.5, amp, bfac)
        call assert_real(15.0, bfac, 0.01,  'stage 1: noiseless samples return their B')
        call assert_real(2.0,  amp,  0.002, 'stage 1: noiseless samples return their amplitude')
        call fit_gauss_width(0. * y, r2, 1.e6, log(20.), 0.1, amp, bfac)
        call assert_real(20.0, bfac, 0.05, 'stage 1: without signal the prior sets B')
        ic = 2.0 * (15.0 / (4. * real(DPI)))**1.5
        call fit_gauss_width(y, r2, 1.e-6, 0., 1., amp, bfac, i_class=ic, tau_class=0.05 * ic)
        call assert_real(15.0, bfac, 0.01, 'stage 2: samples that integrate to the class intensity return their B')
        ! weak data, tied to a class of intensity ic: amplitude times (B / 4 pi)**1.5 is ic
        call gauss_samples(0.5, 15.0, y, r2)
        call fit_gauss_width(y, r2, 1.e4, 0., 1., amp, bfac, i_class=ic, tau_class=0.001 * ic)
        call assert_real(ic, amp * (bfac / (4. * real(DPI)))**1.5, 0.01 * ic, 'stage 2: the tie fixes the integral')
    end subroutine test_gauss_width

    ! a unit impulse becomes the normalised separable Gaussian
    subroutine test_gauss_filter()
        integer, parameter :: N = 21, C = 11
        real, parameter    :: SIG = 1.5
        real    :: a(N,N,N), b(N,N,N), w(-6:6)
        integer :: l
        write(*,'(A)') 'test_gauss_filter'
        a = 0.
        a(C,C,C) = 1.
        call gauss_filter3D(a, SIG, b)
        do l = -6,6
            w(l) = exp(-0.5 * (real(l) / SIG)**2)
        enddo
        w = w / sum(w)
        call assert_real(1., sum(b), 1.e-5, 'gauss_filter3D preserves the sum')
        call assert_real(w(2) * w(-1) * w(0), b(C+2,C-1,C), 1.e-7, 'gauss_filter3D is the separable normalised kernel')
    end subroutine test_gauss_filter

    ! two isolated bumps above threshold and one below; a lone spike lacks neighbours above threshold
    subroutine test_local_maxima()
        integer, parameter :: N = 24
        real,    allocatable :: vals(:)
        integer, allocatable :: ijk(:,:)
        real    :: z(N,N,N)
        logical :: mask(N,N,N)
        integer :: i, j, k
        write(*,'(A)') 'test_local_maxima'
        do k = 1,N
            do j = 1,N
                do i = 1,N
                    z(i,j,k) = 6. * exp(-real((i-6)**2 + (j-6)**2 + (k-6)**2) / 4.) +&
                    &          5.5 * exp(-real((i-18)**2 + (j-12)**2 + (k-8)**2) / 4.) +&
                    &          3. * exp(-real((i-12)**2 + (j-18)**2 + (k-18)**2) / 4.)
                enddo
            enddo
        enddo
        z(20,20,4) = 9.
        mask = .true.
        call local_maxima(z, mask, 4., 3, ijk, vals)
        call assert_int(2, size(vals), 'local_maxima: two maxima above 4 with three neighbours above 4')
        if( size(vals) == 2 )then
            call assert_true(all(ijk(:,1) == [6,6,6]) .and. all(ijk(:,2) == [18,12,8]), 'local_maxima: positions by decreasing value')
        endif
        call local_maxima(z, mask, 4., 0, ijk, vals)
        call assert_int(3, size(vals), 'local_maxima: the lone spike counts without the neighbour condition')
        mask(:,:,1:10) = .false.
        call local_maxima(z, mask, 2., 3, ijk, vals)
        call assert_int(1, size(vals), 'local_maxima: only maxima inside the mask')
    end subroutine test_local_maxima

    subroutine test_robust_spread()
        real :: x(9)
        write(*,'(A)') 'test_robust_spread'
        x = [1., 2., 3., 4., 5., 6., 7., 8., 100.]
        ! median 5, absolute deviations 4 3 2 1 0 1 2 3 95, their median 2
        call assert_real(1.4826 * 2., robust_spread(x), 1.e-5, 'robust_spread is 1.4826 MAD and ignores the outlier')
    end subroutine test_robust_spread

end module simple_nano_species_tester
