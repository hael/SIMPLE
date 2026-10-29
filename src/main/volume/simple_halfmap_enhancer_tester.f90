!@descr: unit tests for local half-map Wiener sharpening and its confidence products
module simple_halfmap_enhancer_tester
use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
use simple_defs,            only: DP, logfhandle
use simple_image,           only: image
use simple_halfmap_enhancer, only: halfmap_enhancer
use simple_test_utils,      only: assert_real, assert_true, set_fixed_seed
implicit none
private

public :: run_all_halfmap_enhancer_tests

integer, parameter :: BOX  = 32
integer, parameter :: SEED = 20261029
real,    parameter :: SMPD = 1.5

contains

    subroutine run_all_halfmap_enhancer_tests()
        write(logfhandle,'(A)') '**** running all half-map enhancer tests ****'
        call test_constant_pair()
        call test_unsupported_pair()
        call test_noisy_blob_pair()
        call test_correlated_background_pair()
    end subroutine run_all_halfmap_enhancer_tests

    !> A constant solvent level contains no half-map-supported structure.  It is
    !! preserved as the estimated background while all confidence products vanish.
    subroutine test_constant_pair()
        real, allocatable :: even(:,:,:), odd(:,:,:), enhanced(:,:,:)
        real, allocatable :: confidence(:,:,:), locres(:,:,:), envelope(:,:,:)
        write(logfhandle,'(A)') '-- test_constant_pair'
        allocate(even(BOX,BOX,BOX), source=2.)
        allocate(odd(BOX,BOX,BOX),  source=2.)
        call run_enhancer(even, odd, enhanced, confidence, locres, envelope)
        call assert_real(2., maxval(enhanced), 1.e-5, 'a constant pair preserves its background level')
        call assert_real(2., minval(enhanced), 1.e-5, 'the preserved constant pair stays spatially uniform')
        call assert_real(0., maxval(confidence), 1.e-6, 'a constant pair has no supported confidence')
        call assert_real(0., maxval(envelope), 1.e-6, 'a constant pair has an empty molecular envelope')
        call assert_real(0., maxval(locres), 1.e-6, 'a constant pair has no local resolution')
    end subroutine test_constant_pair

    !> Equal and opposite half-maps have signal in each input but no common signal.
    subroutine test_unsupported_pair()
        real, allocatable :: even(:,:,:), odd(:,:,:), enhanced(:,:,:), signal(:,:,:)
        real, allocatable :: confidence(:,:,:), locres(:,:,:), envelope(:,:,:)
        write(logfhandle,'(A)') '-- test_unsupported_pair'
        signal = blob_map()
        even = signal
        odd  = -signal
        call run_enhancer(even, odd, enhanced, confidence, locres, envelope)
        call assert_real(0., maxval(abs(enhanced)), 1.e-5, 'opposite half-map coefficients are suppressed')
        call assert_real(0., maxval(confidence), 1.e-6, 'opposite half-maps have zero confidence')
        call assert_real(0., maxval(envelope), 1.e-6, 'opposite half-maps have an empty envelope')
    end subroutine test_unsupported_pair

    !> Independent half-map noise around a common compact phantom is reduced.  The
    !! test expectations come from that construction, not from the implementation.
    subroutine test_noisy_blob_pair()
        real, allocatable :: truth(:,:,:), even(:,:,:), odd(:,:,:), raw(:,:,:)
        real, allocatable :: enhanced(:,:,:), confidence(:,:,:), locres(:,:,:), envelope(:,:,:)
        logical, allocatable :: core(:,:,:), solvent(:,:,:)
        real :: corr_raw, corr_enhanced, core_confidence, solvent_confidence
        write(logfhandle,'(A)') '-- test_noisy_blob_pair'
        call set_fixed_seed(SEED)
        truth = blob_map()
        even  = truth + 0.65 * white_noise()
        odd   = truth + 0.65 * white_noise()
        raw   = 0.5 * (even + odd)
        call run_enhancer(even, odd, enhanced, confidence, locres, envelope)
        core    = radius_sq() < 64.
        solvent = radius_sq() > 169.
        corr_raw      = array_correlation(truth, raw, core)
        corr_enhanced = array_correlation(truth, enhanced, core)
        core_confidence    = sum(confidence, mask=core) / real(count(core))
        solvent_confidence = sum(confidence, mask=solvent) / real(count(solvent))
        write(logfhandle,'(A,4F10.4)') 'correlation raw/enhanced, confidence core/solvent: ', &
            &corr_raw, corr_enhanced, core_confidence, solvent_confidence
        call assert_true(all(ieee_is_finite(enhanced)), 'the enhanced noisy map is finite')
        call assert_true(all(ieee_is_finite(confidence)), 'the confidence map is finite')
        call assert_true(minval(confidence) >= 0. .and. maxval(confidence) <= 1., &
            &'confidence stays in the unit interval')
        call assert_true(minval(envelope) >= 0. .and. maxval(envelope) <= 1., &
            &'the conservative envelope stays in the unit interval')
        call assert_true(corr_enhanced > corr_raw + 0.03, &
            &'local Wiener filtering improves correlation with the independent phantom')
        call assert_true(core_confidence > max(5. * solvent_confidence, solvent_confidence + 0.002), &
            &'confidence is higher in the common-signal core than in solvent')
        call assert_true(any(pack(locres > 0., core)), 'the supported core receives a local resolution')
    end subroutine test_noisy_blob_pair

    !> Shared solvent texture is not independent half-map evidence and must not
    !! expand the molecular envelope across the box.  This emulates archived
    !! half-maps whose outer regions retain correlated processing background.
    subroutine test_correlated_background_pair()
        real, allocatable :: truth(:,:,:), shared(:,:,:), even(:,:,:), odd(:,:,:)
        real, allocatable :: enhanced(:,:,:), confidence(:,:,:), locres(:,:,:), envelope(:,:,:)
        logical, allocatable :: core(:,:,:), solvent(:,:,:)
        real :: core_envelope, solvent_envelope
        write(logfhandle,'(A)') '-- test_correlated_background_pair'
        call set_fixed_seed(SEED + 1)
        truth  = blob_map()
        shared = 0.10 * white_noise()
        even   = truth + shared + 0.20 * white_noise()
        odd    = truth + shared + 0.20 * white_noise()
        call run_enhancer(even, odd, enhanced, confidence, locres, envelope)
        core    = radius_sq() < 64.
        solvent = radius_sq() > 169.
        core_envelope    = sum(envelope, mask=core) / real(count(core))
        solvent_envelope = sum(envelope, mask=solvent) / real(count(solvent))
        write(logfhandle,'(A,2F10.4)') 'correlated-background envelope core/solvent: ', &
            &core_envelope, solvent_envelope
        call assert_true(core_envelope > solvent_envelope + 0.10, &
            &'the molecular core remains distinguishable from correlated solvent background')
        call assert_true(solvent_envelope < 0.40, &
            &'correlated solvent background does not fill the molecular envelope')
    end subroutine test_correlated_background_pair

    subroutine run_enhancer( even_r, odd_r, enhanced_r, confidence_r, locres_r, envelope_r )
        real, intent(in) :: even_r(:,:,:), odd_r(:,:,:)
        real, allocatable, intent(out) :: enhanced_r(:,:,:), confidence_r(:,:,:)
        real, allocatable, intent(out) :: locres_r(:,:,:), envelope_r(:,:,:)
        type(image) :: even, odd, enhanced, confidence, locres, envelope
        class(halfmap_enhancer), allocatable :: enhancer
        allocate(enhancer)
        call even%new([BOX,BOX,BOX], SMPD, wthreads=.false.)
        call odd%new([BOX,BOX,BOX], SMPD, wthreads=.false.)
        call even%set_rmat(even_r, .false.)
        call odd%set_rmat(odd_r, .false.)
        call enhancer%new([BOX,BOX,BOX], SMPD, window=16, stride=8)
        call enhancer%enhance(even, odd, enhanced, confidence, locres, envelope)
        enhanced_r   = enhanced%get_rmat()
        confidence_r = confidence%get_rmat()
        locres_r     = locres%get_rmat()
        envelope_r   = envelope%get_rmat()
        call enhancer%kill()
        deallocate(enhancer)
        call even%kill()
        call odd%kill()
        call enhanced%kill()
        call confidence%kill()
        call locres%kill()
        call envelope%kill()
    end subroutine run_enhancer

    function blob_map() result( rmat )
        real, allocatable :: rmat(:,:,:)
        real, parameter :: CENTRES(3,3) = reshape([-4.,2.,1., 5.,-3.,0., 1.,5.,-4.], [3,3])
        real :: d2
        integer :: i, j, k, b
        allocate(rmat(BOX,BOX,BOX), source=0.)
        do k = 1, BOX
            do j = 1, BOX
                do i = 1, BOX
                    do b = 1, size(CENTRES,2)
                        d2 = (real(i - BOX/2 - 1) - CENTRES(1,b))**2 &
                            &+ (real(j - BOX/2 - 1) - CENTRES(2,b))**2 &
                            &+ (real(k - BOX/2 - 1) - CENTRES(3,b))**2
                        rmat(i,j,k) = rmat(i,j,k) + exp(-0.5 * d2 / 4.)
                    enddo
                enddo
            enddo
        enddo
    end function blob_map

    function white_noise() result( rmat )
        real, allocatable :: rmat(:,:,:)
        allocate(rmat(BOX,BOX,BOX))
        call random_number(rmat)
        rmat = sqrt(12.) * (rmat - 0.5)
    end function white_noise

    function radius_sq() result( r2 )
        real, allocatable :: r2(:,:,:)
        integer :: i, j, k
        allocate(r2(BOX,BOX,BOX))
        do k = 1, BOX
            do j = 1, BOX
                do i = 1, BOX
                    r2(i,j,k) = real((i - BOX/2 - 1)**2 + (j - BOX/2 - 1)**2 + (k - BOX/2 - 1)**2)
                enddo
            enddo
        enddo
    end function radius_sq

    real function array_correlation( a, b, mask ) result( corr )
        real,    intent(in) :: a(:,:,:), b(:,:,:)
        logical, intent(in) :: mask(:,:,:)
        real(DP) :: amean, bmean, numer, aden, bden
        integer :: n
        n = count(mask)
        amean = sum(real(a,DP), mask=mask) / real(n,DP)
        bmean = sum(real(b,DP), mask=mask) / real(n,DP)
        numer = sum((real(a,DP) - amean) * (real(b,DP) - bmean), mask=mask)
        aden  = sum((real(a,DP) - amean)**2, mask=mask)
        bden  = sum((real(b,DP) - bmean)**2, mask=mask)
        corr = real(numer / sqrt(aden * bden))
    end function array_correlation

end module simple_halfmap_enhancer_tester
