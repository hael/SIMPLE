!@descr: unit tests for the neutral Cartesian Fourier layer (simple_cartesian_fourier, simple_kbinterpol)
! The Kaiser-Bessel machinery under the Cartesian pose refiner and the PCG operator:
! the fast KB polynomial and its derivative against the ideal Bessel window, the
! normalised 3-D stencil derivatives and partition of unity, the stencil-switch
! discontinuity, the packed/Friedel gather derivative against finite differences on a
! Gaussian-blob phantom, and the parity of the extracted neutral operations (centred
! embed/crop, crop envelope, packed gathers, native plane extraction) with the
! pre-extraction implementations retained here as oracles.
module simple_cartesian_fourier_tester
use ieee_arithmetic,          only: ieee_is_finite
use simple_defs,              only: dp, sp, KBALPHA, KBWINSZ, KB_BETA_KB15_A2
use simple_cartesian_fourier, only: center_embed_real3d, center_crop_real3d, &
    &extract_native_fourier_plane, gather_packed_window, gather_packed_window_grad
use simple_core_module_api,   only: cyci_1d
use simple_gridding,          only: kb_stencil_envelope_1d, kb_stencil_centered_crop_inv_envelope_1d
use simple_image,             only: image
use simple_kbinterpol,        only: kbinterpol
use simple_ori,               only: ori
use simple_reconstructor_pcg, only: reconstructor_pcg
use simple_test_utils
implicit none
private
public :: run_all_cartesian_fourier_tests

integer, parameter :: TRUTH_VOLUME_BOX = 24
integer, parameter :: NBLOBS = 4
real,    parameter :: CTRS(3,NBLOBS) = reshape([&
    &-5.0, -3.0,  2.0, &
    & 4.0,  5.0, -3.0, &
    & 0.0, -6.0, -5.0, &
    & 3.0, -2.0,  6.0], [3,NBLOBS])
real,    parameter :: SIGMAS(NBLOBS) = [2.0, 2.5, 1.8, 2.2]
real,    parameter :: AMPS(NBLOBS)   = [1.0, 0.8, 0.6, 0.5]
! packed gather
integer, parameter :: N_FD_POINTS = 3
real(sp), parameter :: FD_STEP = 5.e-3_sp
real(dp), parameter :: FORWARD_TOL = 2.e-6_dp
real(dp), parameter :: FRIEDEL_TOL = 2.e-5_dp
real(dp), parameter :: GRADIENT_TOL = 5.e-3_dp
! KB derivative
integer, parameter :: N_IDEAL_SAMPLES = 281
integer, parameter :: N_VALUE_SAMPLES = 321
integer, parameter :: N_STENCIL_CASES = 3
real(sp), parameter :: SWITCH_OFFSET = 1.e-4_sp
real(sp), parameter :: FD_POINTS(7) = [-1.20_sp, -0.75_sp, -0.27_sp, 0._sp, 0.36_sp, 0.81_sp, 1.19_sp]
! neutral extraction (box 8 padded to 16)
integer, parameter :: BOX = 8
integer, parameter :: BOXPD = 16
integer, parameter :: WDIM = 3
real, parameter :: ENVELOPE_TOL = 2.e-6

contains

    subroutine run_all_cartesian_fourier_tests()
        write(*,'(A)') '**** running all Cartesian Fourier tests ****'
        write(*,'(A)') 'test_kb_derivative'
        call run_kb_derivative()
        write(*,'(A)') 'test_packed_gather_derivative'
        call run_packed_gather_derivative()
        write(*,'(A)') 'test_neutral_extract'
        call run_neutral_extract()
    end subroutine run_all_cartesian_fourier_tests

    !> the asymmetric Gaussian-blob phantom shared with the PCG testers
    subroutine build_truth_volume( volume )
        real, allocatable, intent(out) :: volume(:,:,:)
        real    :: ctr, dx, dy, dz
        integer :: b, i, j, k
        allocate(volume(TRUTH_VOLUME_BOX,TRUTH_VOLUME_BOX,TRUTH_VOLUME_BOX), source=0.)
        ctr = real(TRUTH_VOLUME_BOX)/2. + 0.5
        do k = 1, TRUTH_VOLUME_BOX
            do j = 1, TRUTH_VOLUME_BOX
                do i = 1, TRUTH_VOLUME_BOX
                    do b = 1, NBLOBS
                        dx = real(i) - ctr - CTRS(1,b); dy = real(j) - ctr - CTRS(2,b); dz = real(k) - ctr - CTRS(3,b)
                        volume(i,j,k) = volume(i,j,k) + AMPS(b) * exp(-(dx*dx + dy*dy + dz*dz) / (2. * SIGMAS(b)**2))
                    enddo
                enddo
            enddo
        enddo
    end subroutine build_truth_volume


!> Compare the packed/Friedel Fourier gather derivative with fixed-cell finite differences.
subroutine run_packed_gather_derivative()
    real(sp), parameter :: locs(3,N_FD_POINTS) = reshape([&
        & 2.17_sp, -1.23_sp,  0.31_sp, &
        &-3.37_sp,  2.18_sp, -0.29_sp, &
        & 0.21_sp,  1.34_sp, -2.27_sp], [3,N_FD_POINTS])
    type(reconstructor_pcg) :: pcgop
    type(kbinterpol) :: kbwin
    type(image) :: padded_image
    type(ori) :: orientation
    complex, allocatable :: plane(:,:), cmat(:,:,:)
    real, allocatable :: phantom(:,:,:)
    integer, allocatable :: wrap(:)
    complex :: derivative(3), derivative_minus(3), fd_derivative
    complex :: value, value_minus, value_plus
    real(sp) :: loc(3), margin(3), margin_minus(3), margin_plus(3)
    real(dp) :: max_fd_error, max_forward_error, max_friedel_deriv_error
    real(dp) :: max_friedel_value_error
    integer :: axis, h, icase, k, lims2(2,2), lims3(3,2), wlims(2), lo, hi, i

    call build_truth_volume(phantom)
    call pcgop%new(TRUTH_VOLUME_BOX, 1._sp)
    call pcgop%set_volume(phantom)
    call orientation%new(.false.)
    call orientation%set_euler([0._sp,0._sp,0._sp])
    lims2 = pcgop%get_lims2()
    allocate(plane(lims2(1,1):lims2(1,2),lims2(2,1):lims2(2,2)))
    call pcgop%forward_plane(orientation, plane)
    kbwin = kbinterpol(KBWINSZ,KBALPHA)
    call padded_image%new([2*TRUTH_VOLUME_BOX,2*TRUTH_VOLUME_BOX,2*TRUTH_VOLUME_BOX],1._sp)
    call padded_image%set_rmat(center_embed_real3d(phantom,2*TRUTH_VOLUME_BOX),.false.)
    call padded_image%fft()
    cmat = padded_image%get_cmat()
    lims3 = padded_image%loop_lims(3)
    wlims = lims3(2,:)
    lo = wlims(1)-ceiling(kbwin%get_winsz()-0.5)-1
    hi = wlims(2)+ceiling(kbwin%get_winsz()-0.5)+1
    allocate(wrap(lo:hi))
    do i = lo, hi
        wrap(i) = cyci_1d(wlims,i)
    enddo

    max_forward_error = 0._dp
    do icase = 1, 2
        if( icase == 1 )then
            h = 2;  k = -1
        else
            h = -3; k = 2
        endif
        loc = real(KBALPHA,sp) * real([h,k,0],sp)
        call sample_neutral(cmat,lo,wrap,kbwin,loc,value,derivative,margin)
        max_forward_error = max(max_forward_error, relative_complex_error(value, plane(h,k)))
    enddo
    call assert_true(max_forward_error < FORWARD_TOL, &
        &'workspace packed gather value disagrees with forward_plane')

    max_fd_error = 0._dp
    max_friedel_value_error = 0._dp
    max_friedel_deriv_error = 0._dp
    do icase = 1, N_FD_POINTS
        loc = locs(:,icase)
        call sample_neutral(cmat,lo,wrap,kbwin,loc,value,derivative,margin)
        call assert_true(all(margin > FD_STEP), &
            &'packed-gather finite-difference point is too close to a stencil switch')
        call assert_true(complex_is_finite(value) .and. complex_vector_is_finite(derivative), &
            &'packed gather returned a non-finite value or derivative')
        do axis = 1, 3
            loc         = locs(:,icase)
            loc(axis)   = loc(axis) + FD_STEP
            call sample_neutral(cmat,lo,wrap,kbwin,loc,value_plus,derivative_minus,margin_plus)
            loc(axis)   = loc(axis) - 2._sp * FD_STEP
            call sample_neutral(cmat,lo,wrap,kbwin,loc,value_minus,derivative_minus,margin_minus)
            fd_derivative = (value_plus - value_minus) / (2._sp * FD_STEP)
            max_fd_error = max(max_fd_error, relative_complex_error(fd_derivative, derivative(axis)))
        enddo

        ! F(-loc) = conjg(F(loc)) for a real-space volume.
        call sample_neutral(cmat,lo,wrap,kbwin,-locs(:,icase),value_minus,derivative_minus,margin_minus)
        max_friedel_value_error = max(max_friedel_value_error, &
            &relative_complex_error(value_minus, conjg(value)))
        ! grad F(-loc) = -conjg(grad F(loc)).
        do axis = 1, 3
            max_friedel_deriv_error = max(max_friedel_deriv_error, &
                &relative_complex_error(derivative_minus(axis), -conjg(derivative(axis))))
        enddo
    enddo
    call assert_true(max_fd_error < GRADIENT_TOL, &
        &'packed gather derivative disagrees with fixed-cell finite differences')
    call assert_true(max_friedel_value_error < FRIEDEL_TOL, &
        &'packed gather value violates Friedel conjugacy')
    call assert_true(max_friedel_deriv_error < FRIEDEL_TOL, &
        &'packed gather derivative violates differentiated Friedel conjugacy')

    write(*,'(a,4(es14.6,1x))') 'CARTESIAN_FOURIER_KB packed forward/fd/Friedel-value/Friedel-grad: ', &
        &max_forward_error, max_fd_error, max_friedel_value_error, max_friedel_deriv_error
    call padded_image%kill
    call pcgop%kill
    deallocate(plane,phantom,cmat,wrap)
end subroutine run_packed_gather_derivative

subroutine sample_neutral(cmat,wrap_lower,wrap,kbwin,loc,value,gradient,switch_margin)
    complex, intent(in) :: cmat(:,:,:)
    integer, intent(in) :: wrap_lower, wrap(wrap_lower:)
    type(kbinterpol), intent(in) :: kbwin
    real(sp), intent(in) :: loc(3)
    complex, intent(out) :: value, gradient(3)
    real(sp), intent(out) :: switch_margin(3)
    integer :: iwinsz, wdim, i0(3)
    real(sp), allocatable :: weights(:,:,:), derivatives(:,:,:,:)

    iwinsz = ceiling(kbwin%get_winsz()-0.5)
    wdim = 2*iwinsz+1
    allocate(weights(wdim,wdim,wdim),derivatives(wdim,wdim,wdim,3))
    call kbwin%apod_mat_3d_fast_grad(loc,iwinsz,wdim,i0,switch_margin,weights,derivatives)
    call gather_packed_window_grad(cmat,wrap_lower,wrap,i0,weights,derivatives,value,gradient)
    value = real(KBALPHA)**3*value
    gradient = real(KBALPHA)**3*gradient
    deallocate(weights,derivatives)
end subroutine sample_neutral

!> Return relative error between two complex Fourier samples.
pure real(dp) function relative_complex_error(actual, expected) result(error)
    complex, intent(in) :: actual, expected
    error = real(abs(actual-expected),dp) / max(1._dp,real(abs(actual),dp),real(abs(expected),dp))
end function relative_complex_error

!> Report whether all real and imaginary vector components are finite.
pure logical function complex_vector_is_finite(values) result(finite)
    complex, intent(in) :: values(:)
    finite = all(ieee_is_finite(real(values))) .and. all(ieee_is_finite(aimag(values)))
end function complex_vector_is_finite

!> Report whether both components of one complex value are finite.
pure logical function complex_is_finite(value) result(finite)
    complex, intent(in) :: value
    finite = ieee_is_finite(real(value)) .and. ieee_is_finite(aimag(value))
end function complex_is_finite


!> Verify the executed fast KB polynomial and normalized-stencil derivatives.
subroutine run_kb_derivative()
    type(kbinterpol) :: kbwin
    integer :: iwinsz, wdim

    kbwin  = kbinterpol(KBWINSZ, KBALPHA)
    iwinsz = ceiling(KBWINSZ - 0.5_sp)
    wdim   = 2 * iwinsz + 1
    call assert_true(wdim == 3, 'KB derivative test requires the standard three-tap stencil')

    call test_fast_polynomial(kbwin)
    call test_normalized_stencil(kbwin, iwinsz, wdim)
    call test_stencil_switch(kbwin, iwinsz, wdim)
end subroutine run_kb_derivative

!> Compare the fast polynomial derivative with fixed-cell finite differences.
subroutine test_fast_polynomial(kbwin)
    type(kbinterpol), intent(in) :: kbwin
    real(dp) :: ideal_deriv, ideal_value, max_ideal_deriv_error
    real(dp) :: max_ideal_value_error, max_fd_error, xdp
    real(sp) :: derivative, fd_derivative, outside_derivative, outside_value
    real(sp) :: value, value_minus, value_plus, x
    logical :: values_match
    integer :: i

    values_match = .true.
    do i = 0, N_VALUE_SAMPLES - 1
        x = -1.6_sp + 3.2_sp * real(i,sp) / real(N_VALUE_SAMPLES-1,sp)
        call kbwin%apod_fast_value_deriv(x, value, derivative)
        values_match = values_match .and. value == kbwin%apod_fast(x)
        call assert_true(ieee_is_finite(value) .and. ieee_is_finite(derivative), &
            &'joint fast KB evaluator returned a non-finite value')
    enddo
    call assert_true(values_match, 'joint fast KB value is not bit-identical to apod_fast')

    max_fd_error = 0._dp
    do i = 0, 6
        x = FD_POINTS(i+1)
        call kbwin%apod_fast_value_deriv(x, value, derivative)
        value_plus  = kbwin%apod_fast(x + FD_STEP)
        value_minus = kbwin%apod_fast(x - FD_STEP)
        fd_derivative = (value_plus - value_minus) / (2._sp * FD_STEP)
        max_fd_error = max(max_fd_error, abs(real(fd_derivative-derivative,dp)))
    enddo
    call assert_true(max_fd_error < 1.5e-3_dp, &
        &'fast KB analytic derivative disagrees with fixed-support finite differences')

    max_ideal_value_error = 0._dp
    max_ideal_deriv_error = 0._dp
    do i = 0, N_IDEAL_SAMPLES - 1
        xdp = -1.4_dp + 2.8_dp * real(i,dp) / real(N_IDEAL_SAMPLES-1,dp)
        x   = real(xdp,sp)
        call kbwin%apod_fast_value_deriv(x, value, derivative)
        call ideal_kb_reference(xdp, ideal_value, ideal_deriv)
        max_ideal_value_error = max(max_ideal_value_error, abs(real(value,dp)-ideal_value))
        max_ideal_deriv_error = max(max_ideal_deriv_error, abs(real(derivative,dp)-ideal_deriv))
    enddo
    call assert_true(max_ideal_value_error < 1.e-4_dp, &
        &'fast KB value differs excessively from the ideal Bessel reference')
    call assert_true(max_ideal_deriv_error < 2.e-4_dp, &
        &'fast KB derivative differs excessively from the ideal Bessel reference')

    call kbwin%apod_fast_value_deriv(KBWINSZ, value, derivative)
    call kbwin%apod_fast_value_deriv(KBWINSZ+SWITCH_OFFSET, outside_value, outside_derivative)
    call assert_true(value > 0._sp, 'fast KB endpoint should be nonzero inside hard support')
    call assert_true(outside_value == 0._sp .and. outside_derivative == 0._sp, &
        &'fast KB value and derivative should be zero outside hard support')

    write(*,'(a,3(es14.6,1x))') 'CARTESIAN_FOURIER_KB fast fd/ideal-value/ideal-derivative max abs: ', &
        &max_fd_error, max_ideal_value_error, max_ideal_deriv_error
    write(*,'(a,3(es14.6,1x))') 'CARTESIAN_FOURIER_KB support inside-value/inside-derivative/outside-value: ', &
        &real(value,dp), real(derivative,dp), real(outside_value,dp)
end subroutine test_fast_polynomial

!> Verify normalized 3-D stencil derivatives and the partition-of-unity derivative.
subroutine test_normalized_stencil(kbwin, iwinsz, wdim)
    type(kbinterpol), intent(in) :: kbwin
    integer,            intent(in) :: iwinsz, wdim
    real(sp), parameter :: locs(3,N_STENCIL_CASES) = reshape([&
        & 0.17_sp, -0.23_sp,  0.31_sp, &
        &-0.41_sp,  0.08_sp, -0.34_sp, &
        & 0.49_sp, -0.47_sp,  0.02_sp], [3,N_STENCIL_CASES])
    real(sp) :: dw(wdim,wdim,wdim,3), fd(wdim,wdim,wdim)
    real(sp) :: loc(3), loc_minus(3), loc_plus(3), margin(3)
    real(sp) :: w(wdim,wdim,wdim), w_minus(wdim,wdim,wdim)
    real(sp) :: w_plus(wdim,wdim,wdim), w_reference(wdim,wdim,wdim)
    real(dp) :: max_derivative_sum, max_fd_error, max_value_difference
    integer :: axis, icase, i0(3)

    max_derivative_sum = 0._dp
    max_fd_error        = 0._dp
    max_value_difference = 0._dp
    do icase = 1, N_STENCIL_CASES
        loc = locs(:,icase)
        call kbwin%apod_mat_3d_fast_grad(loc, iwinsz, wdim, i0, margin, w, dw)
        call kbwin%apod_mat_3d_fast(loc, iwinsz, wdim, w_reference)
        max_value_difference = max(max_value_difference, &
            &maxval(abs(real(w,dp)-real(w_reference,dp))))
        call assert_true(all(margin > FD_STEP), &
            &'finite-difference stencil case is too close to an nint switch')
        call assert_true(abs(sum(real(w,dp)) - 1._dp) < 5.e-6_dp, &
            &'normalized fast KB stencil does not sum to one')
        do axis = 1, 3
            max_derivative_sum = max(max_derivative_sum, abs(sum(real(dw(:,:,:,axis),dp))))
            loc_plus       = loc
            loc_minus      = loc
            loc_plus(axis) = loc_plus(axis) + FD_STEP
            loc_minus(axis)= loc_minus(axis) - FD_STEP
            call assert_true(all(nint(loc_plus)-iwinsz == i0) .and. &
                &all(nint(loc_minus)-iwinsz == i0), &
                &'stencil finite difference crossed an nint switch')
            call kbwin%apod_mat_3d_fast(loc_plus, iwinsz, wdim, w_plus)
            call kbwin%apod_mat_3d_fast(loc_minus, iwinsz, wdim, w_minus)
            fd = (w_plus - w_minus) / (2._sp * FD_STEP)
            max_fd_error = max(max_fd_error, &
                &maxval(abs(real(fd,dp)-real(dw(:,:,:,axis),dp))))
        enddo
    enddo
    call assert_true(max_value_difference < 2.e-6_dp, &
        &'gradient stencil value differs from apod_mat_3d_fast')
    call assert_true(max_derivative_sum < 2.e-5_dp, &
        &'normalized fast KB derivative stencil does not sum to zero')
    call assert_true(max_fd_error < 5.e-4_dp, &
        &'normalized fast KB derivative disagrees with fixed-cell finite differences')

    write(*,'(a,3(es14.6,1x))') 'CARTESIAN_FOURIER_KB stencil value/derivative-sum/fd max abs: ', &
        &max_value_difference, max_derivative_sum, max_fd_error
end subroutine test_normalized_stencil

!> Measure the discontinuity when nearest-grid stencil ownership changes.
subroutine test_stencil_switch(kbwin, iwinsz, wdim)
    type(kbinterpol), intent(in) :: kbwin
    integer,            intent(in) :: iwinsz, wdim
    real(sp) :: dw_left(wdim,wdim,wdim,3), dw_right(wdim,wdim,wdim,3)
    real(sp) :: loc_left(3), loc_right(3), margin_left(3), margin_right(3)
    real(sp) :: w_left(wdim,wdim,wdim), w_right(wdim,wdim,wdim)
    real(dp) :: deriv_left, deriv_right, jump, value_left, value_right
    integer :: i0_left(3), i0_right(3)

    loc_left  = [0.5_sp-SWITCH_OFFSET, 0.13_sp, -0.21_sp]
    loc_right = [0.5_sp+SWITCH_OFFSET, 0.13_sp, -0.21_sp]
    call kbwin%apod_mat_3d_fast_grad(loc_left, iwinsz, wdim, i0_left, &
        &margin_left, w_left, dw_left)
    call kbwin%apod_mat_3d_fast_grad(loc_right, iwinsz, wdim, i0_right, &
        &margin_right, w_right, dw_right)
    call assert_true(all(i0_right == i0_left + [1,0,0]), &
        &'the two switch probes did not select adjacent x stencils')
    call assert_true(abs(real(margin_left(1),dp)-real(SWITCH_OFFSET,dp)) < 2.e-6_dp .and. &
        &abs(real(margin_right(1),dp)-real(SWITCH_OFFSET,dp)) < 2.e-6_dp, &
        &'reported switch margin does not match the half-integer distance')

    value_left  = gather_synthetic_field(w_left, i0_left)
    value_right = gather_synthetic_field(w_right, i0_right)
    deriv_left  = gather_synthetic_field(dw_left(:,:,:,1), i0_left)
    deriv_right = gather_synthetic_field(dw_right(:,:,:,1), i0_right)
    jump = value_right - value_left
    call assert_true(all(ieee_is_finite([value_left,value_right,deriv_left,deriv_right,jump])), &
        &'stencil-switch measurement produced a non-finite value')
    call assert_true(abs(jump) > 1.e-6_dp, &
        &'asymmetric synthetic field did not expose the expected stencil-switch jump')

    write(*,'(a,3(i0,1x),a,3(i0,1x))') 'CARTESIAN_FOURIER_KB switch i0 left/right: ', &
        &i0_left, '/ ', i0_right
    write(*,'(a,5(es14.6,1x))') 'CARTESIAN_FOURIER_KB switch left/right/jump/dleft/dright: ', &
        &value_left, value_right, jump, deriv_left, deriv_right
end subroutine test_stencil_switch

!> Evaluate the ideal Bessel KB value and derivative only as an accuracy reference.
pure subroutine ideal_kb_reference(x, value, derivative)
    real(dp), intent(in)  :: x
    real(dp), intent(out) :: value, derivative
    integer, parameter :: REFERENCE_DEGREE = 40
    real(dp) :: coeff(0:REFERENCE_DEGREE), dpdu, p, scale, u, width
    integer :: k

    ! Independently generate a double-precision I0 power series rather than
    ! reusing the production's 15 stored single-precision coefficients.
    width    = 2._dp * real(KBWINSZ,dp)
    scale    = real(KB_BETA_KB15_A2,dp)**2 / 4._dp
    u        = 1._dp - (2._dp*x/width)**2
    coeff(0) = 1._dp
    do k = 1, REFERENCE_DEGREE
        coeff(k) = coeff(k-1) * scale / real(k*k,dp)
    enddo
    p    = coeff(REFERENCE_DEGREE)
    dpdu = 0._dp
    do k = REFERENCE_DEGREE-1, 0, -1
        dpdu = dpdu*u + p
        p    = p*u + coeff(k)
    enddo
    value      = p / width
    derivative = dpdu * (-8._dp*x/width**2) / width
end subroutine ideal_kb_reference

!> Gather the analytic synthetic field with one normalized KB stencil.
pure real(dp) function gather_synthetic_field(weights, i0) result(value)
    real(sp), intent(in) :: weights(:,:,:)
    integer,  intent(in) :: i0(3)
    real(dp) :: sample
    integer :: i, ix, iy, iz, j, k

    value = 0._dp
    do k = 1, size(weights,3)
        iz = i0(3) + k - 1
        do j = 1, size(weights,2)
            iy = i0(2) + j - 1
            do i = 1, size(weights,1)
                ix = i0(1) + i - 1
                sample = real(ix*ix,dp) + 0.3_dp*real(ix,dp) + &
                    &0.07_dp*real(iy,dp) - 0.11_dp*real(iz,dp)
                value = value + real(weights(i,j,k),dp) * sample
            enddo
        enddo
    enddo
end function gather_synthetic_field


!> Compare each Phase 2 neutral operation with the pre-extraction algorithm.
!! The hypothesis is behavioral parity after moving state-independent Cartesian
!! operations out of PCG. Exact index operations must match bit for bit; the
!! separable envelope may differ only by declared single-precision roundoff.
subroutine run_neutral_extract()
    call test_centered_embed_crop()
    call test_centered_crop_envelope()
    call test_packed_gathers()
    call test_native_plane_extraction()
    write (*, '(a)') 'CONTINUOUS_3D_NEUTRAL_EXTRACT: PASS'
end subroutine run_neutral_extract

!> Verify exact centered indices and the embed/crop adjoint pair.
subroutine test_centered_embed_crop()
    real :: native(BOX, BOX, BOX), arbitrary_padded(BOXPD, BOXPD, BOXPD)
    real, allocatable :: legacy_crop(:, :, :), legacy_embed(:, :, :)
    real, allocatable :: neutral_crop(:, :, :), neutral_embed(:, :, :)
    integer :: i, j, k

    do k = 1, BOX
        do j = 1, BOX
            do i = 1, BOX
                native(i, j, k) = real(100*i + 10*j + k)
            end do
        end do
    end do
    do k = 1, BOXPD
        do j = 1, BOXPD
            do i = 1, BOXPD
                arbitrary_padded(i, j, k) = real(10000*i + 100*j + k)
            end do
        end do
    end do

    neutral_embed = center_embed_real3d(native, BOXPD)
    legacy_embed = legacy_center_embed(native, BOXPD)
    call assert_true(all(neutral_embed == legacy_embed), &
        &'neutral centered embed changed an index or value')
    neutral_crop = center_crop_real3d(neutral_embed, BOX)
    call assert_true(all(neutral_crop == native), &
        &'neutral crop did not invert the centered embed')
    neutral_crop = center_crop_real3d(arbitrary_padded, BOX)
    legacy_crop = legacy_center_crop(arbitrary_padded, BOX)
    call assert_true(all(neutral_crop == legacy_crop), &
        &'neutral centered crop changed an index or value')

    write (*, '(a)') 'CONTINUOUS_3D_NEUTRAL centered embed/crop: exact'
end subroutine test_centered_embed_crop

!> Verify padded-period native factors and the unchanged PCG 3-D envelope.
subroutine test_centered_crop_envelope()
    type(kbinterpol) :: kbwin
    type(reconstructor_pcg) :: pcgop
    real, allocatable :: env_padded(:), inv1d(:), legacy_inv1d(:)
    real, allocatable :: legacy_env3(:, :, :), legacy_invenv3(:, :, :)
    real, allocatable :: pcg_env(:, :, :), pcg_invenv(:, :, :)
    real :: center_value, max_env_error, max_factor_error, max_invenv_error
    integer :: i, j, k, offset

    kbwin = kbinterpol(KBWINSZ, KBALPHA)
    call kb_stencil_envelope_1d(kbwin, BOXPD, env_padded)
    call kb_stencil_centered_crop_inv_envelope_1d(kbwin, BOXPD, BOX, inv1d)
    offset = (BOXPD - BOX)/2
    allocate (legacy_inv1d(BOX), source=0.0)
    do i = 1, BOX
        if (abs(env_padded(offset + i)) >= 1.e-8) legacy_inv1d(i) = 1.0/env_padded(offset + i)
    end do
    max_factor_error = maxval(abs(inv1d - legacy_inv1d))
    call assert_true(max_factor_error == 0.0, &
        &'centered-crop inverse-envelope factors changed the old calculation')

    allocate (legacy_env3(BOX, BOX, BOX), legacy_invenv3(BOX, BOX, BOX))
    ! The 3-D KB envelope is the separable product E(i,j,k)=e(i)e(j)e(k).
    do k = 1, BOX
        do j = 1, BOX
            do i = 1, BOX
                legacy_env3(i, j, k) = env_padded(offset + i)*env_padded(offset + j)*env_padded(offset + k)
            end do
        end do
    end do
    center_value = legacy_env3(BOX/2 + 1, BOX/2 + 1, BOX/2 + 1)
    if (abs(center_value) < 1.e-8) center_value = 1.0
    legacy_env3 = legacy_env3/center_value
    legacy_invenv3 = 1.0
    ! Guard the reciprocal so unsupported envelope values cannot amplify noise.
    where (abs(legacy_env3) < 1.e-8)
        legacy_invenv3 = 0.0
    elsewhere
        legacy_invenv3 = 1.0/legacy_env3
    end where

    call pcgop%new(BOX, 1.0)
    pcg_env = pcgop%get_env()
    pcg_invenv = pcgop%get_invenv()
    max_env_error = maxval(abs(pcg_env - legacy_env3))
    max_invenv_error = maxval(abs(pcg_invenv - legacy_invenv3)/max(1.0, abs(legacy_invenv3)))
    call assert_true(max_env_error <= ENVELOPE_TOL, &
        &'PCG envelope changed after the neutral extraction')
    call assert_true(max_invenv_error <= ENVELOPE_TOL, &
        &'PCG inverse envelope changed beyond single-precision roundoff')

    write (*, '(a,3(es14.6,1x))') 'CONTINUOUS_3D_NEUTRAL factor/env-abs/invenv-rel max: ', &
        &max_factor_error, max_env_error, max_invenv_error
    call pcgop%kill
end subroutine test_centered_crop_envelope

!> Verify the 27-tap order, packed/Friedel map, and negative wrap-table bound.
subroutine test_packed_gathers()
    complex :: cmat(BOXPD/2 + 1, BOXPD, BOXPD)
    complex :: legacy_gradient(3), legacy_value, neutral_gradient(3), neutral_value
    complex :: neutral_value_only
    real :: dw(WDIM, WDIM, WDIM, 3), w(WDIM, WDIM, WDIM)
    integer :: i, j, k, axis, i0(3), wrap(-BOXPD - 2:BOXPD + 2)

    do k = 1, BOXPD
        do j = 1, BOXPD
            do i = 1, BOXPD/2 + 1
                cmat(i, j, k) = cmplx(real(100*i + 3*j + k), real(-2*i + j - 4*k))
            end do
        end do
    end do
    do i = lbound(wrap, 1), ubound(wrap, 1)
        wrap(i) = modulo(i + BOXPD/2, BOXPD) - BOXPD/2
    end do
    do k = 1, WDIM
        do j = 1, WDIM
            do i = 1, WDIM
                w(i, j, k) = real(i + 2*j + 3*k)/324.0
                do axis = 1, 3
                    dw(i, j, k, axis) = real((i - axis)*(j + axis) - k)/1000.0
                end do
            end do
        end do
    end do
    i0 = [-9, 6, -2]

    ! Compare the identical 27-tap sum and its three spatial derivatives.
    neutral_value_only = gather_packed_window(cmat, lbound(wrap, 1), wrap, i0, w)
    call gather_packed_window_grad(cmat, lbound(wrap, 1), wrap, i0, w, dw, &
        &neutral_value, neutral_gradient)
    call legacy_packed_gather(cmat, lbound(wrap, 1), wrap, i0, w, dw, &
        &legacy_value, legacy_gradient)
    call assert_true(neutral_value_only == legacy_value .and. neutral_value == legacy_value, &
        &'neutral packed value gather changed the old traversal')
    call assert_true(all(neutral_gradient == legacy_gradient), &
        &'neutral packed gradient gather changed the old traversal')

    write (*, '(a,i0,a,i0)') 'CONTINUOUS_3D_NEUTRAL packed gather wrap bounds: ', &
        &lbound(wrap, 1), ':', ubound(wrap, 1)
end subroutine test_packed_gathers

!> Verify the neutral extraction and retained PCG wrapper against old code.
subroutine test_native_plane_extraction()
    type(image) :: img2d
    type(reconstructor_pcg) :: pcgop
    real :: pixels(BOX, BOX, 1)
    complex, allocatable :: legacy_plane(:, :), neutral_plane(:, :), wrapper_plane(:, :)
    integer :: i, j, lims2(2, 2), sqlp

    do j = 1, BOX
        do i = 1, BOX
            pixels(i, j, 1) = real(3*i*i + 5*j + 2*i*j)
        end do
    end do
    call img2d%new([BOX, BOX, 1], 1.0)
    call img2d%set_rmat(pixels, .false.)
    call img2d%fft()
    call pcgop%new(BOX, 1.0)
    lims2 = pcgop%get_lims2()
    sqlp = (BOX/2)**2
    neutral_plane = extract_native_fourier_plane(img2d, lims2, sqlp)
    wrapper_plane = pcgop%extract_native_plane(img2d)
    legacy_plane = legacy_native_plane(img2d, lims2, sqlp)
    call assert_true(all(neutral_plane == legacy_plane), &
        &'neutral native-plane extraction changed the old result')
    call assert_true(all(wrapper_plane == legacy_plane), &
        &'PCG native-plane wrapper changed the old result')

    write (*, '(a)') 'CONTINUOUS_3D_NEUTRAL native plane: exact'
    call pcgop%kill
    call img2d%kill
end subroutine test_native_plane_extraction

!> Pre-extraction centered embedding retained as a direct comparison oracle.
function legacy_center_embed(native, padded_box) result(padded)
    real, intent(in) :: native(:, :, :)
    integer, intent(in) :: padded_box
    real, allocatable :: padded(:, :, :)
    integer :: native_box, offset
    native_box = size(native, 1)
    offset = (padded_box - native_box)/2
    allocate (padded(padded_box, padded_box, padded_box), source=0.0)
    padded(offset + 1:offset + native_box, offset + 1:offset + native_box, &
        &offset + 1:offset + native_box) = native
end function legacy_center_embed

!> Pre-extraction centered crop retained as a direct comparison oracle.
function legacy_center_crop(padded, native_box) result(native)
    real, intent(in) :: padded(:, :, :)
    integer, intent(in) :: native_box
    real, allocatable :: native(:, :, :)
    integer :: offset
    offset = (size(padded, 1) - native_box)/2
    allocate (native(native_box, native_box, native_box), &
        &source=padded(offset + 1:offset + native_box, offset + 1:offset + native_box, &
        &offset + 1:offset + native_box))
end function legacy_center_crop

!> Pre-extraction packed gather retained as a direct comparison oracle.
pure subroutine legacy_packed_gather(cmat, wrap_lower, wrap, i0, w, dw, value, gradient)
    complex, intent(in) :: cmat(:, :, :)
    integer, intent(in) :: wrap_lower, wrap(wrap_lower:), i0(3)
    real, intent(in) :: w(:, :, :), dw(:, :, :, :)
    complex, intent(out) :: value, gradient(3)
    complex :: fcomp
    integer :: di, dj, dk, hh, kk, mm, ph, pk, pm, ny, nz
    ny = size(cmat, 2)
    nz = size(cmat, 3)
    value = cmplx(0., 0.)
    gradient = cmplx(0., 0.)
    do dk = 1, size(w, 3)
        mm = wrap(i0(3) + dk - 1)
        do dj = 1, size(w, 2)
            kk = wrap(i0(2) + dj - 1)
            do di = 1, size(w, 1)
                hh = wrap(i0(1) + di - 1)
                if (hh >= 0) then
                    ph = hh + 1
                    pk = kk + 1; if (kk < 0) pk = pk + ny
                    pm = mm + 1; if (mm < 0) pm = pm + nz
                    fcomp = cmat(ph, pk, pm)
                else
                    ph = -hh + 1
                    pk = -kk + 1; if (-kk < 0) pk = pk + ny
                    pm = -mm + 1; if (-mm < 0) pm = pm + nz
                    fcomp = conjg(cmat(ph, pk, pm))
                end if
                value = value + w(di, dj, dk)*fcomp
                gradient = gradient + dw(di, dj, dk, :)*fcomp
            end do
        end do
    end do
end subroutine legacy_packed_gather

!> Pre-extraction native-plane loop retained as a direct comparison oracle.
function legacy_native_plane(img2d, lims2, sqlp) result(plane)
    class(image), intent(in) :: img2d
    integer, intent(in) :: lims2(2, 2), sqlp
    complex :: plane(lims2(1, 1):lims2(1, 2), lims2(2, 1):lims2(2, 2))
    integer :: h, k, phys(3)
    plane = cmplx(0., 0.)
    do k = lims2(2, 1), lims2(2, 2)
        do h = lims2(1, 1), lims2(1, 2)
            if (h*h + k*k > sqlp) cycle
            phys = img2d%comp_addr_phys(h, k, 0)
            plane(h, k) = img2d%get_fcomp([h, k, 0], phys)
        end do
    end do
end function legacy_native_plane

end module simple_cartesian_fourier_tester
