module continuous_3D_pcg_refinement_kb_test
use continuous_3D_pcg_refinement_kb_gather_test, only: run_packed_gather_derivative
use continuous_3D_pcg_refinement_test_helpers, only: assert_true
use ieee_arithmetic, only: ieee_is_finite
use simple_defs, only: dp, sp, KBALPHA, KBWINSZ, KB_BETA_KB15_A2
use simple_kbinterpol, only: kbinterpol
implicit none
private
public :: run_kb_derivative

integer, parameter :: N_IDEAL_SAMPLES = 281
integer, parameter :: N_VALUE_SAMPLES = 321
integer, parameter :: N_STENCIL_CASES = 3
real(sp), parameter :: FD_STEP = 5.e-3_sp
real(sp), parameter :: SWITCH_OFFSET = 1.e-4_sp
real(sp), parameter :: FD_POINTS(7) = [-1.20_sp, -0.75_sp, -0.27_sp, 0._sp, 0.36_sp, 0.81_sp, 1.19_sp]

contains

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
    call run_packed_gather_derivative()
    write(*,'(a)') 'CONTINUOUS_3D_PCG_KB_DERIVATIVE: PASS'
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

    write(*,'(a,3(es14.6,1x))') 'CONTINUOUS_3D_PCG_KB fast fd/ideal-value/ideal-derivative max abs: ', &
        &max_fd_error, max_ideal_value_error, max_ideal_deriv_error
    write(*,'(a,3(es14.6,1x))') 'CONTINUOUS_3D_PCG_KB support inside-value/inside-derivative/outside-value: ', &
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

    write(*,'(a,3(es14.6,1x))') 'CONTINUOUS_3D_PCG_KB stencil value/derivative-sum/fd max abs: ', &
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

    write(*,'(a,3(i0,1x),a,3(i0,1x))') 'CONTINUOUS_3D_PCG_KB switch i0 left/right: ', &
        &i0_left, '/ ', i0_right
    write(*,'(a,5(es14.6,1x))') 'CONTINUOUS_3D_PCG_KB switch left/right/jump/dleft/dright: ', &
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

end module continuous_3D_pcg_refinement_kb_test
