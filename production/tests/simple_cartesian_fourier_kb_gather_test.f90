module cartesian_fourier_kb_gather_test
use cartesian_fourier_test_helpers, only: assert_true, build_truth_volume, TRUTH_VOLUME_BOX
use ieee_arithmetic, only: ieee_is_finite
use simple_defs, only: dp, sp, KBALPHA, KBWINSZ
use simple_cartesian_fourier, only: center_embed_real3d, gather_packed_window_grad
use simple_core_module_api, only: kbinterpol, cyci_1d
use simple_image, only: image
use simple_ori, only: ori
use simple_reconstructor_pcg, only: reconstructor_pcg
implicit none
private
public :: run_packed_gather_derivative

integer, parameter :: N_FD_POINTS = 3
real(sp), parameter :: FD_STEP = 5.e-3_sp
real(dp), parameter :: FORWARD_TOL = 2.e-6_dp
real(dp), parameter :: FRIEDEL_TOL = 2.e-5_dp
real(dp), parameter :: GRADIENT_TOL = 5.e-3_dp

contains

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

end module cartesian_fourier_kb_gather_test
