module continuous_3D_pcg_refinement_kb_gather_test
use continuous_3D_pcg_refinement_test_helpers, only: assert_true, build_truth_volume, TRUTH_VOLUME_BOX
use ieee_arithmetic, only: ieee_is_finite
use simple_defs, only: dp, sp, KBALPHA
use simple_ori, only: ori
use simple_reconstructor_pcg, only: pcg_fourier_workspace, reconstructor_pcg
implicit none
private
public :: run_packed_gather_derivative

integer, parameter :: N_FD_POINTS = 3
real(sp), parameter :: FD_STEP = 5.e-3_sp
real(dp), parameter :: FORWARD_TOL = 2.e-6_dp
real(dp), parameter :: FRIEDEL_TOL = 2.e-5_dp
real(dp), parameter :: GRADIENT_TOL = 5.e-3_dp

contains

subroutine run_packed_gather_derivative()
    real(sp), parameter :: locs(3,N_FD_POINTS) = reshape([&
        & 2.17_sp, -1.23_sp,  0.31_sp, &
        &-3.37_sp,  2.18_sp, -0.29_sp, &
        & 0.21_sp,  1.34_sp, -2.27_sp], [3,N_FD_POINTS])
    type(pcg_fourier_workspace) :: workspace
    type(reconstructor_pcg) :: pcgop
    type(ori) :: orientation
    complex, allocatable :: plane(:,:)
    real, allocatable :: phantom(:,:,:)
    complex :: derivative(3), derivative_minus(3), fd_derivative
    complex :: value, value_minus, value_plus
    real(sp) :: loc(3), margin(3), margin_minus(3), margin_plus(3)
    real(dp) :: max_fd_error, max_forward_error, max_friedel_deriv_error
    real(dp) :: max_friedel_value_error
    integer :: axis, h, icase, k, lims2(2,2)

    call build_truth_volume(phantom)
    call pcgop%new(TRUTH_VOLUME_BOX, 1._sp)
    call pcgop%set_volume(phantom)
    call orientation%new(.false.)
    call orientation%set_euler([0._sp,0._sp,0._sp])
    lims2 = pcgop%get_lims2()
    allocate(plane(lims2(1,1):lims2(1,2),lims2(2,1):lims2(2,2)))
    call pcgop%forward_plane(orientation, plane)
    ! Snapshot the padded Fourier lattice once for every sample below.
    call pcgop%begin_fourier_workspace(workspace)

    max_forward_error = 0._dp
    do icase = 1, 2
        if( icase == 1 )then
            h = 2;  k = -1
        else
            h = -3; k = 2
        endif
        loc = real(KBALPHA,sp) * real([h,k,0],sp)
        call workspace%sample_with_grad(loc, value, derivative, margin)
        max_forward_error = max(max_forward_error, relative_complex_error(value, plane(h,k)))
    enddo
    call assert_true(max_forward_error < FORWARD_TOL, &
        &'workspace packed gather value disagrees with forward_plane')

    max_fd_error = 0._dp
    max_friedel_value_error = 0._dp
    max_friedel_deriv_error = 0._dp
    do icase = 1, N_FD_POINTS
        loc = locs(:,icase)
        call workspace%sample_with_grad(loc, value, derivative, margin)
        call assert_true(all(margin > FD_STEP), &
            &'packed-gather finite-difference point is too close to a stencil switch')
        call assert_true(complex_is_finite(value) .and. complex_vector_is_finite(derivative), &
            &'packed gather returned a non-finite value or derivative')
        do axis = 1, 3
            loc         = locs(:,icase)
            loc(axis)   = loc(axis) + FD_STEP
            call workspace%sample_with_grad(loc, value_plus, derivative_minus, margin_plus)
            loc(axis)   = loc(axis) - 2._sp * FD_STEP
            call workspace%sample_with_grad(loc, value_minus, derivative_minus, margin_minus)
            fd_derivative = (value_plus - value_minus) / (2._sp * FD_STEP)
            max_fd_error = max(max_fd_error, relative_complex_error(fd_derivative, derivative(axis)))
        enddo

        ! F(-loc) = conjg(F(loc)) for a real-space volume.
        call workspace%sample_with_grad(-locs(:,icase), value_minus, derivative_minus, margin_minus)
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

    write(*,'(a,4(es14.6,1x))') 'CONTINUOUS_3D_PCG_KB packed forward/fd/Friedel-value/Friedel-grad: ', &
        &max_forward_error, max_fd_error, max_friedel_value_error, max_friedel_deriv_error
    call workspace%kill
    call pcgop%kill
    deallocate(plane, phantom)
end subroutine run_packed_gather_derivative

pure real(dp) function relative_complex_error(actual, expected) result(error)
    complex, intent(in) :: actual, expected
    error = real(abs(actual-expected),dp) / max(1._dp,real(abs(actual),dp),real(abs(expected),dp))
end function relative_complex_error

pure logical function complex_vector_is_finite(values) result(finite)
    complex, intent(in) :: values(:)
    finite = all(ieee_is_finite(real(values))) .and. all(ieee_is_finite(aimag(values)))
end function complex_vector_is_finite

pure logical function complex_is_finite(value) result(finite)
    complex, intent(in) :: value
    finite = ieee_is_finite(real(value)) .and. ieee_is_finite(aimag(value))
end function complex_is_finite

end module continuous_3D_pcg_refinement_kb_gather_test
