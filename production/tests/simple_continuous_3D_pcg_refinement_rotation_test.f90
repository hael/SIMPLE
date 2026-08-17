module continuous_3D_pcg_refinement_rotation_test
use continuous_3D_pcg_refinement_test_helpers, only: assert_true, build_truth_volume, TRUTH_VOLUME_BOX
use ieee_arithmetic, only: ieee_is_finite
use simple_defs, only: dp, sp
use simple_core_module_api, only: euler2m
use simple_reconstructor_pcg, only: pcg_fourier_workspace, reconstructor_pcg, right_increment_rotation
use simple_type_defs, only: ctfparams, CTFFLAG_YES
implicit none
private
public :: run_rotation_gradient

integer, parameter :: N_FD_STEPS = 3
integer, parameter :: MAX_FD_SHELL = 4
real(dp), parameter :: FD_STEPS(N_FD_STEPS) = [2.e-4_dp,1.e-4_dp,5.e-5_dp]
real(dp), parameter :: ROTATION_JVP_TOL = 8.e-3_dp
real(dp), parameter :: OBJECTIVE_TOL = 8.e-3_dp
real(dp), parameter :: ROTATION_UPDATE_TOL = 2.e-12_dp

contains

!> Compare analytic right-tangent rotation derivatives with independent finite differences.
subroutine run_rotation_gradient()
    type(pcg_fourier_workspace) :: workspace
    type(reconstructor_pcg) :: pcgop
    type(ctfparams) :: ctfparms
    real, allocatable :: phantom(:,:,:), sig2(:)
    complex, allocatable :: observed(:,:), zero_plane(:,:), residual_minus(:,:), residual_plus(:,:)
    complex, allocatable :: jv(:,:), transfer(:,:), weighted_observed(:,:)
    real(dp) :: rotmat(3,3), true_rotmat(3,3), minus_rotmat(3,3), plus_rotmat(3,3), updated_rotmat(3,3)
    real(dp) :: independent_rotmat(3,3), shift(2), true_shift(2), rotation_direction(3), pose_direction(5)
    real(dp) :: gradient(5), ignored_objective, objective_minus, objective_plus, directional_derivative
    real(dp) :: jvp_errors(N_FD_STEPS), objective_errors(N_FD_STEPS), weighted_errors(N_FD_STEPS)
    real(dp) :: update_error, orthogonality_error, determinant_error
    real(dp) :: input_orthogonality_error, input_determinant_error
    integer :: fixed_cell_count, h, istep, lims2(2,2), minus_switches, plus_switches, measured_switches

    call build_truth_volume(phantom)
    call pcgop%new(TRUTH_VOLUME_BOX,1._sp)
    call pcgop%set_volume(phantom)
    call pcgop%begin_fourier_workspace(workspace)
    ! A small shell set permits a true whole-plane fixed-cell finite difference.
    call workspace%set_shell_range([2,MAX_FD_SHELL])
    lims2 = pcgop%get_lims2()
    allocate(observed(lims2(1,1):lims2(1,2),lims2(2,1):lims2(2,2)))
    allocate(zero_plane,mold=observed)
    allocate(residual_minus,mold=observed)
    allocate(residual_plus,mold=observed)
    allocate(jv,mold=observed)
    zero_plane = cmplx(0.,0.)

    rotmat = real(euler2m([17.,31.,23.]),dp)
    shift = [0.21_dp,-0.16_dp]
    true_rotmat = independent_right_increment(rotmat,[0.010_dp,-0.008_dp,0.006_dp])
    true_shift = [0.25_dp,-0.18_dp]
    rotation_direction = [0.36_dp,-0.48_dp,0.80_dp]
    pose_direction = [rotation_direction,0.30_dp,-0.20_dp]
    call workspace%shift_residual(true_rotmat,true_shift,zero_plane,observed,ignored_objective)
    call workspace%rotation_jvp(rotmat,shift,rotation_direction,jv)

    jvp_errors = huge(0._dp)
    objective_errors = huge(0._dp)
    fixed_cell_count = 0
    do istep = 1, N_FD_STEPS
        minus_rotmat = independent_right_increment(rotmat,-FD_STEPS(istep)*rotation_direction)
        plus_rotmat = independent_right_increment(rotmat,FD_STEPS(istep)*rotation_direction)
        minus_switches = workspace%count_stencil_switches(rotmat,minus_rotmat)
        plus_switches = workspace%count_stencil_switches(rotmat,plus_rotmat)
        if( minus_switches+plus_switches /= 0 ) cycle
        fixed_cell_count = fixed_cell_count+1
        call workspace%shift_residual(minus_rotmat,shift,observed,residual_minus,objective_minus)
        call workspace%shift_residual(plus_rotmat,shift,observed,residual_plus,objective_plus)
        jvp_errors(istep) = plane_relative_error( &
            &(residual_plus-residual_minus)/real(2._dp*FD_STEPS(istep),sp),jv,lims2)
        call workspace%pose_objective_gradient(rotmat,shift,observed,ignored_objective,gradient)
        minus_rotmat = independent_right_increment(rotmat,-FD_STEPS(istep)*pose_direction(1:3))
        plus_rotmat = independent_right_increment(rotmat,FD_STEPS(istep)*pose_direction(1:3))
        call workspace%pose_objective_gradient(minus_rotmat,shift-FD_STEPS(istep)*pose_direction(4:5), &
            &observed,objective_minus,gradient)
        call workspace%pose_objective_gradient(plus_rotmat,shift+FD_STEPS(istep)*pose_direction(4:5), &
            &observed,objective_plus,gradient)
        call workspace%pose_objective_gradient(rotmat,shift,observed,ignored_objective,gradient)
        directional_derivative = dot_product(gradient,pose_direction)
        objective_errors(istep) = abs((objective_plus-objective_minus)/(2._dp*FD_STEPS(istep))- &
            &directional_derivative)/max(1._dp,abs(directional_derivative))
    enddo
    call assert_true(fixed_cell_count > 0, &
        &'rotation finite differences could not find a fixed interpolation cell')
    call assert_true(minval(jvp_errors) < ROTATION_JVP_TOL, &
        &'right-rotation Jacobian disagrees with fixed-cell centred differences')
    call assert_true(minval(objective_errors) < OBJECTIVE_TOL, &
        &'joint pose objective gradient disagrees with centred differences')

    ! This independent exponential detects a left/right or skew-sign mismatch.
    input_orthogonality_error = sqrt(sum((matmul(transpose(rotmat),rotmat)-identity3())**2))
    input_determinant_error = abs(determinant3(rotmat)-1._dp)
    updated_rotmat = right_increment_rotation(rotmat,[0.013_dp,-0.017_dp,0.011_dp])
    independent_rotmat = independent_right_increment(rotmat,[0.013_dp,-0.017_dp,0.011_dp])
    update_error = sqrt(sum((updated_rotmat-independent_rotmat)**2))
    orthogonality_error = sqrt(sum((matmul(transpose(updated_rotmat),updated_rotmat)-identity3())**2))
    determinant_error = abs(determinant3(updated_rotmat)-1._dp)
    call assert_true(update_error < ROTATION_UPDATE_TOL, &
        &'production right-rotation update disagrees with the independent exponential')
    call assert_true(orthogonality_error <= input_orthogonality_error+ROTATION_UPDATE_TOL .and. &
        &determinant_error <= input_determinant_error+ROTATION_UPDATE_TOL, &
        &'right-rotation update increased the input SO(3) storage error')
    measured_switches = workspace%count_stencil_switches(rotmat, &
        &independent_right_increment(rotmat,[0.05_dp,-0.03_dp,0.04_dp]))
    call assert_true(measured_switches > 0, &
        &'rotation stencil-switch telemetry did not detect a deliberate cell crossing')

    ! Repeat the five-vector gradient check through production CTF/noise weighting.
    allocate(sig2(0:lims2(1,2)))
    sig2 = [(1.0+0.04*real(h),h=0,lims2(1,2))]
    allocate(transfer,mold=observed)
    allocate(weighted_observed,mold=observed)
    ctfparms%smpd = 1.0
    ctfparms%kv = 300.0
    ctfparms%cs = 2.7
    ctfparms%fraca = 0.1
    ctfparms%dfx = 1.4
    ctfparms%dfy = 1.55
    ctfparms%angast = 23.0
    ctfparms%phshift = 0.0
    ctfparms%ctfflag = CTFFLAG_YES
    transfer = pcgop%build_transfer(ctfparms,[0.,0.],sig2)
    weighted_observed = transfer*observed
    weighted_errors = huge(0._dp)
    do istep = 1, N_FD_STEPS
        minus_rotmat = independent_right_increment(rotmat,-FD_STEPS(istep)*pose_direction(1:3))
        plus_rotmat = independent_right_increment(rotmat,FD_STEPS(istep)*pose_direction(1:3))
        if( workspace%count_stencil_switches(rotmat,minus_rotmat)+ &
            &workspace%count_stencil_switches(rotmat,plus_rotmat) /= 0 ) cycle
        call workspace%pose_objective_gradient(rotmat,shift,weighted_observed, &
            &ignored_objective,gradient,transfer)
        call workspace%pose_objective_gradient(minus_rotmat,shift-FD_STEPS(istep)*pose_direction(4:5), &
            &weighted_observed,objective_minus,gradient,transfer)
        call workspace%pose_objective_gradient(plus_rotmat,shift+FD_STEPS(istep)*pose_direction(4:5), &
            &weighted_observed,objective_plus,gradient,transfer)
        call workspace%pose_objective_gradient(rotmat,shift,weighted_observed, &
            &ignored_objective,gradient,transfer)
        directional_derivative = dot_product(gradient,pose_direction)
        weighted_errors(istep) = abs((objective_plus-objective_minus)/(2._dp*FD_STEPS(istep))- &
            &directional_derivative)/max(1._dp,abs(directional_derivative))
    enddo
    call assert_true(minval(weighted_errors) < OBJECTIVE_TOL, &
        &'CTF/sigma-weighted joint pose gradient disagrees with centred differences')
    call assert_true(all(ieee_is_finite([minval(jvp_errors),minval(objective_errors), &
        &minval(weighted_errors),update_error,orthogonality_error,determinant_error])), &
        &'rotation derivative test produced non-finite diagnostics')

    write(*,'(a,3(es14.6,1x))') 'CONTINUOUS_3D_PCG_ROTATION Jv fixed-cell errors: ', jvp_errors
    write(*,'(a,3(es14.6,1x))') 'CONTINUOUS_3D_PCG_ROTATION pose-gradient errors: ', objective_errors
    write(*,'(a,3(es14.6,1x))') 'CONTINUOUS_3D_PCG_ROTATION weighted-gradient errors: ', weighted_errors
    write(*,'(a,3(es14.6,1x),2(i0,1x))') 'CONTINUOUS_3D_PCG_ROTATION update/orthogonal/determinant/fixed/switch: ', &
        &update_error,orthogonality_error,determinant_error,fixed_cell_count,measured_switches
    write(*,'(a,2(es14.6,1x))') 'CONTINUOUS_3D_PCG_ROTATION input orthogonal/determinant errors: ', &
        &input_orthogonality_error,input_determinant_error
    write(*,'(a)') 'CONTINUOUS_3D_PCG_ROTATION_GRADIENT: PASS'
    deallocate(phantom,observed,zero_plane,residual_minus,residual_plus,jv,transfer,weighted_observed,sig2)
    call workspace%kill
    call pcgop%kill
end subroutine run_rotation_gradient

!> Apply an independent Rodrigues implementation for derivative-sign verification.
pure function independent_right_increment(rotmat,omega) result(updated_rotmat)
    real(dp), intent(in) :: rotmat(3,3), omega(3)
    real(dp) :: updated_rotmat(3,3), theta, axis(3), cross_matrix(3,3), exponential(3,3)

    theta = sqrt(dot_product(omega,omega))
    if( theta <= tiny(1._dp) )then
        updated_rotmat = rotmat
        return
    endif
    axis = omega/theta
    cross_matrix = 0._dp
    cross_matrix(1,2) = -axis(3)
    cross_matrix(1,3) = axis(2)
    cross_matrix(2,1) = axis(3)
    cross_matrix(2,3) = -axis(1)
    cross_matrix(3,1) = -axis(2)
    cross_matrix(3,2) = axis(1)
    exponential = cos(theta)*identity3()+(1._dp-cos(theta))*outer3(axis)+sin(theta)*cross_matrix
    updated_rotmat = matmul(rotmat,exponential)
end function independent_right_increment

!> Return a 3-by-3 identity matrix for the independent exponential map.
pure function identity3() result(identity)
    real(dp) :: identity(3,3)
    identity = 0._dp
    identity(1,1) = 1._dp
    identity(2,2) = 1._dp
    identity(3,3) = 1._dp
end function identity3

!> Return the outer product used by the independent Rodrigues formula.
pure function outer3(vector) result(product)
    real(dp), intent(in) :: vector(3)
    real(dp) :: product(3,3)
    integer :: i, j
    do j = 1, 3
        do i = 1, 3
            product(i,j) = vector(i)*vector(j)
        enddo
    enddo
end function outer3

!> Return a 3-by-3 determinant for rotation-validity checks.
pure function determinant3(matrix) result(determinant)
    real(dp), intent(in) :: matrix(3,3)
    real(dp) :: determinant
    determinant = matrix(1,1)*(matrix(2,2)*matrix(3,3)-matrix(2,3)*matrix(3,2))- &
        &matrix(1,2)*(matrix(2,1)*matrix(3,3)-matrix(2,3)*matrix(3,1))+ &
        &matrix(1,3)*(matrix(2,1)*matrix(3,2)-matrix(2,2)*matrix(3,1))
end function determinant3

!> Return the relative L2 error between two packed Fourier planes.
function plane_relative_error(actual,expected,lims2) result(error)
    integer, intent(in) :: lims2(2,2)
    complex, intent(in) :: actual(lims2(1,1):lims2(1,2),lims2(2,1):lims2(2,2))
    complex, intent(in) :: expected(lims2(1,1):lims2(1,2),lims2(2,1):lims2(2,2))
    real(dp) :: error, numerator, denominator
    integer :: h, k
    numerator = 0._dp
    denominator = 0._dp
    do k = lims2(2,1), lims2(2,2)
        do h = lims2(1,1), lims2(1,2)
            if( h*h+k*k < 4 .or. h*h+k*k > MAX_FD_SHELL**2 ) cycle
            numerator = numerator+real(abs(actual(h,k)-expected(h,k)),dp)**2
            denominator = denominator+real(abs(expected(h,k)),dp)**2
        enddo
    enddo
    error = sqrt(numerator/max(denominator,epsilon(1._dp)))
end function plane_relative_error

end module continuous_3D_pcg_refinement_rotation_test
