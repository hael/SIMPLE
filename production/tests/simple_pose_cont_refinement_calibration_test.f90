module pose_cont_refinement_calibration_test
use pose_cont_refinement_calibration_helpers
use pose_cont_refinement_test_helpers, only: assert_true
use simple_defs, only: dp, sp, KBALPHA, KBWINSZ
use simple_cartesian_pose_refiner, only: cartesian_pose_refiner
use simple_kbinterpol, only: kbinterpol
implicit none
private
public :: run_tolerance_calibration

contains

!> Calibrate every tolerance family without constructing an acceptance fixture.
subroutine run_tolerance_calibration()
    type(calibration_fixture) :: first(2), repeated(2)
    type(tolerance_record) :: records(N_TOLERANCE_FAMILIES), repeated_records(N_TOLERANCE_FAMILIES)
    character(len=:), allocatable :: evidence_directory
    integer :: fixture_index

    evidence_directory = required_argument('evidence_dir')
    call initialize_tolerances(records)
    call initialize_tolerances(repeated_records)
    do fixture_index = 1, size(CALIBRATION_BOXES)
        call build_calibration_fixture(CALIBRATION_BOXES(fixture_index),first(fixture_index))
        call build_calibration_fixture(CALIBRATION_BOXES(fixture_index),repeated(fixture_index))
        call assert_true(fixtures_are_identical(first(fixture_index),repeated(fixture_index)), &
            &'calibration fixture is not deterministic')
        call collect_fixture_observations(first(fixture_index),records)
        call collect_fixture_observations(repeated(fixture_index),repeated_records)
    enddo
    call finalize_tolerances(records)
    call finalize_tolerances(repeated_records)
    call require_repeatable_records(records,repeated_records)
    call check_combined_metric()
    call assert_true(.not. any(CALIBRATION_BOXES(1) == ACCEPTANCE_BOXES) .and. &
        &.not. any(CALIBRATION_BOXES(2) == ACCEPTANCE_BOXES), &
        &'calibration and acceptance box sets overlap')
    call write_calibration_artifacts(evidence_directory,records,first)

    write(*,'(a)') 'POSE_CONT_REFINEMENT_CALIBRATION boxes sampled: 8 12'
    write(*,'(a)') 'POSE_CONT_REFINEMENT_CALIBRATION acceptance boxes not sampled: 10 16'
    write(*,'(a,es14.6)') 'POSE_CONT_REFINEMENT_CALIBRATION safety factor: ',CALIBRATION_SAFETY_FACTOR
    write(*,'(a)') 'POSE_CONT_REFINEMENT_TOLERANCE_CALIBRATION: PASS'
end subroutine run_tolerance_calibration

!> Verify the componentwise absolute-plus-relative comparison boundary.
subroutine check_combined_metric()
    real(dp), parameter :: absolute_tolerance = 1.e-4_dp, relative_tolerance = 2.e-3_dp

    call assert_true(combined_real_passes([1._dp,0._dp],[1.001_dp,5.e-5_dp], &
        &absolute_tolerance,relative_tolerance), &
        'combined real metric rejected an in-bound vector')
    call assert_true(.not. combined_real_passes([1._dp],[1.01_dp],absolute_tolerance,relative_tolerance), &
        &'combined real metric accepted an out-of-bound component')
    call assert_true(combined_complex_passes([cmplx(1._dp,1._dp,kind=dp)], &
        &[cmplx(1.001_dp,1._dp,kind=dp)],absolute_tolerance,relative_tolerance), &
        &'combined complex metric rejected an in-bound value')
end subroutine check_combined_metric

!> Add deterministic numerical probes for all predeclared tolerance families.
subroutine collect_fixture_observations(fixture,records)
    type(calibration_fixture), intent(in) :: fixture
    type(tolerance_record), intent(inout) :: records(N_TOLERANCE_FAMILIES)
    type(cartesian_pose_refiner) :: workspace
    complex(dp) :: analytic_values(4), brute_values(4), direct_values(4), executed_values(4)
    complex(dp) :: slow_values(4), fast_values(4), projected_values(4), projection_oracle(4)
    complex, allocatable :: zero_plane(:,:), executed_plane(:,:)
    real, allocatable :: physical_volume(:,:,:)
    real(dp) :: exact_values(5), rounded_real(5), derivative_exact(5), derivative_fd(5)
    real(dp) :: lm_matrix(5,5), lm_rhs(5), lm_solution(5), lm_residual(5)
    real(dp) :: identity_rotation(3,3), ignored_objective, volume_sum_dp, volume_sum_sp
    integer, parameter :: h_modes(4) = [0,1,-2,1]
    integer, parameter :: k_modes(4) = [0,-1,1,2]
    integer :: axis, lims2(2,2), mode

    volume_sum_dp = sum(fixture%volume)
    volume_sum_sp = real(sum(real(fixture%volume,sp)),dp)
    call observe_real_comparison(records(1),[volume_sum_sp],[volume_sum_dp])

    do axis = 1, 5
        exact_values(axis) = fixture%nonstationary_pose(axis)+0.03_dp*real(axis,dp)
        rounded_real(axis) = real(real(exact_values(axis),sp),dp)
    enddo
    lm_matrix = 0._dp
    do axis = 1, 5
        lm_matrix(axis,axis) = 1._dp+0.2_dp*real(axis,dp)
    enddo
    lm_rhs = matmul(lm_matrix,exact_values)
    lm_solution = rounded_real
    lm_residual = matmul(lm_matrix,lm_solution)-lm_rhs
    call observe_real_comparison(records(2),lm_solution,exact_values)
    call observe_real_comparison(records(2),lm_residual,[0._dp,0._dp,0._dp,0._dp,0._dp])

    do axis = 1, 5
        derivative_exact(axis) = cos(fixture%nonstationary_pose(axis))+ &
            &0.2_dp*fixture%nonstationary_pose(axis)
        derivative_fd(axis) = centred_probe_derivative(fixture%nonstationary_pose(axis))
    enddo
    call observe_real_comparison(records(3),derivative_fd,derivative_exact)

    allocate(physical_volume,source=real(fixture%volume,sp))
    call workspace%new(physical_volume)
    call workspace%set_shell_range([0,fixture%box/2])
    lims2 = workspace%get_lims2()
    allocate(zero_plane(lims2(1,1):lims2(1,2),lims2(2,1):lims2(2,2)),source=cmplx(0.,0.))
    allocate(executed_plane,mold=zero_plane)
    identity_rotation = 0._dp
    identity_rotation(1,1) = 1._dp
    identity_rotation(2,2) = 1._dp
    identity_rotation(3,3) = 1._dp
    call workspace%shift_residual(identity_rotation,[0._dp,0._dp],zero_plane,executed_plane,ignored_objective)
    do mode = 1, 4
        analytic_values(mode) = analytic_point_dft(fixture%box,mode)
        brute_values(mode) = brute_point_dft(fixture%box,mode)
        direct_values(mode) = direct_fourier_plane_probe(fixture,mode)
        executed_values(mode) = cmplx(executed_plane(h_modes(mode),k_modes(mode)),kind=dp)
        call kb_gather_probe(fixture,mode,slow_values(mode),fast_values(mode))
        projection_oracle(mode) = direct_fourier_plane_probe(fixture,mode)
        projected_values(mode) = identity_projection_probe(fixture,mode)
    enddo
    call observe_complex_comparison(records(4),brute_values,analytic_values)
    call observe_complex_comparison(records(5),executed_values,direct_values)
    call observe_complex_comparison(records(6),fast_values,slow_values)
    call observe_complex_comparison(records(7),projected_values,projection_oracle)
    deallocate(physical_volume,zero_plane,executed_plane)
    call workspace%kill
end subroutine collect_fixture_observations

!> Return the finite-difference derivative of the independent scalar probe.
pure real(dp) function centred_probe_derivative(value) result(derivative)
    real(dp), intent(in) :: value
    real(dp), parameter :: step = 2.e-5_dp
    derivative = (probe_value(value+step)-probe_value(value-step))/(2._dp*step)
end function centred_probe_derivative

!> Return a smooth independent scalar probe for derivative calibration.
pure real(dp) function probe_value(value) result(probe)
    real(dp), intent(in) :: value
    probe = sin(value)+0.1_dp*value*value
end function probe_value

!> Compare the executed normalized fast stencil with an independent slow traversal.
subroutine kb_gather_probe(fixture,mode,slow_value,fast_value)
    type(calibration_fixture), intent(in) :: fixture
    integer, intent(in) :: mode
    complex(dp), intent(out) :: slow_value, fast_value
    type(kbinterpol) :: kbwin
    complex(dp) :: samples(3,3,3)
    real(sp) :: fast_weights(3,3,3), loc(3), one_dimensional(3,3), weight_sum
    integer :: coordinate, i, j, k, stencil_origin, tap

    kbwin = kbinterpol(KBWINSZ,KBALPHA)
    loc = real([0.13_dp*mode,-0.17_dp*mode,0.09_dp*(mode-2)],sp)
    do k = 1, 3
        do j = 1, 3
            do i = 1, 3
                samples(i,j,k) = cmplx(0.03_dp*fixture%fixture_id+real(2*i-j+3*k,dp), &
                    &real(i+4*j-2*k,dp)/real(fixture%box,dp),kind=dp)
            enddo
        enddo
    enddo
    call kbwin%apod_mat_3d_fast(loc,1,3,fast_weights)
    fast_value = sum(samples*cmplx(fast_weights,kind=dp))
    do coordinate = 1, 3
        stencil_origin = nint(loc(coordinate))-1
        do tap = 1, 3
            one_dimensional(tap,coordinate) = kbwin%apod_fast( &
                &real(stencil_origin+tap-1,sp)-loc(coordinate))
        enddo
        weight_sum = sum(one_dimensional(:,coordinate))
        one_dimensional(:,coordinate) = one_dimensional(:,coordinate)/weight_sum
    enddo
    slow_value = cmplx(0._dp,0._dp,kind=dp)
    do k = 1, 3
        do j = 1, 3
            do i = 1, 3
                slow_value = slow_value+samples(i,j,k)*real(one_dimensional(i,1)* &
                    &one_dimensional(j,2)*one_dimensional(k,3),dp)
            enddo
        enddo
    enddo
end subroutine kb_gather_probe

!> Evaluate the direct three-dimensional DFT on the identity projection plane.
function direct_fourier_plane_probe(fixture,mode) result(value)
    type(calibration_fixture), intent(in) :: fixture
    integer, intent(in) :: mode
    complex(dp) :: value
    real(dp) :: angle, origin, twopi
    integer, parameter :: h_modes(4) = [0,1,-2,1]
    integer, parameter :: k_modes(4) = [0,-1,1,2]
    integer :: h, k, i, j, iz

    h = h_modes(mode)
    k = k_modes(mode)
    twopi = 2._dp*acos(-1._dp)
    origin = real(fixture%box/2+1,dp)
    value = cmplx(0._dp,0._dp,kind=dp)
    do iz = 1, fixture%box
        do j = 1, fixture%box
            do i = 1, fixture%box
                angle = -twopi*(real(h,dp)*(real(i,dp)-origin)+ &
                    &real(k,dp)*(real(j,dp)-origin))/real(fixture%box,dp)
                value = value+fixture%volume(i,j,iz)*cmplx(cos(angle),sin(angle),kind=dp)
            enddo
        enddo
    enddo
    value = value/real(fixture%box**3,dp)
end function direct_fourier_plane_probe

!> Evaluate the identity line-sum projection followed by a normalized 2-D DFT.
function identity_projection_probe(fixture,mode) result(value)
    type(calibration_fixture), intent(in) :: fixture
    integer, intent(in) :: mode
    complex(dp) :: value
    real(dp) :: angle, origin, line_sum, twopi
    integer, parameter :: h_modes(4) = [0,1,-2,1]
    integer, parameter :: k_modes(4) = [0,-1,1,2]
    integer :: h, k, i, j

    h = h_modes(mode)
    k = k_modes(mode)
    twopi = 2._dp*acos(-1._dp)
    origin = real(fixture%box/2+1,dp)
    value = cmplx(0._dp,0._dp,kind=dp)
    do j = 1, fixture%box
        do i = 1, fixture%box
            line_sum = sum(fixture%volume(i,j,:))
            angle = -twopi*(real(h,dp)*(real(i,dp)-origin)+ &
                &real(k,dp)*(real(j,dp)-origin))/real(fixture%box,dp)
            value = value+line_sum*cmplx(cos(angle),sin(angle),kind=dp)
        enddo
    enddo
    value = value/real(fixture%box**3,dp)
end function identity_projection_probe

!> Evaluate the closed-form DFT of a non-collinear unequal-amplitude point set.
function analytic_point_dft(box,mode) result(value)
    integer, intent(in) :: box, mode
    complex(dp) :: value
    integer, parameter :: points(3,4) = reshape([0,0,0, 1,-2,1, -2,1,-1, 2,2,-2],[3,4])
    integer, parameter :: h_modes(4) = [0,1,-2,1]
    integer, parameter :: k_modes(4) = [0,-1,1,2]
    integer, parameter :: l_modes(4) = [0,1,1,-1]
    real(dp), parameter :: amplitudes(4) = [1._dp,0.61_dp,-0.37_dp,0.23_dp]
    real(dp) :: angle, twopi
    integer :: point

    twopi = 2._dp*acos(-1._dp)
    value = cmplx(0._dp,0._dp,kind=dp)
    do point = 1, 4
        angle = -twopi*real(h_modes(mode)*points(1,point)+k_modes(mode)*points(2,point)+ &
            &l_modes(mode)*points(3,point),dp)/real(box,dp)
        value = value+amplitudes(point)*cmplx(cos(angle),sin(angle),kind=dp)
    enddo
    value = value/real(box**3,dp)
end function analytic_point_dft

!> Evaluate the same point set through a full voxel-loop brute-force DFT.
function brute_point_dft(box,mode) result(value)
    integer, intent(in) :: box, mode
    complex(dp) :: value
    integer, parameter :: points(3,4) = reshape([0,0,0, 1,-2,1, -2,1,-1, 2,2,-2],[3,4])
    integer, parameter :: h_modes(4) = [0,1,-2,1]
    integer, parameter :: k_modes(4) = [0,-1,1,2]
    integer, parameter :: l_modes(4) = [0,1,1,-1]
    real(dp), parameter :: amplitudes(4) = [1._dp,0.61_dp,-0.37_dp,0.23_dp]
    real(dp) :: angle, sample, twopi
    integer :: h, k, l, point, x, y, z

    h = h_modes(mode)
    k = k_modes(mode)
    l = l_modes(mode)
    twopi = 2._dp*acos(-1._dp)
    value = cmplx(0._dp,0._dp,kind=dp)
    do z = -box/2, box/2-1
        do y = -box/2, box/2-1
            do x = -box/2, box/2-1
                sample = 0._dp
                do point = 1, 4
                    if( all([x,y,z] == points(:,point)) ) sample = sample+amplitudes(point)
                enddo
                angle = -twopi*real(h*x+k*y+l*z,dp)/real(box,dp)
                value = value+sample*cmplx(cos(angle),sin(angle),kind=dp)
            enddo
        enddo
    enddo
    value = value/real(box**3,dp)
end function brute_point_dft

!> Require bitwise-identical calibration records from the repeated run.
subroutine require_repeatable_records(first,repeated)
    type(tolerance_record), intent(in) :: first(N_TOLERANCE_FAMILIES)
    type(tolerance_record), intent(in) :: repeated(N_TOLERANCE_FAMILIES)
    integer :: family

    do family = 1, N_TOLERANCE_FAMILIES
        call assert_true(first(family)%name == repeated(family)%name, &
            &'repeated calibration changed a family name')
        call assert_true(first(family)%observations == repeated(family)%observations, &
            &'repeated calibration changed an observation count')
        call assert_true(first(family)%maximum_absolute_error == &
            &repeated(family)%maximum_absolute_error, &
            'repeated calibration changed an absolute maximum')
        call assert_true(first(family)%maximum_scaled_relative_error == &
            &repeated(family)%maximum_scaled_relative_error, &
            &'repeated calibration changed a relative maximum')
        call assert_true(all(first(family)%absolute_errors == repeated(family)%absolute_errors) .and. &
            &all(first(family)%scaled_relative_errors == repeated(family)%scaled_relative_errors), &
            &'repeated calibration changed a raw error table')
        call assert_true(first(family)%absolute_tolerance == repeated(family)%absolute_tolerance .and. &
            &first(family)%relative_tolerance == repeated(family)%relative_tolerance, &
            &'repeated calibration changed a derived tolerance')
    enddo
end subroutine require_repeatable_records

!> Read one required key=value command argument.
function required_argument(key) result(value)
    character(len=*), intent(in) :: key
    character(len=:), allocatable :: value
    character(len=4096) :: argument
    integer :: argument_index, argument_status, separator

    value = ''
    do argument_index = 1, command_argument_count()
        call get_command_argument(argument_index,argument,status=argument_status)
        if( argument_status /= 0 ) error stop 'could not read calibration argument'
        separator = index(argument,'=')
        if( separator <= 0 ) cycle
        if( trim(argument(:separator-1)) /= key ) cycle
        value = trim(argument(separator+1:))
    enddo
    if( len(value) == 0 ) error stop 'tolerance calibration requires evidence_dir=<existing-directory>'
end function required_argument

end module pose_cont_refinement_calibration_test
