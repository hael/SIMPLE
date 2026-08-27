! Acceptance test for the prepared Cartesian pose objective and normal equations.
module pose_cont_refinement_objective_normals_test
use pose_cont_refinement_calibration_helpers, only: calibration_fixture, ACCEPTANCE_BOXES, &
    &FROZEN_ABSOLUTE_TOLERANCES, FROZEN_RELATIVE_TOLERANCES, build_acceptance_fixture, &
    &combined_real_passes, combined_complex_passes
use pose_cont_refinement_test_helpers, only: assert_true
use simple_defs, only: dp, sp, PI, OSMPL_PAD_FAC
use simple_cartesian_pose_refiner, only: cartesian_pose_refiner, cartesian_pose_data, &
    &right_increment_rotation, POSE_DATA_VALID
use simple_type_defs, only: ctfparams, CTFFLAG_NO
implicit none
private
public :: run_objective_normals

integer, parameter :: ALGEBRAIC_TOLERANCE = 1
integer, parameter :: DERIVATIVE_TOLERANCE = 3
integer, parameter :: N_VARIANCE_PROFILES = 3
integer, parameter :: N_POSE_STATES = 2
character(len=16), parameter :: VARIANCE_NAMES(N_VARIANCE_PROFILES) = [character(len=16) :: &
    &'unit', 'constant', 'varying']
character(len=16), parameter :: POSE_NAMES(N_POSE_STATES) = [character(len=16) :: &
    &'exact', 'nonstationary']

contains

!> Compare the fused implementation with an independently accumulated scalar oracle.
subroutine run_objective_normals()
    type(calibration_fixture) :: fixture
    character(len=:), allocatable :: evidence_directory
    character(len=:), allocatable :: component_path, normal_path, summary_path
    integer :: box_index, component_unit, normal_unit, summary_unit

    evidence_directory = required_argument('evidence_dir')
    component_path = evidence_directory//'/residual_jacobian_components.tsv'
    normal_path = evidence_directory//'/normal_equation_components.tsv'
    summary_path = evidence_directory//'/objective_normal_summary.tsv'
    open(newunit=component_unit,file=component_path,status='replace',action='write')
    open(newunit=normal_unit,file=normal_path,status='replace',action='write')
    open(newunit=summary_unit,file=summary_path,status='replace',action='write')
    write(component_unit,'(a)') 'box'//achar(9)//'variance'//achar(9)//'pose'//achar(9)// &
        &'h'//achar(9)//'k'//achar(9)//'component'//achar(9)//'actual_real'//achar(9)// &
        &'actual_imag'//achar(9)//'expected_real'//achar(9)//'expected_imag'//achar(9)// &
        &'absolute_error'
    write(normal_unit,'(a)') 'box'//achar(9)//'variance'//achar(9)//'pose'//achar(9)// &
        &'component'//achar(9)//'actual'//achar(9)//'expected'//achar(9)//'absolute_error'
    write(summary_unit,'(a)') 'box'//achar(9)//'variance'//achar(9)//'pose'//achar(9)// &
        &'active_samples'//achar(9)//'actual_objective'//achar(9)//'expected_objective'// &
        &achar(9)//'max_residual_error'//achar(9)//'max_jacobian_error'//achar(9)// &
        &'max_gradient_error'//achar(9)//'max_hessian_error'//achar(9)//'hessian_asymmetry'// &
        &achar(9)//'minimum_switch_margin'

    do box_index = 1, size(ACCEPTANCE_BOXES)
        call build_acceptance_fixture(ACCEPTANCE_BOXES(box_index),fixture)
        call test_fixture(fixture,component_unit,normal_unit,summary_unit)
    enddo
    close(component_unit)
    close(normal_unit)
    close(summary_unit)

    write(*,'(a)') 'POSE_CONT_REFINEMENT_OBJECTIVE_NORMALS acceptance boxes: 10 16'
    write(*,'(a)') 'POSE_CONT_REFINEMENT_OBJECTIVE_NORMALS variance profiles: unit constant varying'
    write(*,'(a)') 'POSE_CONT_REFINEMENT_OBJECTIVE_NORMALS pose states: exact nonstationary'
    write(*,'(a)') 'POSE_CONT_REFINEMENT_OBJECTIVE_NORMALS: PASS'
end subroutine run_objective_normals

!> Exercise all variance and pose combinations for one untouched acceptance fixture.
subroutine test_fixture(fixture,component_unit,normal_unit,summary_unit)
    type(calibration_fixture), intent(in) :: fixture
    integer, intent(in) :: component_unit, normal_unit, summary_unit
    type(cartesian_pose_refiner) :: workspace
    complex, allocatable :: raw_observed(:,:)
    real, allocatable :: physical_volume(:,:,:), sigma2(:)
    real(dp) :: identity(3,3), truth_rotmat(3,3)
    integer :: pose_state, variance_profile

    identity = reshape([1._dp,0._dp,0._dp,0._dp,1._dp,0._dp,0._dp,0._dp,1._dp],[3,3])
    physical_volume = real(fixture%volume,sp)
    call workspace%new(physical_volume)
    call workspace%set_shell_range([0,fixture%box/2])
    truth_rotmat = right_increment_rotation(identity,fixture%exact_pose(1:3))
    call build_raw_observation(workspace,truth_rotmat,fixture%exact_pose(4:5),raw_observed)

    do variance_profile = 1, N_VARIANCE_PROFILES
        call select_variance(fixture,variance_profile,sigma2)
        do pose_state = 1, N_POSE_STATES
            call test_case(workspace,fixture,raw_observed,sigma2,variance_profile,pose_state, &
                &component_unit,normal_unit,summary_unit)
        enddo
    enddo
    call workspace%kill
end subroutine test_fixture

!> Construct the unweighted observation at the known exact pose.
subroutine build_raw_observation(workspace,rotmat,shift,raw_observed)
    type(cartesian_pose_refiner), intent(in) :: workspace
    real(dp), intent(in) :: rotmat(3,3), shift(2)
    complex, allocatable, intent(out) :: raw_observed(:,:)
    complex :: value, derivative(3), phase
    real :: margin(3)
    real(sp) :: loc(3)
    real(dp) :: argument
    integer :: box, h, k, lims2(2,2)

    lims2 = workspace%get_lims2()
    box = lims2(1,2)-lims2(1,1)
    allocate(raw_observed(lims2(1,1):lims2(1,2),lims2(2,1):lims2(2,2)),source=cmplx(0.,0.))
    do k = lims2(2,1), lims2(2,2)
        do h = lims2(1,1), lims2(1,2)
            if( h*h+k*k > (box/2)**2 ) cycle
            loc = real(OSMPL_PAD_FAC,sp)*real(matmul(real([h,k,0],dp),rotmat),sp)
            call workspace%sample_with_grad(loc,value,derivative,margin)
            argument = 2._dp*real(PI,dp)*(real(h,dp)*shift(1)+real(k,dp)*shift(2))/real(box,dp)
            phase = cmplx(cos(argument),sin(argument),kind=sp)
            raw_observed(h,k) = phase*value
        enddo
    enddo
end subroutine build_raw_observation

!> Run one weighted pose case and retain each compared component.
subroutine test_case(workspace,fixture,raw_observed,sigma2,variance_profile,pose_state, &
    &component_unit,normal_unit,summary_unit)
    type(cartesian_pose_refiner), intent(in) :: workspace
    type(calibration_fixture), intent(in) :: fixture
    complex, intent(in) :: raw_observed(-fixture%box/2:,-fixture%box/2:)
    real, intent(in) :: sigma2(0:)
    integer, intent(in) :: variance_profile, pose_state, component_unit, normal_unit, summary_unit
    type(cartesian_pose_data) :: data
    type(ctfparams) :: no_ctf
    complex(dp), allocatable :: actual_residual(:,:), actual_jacobian(:,:,:)
    complex(dp), allocatable :: expected_residual(:,:), expected_jacobian(:,:,:)
    real(dp) :: pose(5), rotmat(3,3), identity(3,3)
    real(dp) :: actual_objective, expected_objective
    real(dp) :: actual_gradient(5), expected_gradient(5)
    real(dp) :: actual_hessian(5,5), expected_hessian(5,5)
    real(dp) :: actual_margin, expected_margin
    real(dp) :: max_residual_error, max_jacobian_error, max_gradient_error, max_hessian_error
    real(dp) :: hessian_asymmetry
    integer :: active_samples, data_status

    identity = reshape([1._dp,0._dp,0._dp,0._dp,1._dp,0._dp,0._dp,0._dp,1._dp],[3,3])
    if( pose_state == 1 )then
        pose = fixture%exact_pose
    else
        pose = fixture%nonstationary_pose
    endif
    rotmat = right_increment_rotation(identity,pose(1:3))
    no_ctf%ctfflag = CTFFLAG_NO
    call workspace%prepare_particle(raw_observed,no_ctf,sigma2,[0,fixture%box/2],data,data_status)
    call assert_true(data_status == POSE_DATA_VALID .and. data%is_valid(), &
        &'objective-normal fixture did not produce valid prepared data')
    call workspace%prepared_normal_terms_test(rotmat,pose(4:5),data,actual_objective, &
        &actual_gradient,actual_hessian,actual_margin)
    call workspace%prepared_residual_jacobian_test(rotmat,pose(4:5),data,actual_residual, &
        &actual_jacobian,expected_margin)
    call independent_normal_terms(workspace,fixture%box,rotmat,pose(4:5),raw_observed,sigma2, &
        &expected_residual,expected_jacobian,expected_objective,expected_gradient, &
        &expected_hessian,expected_margin,active_samples)

    call compare_components(fixture%box,variance_profile,pose_state,actual_residual, &
        &expected_residual,actual_jacobian,expected_jacobian,component_unit, &
        &max_residual_error,max_jacobian_error)
    call compare_normal_terms(fixture%box,variance_profile,pose_state,actual_objective, &
        &expected_objective,actual_gradient,expected_gradient,actual_hessian,expected_hessian, &
        &normal_unit,max_gradient_error,max_hessian_error)
    hessian_asymmetry = maxval(abs(actual_hessian-transpose(actual_hessian)))
    call assert_true(combined_real_passes([actual_objective],[expected_objective], &
        &FROZEN_ABSOLUTE_TOLERANCES(ALGEBRAIC_TOLERANCE), &
        &FROZEN_RELATIVE_TOLERANCES(ALGEBRAIC_TOLERANCE)), &
        &'prepared objective differs from independent scalar accumulation')
    call assert_true(combined_real_passes(reshape(actual_hessian,[25]), &
        &reshape(expected_hessian,[25]),FROZEN_ABSOLUTE_TOLERANCES(ALGEBRAIC_TOLERANCE), &
        &FROZEN_RELATIVE_TOLERANCES(ALGEBRAIC_TOLERANCE)), &
        &'prepared Gauss-Newton matrix differs from independent accumulation')
    call assert_true(combined_real_passes(actual_gradient,expected_gradient, &
        &FROZEN_ABSOLUTE_TOLERANCES(ALGEBRAIC_TOLERANCE), &
        &FROZEN_RELATIVE_TOLERANCES(ALGEBRAIC_TOLERANCE)), &
        &'prepared JHr vector differs from independent accumulation')
    call assert_true(hessian_asymmetry <= FROZEN_ABSOLUTE_TOLERANCES(ALGEBRAIC_TOLERANCE), &
        &'prepared Gauss-Newton matrix is not symmetric')
    call write_summary(summary_unit,fixture%box,variance_profile,pose_state,active_samples, &
        &actual_objective,expected_objective,max_residual_error,max_jacobian_error, &
        &max_gradient_error,max_hessian_error,hessian_asymmetry,actual_margin)
end subroutine test_case

!> Recompute residuals, Jacobians, JHr, and JHJ without calling fused normal terms.
subroutine independent_normal_terms(workspace,box,rotmat,shift,raw_observed,sigma2,residual, &
    &jacobian,objective,gradient,hessian,min_switch_margin,active_samples)
    type(cartesian_pose_refiner), intent(in) :: workspace
    integer, intent(in) :: box
    real(dp), intent(in) :: rotmat(3,3), shift(2)
    complex, intent(in) :: raw_observed(-box/2:,-box/2:)
    real, intent(in) :: sigma2(0:)
    complex(dp), allocatable, intent(out) :: residual(:,:), jacobian(:,:,:)
    real(dp), intent(out) :: objective, gradient(5), hessian(5,5), min_switch_margin
    integer, intent(out) :: active_samples
    complex :: value, derivative(3), phase
    complex(dp) :: model, weighted_phase
    real :: margin(3), inverse_sigma
    real(sp) :: loc(3)
    real(dp) :: argument, dloc(3,3), frequency(2)
    integer :: axis, h, jaxis, k, lims2(2,2), shell

    lims2 = workspace%get_lims2()
    allocate(residual(lims2(1,1):lims2(1,2),lims2(2,1):lims2(2,2)), &
        &source=cmplx(0._dp,0._dp,kind=dp))
    allocate(jacobian(lims2(1,1):lims2(1,2),lims2(2,1):lims2(2,2),5), &
        &source=cmplx(0._dp,0._dp,kind=dp))
    objective = 0._dp
    gradient = 0._dp
    hessian = 0._dp
    min_switch_margin = huge(0._dp)
    active_samples = 0
    do k = lims2(2,1), lims2(2,2)
        do h = lims2(1,1), lims2(1,2)
            if( h*h+k*k > (box/2)**2 ) cycle
            shell = nint(sqrt(real(h*h+k*k)))
            inverse_sigma = 1.0/sqrt(sigma2(shell))
            loc = real(OSMPL_PAD_FAC,sp)*real(matmul(real([h,k,0],dp),rotmat),sp)
            call workspace%sample_with_grad(loc,value,derivative,margin)
            min_switch_margin = min(min_switch_margin,real(minval(margin),dp))
            dloc(:,1) = [0._dp,real(loc(3),dp),-real(loc(2),dp)]
            dloc(:,2) = [-real(loc(3),dp),0._dp,real(loc(1),dp)]
            dloc(:,3) = [real(loc(2),dp),-real(loc(1),dp),0._dp]
            argument = 2._dp*real(PI,dp)*(real(h,dp)*shift(1)+real(k,dp)*shift(2))/real(box,dp)
            phase = cmplx(cos(argument),sin(argument),kind=sp)
            weighted_phase = cmplx(phase,kind=dp)*real(inverse_sigma,dp)
            model = weighted_phase*cmplx(value,kind=dp)
            residual(h,k) = model-cmplx(raw_observed(h,k)*inverse_sigma,kind=dp)
            do axis = 1, 3
                jacobian(h,k,axis) = weighted_phase* &
                    &sum(cmplx(derivative,kind=dp)*dloc(:,axis))
            enddo
            frequency = 2._dp*real(PI,dp)*real([h,k],dp)/real(box,dp)
            jacobian(h,k,4:5) = cmplx(0._dp,frequency,kind=dp)*model
            objective = objective+0.5_dp*real(conjg(residual(h,k))*residual(h,k),dp)
            do axis = 1, 5
                gradient(axis) = gradient(axis)+ &
                    &real(conjg(jacobian(h,k,axis))*residual(h,k),dp)
                do jaxis = 1, 5
                    hessian(axis,jaxis) = hessian(axis,jaxis)+ &
                        &real(conjg(jacobian(h,k,axis))*jacobian(h,k,jaxis),dp)
                enddo
            enddo
            active_samples = active_samples+1
        enddo
    enddo
end subroutine independent_normal_terms

!> Check and retain all active complex residual and Jacobian components.
subroutine compare_components(box,variance_profile,pose_state,actual_residual,expected_residual, &
    &actual_jacobian,expected_jacobian,unit,max_residual_error,max_jacobian_error)
    integer, intent(in) :: box, variance_profile, pose_state, unit
    complex(dp), intent(in) :: actual_residual(-box/2:,-box/2:), expected_residual(-box/2:,-box/2:)
    complex(dp), intent(in) :: actual_jacobian(-box/2:,-box/2:,:)
    complex(dp), intent(in) :: expected_jacobian(-box/2:,-box/2:,:)
    real(dp), intent(out) :: max_residual_error, max_jacobian_error
    character(len=24) :: component
    integer :: axis, h, k

    max_residual_error = 0._dp
    max_jacobian_error = 0._dp
    do k = -box/2, box/2
        do h = -box/2, box/2
            if( h*h+k*k > (box/2)**2 ) cycle
            max_residual_error = max(max_residual_error,abs(actual_residual(h,k)-expected_residual(h,k)))
            call assert_true(combined_complex_passes([actual_residual(h,k)],[expected_residual(h,k)], &
                &FROZEN_ABSOLUTE_TOLERANCES(ALGEBRAIC_TOLERANCE), &
                &FROZEN_RELATIVE_TOLERANCES(ALGEBRAIC_TOLERANCE)), &
                &'prepared complex residual differs from its independent component')
            call write_complex_component(unit,box,variance_profile,pose_state,h,k,'residual', &
                &actual_residual(h,k),expected_residual(h,k))
            do axis = 1, 5
                max_jacobian_error = max(max_jacobian_error, &
                    &abs(actual_jacobian(h,k,axis)-expected_jacobian(h,k,axis)))
                call assert_true(combined_complex_passes([actual_jacobian(h,k,axis)], &
                    &[expected_jacobian(h,k,axis)], &
                    &FROZEN_ABSOLUTE_TOLERANCES(DERIVATIVE_TOLERANCE), &
                    &FROZEN_RELATIVE_TOLERANCES(DERIVATIVE_TOLERANCE)), &
                    &'prepared complex Jacobian column differs from its independent component')
                write(component,'(a,i0)') 'jacobian_',axis
                call write_complex_component(unit,box,variance_profile,pose_state,h,k,component, &
                    &actual_jacobian(h,k,axis),expected_jacobian(h,k,axis))
            enddo
        enddo
    enddo
end subroutine compare_components

!> Retain scalar, vector, and matrix components and their absolute errors.
subroutine compare_normal_terms(box,variance_profile,pose_state,actual_objective, &
    &expected_objective,actual_gradient,expected_gradient,actual_hessian,expected_hessian, &
    &unit,max_gradient_error,max_hessian_error)
    integer, intent(in) :: box, variance_profile, pose_state, unit
    real(dp), intent(in) :: actual_objective, expected_objective
    real(dp), intent(in) :: actual_gradient(5), expected_gradient(5)
    real(dp), intent(in) :: actual_hessian(5,5), expected_hessian(5,5)
    real(dp), intent(out) :: max_gradient_error, max_hessian_error
    character(len=24) :: component
    integer :: axis, jaxis

    call write_real_component(unit,box,variance_profile,pose_state,'objective', &
        &actual_objective,expected_objective)
    max_gradient_error = maxval(abs(actual_gradient-expected_gradient))
    max_hessian_error = maxval(abs(actual_hessian-expected_hessian))
    do axis = 1, 5
        write(component,'(a,i0)') 'gradient_',axis
        call write_real_component(unit,box,variance_profile,pose_state,component, &
            &actual_gradient(axis),expected_gradient(axis))
        do jaxis = 1, 5
            write(component,'(a,i0,a,i0)') 'hessian_',axis,'_',jaxis
            call write_real_component(unit,box,variance_profile,pose_state,component, &
                &actual_hessian(axis,jaxis),expected_hessian(axis,jaxis))
        enddo
    enddo
end subroutine compare_normal_terms

!> Select one positive native-shell variance profile.
subroutine select_variance(fixture,profile,sigma2)
    type(calibration_fixture), intent(in) :: fixture
    integer, intent(in) :: profile
    real, allocatable, intent(out) :: sigma2(:)

    select case(profile)
    case(1)
        allocate(sigma2(0:fixture%box/2),source=1.0)
    case(2)
        sigma2 = real(fixture%constant_sigma2,sp)
    case(3)
        sigma2 = real(fixture%varying_sigma2,sp)
    case default
        error stop 'unknown objective-normal variance profile'
    end select
end subroutine select_variance

!> Write one complex component as a stable tab-separated row.
subroutine write_complex_component(unit,box,variance_profile,pose_state,h,k,component,actual,expected)
    integer, intent(in) :: unit, box, variance_profile, pose_state, h, k
    character(len=*), intent(in) :: component
    complex(dp), intent(in) :: actual, expected

    write(unit,'(i0,a,a,a,a,a,i0,a,i0,a,a)',advance='no') box,achar(9), &
        &trim(VARIANCE_NAMES(variance_profile)),achar(9),trim(POSE_NAMES(pose_state)), &
        &achar(9),h,achar(9),k,achar(9),trim(component)
    write(unit,'(5(a,es24.16))') achar(9),real(actual,dp),achar(9),aimag(actual), &
        &achar(9),real(expected,dp),achar(9),aimag(expected),achar(9),abs(actual-expected)
end subroutine write_complex_component

!> Write one real normal-equation component as a stable tab-separated row.
subroutine write_real_component(unit,box,variance_profile,pose_state,component,actual,expected)
    integer, intent(in) :: unit, box, variance_profile, pose_state
    character(len=*), intent(in) :: component
    real(dp), intent(in) :: actual, expected

    write(unit,'(i0,a,a,a,a,a,a)',advance='no') box,achar(9), &
        &trim(VARIANCE_NAMES(variance_profile)),achar(9),trim(POSE_NAMES(pose_state)), &
        &achar(9),trim(component)
    write(unit,'(3(a,es24.16))') achar(9),actual,achar(9),expected,achar(9),abs(actual-expected)
end subroutine write_real_component

!> Write one case summary with all maximum component errors.
subroutine write_summary(unit,box,variance_profile,pose_state,active_samples,actual_objective, &
    &expected_objective,max_residual_error,max_jacobian_error,max_gradient_error, &
    &max_hessian_error,hessian_asymmetry,min_switch_margin)
    integer, intent(in) :: unit, box, variance_profile, pose_state, active_samples
    real(dp), intent(in) :: actual_objective, expected_objective, max_residual_error
    real(dp), intent(in) :: max_jacobian_error, max_gradient_error, max_hessian_error
    real(dp), intent(in) :: hessian_asymmetry, min_switch_margin

    write(unit,'(i0,a,a,a,a,a,i0)',advance='no') box,achar(9), &
        &trim(VARIANCE_NAMES(variance_profile)),achar(9),trim(POSE_NAMES(pose_state)), &
        &achar(9),active_samples
    write(unit,'(8(a,es24.16))') achar(9),actual_objective,achar(9),expected_objective, &
        &achar(9),max_residual_error,achar(9),max_jacobian_error,achar(9),max_gradient_error, &
        &achar(9),max_hessian_error,achar(9),hessian_asymmetry,achar(9),min_switch_margin
end subroutine write_summary

!> Read one required key=value command argument.
function required_argument(key) result(value)
    character(len=*), intent(in) :: key
    character(len=:), allocatable :: value
    character(len=4096) :: argument
    integer :: argument_index, argument_status, separator

    value = ''
    do argument_index = 1, command_argument_count()
        call get_command_argument(argument_index,argument,status=argument_status)
        if( argument_status /= 0 ) error stop 'could not read objective-normal argument'
        separator = index(argument,'=')
        if( separator <= 0 ) cycle
        if( trim(argument(:separator-1)) /= key ) cycle
        value = trim(argument(separator+1:))
    enddo
    if( len(value) == 0 ) &
        &error stop 'objective-normal test requires evidence_dir=<existing-directory>'
end function required_argument

end module pose_cont_refinement_objective_normals_test
