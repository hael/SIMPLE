! Complete Phase 9 configuration matrix for the prepared Cartesian pose owner.
module pose_cont_refinement_configuration_matrix_test
use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
use pose_cont_refinement_calibration_helpers, only: FROZEN_ABSOLUTE_TOLERANCES, &
    &FROZEN_RELATIVE_TOLERANCES, combined_real_passes, combined_complex_passes
use pose_cont_refinement_test_helpers, only: assert_true
use simple_defs, only: dp, sp, PI, OSMPL_PAD_FAC
use simple_cartesian_pose_refiner, only: cartesian_pose_refiner, cartesian_pose_data, &
    &right_increment_rotation, POSE_DATA_VALID, POSE_LM_INVALID_NUMERICS
use simple_ctf, only: ctf
use simple_type_defs, only: ctfparams, CTFFLAG_NO, CTFFLAG_YES, CTFFLAG_FLIP
implicit none
private
public :: run_configuration_matrix

integer, parameter :: MATRIX_BOXES(4) = [10,14,16,18]
integer, parameter :: N_VOLUME_VARIANTS = 3
integer, parameter :: N_ROTATION_AXES = 3
integer, parameter :: N_SHIFT_SIGNS = 2
integer, parameter :: N_CTF_MODES = 3
integer, parameter :: N_VARIANCE_PROFILES = 3
integer, parameter :: N_FOURIER_RANGES = 2
integer, parameter :: ALGEBRAIC_TOLERANCE = 1
integer, parameter :: SLOW_GATHER_TOLERANCE = 6
integer, parameter :: MAX_LM_ITERATIONS = 2
integer, parameter :: EXPECTED_CASES = size(MATRIX_BOXES)*N_VOLUME_VARIANTS* &
    &N_ROTATION_AXES*N_SHIFT_SIGNS*N_CTF_MODES*N_VARIANCE_PROFILES*N_FOURIER_RANGES
character(len=12), parameter :: CTF_NAMES(N_CTF_MODES) = [character(len=12) :: &
    &'disabled', 'signed', 'phase_flip']
character(len=12), parameter :: VARIANCE_NAMES(N_VARIANCE_PROFILES) = [character(len=12) :: &
    &'constant', 'varying', 'short']
character(len=8), parameter :: RANGE_NAMES(N_FOURIER_RANGES) = [character(len=8) :: &
    &'full', 'inner']

contains

!> Run the frozen complete matrix and retain one summary row per configuration.
subroutine run_configuration_matrix()
    type(cartesian_pose_refiner) :: workspace
    type(cartesian_pose_data) :: data
    type(ctfparams) :: ctfparms
    complex, allocatable :: raw_observed(:,:)
    real, allocatable :: volume(:,:,:), sigma2(:)
    real(dp) :: accepted_objectives(0:MAX_LM_ITERATIONS)
    real(dp) :: exact_gradient(5), start_gradient(5), exact_objective, start_objective
    real(dp) :: identity(3,3), truth_rotmat(3,3), start_rotmat(3,3), refined_rotmat(3,3)
    real(dp) :: truth_shift(2), start_shift(2), refined_shift(2), omega(3), perturbation(3)
    real(dp) :: max_rotation_step, max_shift_step
    real(dp) :: max_slow_gather_error
    character(len=:), allocatable :: evidence_directory
    integer :: axis, box, box_index, case_count, ctf_mode, data_status, matrix_unit
    integer :: lm_status, nattempted, naccepted, nstencil_switches, range_index
    integer :: requested_range(2), shift_sign, variance_profile, variant

    evidence_directory = required_argument('evidence_dir')
    open(newunit=matrix_unit,file=evidence_directory//'/configuration_matrix.tsv', &
        &status='replace',action='write')
    write(matrix_unit,'(a)') 'case'//achar(9)//'box'//achar(9)//'variant'//achar(9)// &
        &'axis'//achar(9)//'shift_sign'//achar(9)//'ctf_index'//achar(9)//'ctf'//achar(9)// &
        &'variance'//achar(9)//'range'//achar(9)//'effective_low'//achar(9)// &
        &'effective_high'//achar(9)// &
        &'active_samples'//achar(9)//'exact_objective'//achar(9)//'nonstationary_objective'// &
        &achar(9)//'gradient_norm'//achar(9)//'max_slow_gather_error'//achar(9)// &
        &'lm_status'//achar(9)//'attempted'//achar(9)// &
        &'accepted'//achar(9)//'final_objective'//achar(9)//'max_rotation_step'//achar(9)// &
        &'max_shift_step'//achar(9)//'stencil_switches'

    identity = reshape([1._dp,0._dp,0._dp,0._dp,1._dp,0._dp,0._dp,0._dp,1._dp],[3,3])
    case_count = 0
    do box_index = 1, size(MATRIX_BOXES)
        box = MATRIX_BOXES(box_index)
        do variant = 1, N_VOLUME_VARIANTS
            call build_matrix_volume(box,variant,volume)
            call workspace%new(volume)
            do axis = 1, N_ROTATION_AXES
                omega = 0._dp
                omega(axis) = (0.009_dp+0.002_dp*real(variant,dp))*merge(1._dp,-1._dp,mod(box,4) == 2)
                truth_rotmat = right_increment_rotation(identity,omega)
                perturbation = 0._dp
                perturbation(axis) = 0.018_dp+0.003_dp*real(variant,dp)
                do shift_sign = -1, 1, 2
                    truth_shift = real(shift_sign,dp)*[0.13_dp,-0.09_dp]
                    start_rotmat = right_increment_rotation(truth_rotmat,perturbation)
                    start_shift = truth_shift+real(shift_sign,dp)*[0.08_dp,-0.06_dp]
                    do ctf_mode = 1, N_CTF_MODES
                        call make_ctf_parameters(ctf_mode,ctfparms)
                        call build_raw_observation(workspace,truth_rotmat,truth_shift,ctfparms, &
                            &raw_observed,max_slow_gather_error)
                        do variance_profile = 1, N_VARIANCE_PROFILES
                            call build_variance(box,variant,variance_profile,sigma2)
                            do range_index = 1, N_FOURIER_RANGES
                                call select_range(box,range_index,requested_range)
                                call workspace%prepare_particle(raw_observed,ctfparms,sigma2, &
                                    &requested_range,data,data_status)
                                call assert_true(data_status == POSE_DATA_VALID .and. data%is_valid(), &
                                    &'configuration matrix produced invalid prepared data')
                                call workspace%prepared_objective_gradient(truth_rotmat,truth_shift,data, &
                                    &exact_objective,exact_gradient)
                                call workspace%prepared_objective_gradient(start_rotmat,start_shift,data, &
                                    &start_objective,start_gradient)
                                call assert_case_state(exact_objective,exact_gradient,start_objective, &
                                    &start_gradient)
                                refined_rotmat = start_rotmat
                                refined_shift = start_shift
                                call workspace%refine_prepared_pose_lm(refined_rotmat,refined_shift,data, &
                                    &real(box,dp)/2._dp,MAX_LM_ITERATIONS,accepted_objectives,naccepted, &
                                    &lm_status,nattempted,max_rotation_step,max_shift_step,nstencil_switches)
                                call assert_lm_trace(start_objective,accepted_objectives,naccepted,lm_status)
                                case_count = case_count+1
                                call write_case(matrix_unit,case_count,box,variant,axis,shift_sign,ctf_mode, &
                                    &variance_profile,range_index,data,exact_objective,start_objective, &
                                    &start_gradient,max_slow_gather_error,lm_status,nattempted,naccepted, &
                                    &accepted_objectives, &
                                    &max_rotation_step,max_shift_step,nstencil_switches)
                            enddo
                        enddo
                    enddo
                enddo
            enddo
            call workspace%kill
        enddo
    enddo
    close(matrix_unit)
    call assert_true(case_count == EXPECTED_CASES, &
        &'configuration matrix did not execute every declared combination')
    call write_manifest(evidence_directory,case_count)
    write(*,'(a,i0)') 'POSE_CONT_REFINEMENT_CONFIGURATION_MATRIX cases: ',case_count
    write(*,'(a)') 'POSE_CONT_REFINEMENT_CONFIGURATION_MATRIX boxes: 10 14 16 18'
    write(*,'(a)') 'POSE_CONT_REFINEMENT_CONFIGURATION_MATRIX: PASS'
end subroutine run_configuration_matrix

!> Build one of three deterministic asymmetric physical volumes.
subroutine build_matrix_volume(box,variant,volume)
    integer, intent(in) :: box, variant
    real, allocatable, intent(out) :: volume(:,:,:)
    real(dp) :: amplitude(4), centre, centres(3,4), dx, dy, dz, radius2, width(4)
    integer :: blob, i, j, k

    if( variant < 1 .or. variant > N_VOLUME_VARIANTS ) &
        &error stop 'configuration matrix requested an invalid volume variant'
    allocate(volume(box,box,box),source=0.)
    centre = real(box,dp)/2._dp+0.5_dp
    centres = reshape([ &
        &-0.24_dp, 0.11_dp, 0.18_dp, 0.19_dp,-0.27_dp,-0.09_dp, &
        & 0.07_dp, 0.29_dp,-0.23_dp,-0.30_dp,-0.06_dp, 0.26_dp],[3,4])
    centres = real(box,dp)*cshift(centres,shift=variant-1,dim=2)
    amplitude = cshift([0.93_dp,-0.58_dp,0.36_dp,0.21_dp],shift=variant-1)
    width = real(box,dp)*cshift([0.091_dp,0.126_dp,0.149_dp,0.078_dp],shift=variant-1)
    do k = 1, box
        do j = 1, box
            do i = 1, box
                do blob = 1, 4
                    dx = real(i,dp)-centre-centres(1,blob)
                    dy = real(j,dp)-centre-centres(2,blob)
                    dz = real(k,dp)-centre-centres(3,blob)
                    radius2 = (0.84_dp+0.04_dp*variant)*dx*dx+ &
                        &(1.16_dp-0.03_dp*variant)*dy*dy+dz*dz+0.05_dp*dx*dz
                    volume(i,j,k) = volume(i,j,k)+real(amplitude(blob)* &
                        &exp(-radius2/(2._dp*width(blob)**2)),sp)
                enddo
            enddo
        enddo
    enddo
end subroutine build_matrix_volume

!> Build raw data at the exact pose under the selected signed CTF convention.
subroutine build_raw_observation(workspace,rotmat,shift,ctfparms,raw_observed,max_slow_gather_error)
    type(cartesian_pose_refiner), intent(in) :: workspace
    real(dp), intent(in) :: rotmat(3,3), shift(2)
    type(ctfparams), intent(in) :: ctfparms
    complex, allocatable, intent(out) :: raw_observed(:,:)
    real(dp), intent(out) :: max_slow_gather_error
    type(ctf) :: tfun
    complex :: derivative(3), phase, value
    complex(dp) :: slow_value
    real :: cval, margin(3)
    real(sp) :: loc(3)
    real(dp) :: argument
    integer :: box, h, k, lims2(2,2)
    logical :: use_ctf

    lims2 = workspace%get_lims2()
    box = lims2(1,2)-lims2(1,1)
    allocate(raw_observed(lims2(1,1):lims2(1,2),lims2(2,1):lims2(2,2)), &
        &source=cmplx(0.,0.))
    use_ctf = ctfparms%ctfflag /= CTFFLAG_NO
    if( use_ctf )then
        tfun = ctf(ctfparms%smpd,ctfparms%kv,ctfparms%cs,ctfparms%fraca)
        call tfun%init(ctfparms%dfx,ctfparms%dfy,ctfparms%angast)
    endif
    max_slow_gather_error = 0._dp
    do k = lims2(2,1), lims2(2,2)
        do h = lims2(1,1), lims2(1,2)
            if( h*h+k*k > (box/2)**2 ) cycle
            loc = real(OSMPL_PAD_FAC,sp)*real(matmul(real([h,k,0],dp),rotmat),sp)
            call workspace%sample_with_grad(loc,value,derivative,margin)
            call workspace%sample_slow_test(loc,slow_value)
            max_slow_gather_error = max(max_slow_gather_error, &
                &abs(cmplx(value,kind=dp)-cmplx(slow_value,kind=dp)))
            call assert_true(combined_complex_passes([cmplx(value,kind=dp)], &
                &[cmplx(slow_value,kind=dp)],FROZEN_ABSOLUTE_TOLERANCES(SLOW_GATHER_TOLERANCE), &
                &FROZEN_RELATIVE_TOLERANCES(SLOW_GATHER_TOLERANCE)), &
                &'configuration matrix fast and slow gathers differ')
            argument = 2._dp*real(PI,dp)*(real(h,dp)*shift(1)+real(k,dp)*shift(2))/real(box,dp)
            phase = cmplx(cos(argument),sin(argument),kind=sp)
            cval = 1.
            if( use_ctf ) cval = scalar_ctf_value(tfun,ctfparms,box,h,k)
            raw_observed(h,k) = cval*phase*value
        enddo
    enddo
end subroutine build_raw_observation

!> Set one supported CTF application mode with a nonzero stored phase.
subroutine make_ctf_parameters(mode,ctfparms)
    integer, intent(in) :: mode
    type(ctfparams), intent(out) :: ctfparms

    ctfparms%smpd = 1.0
    ctfparms%kv = 300.0
    ctfparms%cs = 2.7
    ctfparms%fraca = 0.1
    ctfparms%dfx = 1.4
    ctfparms%dfy = 1.65
    ctfparms%angast = 23.0
    ctfparms%phshift = 0.37
    select case(mode)
    case(1)
        ctfparms%ctfflag = CTFFLAG_NO
    case(2)
        ctfparms%ctfflag = CTFFLAG_YES
    case(3)
        ctfparms%ctfflag = CTFFLAG_FLIP
    case default
        error stop 'configuration matrix requested an invalid CTF mode'
    end select
end subroutine make_ctf_parameters

!> Evaluate the scalar signed CTF without calling the pose transfer routine.
real function scalar_ctf_value(tfun,ctfparms,box,h,k) result(cval)
    type(ctf), intent(in) :: tfun
    type(ctfparams), intent(in) :: ctfparms
    integer, intent(in) :: box, h, k
    real :: angle, spatial_frequency_squared

    spatial_frequency_squared = real(h*h+k*k)/(real(box)*ctfparms%smpd)**2
    angle = atan2(real(k),real(h))
    cval = tfun%eval(spatial_frequency_squared,angle,ctfparms%phshift)
    if( ctfparms%ctfflag == CTFFLAG_FLIP ) cval = abs(cval)
end function scalar_ctf_value

!> Build constant, nonstationary, or intentionally short positive variance.
subroutine build_variance(box,variant,profile,sigma2)
    integer, intent(in) :: box, variant, profile
    real, allocatable, intent(out) :: sigma2(:)
    integer :: high, shell

    high = box/2
    if( profile == 3 ) high = max(2,box/2-2)
    allocate(sigma2(0:high))
    select case(profile)
    case(1)
        sigma2 = 0.91+0.07*real(variant)
    case(2,3)
        sigma2 = [(0.79+0.041*real(shell)+0.009*real(shell*shell)+ &
            &0.025*real(variant),shell=0,high)]
    case default
        error stop 'configuration matrix requested an invalid variance profile'
    end select
end subroutine build_variance

!> Select the complete disk or an interior annular Fourier range.
subroutine select_range(box,range_index,shell_range)
    integer, intent(in) :: box, range_index
    integer, intent(out) :: shell_range(2)

    select case(range_index)
    case(1)
        shell_range = [0,box/2]
    case(2)
        shell_range = [1,max(2,box/3)]
    case default
        error stop 'configuration matrix requested an invalid Fourier range'
    end select
end subroutine select_range

!> Require finite exact and nonstationary objective states under frozen algebraic error.
subroutine assert_case_state(exact_objective,exact_gradient,start_objective,start_gradient)
    real(dp), intent(in) :: exact_objective, exact_gradient(5), start_objective, start_gradient(5)

    call assert_true(ieee_is_finite(exact_objective) .and. &
        &all(ieee_is_finite(exact_gradient)), 'exact matrix pose produced nonfinite terms')
    call assert_true(ieee_is_finite(start_objective) .and. &
        &all(ieee_is_finite(start_gradient)), 'nonstationary matrix pose produced nonfinite terms')
    call assert_true(combined_real_passes([exact_objective],[0._dp], &
        &FROZEN_ABSOLUTE_TOLERANCES(ALGEBRAIC_TOLERANCE), &
        &FROZEN_RELATIVE_TOLERANCES(ALGEBRAIC_TOLERANCE)), &
        &'exact matrix pose did not reproduce its prepared observation')
    call assert_true(start_objective+FROZEN_ABSOLUTE_TOLERANCES(ALGEBRAIC_TOLERANCE) >= &
        &exact_objective, 'nonstationary matrix objective was below its exact-pose objective')
end subroutine assert_case_state

!> Require a finite, monotone accepted LM trace and no invalid-numerics outcome.
subroutine assert_lm_trace(start_objective,trace,naccepted,status)
    real(dp), intent(in) :: start_objective, trace(0:)
    integer, intent(in) :: naccepted, status
    integer :: accepted

    call assert_true(status /= POSE_LM_INVALID_NUMERICS, &
        &'configuration matrix LM returned invalid numerics')
    call assert_true(naccepted >= 0 .and. naccepted <= ubound(trace,1), &
        &'configuration matrix LM returned an invalid accepted-step count')
    call assert_true(ieee_is_finite(trace(0)) .and. &
        &trace(0) <= start_objective+FROZEN_ABSOLUTE_TOLERANCES(ALGEBRAIC_TOLERANCE), &
        &'configuration matrix LM initial objective changed unexpectedly')
    do accepted = 1, naccepted
        call assert_true(ieee_is_finite(trace(accepted)) .and. &
            &trace(accepted) <= trace(accepted-1), &
            &'configuration matrix LM accepted a nonmonotone objective')
    enddo
end subroutine assert_lm_trace

!> Write one complete matrix outcome.
subroutine write_case(unit,case_index,box,variant,axis,shift_sign,ctf_mode,variance_profile, &
    &range_index,data,exact_objective,start_objective,start_gradient,max_slow_gather_error, &
    &lm_status,nattempted, &
    &naccepted,trace,max_rotation_step,max_shift_step,nstencil_switches)
    integer, intent(in) :: unit, case_index, box, variant, axis, shift_sign, ctf_mode
    integer, intent(in) :: variance_profile, range_index, lm_status, nattempted, naccepted
    integer, intent(in) :: nstencil_switches
    type(cartesian_pose_data), intent(in) :: data
    real(dp), intent(in) :: exact_objective, start_objective, start_gradient(5), trace(0:)
    real(dp), intent(in) :: max_slow_gather_error
    real(dp), intent(in) :: max_rotation_step, max_shift_step
    integer :: effective_range(2)

    effective_range = data%get_shell_range()
    write(unit,'(i0,5(a,i0),3(a,a),3(a,i0),4(a,es24.16),3(a,i0), &
        &3(a,es24.16),a,i0)') case_index,achar(9),box,achar(9),variant,achar(9),axis, &
        &achar(9),shift_sign,achar(9),ctf_mode,achar(9),trim(CTF_NAMES(ctf_mode)),achar(9), &
        &trim(VARIANCE_NAMES(variance_profile)),achar(9),trim(RANGE_NAMES(range_index)), &
        &achar(9),effective_range(1),achar(9),effective_range(2),achar(9), &
        &data%get_active_sample_count(),achar(9),exact_objective,achar(9),start_objective, &
        &achar(9),sqrt(sum(start_gradient*start_gradient)),achar(9),max_slow_gather_error, &
        &achar(9),lm_status,achar(9), &
        &nattempted,achar(9),naccepted,achar(9),trace(naccepted),achar(9),max_rotation_step, &
        &achar(9),max_shift_step,achar(9),nstencil_switches
end subroutine write_case

!> Retain the frozen matrix dimensions and expected row count.
subroutine write_manifest(evidence_directory,case_count)
    character(len=*), intent(in) :: evidence_directory
    integer, intent(in) :: case_count
    integer :: unit

    open(newunit=unit,file=evidence_directory//'/configuration_matrix_manifest.tsv', &
        &status='replace',action='write')
    write(unit,'(a)') 'field'//achar(9)//'value'
    write(unit,'(a,a,i0)') 'boxes',achar(9),size(MATRIX_BOXES)
    write(unit,'(a,a,i0)') 'volume_variants',achar(9),N_VOLUME_VARIANTS
    write(unit,'(a,a,i0)') 'rotation_axes',achar(9),N_ROTATION_AXES
    write(unit,'(a,a,i0)') 'shift_signs',achar(9),N_SHIFT_SIGNS
    write(unit,'(a,a,i0)') 'ctf_modes',achar(9),N_CTF_MODES
    write(unit,'(a,a,i0)') 'variance_profiles',achar(9),N_VARIANCE_PROFILES
    write(unit,'(a,a,i0)') 'fourier_ranges',achar(9),N_FOURIER_RANGES
    write(unit,'(a,a,i0)') 'expected_cases',achar(9),EXPECTED_CASES
    write(unit,'(a,a,i0)') 'executed_cases',achar(9),case_count
    close(unit)
end subroutine write_manifest

!> Read one required old-style key=value argument.
function required_argument(key) result(value)
    character(len=*), intent(in) :: key
    character(len=:), allocatable :: value
    character(len=4096) :: argument
    integer :: argument_status, iarg, separator

    value = ''
    do iarg = 1, command_argument_count()
        call get_command_argument(iarg,argument,status=argument_status)
        if( argument_status /= 0 ) error stop 'could not read configuration-matrix argument'
        separator = index(argument,'=')
        if( separator <= 0 ) cycle
        if( trim(argument(:separator-1)) /= key ) cycle
        value = trim(argument(separator+1:))
    enddo
    if( len(value) == 0 ) &
        &error stop 'configuration matrix requires evidence_dir=<existing-directory>'
end function required_argument

end module pose_cont_refinement_configuration_matrix_test
