! Test-only calibration infrastructure discovered by the pose-suite dependency glob.
module pose_cont_refinement_calibration_helpers
use simple_defs, only: dp, sp
implicit none
private

integer, parameter, public :: CALIBRATION_BOXES(2) = [8,12]
integer, parameter, public :: ACCEPTANCE_BOXES(2) = [10,16]
integer, parameter, public :: N_TOLERANCE_FAMILIES = 7
real(dp), parameter, public :: CALIBRATION_SAFETY_FACTOR = 8._dp
! Phase 4 Oracle artifacts froze these values before acceptance-fixture tests.
real(dp), parameter, public :: FROZEN_ABSOLUTE_TOLERANCES(N_TOLERANCE_FAMILIES) = [ &
    &4.9353318740941177e-4_dp, 2.0599365235796085e-7_dp, 1.e-8_dp, 1.e-10_dp, &
    &1.e-5_dp, 1.1541020711748716e-4_dp, 1.e-5_dp]
real(dp), parameter, public :: FROZEN_RELATIVE_TOLERANCES(N_TOLERANCE_FAMILIES) = [ &
    &5.7686781070864920e-6_dp, 2.e-5_dp, 5.e-3_dp, 2.e-5_dp, 3.e-2_dp, &
    &2.e-5_dp, 5.e-2_dp]
integer, parameter :: MAX_CALIBRATION_OBSERVATIONS = 64

character(len=32), parameter, public :: TOLERANCE_NAMES(N_TOLERANCE_FAMILIES) = [character(len=32) :: &
    &'algebraic', 'lm_system', 'derivative', 'analytic_dft', &
    &'executed_dft', 'slow_gather', 'finite_projection']
real(dp), parameter :: ABSOLUTE_FLOORS(N_TOLERANCE_FAMILIES) = [ &
    &1.e-10_dp, 1.e-10_dp, 1.e-8_dp, 1.e-10_dp, 1.e-5_dp, 1.e-7_dp, 1.e-5_dp]
real(dp), parameter :: RELATIVE_FLOORS(N_TOLERANCE_FAMILIES) = [ &
    &5.e-6_dp, 2.e-5_dp, 5.e-3_dp, 2.e-5_dp, 3.e-2_dp, 2.e-5_dp, 5.e-2_dp]
real(dp), parameter :: RELATIVE_SCALE_FLOORS(N_TOLERANCE_FAMILIES) = [ &
    &1.e-12_dp, 1._dp, 1.e-8_dp, 1.e-12_dp, 1.e-3_dp, 1.e-10_dp, 1.e-3_dp]
real(dp), parameter :: LOOSE_ABSOLUTE_LIMITS(N_TOLERANCE_FAMILIES) = [ &
    &1.e-3_dp, 1.e-3_dp, 5.e-2_dp, 1.e-3_dp, 1.e-1_dp, 1.e-2_dp, 2.e-1_dp]
real(dp), parameter :: LOOSE_RELATIVE_LIMITS(N_TOLERANCE_FAMILIES) = [ &
    &1.e-3_dp, 1.e-3_dp, 3.e-2_dp, 1.e-3_dp, 8.e-2_dp, 1.e-2_dp, 1.e-1_dp]

type, public :: calibration_fixture
    integer :: box = 0
    integer :: fixture_id = 0
    real(dp), allocatable :: volume(:,:,:)
    real(dp) :: exact_pose(5) = 0._dp
    real(dp) :: nonstationary_pose(5) = 0._dp
    real(dp), allocatable :: constant_sigma2(:)
    real(dp), allocatable :: varying_sigma2(:)
    real(dp) :: ordinary_ctf(4) = 0._dp
    real(dp) :: attenuated_ctf(4) = 0._dp
end type calibration_fixture

type, public :: tolerance_record
    character(len=32) :: name = ''
    real(dp) :: absolute_floor = 0._dp
    real(dp) :: relative_floor = 0._dp
    real(dp) :: relative_scale_floor = 0._dp
    real(dp) :: maximum_absolute_error = 0._dp
    real(dp) :: maximum_scaled_relative_error = 0._dp
    real(dp) :: absolute_tolerance = 0._dp
    real(dp) :: relative_tolerance = 0._dp
    integer :: observations = 0
    real(dp) :: absolute_errors(MAX_CALIBRATION_OBSERVATIONS) = 0._dp
    real(dp) :: scaled_relative_errors(MAX_CALIBRATION_OBSERVATIONS) = 0._dp
end type tolerance_record

public :: initialize_tolerances
public :: build_calibration_fixture
public :: build_acceptance_fixture
public :: observe_real_comparison
public :: observe_complex_comparison
public :: combined_real_passes
public :: combined_complex_passes
public :: finalize_tolerances
public :: write_calibration_artifacts
public :: fixtures_are_identical

contains

!> Initialize all tolerance families with their predeclared floors.
subroutine initialize_tolerances(records)
    type(tolerance_record), intent(out) :: records(N_TOLERANCE_FAMILIES)
    integer :: family

    do family = 1, N_TOLERANCE_FAMILIES
        records(family)%name = TOLERANCE_NAMES(family)
        records(family)%absolute_floor = ABSOLUTE_FLOORS(family)
        records(family)%relative_floor = RELATIVE_FLOORS(family)
        records(family)%relative_scale_floor = RELATIVE_SCALE_FLOORS(family)
    enddo
end subroutine initialize_tolerances

!> Build one deterministic asymmetric calibration fixture.
subroutine build_calibration_fixture(box,fixture)
    integer, intent(in) :: box
    type(calibration_fixture), intent(out) :: fixture
    real(dp) :: centre, dx, dy, dz, radius2
    real(dp) :: centres(3,4), amplitudes(4), widths(4)
    integer :: blob, i, j, k

    if( .not. any(box == CALIBRATION_BOXES) ) &
        &error stop 'calibration requested a noncalibration box'
    fixture%box = box
    fixture%fixture_id = merge(4081,4127,box == CALIBRATION_BOXES(1))
    allocate(fixture%volume(box,box,box),source=0._dp)
    allocate(fixture%constant_sigma2(0:box/2),source=1.17_dp)
    allocate(fixture%varying_sigma2(0:box/2))
    centre = real(box,dp)/2._dp+0.5_dp
    centres = reshape([ &
        &-0.23_dp*box,-0.11_dp*box, 0.17_dp*box, &
        & 0.19_dp*box, 0.24_dp*box,-0.16_dp*box, &
        &-0.04_dp*box, 0.27_dp*box, 0.08_dp*box, &
        & 0.28_dp*box,-0.18_dp*box, 0.25_dp*box],[3,4])
    amplitudes = [1._dp,0.73_dp,-0.41_dp,0.29_dp]
    widths = real(box,dp)*[0.105_dp,0.137_dp,0.083_dp,0.119_dp]
    do k = 1, box
        do j = 1, box
            do i = 1, box
                do blob = 1, 4
                    dx = real(i,dp)-centre-centres(1,blob)
                    dy = real(j,dp)-centre-centres(2,blob)
                    dz = real(k,dp)-centre-centres(3,blob)
                    radius2 = dx*dx+1.13_dp*dy*dy+0.91_dp*dz*dz+0.07_dp*dx*dz
                    fixture%volume(i,j,k) = fixture%volume(i,j,k)+ &
                        &amplitudes(blob)*exp(-radius2/(2._dp*widths(blob)**2))
                enddo
            enddo
        enddo
    enddo
    fixture%exact_pose = merge([0.011_dp,-0.008_dp,0.006_dp,0.13_dp,-0.09_dp], &
        &[-0.009_dp,0.013_dp,-0.007_dp,-0.15_dp,0.11_dp],box == CALIBRATION_BOXES(1))
    fixture%nonstationary_pose = merge([0.031_dp,-0.022_dp,0.017_dp,0.27_dp,-0.19_dp], &
        &[-0.024_dp,0.035_dp,-0.019_dp,-0.31_dp,0.23_dp],box == CALIBRATION_BOXES(1))
    do i = 0, box/2
        fixture%varying_sigma2(i) = 0.82_dp+0.055_dp*real(i,dp)+0.013_dp*real(i*i,dp)
    enddo
    fixture%ordinary_ctf = merge([1._dp,0.18_dp,1.37_dp,0.11_dp], &
        &[0.74_dp,-0.27_dp,1.61_dp,0.16_dp],box == CALIBRATION_BOXES(1))
    fixture%attenuated_ctf = 0.43_dp*fixture%ordinary_ctf
end subroutine build_calibration_fixture

!> Build one deterministic acceptance fixture for later phases only.
subroutine build_acceptance_fixture(box,fixture)
    integer, intent(in) :: box
    type(calibration_fixture), intent(out) :: fixture
    real(dp) :: centre, dx, dy, dz, radius2
    real(dp) :: centres(3,5), amplitudes(5), widths(5)
    integer :: blob, i, j, k

    if( .not. any(box == ACCEPTANCE_BOXES) ) &
        &error stop 'acceptance requested a nonacceptance box'
    fixture%box = box
    fixture%fixture_id = 4900+box
    allocate(fixture%volume(box,box,box),source=0._dp)
    allocate(fixture%constant_sigma2(0:box/2),source=0.93_dp)
    allocate(fixture%varying_sigma2(0:box/2))
    centre = real(box,dp)/2._dp+0.5_dp
    centres = reshape([ &
        & 0.21_dp*box,-0.26_dp*box, 0.09_dp*box, &
        &-0.29_dp*box, 0.14_dp*box,-0.19_dp*box, &
        & 0.07_dp*box, 0.31_dp*box, 0.22_dp*box, &
        &-0.16_dp*box,-0.08_dp*box, 0.28_dp*box, &
        & 0.30_dp*box, 0.04_dp*box,-0.27_dp*box],[3,5])
    amplitudes = [0.91_dp,-0.52_dp,0.37_dp,0.24_dp,-0.18_dp]
    widths = real(box,dp)*[0.097_dp,0.128_dp,0.076_dp,0.143_dp,0.089_dp]
    do k = 1, box
        do j = 1, box
            do i = 1, box
                do blob = 1, 5
                    dx = real(i,dp)-centre-centres(1,blob)
                    dy = real(j,dp)-centre-centres(2,blob)
                    dz = real(k,dp)-centre-centres(3,blob)
                    radius2 = 0.87_dp*dx*dx+1.19_dp*dy*dy+dz*dz-0.06_dp*dy*dz
                    fixture%volume(i,j,k) = fixture%volume(i,j,k)+ &
                        &amplitudes(blob)*exp(-radius2/(2._dp*widths(blob)**2))
                enddo
            enddo
        enddo
    enddo
    fixture%exact_pose = merge([-0.008_dp,0.012_dp,0.009_dp,-0.12_dp,0.16_dp], &
        &[0.014_dp,-0.011_dp,0.005_dp,0.17_dp,-0.13_dp],box == ACCEPTANCE_BOXES(1))
    fixture%nonstationary_pose = merge([-0.027_dp,0.019_dp,0.033_dp,-0.22_dp,0.34_dp], &
        &[0.038_dp,-0.029_dp,0.014_dp,0.36_dp,-0.25_dp],box == ACCEPTANCE_BOXES(1))
    do i = 0, box/2
        fixture%varying_sigma2(i) = 0.91_dp+0.037_dp*real(i,dp)+0.009_dp*real(i*i,dp)
    enddo
    fixture%ordinary_ctf = merge([0.83_dp,0.22_dp,1.49_dp,-0.09_dp], &
        &[0.68_dp,-0.31_dp,1.72_dp,0.14_dp],box == ACCEPTANCE_BOXES(1))
    fixture%attenuated_ctf = 0.39_dp*fixture%ordinary_ctf
end subroutine build_acceptance_fixture

!> Record componentwise absolute and scaled-relative errors for real values.
subroutine observe_real_comparison(record,actual,expected)
    type(tolerance_record), intent(inout) :: record
    real(dp), intent(in) :: actual(:), expected(:)
    real(dp) :: absolute_error, relative_error, scale
    integer :: component

    if( size(actual) /= size(expected) ) error stop 'calibration comparison shapes differ'
    do component = 1, size(actual)
        absolute_error = abs(actual(component)-expected(component))
        scale = max(record%relative_scale_floor,abs(actual(component)),abs(expected(component)))
        relative_error = absolute_error/scale
        if( record%observations == MAX_CALIBRATION_OBSERVATIONS ) &
            &error stop 'real calibration observation capacity exceeded'
        record%observations = record%observations+1
        record%absolute_errors(record%observations) = absolute_error
        record%scaled_relative_errors(record%observations) = relative_error
        record%maximum_absolute_error = max(record%maximum_absolute_error,absolute_error)
        record%maximum_scaled_relative_error = max(record%maximum_scaled_relative_error,relative_error)
    enddo
end subroutine observe_real_comparison

!> Record componentwise absolute and scaled-relative errors for complex values.
subroutine observe_complex_comparison(record,actual,expected)
    type(tolerance_record), intent(inout) :: record
    complex(dp), intent(in) :: actual(:), expected(:)
    real(dp) :: absolute_error, relative_error, scale
    integer :: component

    if( size(actual) /= size(expected) ) error stop 'complex calibration comparison shapes differ'
    do component = 1, size(actual)
        absolute_error = abs(actual(component)-expected(component))
        scale = max(record%relative_scale_floor,abs(actual(component)),abs(expected(component)))
        relative_error = absolute_error/scale
        if( record%observations == MAX_CALIBRATION_OBSERVATIONS ) &
            &error stop 'complex calibration observation capacity exceeded'
        record%observations = record%observations+1
        record%absolute_errors(record%observations) = absolute_error
        record%scaled_relative_errors(record%observations) = relative_error
        record%maximum_absolute_error = max(record%maximum_absolute_error,absolute_error)
        record%maximum_scaled_relative_error = max(record%maximum_scaled_relative_error,relative_error)
    enddo
end subroutine observe_complex_comparison

!> Apply the frozen componentwise combined rule to real values.
pure logical function combined_real_passes(actual,expected,absolute_tolerance,relative_tolerance) result(passes)
    real(dp), intent(in) :: actual(:), expected(:), absolute_tolerance, relative_tolerance

    passes = size(actual) == size(expected) .and. absolute_tolerance >= 0._dp .and. &
        &relative_tolerance >= 0._dp
    if( .not. passes ) return
    passes = all(abs(actual-expected) <= absolute_tolerance+ &
        &relative_tolerance*max(abs(actual),abs(expected)))
end function combined_real_passes

!> Apply the frozen componentwise combined rule to complex values.
pure logical function combined_complex_passes(actual,expected,absolute_tolerance,relative_tolerance) result(passes)
    complex(dp), intent(in) :: actual(:), expected(:)
    real(dp), intent(in) :: absolute_tolerance, relative_tolerance

    passes = size(actual) == size(expected) .and. absolute_tolerance >= 0._dp .and. &
        &relative_tolerance >= 0._dp
    if( .not. passes ) return
    passes = all(abs(actual-expected) <= absolute_tolerance+ &
        &relative_tolerance*max(abs(actual),abs(expected)))
end function combined_complex_passes

!> Derive frozen tolerances and reject scientifically loose calibration results.
subroutine finalize_tolerances(records)
    type(tolerance_record), intent(inout) :: records(N_TOLERANCE_FAMILIES)
    integer :: family

    do family = 1, N_TOLERANCE_FAMILIES
        if( records(family)%observations < 1 ) error stop 'tolerance family has no calibration observation'
        records(family)%absolute_tolerance = max(records(family)%absolute_floor, &
            &CALIBRATION_SAFETY_FACTOR*records(family)%maximum_absolute_error)
        records(family)%relative_tolerance = max(records(family)%relative_floor, &
            &CALIBRATION_SAFETY_FACTOR*records(family)%maximum_scaled_relative_error)
        if( records(family)%absolute_tolerance > LOOSE_ABSOLUTE_LIMITS(family) .or. &
            &records(family)%relative_tolerance > LOOSE_RELATIVE_LIMITS(family) ) &
            &error stop 'calibration produced a scientifically loose tolerance'
    enddo
end subroutine finalize_tolerances

!> Compare all fixture fields exactly for deterministic-repeatability evidence.
logical function fixtures_are_identical(left,right) result(identical)
    type(calibration_fixture), intent(in) :: left, right

    identical = left%box == right%box .and. left%fixture_id == right%fixture_id
    identical = identical .and. all(left%exact_pose == right%exact_pose)
    identical = identical .and. all(left%nonstationary_pose == right%nonstationary_pose)
    identical = identical .and. all(left%ordinary_ctf == right%ordinary_ctf)
    identical = identical .and. all(left%attenuated_ctf == right%attenuated_ctf)
    identical = identical .and. all(left%constant_sigma2 == right%constant_sigma2)
    identical = identical .and. all(left%varying_sigma2 == right%varying_sigma2)
    identical = identical .and. all(left%volume == right%volume)
end function fixtures_are_identical

!> Write the complete pre-acceptance calibration package.
subroutine write_calibration_artifacts(output_directory,records,fixtures)
    character(len=*), intent(in) :: output_directory
    type(tolerance_record), intent(in) :: records(N_TOLERANCE_FAMILIES)
    type(calibration_fixture), intent(in) :: fixtures(:)
    integer :: family, fixture_index, observation, unit

    open(newunit=unit,file=trim(output_directory)//'/calibration_raw_errors.tsv', &
        &status='replace',action='write')
    write(unit,'(7a)') 'family',achar(9),'observation',achar(9),'absolute_error',achar(9), &
        &'scaled_relative_error'
    do family = 1, N_TOLERANCE_FAMILIES
        do observation = 1, records(family)%observations
            write(unit,'(a,a,i0,a,es24.16,a,es24.16)') trim(records(family)%name),achar(9), &
                &observation,achar(9),records(family)%absolute_errors(observation),achar(9), &
                &records(family)%scaled_relative_errors(observation)
        enddo
    enddo
    close(unit)

    open(newunit=unit,file=trim(output_directory)//'/frozen_tolerances.tsv', &
        &status='replace',action='write')
    write(unit,'(19a)') 'family',achar(9),'abs_floor',achar(9),'rel_floor',achar(9), &
        &'rel_scale_floor',achar(9),'max_absolute',achar(9),'max_scaled_relative',achar(9), &
        &'abs_tol',achar(9),'rel_tol',achar(9),'safety_factor',achar(9),'justification'
    do family = 1, N_TOLERANCE_FAMILIES
        write(unit,'(a,8(a,es24.16),a)') trim(records(family)%name),achar(9), &
            &records(family)%absolute_floor,achar(9),records(family)%relative_floor,achar(9), &
            &records(family)%relative_scale_floor,achar(9),records(family)%maximum_absolute_error, &
            &achar(9),records(family)%maximum_scaled_relative_error, &
            &achar(9),records(family)%absolute_tolerance,achar(9),records(family)%relative_tolerance, &
            &achar(9),CALIBRATION_SAFETY_FACTOR, &
            &achar(9)//trim(tolerance_justification(family))
    enddo
    close(unit)

    open(newunit=unit,file=trim(output_directory)//'/fixture_manifest.tsv', &
        &status='replace',action='write')
    write(unit,'(15a)') 'role',achar(9),'box',achar(9),'fixture_id',achar(9),'volume_sum',achar(9), &
        &'constant_sigma_sum',achar(9),'varying_sigma_sum',achar(9),'pose_l1',achar(9),'ctf_l1'
    do fixture_index = 1, size(fixtures)
        write(unit,'(a,a,i0,a,i0,5(a,es24.16))') 'calibration',achar(9),fixtures(fixture_index)%box, &
            &achar(9),fixtures(fixture_index)%fixture_id,achar(9),sum(fixtures(fixture_index)%volume), &
            &achar(9),sum(fixtures(fixture_index)%constant_sigma2), &
            &achar(9),sum(fixtures(fixture_index)%varying_sigma2), &
            &achar(9),sum(abs(fixtures(fixture_index)%exact_pose))+ &
            &sum(abs(fixtures(fixture_index)%nonstationary_pose)), &
            &achar(9),sum(abs(fixtures(fixture_index)%ordinary_ctf))+ &
            &sum(abs(fixtures(fixture_index)%attenuated_ctf))
    enddo
    do fixture_index = 1, size(ACCEPTANCE_BOXES)
        write(unit,'(a,a,i0,a,i0,10a)') 'reserved_acceptance_not_sampled',achar(9), &
            &ACCEPTANCE_BOXES(fixture_index),achar(9),4900+ACCEPTANCE_BOXES(fixture_index), &
            &achar(9),'NA',achar(9),'NA',achar(9),'NA',achar(9),'NA',achar(9),'NA'
    enddo
    close(unit)

    open(newunit=unit,file=trim(output_directory)//'/PRE_ACCEPTANCE_TOLERANCES.FROZEN', &
        &status='new',action='write')
    write(unit,'(a)') 'status=FROZEN_PRE_ACCEPTANCE'
    write(unit,'(a,es24.16)') 'safety_factor=',CALIBRATION_SAFETY_FACTOR
    write(unit,'(a)') 'calibration_boxes=8,12'
    write(unit,'(a)') 'acceptance_boxes_reserved_not_sampled=10,16'
    write(unit,'(a)') 'rule=acceptance_failure_reopens_protocol_or_implementation_never_threshold_in_place'
    close(unit)
end subroutine write_calibration_artifacts

!> Explain the numerical contract represented by one tolerance family.
pure function tolerance_justification(family) result(justification)
    integer, intent(in) :: family
    character(len=96) :: justification

    select case(family)
    case(1)
        justification = 'componentwise scalar vector matrix and complex accumulation'
    case(2)
        justification = 'scaled five-by-five solve step reduction ratio and backward error'
    case(3)
        justification = 'fixed-cell residual Jacobian and objective-gradient differences'
    case(4)
        justification = 'analytic object versus normalized brute-force DFT'
    case(5)
        justification = 'executed deapodized gather versus brute-force DFT'
    case(6)
        justification = 'independent slow KB gather versus executed packed gather'
    case(7)
        justification = 'line-sum projection and normalized 2-D FFT versus direct DFT'
    case default
        justification = 'invalid tolerance family'
    end select
end function tolerance_justification

end module pose_cont_refinement_calibration_helpers
