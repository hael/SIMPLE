module pose_cont_refinement_forward_hierarchy_test
use pose_cont_refinement_calibration_helpers, only: combined_complex_passes, FROZEN_ABSOLUTE_TOLERANCES, &
    &FROZEN_RELATIVE_TOLERANCES
use pose_cont_refinement_test_helpers, only: assert_true
use simple_defs, only: dp, sp
use simple_image, only: image
use simple_cartesian_fourier, only: extract_native_fourier_plane
use simple_cartesian_pose_refiner, only: cartesian_pose_refiner
use simple_ori_utils, only: euler2m
implicit none
private
public :: run_forward_hierarchy

integer, parameter :: ANALYTIC_DFT_FAMILY = 4
integer, parameter :: EXECUTED_DFT_FAMILY = 5
integer, parameter :: SLOW_GATHER_FAMILY = 6
integer, parameter :: FINITE_PROJECTION_FAMILY = 7
integer, parameter :: N_ANALYTIC_OBJECTS = 4
integer, parameter :: N_ROTATIONS = 4
integer, parameter :: N_HOLDOUT_VARIANTS = 2
integer, parameter :: HOLDOUT_BOXES(2) = [14,18]
real(dp), parameter :: RELATIVE_SCALE_FLOORS(7) = [1.e-12_dp,1._dp,1.e-8_dp,1.e-12_dp, &
    &1.e-3_dp,1.e-10_dp,1.e-3_dp]
integer, parameter :: N_COMPARISONS = size(HOLDOUT_BOXES)*(N_ANALYTIC_OBJECTS*N_ROTATIONS+ &
    &N_HOLDOUT_VARIANTS*N_ROTATIONS*3)

type :: forward_holdout_fixture
    integer :: box = 0
    integer :: variant = 0
    integer :: fixture_id = 0
    real(dp), allocatable :: volume(:,:,:)
end type forward_holdout_fixture

type :: fixture_summary
    integer :: box = 0
    integer :: variant = 0
    integer :: fixture_id = 0
    real(dp) :: volume_sum = 0._dp
    real(dp) :: volume_l1 = 0._dp
end type fixture_summary

type :: comparison_metrics
    character(len=40) :: name = ''
    integer :: family = 0
    integer :: box = 0
    integer :: variant = 0
    integer :: rotation = 0
    integer :: samples = 0
    real(dp) :: maximum_absolute_error = 0._dp
    real(dp) :: maximum_scaled_relative_error = 0._dp
    real(dp) :: clipping_fraction = 0._dp
    real(dp) :: original_volume_error = 0._dp
    real(dp) :: interpolation_error = 0._dp
    real(dp) :: support_margin = 0._dp
    integer :: stencil_switches = 0
    logical :: passed = .false.
    complex(dp), allocatable :: actual(:)
    complex(dp), allocatable :: expected(:)
end type comparison_metrics

contains

!> Verify all four independent forward contracts on frozen fresh holdouts.
subroutine run_forward_hierarchy()
    type(forward_holdout_fixture) :: fixture
    type(fixture_summary) :: fixtures(size(HOLDOUT_BOXES)*N_HOLDOUT_VARIANTS)
    type(comparison_metrics) :: metrics(N_COMPARISONS)
    character(len=:), allocatable :: evidence_directory
    integer :: box_index, fixture_count, metric_count, variant

    evidence_directory = required_argument('evidence_dir')
    fixture_count = 0
    metric_count = 0
    do box_index = 1, size(HOLDOUT_BOXES)
        call verify_analytic_objects(HOLDOUT_BOXES(box_index),metrics,metric_count)
        do variant = 1, N_HOLDOUT_VARIANTS
            call build_forward_holdout_fixture(HOLDOUT_BOXES(box_index),variant,fixture)
            fixture_count = fixture_count+1
            fixtures(fixture_count)%box = fixture%box
            fixtures(fixture_count)%variant = fixture%variant
            fixtures(fixture_count)%fixture_id = fixture%fixture_id
            fixtures(fixture_count)%volume_sum = sum(fixture%volume)
            fixtures(fixture_count)%volume_l1 = sum(abs(fixture%volume))
            call verify_volume_hierarchy(fixture,metrics,metric_count)
        enddo
    enddo
    call assert_true(metric_count == N_COMPARISONS,'forward hierarchy comparison count changed')
    call write_forward_artifacts(evidence_directory,metrics(:metric_count),fixtures(:fixture_count))
    call print_forward_metrics(metrics(:metric_count))
    call assert_true(all(metrics(:metric_count)%passed), &
        &'forward hierarchy comparison exceeded its frozen tolerance')
    write(*,'(a)') 'POSE_CONT_REFINEMENT_FORWARD_HIERARCHY domain: full_redundant_disk'
    write(*,'(a)') 'POSE_CONT_REFINEMENT_FORWARD_HIERARCHY holdout boxes: 14 18'
    write(*,'(a)') 'POSE_CONT_REFINEMENT_FORWARD_HIERARCHY variants: 2 rotations: 4'
    write(*,'(a)') 'POSE_CONT_REFINEMENT_FORWARD_HIERARCHY: PASS'
end subroutine run_forward_hierarchy

!> Compare analytic point-object transforms with a full voxel-loop DFT.
subroutine verify_analytic_objects(box,metrics,metric_count)
    integer, intent(in) :: box
    type(comparison_metrics), intent(inout) :: metrics(:)
    integer, intent(inout) :: metric_count
    complex(dp), allocatable :: analytic(:), brute(:)
    real(dp) :: rotations(3,3,N_ROTATIONS)
    integer :: active_samples, object_index, rotation_index

    call build_rotations(rotations)
    active_samples = full_disk_sample_count(box)
    allocate(analytic(active_samples),brute(active_samples))
    do object_index = 1, N_ANALYTIC_OBJECTS
        do rotation_index = 1, N_ROTATIONS
            call point_object_planes(box,object_index,rotations(:,:,rotation_index),analytic,brute)
            metric_count = metric_count+1
            metrics(metric_count)%name = analytic_label(object_index,rotation_index)
            metrics(metric_count)%box = box
            metrics(metric_count)%rotation = rotation_index
            call require_comparison(ANALYTIC_DFT_FAMILY,brute,analytic,metrics(metric_count))
        enddo
    enddo
end subroutine verify_analytic_objects

!> Compare executed, slow, direct-DFT, and finite-projection paths.
subroutine verify_volume_hierarchy(fixture,metrics,metric_count)
    type(forward_holdout_fixture), intent(in) :: fixture
    type(comparison_metrics), intent(inout) :: metrics(:)
    integer, intent(inout) :: metric_count
    type(cartesian_pose_refiner) :: workspace
    complex(dp), allocatable :: direct(:), executed(:), fast(:), slow(:), projected(:)
    complex(dp), allocatable :: matched_direct(:), original_direct(:)
    complex, allocatable :: zero_plane(:,:), executed_plane(:,:), projection_plane(:,:)
    real, allocatable :: physical_volume(:,:,:)
    real(dp), allocatable :: rotated_volume(:,:,:)
    real(dp) :: rotations(3,3,N_ROTATIONS), ignored_objective, clipping_fraction
    real(dp) :: interpolation_error, support_margin
    integer :: active_samples, lims2(2,2), rotation_index, stencil_switches

    call build_rotations(rotations)
    allocate(physical_volume,source=real(fixture%volume,sp))
    call workspace%new(physical_volume)
    call workspace%set_shell_range([0,fixture%box/2])
    lims2 = workspace%get_lims2()
    allocate(zero_plane(lims2(1,1):lims2(1,2),lims2(2,1):lims2(2,2)),source=cmplx(0.,0.))
    allocate(executed_plane,mold=zero_plane)
    active_samples = full_disk_sample_count(fixture%box)
    allocate(direct(active_samples),executed(active_samples),fast(active_samples), &
        &slow(active_samples),projected(active_samples),matched_direct(active_samples), &
        &original_direct(active_samples))
    allocate(rotated_volume(fixture%box,fixture%box,fixture%box))

    do rotation_index = 1, N_ROTATIONS
        call workspace%shift_residual(rotations(:,:,rotation_index),[0._dp,0._dp], &
            &zero_plane,executed_plane,ignored_objective)
        call collect_direct_and_gather_planes(fixture%volume,workspace, &
            &rotations(:,:,rotation_index),executed_plane,direct,executed,fast,slow,support_margin)

        metric_count = metric_count+1
        metrics(metric_count)%name = rotation_label('executed_dft',rotation_index)
        metrics(metric_count)%box = fixture%box
        metrics(metric_count)%variant = fixture%variant
        metrics(metric_count)%rotation = rotation_index
        metrics(metric_count)%support_margin = support_margin
        call require_comparison(EXECUTED_DFT_FAMILY,executed,direct,metrics(metric_count))

        metric_count = metric_count+1
        metrics(metric_count)%name = rotation_label('slow_gather',rotation_index)
        metrics(metric_count)%box = fixture%box
        metrics(metric_count)%variant = fixture%variant
        metrics(metric_count)%rotation = rotation_index
        metrics(metric_count)%support_margin = support_margin
        call require_comparison(SLOW_GATHER_FAMILY,fast,slow,metrics(metric_count))
        stencil_switches = workspace%count_stencil_switches(rotations(:,:,1), &
            &rotations(:,:,rotation_index),[0,fixture%box/2])
        metrics(metric_count)%stencil_switches = stencil_switches

        call rotate_finite_volume(fixture,rotations(:,:,rotation_index),rotated_volume, &
            &clipping_fraction,interpolation_error)
        call finite_projection_fft(rotated_volume,lims2,projection_plane)
        call flatten_full_disk(projection_plane,fixture%box,projected)
        call collect_direct_planes(fixture%volume,rotated_volume,rotations(:,:,rotation_index), &
            &original_direct,matched_direct)
        metric_count = metric_count+1
        metrics(metric_count)%name = rotation_label('finite_projection',rotation_index)
        metrics(metric_count)%box = fixture%box
        metrics(metric_count)%variant = fixture%variant
        metrics(metric_count)%rotation = rotation_index
        metrics(metric_count)%clipping_fraction = clipping_fraction
        metrics(metric_count)%original_volume_error = maxval(abs(matched_direct-original_direct))
        metrics(metric_count)%interpolation_error = interpolation_error
        metrics(metric_count)%support_margin = support_margin
        metrics(metric_count)%stencil_switches = stencil_switches
        call require_comparison(FINITE_PROJECTION_FAMILY,projected,matched_direct,metrics(metric_count))
        deallocate(projection_plane)
    enddo

    deallocate(physical_volume,zero_plane,executed_plane,direct,executed,fast,slow,projected)
    deallocate(matched_direct,original_direct,rotated_volume)
    call workspace%kill
end subroutine verify_volume_hierarchy

!> Collect one full-disk direct DFT and both prepared-reference gathers.
subroutine collect_direct_and_gather_planes(volume,workspace,rotation,executed_plane, &
    &direct,executed,fast,slow,support_margin)
    real(dp), intent(in) :: volume(:,:,:), rotation(3,3)
    type(cartesian_pose_refiner), intent(in) :: workspace
    complex, intent(in) :: executed_plane(-size(volume,1)/2:,-size(volume,1)/2:)
    complex(dp), intent(out) :: direct(:), executed(:), fast(:), slow(:)
    real(dp), intent(out) :: support_margin
    complex :: fast_value, ignored_gradient(3)
    real(sp) :: loc(3), margin(3)
    integer :: box, h, k, sample

    box = size(volume,1)
    sample = 0
    support_margin = huge(0._dp)
    do k = -box/2, box/2
        do h = -box/2, box/2
            if( h*h+k*k > (box/2)**2 ) cycle
            sample = sample+1
            loc = 2._sp*real(matmul(real([h,k,0],dp),rotation),sp)
            call workspace%sample_with_grad(loc,fast_value,ignored_gradient,margin)
            call workspace%sample_slow_test(loc,slow(sample))
            fast(sample) = cmplx(fast_value,kind=dp)
            executed(sample) = cmplx(executed_plane(h,k),kind=dp)
            direct(sample) = direct_volume_dft(volume,rotation,h,k)
            support_margin = min(support_margin,real(minval(margin),dp))
        enddo
    enddo
    call assert_true(sample == size(direct),'forward hierarchy did not traverse the full disk')
end subroutine collect_direct_and_gather_planes

!> Build a Phase-8-only physical volume that calibration cannot construct.
subroutine build_forward_holdout_fixture(box,variant,fixture)
    integer, intent(in) :: box, variant
    type(forward_holdout_fixture), intent(out) :: fixture
    real(dp) :: coordinate(3)
    integer :: i, j, k

    if( .not. any(box == HOLDOUT_BOXES) ) error stop 'forward holdout requested an unfrozen box'
    if( variant < 1 .or. variant > N_HOLDOUT_VARIANTS ) &
        &error stop 'forward holdout requested an invalid volume variant'
    fixture%box = box
    fixture%variant = variant
    fixture%fixture_id = 8000+100*variant+box
    allocate(fixture%volume(box,box,box))
    do k = 1, box
        do j = 1, box
            do i = 1, box
                coordinate = real([i,j,k]-(box/2+1),dp)
                fixture%volume(i,j,k) = holdout_field(box,variant,coordinate)
            enddo
        enddo
    enddo
end subroutine build_forward_holdout_fixture

!> Evaluate one predeclared asymmetric analytic holdout field.
function holdout_field(box,variant,coordinate) result(value)
    integer, intent(in) :: box, variant
    real(dp), intent(in) :: coordinate(3)
    real(dp) :: value
    real(dp) :: centres(3,5), amplitudes(5), widths(5), delta(3), radius2
    integer :: blob

    if( variant == 1 )then
        centres = reshape([ &
            &-0.27_dp*box, 0.08_dp*box, 0.19_dp*box, &
            & 0.16_dp*box,-0.29_dp*box,-0.11_dp*box, &
            & 0.31_dp*box, 0.21_dp*box,-0.24_dp*box, &
            &-0.09_dp*box,-0.18_dp*box, 0.30_dp*box, &
            & 0.04_dp*box, 0.33_dp*box, 0.06_dp*box],[3,5])
        amplitudes = [0.94_dp,-0.47_dp,0.36_dp,0.22_dp,-0.15_dp]
        widths = real(box,dp)*[0.086_dp,0.121_dp,0.103_dp,0.147_dp,0.074_dp]
    else
        centres = reshape([ &
            & 0.24_dp*box, 0.15_dp*box,-0.28_dp*box, &
            &-0.21_dp*box,-0.32_dp*box, 0.07_dp*box, &
            & 0.09_dp*box, 0.27_dp*box, 0.25_dp*box, &
            &-0.30_dp*box, 0.03_dp*box,-0.17_dp*box, &
            & 0.29_dp*box,-0.14_dp*box, 0.12_dp*box],[3,5])
        amplitudes = [0.82_dp,0.59_dp,-0.43_dp,0.31_dp,-0.19_dp]
        widths = real(box,dp)*[0.112_dp,0.079_dp,0.138_dp,0.091_dp,0.126_dp]
    endif
    value = 0._dp
    do blob = 1, 5
        delta = coordinate-centres(:,blob)
        radius2 = 1.07_dp*delta(1)**2+0.89_dp*delta(2)**2+1.16_dp*delta(3)**2+ &
            &0.05_dp*delta(1)*delta(2)-0.04_dp*delta(2)*delta(3)
        value = value+amplitudes(blob)*exp(-radius2/(2._dp*widths(blob)**2))
    enddo
end function holdout_field

!> Build analytic and brute-force full-disk planes for one point fixture.
subroutine point_object_planes(box,object_index,rotation,analytic,brute)
    integer, intent(in) :: box, object_index
    real(dp), intent(in) :: rotation(3,3)
    complex(dp), intent(out) :: analytic(:), brute(:)
    real(dp) :: amplitudes(4), angle, loc(3), twopi, voxel
    integer :: points(3,4), npoints, h, k, point, sample, x, y, z

    call point_object_definition(object_index,points,amplitudes,npoints)
    twopi = 2._dp*acos(-1._dp)
    sample = 0
    do k = -box/2, box/2
        do h = -box/2, box/2
            if( h*h+k*k > (box/2)**2 ) cycle
            sample = sample+1
            loc = matmul(real([h,k,0],dp),rotation)
            analytic(sample) = cmplx(0._dp,0._dp,kind=dp)
            do point = 1, npoints
                angle = -twopi*dot_product(loc,real(points(:,point),dp))/real(box,dp)
                analytic(sample) = analytic(sample)+amplitudes(point)* &
                    &cmplx(cos(angle),sin(angle),kind=dp)
            enddo
            brute(sample) = cmplx(0._dp,0._dp,kind=dp)
            do z = -box/2, box/2-1
                do y = -box/2, box/2-1
                    do x = -box/2, box/2-1
                        voxel = 0._dp
                        do point = 1, npoints
                            if( all([x,y,z] == points(:,point)) ) voxel = voxel+amplitudes(point)
                        enddo
                        angle = -twopi*dot_product(loc,real([x,y,z],dp))/real(box,dp)
                        brute(sample) = brute(sample)+voxel*cmplx(cos(angle),sin(angle),kind=dp)
                    enddo
                enddo
            enddo
            analytic(sample) = analytic(sample)/real(box**3,dp)
            brute(sample) = brute(sample)/real(box**3,dp)
        enddo
    enddo
end subroutine point_object_planes

!> Define centered, shifted, separated, and non-collinear point fixtures.
subroutine point_object_definition(object_index,points,amplitudes,npoints)
    integer, intent(in) :: object_index
    integer, intent(out) :: points(3,4), npoints
    real(dp), intent(out) :: amplitudes(4)

    points = 0
    amplitudes = 0._dp
    select case(object_index)
    case(1)
        npoints = 1; amplitudes(1) = 1._dp
    case(2)
        npoints = 1; points(:,1) = [2,-1,1]; amplitudes(1) = 0.73_dp
    case(3)
        npoints = 2; points(:,1) = [-2,1,-1]; points(:,2) = [2,-2,1]
        amplitudes(1:2) = [1._dp,-0.41_dp]
    case(4)
        npoints = 4
        points = reshape([0,0,0, 1,-2,1, -2,1,-1, 2,2,-2],[3,4])
        amplitudes = [1._dp,0.61_dp,-0.37_dp,0.23_dp]
    case default
        error stop 'invalid analytic point fixture'
    end select
end subroutine point_object_definition

!> Evaluate the normalized direct 3-D DFT at a rotated plane coordinate.
function direct_volume_dft(volume,rotation,h,k) result(value)
    real(dp), intent(in) :: volume(:,:,:), rotation(3,3)
    integer, intent(in) :: h, k
    complex(dp) :: value
    real(dp) :: angle, loc(3), twopi
    integer :: box, i, j, l, x, y, z

    box = size(volume,1)
    twopi = 2._dp*acos(-1._dp)
    loc = matmul(real([h,k,0],dp),rotation)
    value = cmplx(0._dp,0._dp,kind=dp)
    do l = 1, box
        z = l-(box/2+1)
        do j = 1, box
            y = j-(box/2+1)
            do i = 1, box
                x = i-(box/2+1)
                angle = -twopi*dot_product(loc,real([x,y,z],dp))/real(box,dp)
                value = value+volume(i,j,l)*cmplx(cos(angle),sin(angle),kind=dp)
            enddo
        enddo
    enddo
    value = value/real(box**3,dp)
end function direct_volume_dft

!> Evaluate original-frequency and matched finite-object DFT planes.
subroutine collect_direct_planes(original,rotated,rotation,original_direct,matched_direct)
    real(dp), intent(in) :: original(:,:,:), rotated(:,:,:), rotation(3,3)
    complex(dp), intent(out) :: original_direct(:), matched_direct(:)
    real(dp) :: identity(3,3)
    integer :: box, h, k, sample

    identity = 0._dp
    identity(1,1) = 1._dp
    identity(2,2) = 1._dp
    identity(3,3) = 1._dp
    box = size(original,1)
    sample = 0
    do k = -box/2, box/2
        do h = -box/2, box/2
            if( h*h+k*k > (box/2)**2 ) cycle
            sample = sample+1
            original_direct(sample) = direct_volume_dft(original,rotation,h,k)
            matched_direct(sample) = direct_volume_dft(rotated,identity,h,k)
        enddo
    enddo
    call assert_true(sample == size(original_direct),'direct holdout plane is incomplete')
end subroutine collect_direct_planes

!> Resample one finite holdout object and retain interpolation diagnostics.
subroutine rotate_finite_volume(fixture,rotation,rotated,clipping_fraction,interpolation_error)
    type(forward_holdout_fixture), intent(in) :: fixture
    real(dp), intent(in) :: rotation(3,3)
    real(dp), intent(out) :: rotated(:,:,:), clipping_fraction, interpolation_error
    real(dp) :: coordinate(3), exact_value
    logical :: clipped_sample
    integer :: box, clipped, i, j, k, u, v, w

    box = fixture%box
    clipped = 0
    interpolation_error = 0._dp
    do k = 1, box
        w = k-(box/2+1)
        do j = 1, box
            v = j-(box/2+1)
            do i = 1, box
                u = i-(box/2+1)
                coordinate = matmul(real([u,v,w],dp),rotation)
                call trilinear_finite_sample(fixture%volume,coordinate,rotated(i,j,k), &
                    &clipped,clipped_sample)
                if( .not. clipped_sample )then
                    exact_value = holdout_field(box,fixture%variant,coordinate)
                    interpolation_error = max(interpolation_error,abs(rotated(i,j,k)-exact_value))
                endif
            enddo
        enddo
    enddo
    clipping_fraction = real(clipped,dp)/real(box**3,dp)
end subroutine rotate_finite_volume

!> Form the line-sum projection and normalized FFT of one finite object.
subroutine finite_projection_fft(volume,lims2,plane)
    real(dp), intent(in) :: volume(:,:,:)
    integer, intent(in) :: lims2(2,2)
    complex, allocatable, intent(out) :: plane(:,:)
    type(image) :: projection_image
    real :: projection(size(volume,1),size(volume,1),1)
    integer :: box, i, j

    box = size(volume,1)
    do j = 1, box
        do i = 1, box
            projection(i,j,1) = real(sum(volume(i,j,:)),sp)
        enddo
    enddo
    call projection_image%new([box,box,1],1.0)
    call projection_image%set_rmat(projection,.false.)
    call projection_image%fft()
    plane = extract_native_fourier_plane(projection_image,lims2,(box/2)**2)/real(box,sp)
    call projection_image%kill
end subroutine finite_projection_fft

!> Sample a centered finite box with zero outside and trilinear interpolation.
subroutine trilinear_finite_sample(volume,coordinate,value,clipped,clipped_sample)
    real(dp), intent(in) :: volume(:,:,:), coordinate(3)
    real(dp), intent(out) :: value
    integer, intent(inout) :: clipped
    logical, intent(out) :: clipped_sample
    real(dp) :: fraction(3), weight
    integer :: box, lower(3), upper(3), ix, iy, iz, xindex, yindex, zindex

    box = size(volume,1)
    lower = floor(coordinate)
    fraction = coordinate-real(lower,dp)
    upper = lower+1
    where(abs(fraction) <= 8._dp*epsilon(1._dp)) upper = lower
    if( any(lower < -box/2) .or. any(upper > box/2-1) )then
        value = 0._dp
        clipped = clipped+1
        clipped_sample = .true.
        return
    endif
    clipped_sample = .false.
    value = 0._dp
    do iz = 0, 1
        zindex = merge(lower(3),upper(3),iz == 0)+box/2+1
        do iy = 0, 1
            yindex = merge(lower(2),upper(2),iy == 0)+box/2+1
            do ix = 0, 1
                xindex = merge(lower(1),upper(1),ix == 0)+box/2+1
                weight = merge(1._dp-fraction(1),fraction(1),ix == 0)* &
                    &merge(1._dp-fraction(2),fraction(2),iy == 0)* &
                    &merge(1._dp-fraction(3),fraction(3),iz == 0)
                value = value+weight*volume(xindex,yindex,zindex)
            enddo
        enddo
    enddo
end subroutine trilinear_finite_sample

!> Flatten the explicitly redundant full disk without removing Nyquist mates.
subroutine flatten_full_disk(plane,box,values)
    integer, intent(in) :: box
    complex, intent(in) :: plane(-box/2:,-box/2:)
    complex(dp), intent(out) :: values(:)
    integer :: h, k, sample

    sample = 0
    do k = -box/2, box/2
        do h = -box/2, box/2
            if( h*h+k*k > (box/2)**2 ) cycle
            sample = sample+1
            values(sample) = cmplx(plane(h,k),kind=dp)
        enddo
    enddo
    call assert_true(sample == size(values),'plane flattening omitted full-disk samples')
end subroutine flatten_full_disk

!> Apply one frozen comparison family and retain maximum errors.
subroutine require_comparison(family,actual,expected,metrics)
    integer, intent(in) :: family
    complex(dp), intent(in) :: actual(:), expected(:)
    type(comparison_metrics), intent(inout) :: metrics
    real(dp) :: absolute_error, scale
    integer :: sample

    metrics%family = family
    metrics%samples = size(actual)
    allocate(metrics%actual,source=actual)
    allocate(metrics%expected,source=expected)
    do sample = 1, size(actual)
        absolute_error = abs(actual(sample)-expected(sample))
        scale = max(RELATIVE_SCALE_FLOORS(family),abs(actual(sample)),abs(expected(sample)))
        metrics%maximum_absolute_error = max(metrics%maximum_absolute_error,absolute_error)
        metrics%maximum_scaled_relative_error = max( &
            &metrics%maximum_scaled_relative_error,absolute_error/scale)
    enddo
    metrics%passed = combined_complex_passes(actual,expected, &
        &FROZEN_ABSOLUTE_TOLERANCES(family),FROZEN_RELATIVE_TOLERANCES(family))
    if( .not. metrics%passed )then
        write(*,'(a,1x,a,1x,i0,1x,i0,5(1x,es14.6))') &
            &'POSE_CONT_REFINEMENT_FORWARD_FAILURE',trim(metrics%name),family, &
            &metrics%samples,metrics%maximum_absolute_error, &
            &metrics%maximum_scaled_relative_error,FROZEN_ABSOLUTE_TOLERANCES(family), &
            &FROZEN_RELATIVE_TOLERANCES(family),metrics%clipping_fraction
    endif
end subroutine require_comparison

!> Write aggregate, raw, fixture, and frozen-contract evidence.
subroutine write_forward_artifacts(output_directory,metrics,fixtures)
    character(len=*), intent(in) :: output_directory
    type(comparison_metrics), intent(in) :: metrics(:)
    type(fixture_summary), intent(in) :: fixtures(:)
    real(dp) :: absolute_error, allowed_error, scale
    integer :: h, item, k, sample, unit

    open(newunit=unit,file=trim(output_directory)//'/forward_hierarchy.tsv', &
        &status='replace',action='write')
    write(unit,'(a)') 'comparison'//achar(9)//'family'//achar(9)//'box'//achar(9)// &
        &'variant'//achar(9)//'rotation'//achar(9)//'samples'//achar(9)//'max_absolute'// &
        &achar(9)//'max_scaled_relative'//achar(9)//'frozen_absolute'//achar(9)// &
        &'frozen_relative'//achar(9)//'clipping_fraction'//achar(9)// &
        &'original_volume_error'//achar(9)//'interpolation_error'//achar(9)// &
        &'support_margin'//achar(9)//'stencil_switches'//achar(9)//'passed'
    do item = 1, size(metrics)
        write(unit,'(a,5(a,i0),8(a,es24.16),a,i0,a,l1)') trim(metrics(item)%name),achar(9), &
            &metrics(item)%family,achar(9),metrics(item)%box,achar(9),metrics(item)%variant, &
            &achar(9),metrics(item)%rotation,achar(9),metrics(item)%samples,achar(9), &
            &metrics(item)%maximum_absolute_error,achar(9), &
            &metrics(item)%maximum_scaled_relative_error,achar(9), &
            &FROZEN_ABSOLUTE_TOLERANCES(metrics(item)%family),achar(9), &
            &FROZEN_RELATIVE_TOLERANCES(metrics(item)%family),achar(9), &
            &metrics(item)%clipping_fraction,achar(9),metrics(item)%original_volume_error, &
            &achar(9),metrics(item)%interpolation_error,achar(9),metrics(item)%support_margin, &
            &achar(9),metrics(item)%stencil_switches, &
            &achar(9),metrics(item)%passed
    enddo
    close(unit)

    open(newunit=unit,file=trim(output_directory)//'/forward_hierarchy_components.tsv', &
        &status='replace',action='write')
    write(unit,'(a)') 'comparison'//achar(9)//'box'//achar(9)//'variant'//achar(9)// &
        &'rotation'//achar(9)//'h'//achar(9)//'k'//achar(9)// &
        &'actual_real'//achar(9)//'actual_imag'//achar(9)//'expected_real'//achar(9)// &
        &'expected_imag'//achar(9)//'absolute_error'//achar(9)//'scale'//achar(9)// &
        &'allowed_error'//achar(9)//'passed'
    do item = 1, size(metrics)
        sample = 0
        do k = -metrics(item)%box/2, metrics(item)%box/2
            do h = -metrics(item)%box/2, metrics(item)%box/2
                if( h*h+k*k > (metrics(item)%box/2)**2 ) cycle
                sample = sample+1
                absolute_error = abs(metrics(item)%actual(sample)-metrics(item)%expected(sample))
                scale = max(abs(metrics(item)%actual(sample)),abs(metrics(item)%expected(sample)))
                allowed_error = FROZEN_ABSOLUTE_TOLERANCES(metrics(item)%family)+ &
                    &FROZEN_RELATIVE_TOLERANCES(metrics(item)%family)*scale
                write(unit,'(a,5(a,i0),7(a,es24.16),a,l1)') trim(metrics(item)%name),achar(9), &
                    &metrics(item)%box,achar(9),metrics(item)%variant,achar(9), &
                    &metrics(item)%rotation,achar(9),h,achar(9),k, &
                    &achar(9),real(metrics(item)%actual(sample),dp),achar(9), &
                    &aimag(metrics(item)%actual(sample)),achar(9), &
                    &real(metrics(item)%expected(sample),dp),achar(9), &
                    &aimag(metrics(item)%expected(sample)),achar(9),absolute_error, &
                    &achar(9),scale,achar(9),allowed_error,achar(9),absolute_error <= allowed_error
            enddo
        enddo
        call assert_true(sample == metrics(item)%samples,'raw forward evidence sample count changed')
    enddo
    close(unit)

    open(newunit=unit,file=trim(output_directory)//'/forward_hierarchy_manifest.tsv', &
        &status='replace',action='write')
    write(unit,'(a,a,a)') 'contract',achar(9),'value'
    write(unit,'(a,a,a)') 'domain',achar(9),'full_redundant_disk'
    write(unit,'(a,a,a)') 'fresh_holdout_boxes',achar(9),'14,18'
    write(unit,'(a,a,a)') 'volume_variants',achar(9),'2'
    write(unit,'(a,a,a)') 'rotations',achar(9),'identity;11,-37,26;-34,22,-15;43,7,32'
    write(unit,'(a,a,a)') 'projection_rays',achar(9),'line_sum'
    write(unit,'(a,a,a)') 'projection_fft_scale',achar(9),'explicit_1_over_N'
    write(unit,'(a,a,a)') 'finite_projection_object',achar(9),'same_rotated_resampled_X_R'
    write(unit,'(a,a,a)') 'finite_box_boundary',achar(9),'zero_clipped_trilinear'
    write(unit,'(a,a,a)') 'acceptance_tolerances',achar(9),'phase4_forward_amendment_frozen'
    write(unit,'(a,a,a)') 'observed_diagnostic_boxes_excluded',achar(9),'10,16'
    close(unit)

    open(newunit=unit,file=trim(output_directory)//'/forward_holdout_fixtures.tsv', &
        &status='replace',action='write')
    write(unit,'(a)') 'box'//achar(9)//'variant'//achar(9)//'fixture_id'//achar(9)// &
        &'volume_sum'//achar(9)//'volume_l1'
    do item = 1, size(fixtures)
        write(unit,'(3(i0,a),es24.16,a,es24.16)') fixtures(item)%box,achar(9), &
            &fixtures(item)%variant,achar(9),fixtures(item)%fixture_id,achar(9), &
            &fixtures(item)%volume_sum,achar(9),fixtures(item)%volume_l1
    enddo
    close(unit)
end subroutine write_forward_artifacts

!> Print compact evidence for retained Oracle logs.
subroutine print_forward_metrics(metrics)
    type(comparison_metrics), intent(in) :: metrics(:)
    integer :: item

    do item = 1, size(metrics)
        write(*,'(a,1x,a,5(1x,i0),8(1x,es14.6),1x,i0,1x,l1)') &
            &'POSE_CONT_REFINEMENT_FORWARD_METRIC',trim(metrics(item)%name), &
            &metrics(item)%family,metrics(item)%box,metrics(item)%variant,metrics(item)%rotation, &
            &metrics(item)%samples, &
            &metrics(item)%maximum_absolute_error,metrics(item)%maximum_scaled_relative_error, &
            &FROZEN_ABSOLUTE_TOLERANCES(metrics(item)%family), &
            &FROZEN_RELATIVE_TOLERANCES(metrics(item)%family),metrics(item)%clipping_fraction, &
            &metrics(item)%original_volume_error,metrics(item)%interpolation_error, &
            &metrics(item)%support_margin, &
            &metrics(item)%stencil_switches,metrics(item)%passed
    enddo
end subroutine print_forward_metrics

!> Construct identity and three holdout rotations distinct from calibration.
subroutine build_rotations(rotations)
    real(dp), intent(out) :: rotations(3,3,N_ROTATIONS)
    rotations = 0._dp
    rotations(1,1,1) = 1._dp
    rotations(2,2,1) = 1._dp
    rotations(3,3,1) = 1._dp
    rotations(:,:,2) = real(euler2m([11._sp,-37._sp,26._sp]),dp)
    rotations(:,:,3) = real(euler2m([-34._sp,22._sp,-15._sp]),dp)
    rotations(:,:,4) = real(euler2m([43._sp,7._sp,32._sp]),dp)
end subroutine build_rotations

!> Return the full redundant-disk population including both Nyquist endpoints.
integer function full_disk_sample_count(box) result(count)
    integer, intent(in) :: box
    integer :: h, k
    count = 0
    do k = -box/2, box/2
        do h = -box/2, box/2
            if( h*h+k*k <= (box/2)**2 ) count = count+1
        enddo
    enddo
end function full_disk_sample_count

function analytic_label(object_index,rotation_index) result(label)
    integer, intent(in) :: object_index, rotation_index
    character(len=40) :: label
    character(len=12), parameter :: object_names(4) = [character(len=12) :: &
        &'center_delta','offset_delta','two_points','unequal_set']
    write(label,'(a,"_r",i0)') trim(object_names(object_index)),rotation_index
end function analytic_label

function rotation_label(prefix,rotation_index) result(label)
    character(len=*), intent(in) :: prefix
    integer, intent(in) :: rotation_index
    character(len=40) :: label
    write(label,'(a,"_r",i0)') trim(prefix),rotation_index
end function rotation_label

!> Read one required key=value command argument.
function required_argument(key) result(value)
    character(len=*), intent(in) :: key
    character(len=:), allocatable :: value
    character(len=4096) :: argument
    integer :: argument_index, argument_status, separator

    value = ''
    do argument_index = 1, command_argument_count()
        call get_command_argument(argument_index,argument,status=argument_status)
        if( argument_status /= 0 ) error stop 'could not read forward hierarchy argument'
        separator = index(argument,'=')
        if( separator <= 0 ) cycle
        if( trim(argument(:separator-1)) /= key ) cycle
        value = trim(argument(separator+1:))
    enddo
    if( len(value) == 0 ) error stop 'forward hierarchy requires evidence_dir=<existing-directory>'
end function required_argument

end module pose_cont_refinement_forward_hierarchy_test
