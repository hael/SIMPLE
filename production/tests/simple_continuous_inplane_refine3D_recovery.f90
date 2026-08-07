module continuous_inplane_refine3D_recovery_test
use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
implicit none
private
public :: run_synthetic_recovery_contract

contains

subroutine run_synthetic_recovery_contract()
use simple_pftc_srch_api
use simple_builder, only: builder
use simple_matcher_smpl_and_lplims, only: set_bp_range3D
use simple_projector_pft, only: fproject_polar
use simple_pftc_shsrch_grad, only: pftc_shsrch_grad
implicit none

type(parameters) :: params
type(cmdline) :: cline
type(builder) :: build
type(ori) :: reference_ori, particle_ori
type(pftc_shsrch_grad) :: seed_search, joint_search
real(sp), allocatable, target :: sigma2_noise(:,:)
real(sp), allocatable :: raw_losses(:)
character(len=4096) :: vol_file
real :: mskdiam, smpd, lp, truth_angle, truth_shift(2)
real :: seed_limits(2,2), seed_limits_init(2,2), seed_result(3)
real :: joint_limits(3,2), joint_result(3)
real(dp) :: grid_coord, joint_coord, grid_loss, joint_loss, initial_loss
real(dp) :: grid_angle, joint_angle, grid_angle_error, joint_angle_error
real(dp) :: seed_shift(2), expected_shift(2), recovered_shift(2)
real(dp) :: shift_rms, objective_tolerance
real(dp) :: theta_rad, unused_gradient(3)
integer :: projection_index, grid_index, seed_index, joint_index, nrots
logical :: evaluation_valid, improved, volume_exists

call parse_recovery_arguments(vol_file, mskdiam, smpd, lp, truth_angle, truth_shift)
inquire(file=trim(vol_file), exist=volume_exists)
if( .not. volume_exists ) error stop 'synthetic refine3D recovery volume does not exist'

call cline%set('vol1', trim(vol_file))
call cline%set('mskdiam', mskdiam)
call cline%set('smpd', smpd)
call cline%set('lp', lp)
call cline%set('nptcls', 1.)
call cline%set('trs', 5.)
call cline%set('ctf', 'no')
call cline%set('objfun', 'euclid')
call cline%check
call build%init_params_and_build_strategy3D_tbox(cline, params)
call set_bp_range3D(params, build, cline)
call build%pftc%new(params, 1, [1,1], params%kfromto)
call build%vol%read(params%vols(1))
call build%vol%mask3D_soft(params%msk)
call build%vol%fft()
call build%vol%expand_cmat(params%box)

! Hold the 3D state and projection direction fixed. Only the two image-plane
! shifts and the in-plane Euler angle differ between reference and particle.
projection_index = 1
call build%eulspace%get_ori(projection_index, reference_ori)
call reference_ori%e3set(0.)
particle_ori = reference_ori
call particle_ori%e3set(truth_angle)
call fproject_polar(build%vol, 1, particle_ori, build%pftc, iseven=.true.)
call build%pftc%cp_even_ref2ptcl(1, 1)
call fproject_polar(build%vol, 1, reference_ori, build%pftc, iseven=.true.)
call build%pftc%set_eo(1, .true.)
call build%pftc%shift_ptcl(1, truth_shift)
call build%pftc%memoize_refs
call build%pftc%memoize_ptcls
allocate(sigma2_noise(params%kfromto(1):params%kfromto(2),1), source=1.0_sp)
call build%pftc%assign_sigma2_noise(sigma2_noise)
call build%pftc%memoize_sqsum_ptcl(1)

nrots = build%pftc%get_nrots()
allocate(raw_losses(nrots))
call build%pftc%gen_raw_euclid_vals(1, 1, [0._sp,0._sp], raw_losses)
grid_index = minloc(raw_losses, dim=1)

! Reproduce refine3D's production pre-selection seed. The selected reference
! first receives an alternating discrete-angle/continuous-shift search, and the
! resulting native-frame xy_first is then used to score all grid candidates.
seed_limits(:,1) = -params%trs
seed_limits(:,2) =  params%trs
seed_limits_init(:,1) = -SHC_INPL_TRSHWDTH
seed_limits_init(:,2) =  SHC_INPL_TRSHWDTH
call seed_search%new(build, seed_limits, lims_init=seed_limits_init, &
    &maxits=params%maxits_sh, opt_angle=.true., coarse_init=.true.)
call seed_search%set_indices(1, 1)
seed_index = grid_index
seed_result = seed_search%minimize(irot=seed_index, sh_rot=.false.)
if( seed_index < 1 ) &
    &error stop 'synthetic refine3D pre-selection shift search found no improving seed'
seed_shift = real(seed_result(2:3),dp)
call build%pftc%gen_raw_euclid_vals(1, 1, real(seed_shift,sp), raw_losses)
grid_index = minloc(raw_losses, dim=1)
grid_coord = real(grid_index,dp)
call build%pftc%gen_raw_euclid_grad_at_angle(1, 1, seed_shift, &
    &grid_coord, grid_loss, unused_gradient)

joint_limits(1:2,1) = -5.
joint_limits(1:2,2) =  5.
joint_limits(3,:) = [real(grid_index)-2., real(grid_index)+2.]
call joint_search%new_joint(build, joint_limits, 100)
call joint_search%set_indices(1, 1)
joint_index = grid_index
joint_result = joint_search%minimize_joint(joint_index, real(grid_coord), real(seed_shift), &
    &sh_rot=.false., rotind_frac=joint_coord, evaluation_valid=evaluation_valid, &
    &improved=improved, initial_cost_out=initial_loss)
if( .not. evaluation_valid ) error stop 'synthetic refine3D joint evaluation was invalid'
if( .not. improved .or. joint_index < 1 ) &
    &error stop 'synthetic refine3D joint refinement did not improve the grid seed'
recovered_shift = real(joint_result(2:3),dp)
call build%pftc%gen_raw_euclid_grad_at_angle(1, 1, recovered_shift, &
    &joint_coord, joint_loss, unused_gradient)

theta_rad = real(truth_angle,dp) * acos(-1.d0) / 180.d0
expected_shift = [cos(theta_rad)*real(truth_shift(1),dp) - &
    &sin(theta_rad)*real(truth_shift(2),dp), &
    &sin(theta_rad)*real(truth_shift(1),dp) + &
    &cos(theta_rad)*real(truth_shift(2),dp)]
grid_angle = (grid_coord - 1.d0) * real(build%pftc%get_dang(),dp)
joint_angle = (joint_coord - 1.d0) * real(build%pftc%get_dang(),dp)
grid_angle_error = periodic_angle_error(grid_angle, -real(truth_angle,dp))
joint_angle_error = periodic_angle_error(joint_angle, -real(truth_angle,dp))
shift_rms = sqrt(sum((recovered_shift-expected_shift)**2)/2.d0)
objective_tolerance = 1.d-8 * max(1.d0, abs(grid_loss), abs(joint_loss))

write(*,'(a,i0)') 'REFINE3D_INPL_CONT_RECOVERY PROJECTION_INDEX: ', projection_index
write(*,'(a,3f12.6)') 'REFINE3D_INPL_CONT_RECOVERY TRUTH/GRID/JOINT_ANGLE_DEG: ', &
    truth_angle, grid_angle, joint_angle
write(*,'(a,2es16.8)') 'REFINE3D_INPL_CONT_RECOVERY GRID/JOINT_ANGLE_ERROR_DEG: ', &
    grid_angle_error, joint_angle_error
write(*,'(a,2es16.8)') 'REFINE3D_INPL_CONT_RECOVERY GRID/JOINT_RAW_LOSS: ', &
    grid_loss, joint_loss
write(*,'(a,2f12.6)') 'REFINE3D_INPL_CONT_RECOVERY EXPECTED_SHIFT_XY: ', expected_shift
write(*,'(a,2f12.6)') 'REFINE3D_INPL_CONT_RECOVERY GRID_SEED_SHIFT_XY: ', seed_shift
write(*,'(a,2f12.6)') 'REFINE3D_INPL_CONT_RECOVERY RECOVERED_SHIFT_XY: ', recovered_shift
write(*,'(a,es16.8)') 'REFINE3D_INPL_CONT_RECOVERY SHIFT_RMS_PIX: ', shift_rms

if( .not. all(ieee_is_finite([grid_loss, joint_loss, initial_loss, &
        &grid_angle_error, joint_angle_error, shift_rms])) ) &
    &error stop 'synthetic refine3D recovery produced non-finite results'
if( abs(initial_loss-grid_loss) > 16.d0*objective_tolerance ) &
    &error stop 'synthetic refine3D joint seed does not match the grid objective'
if( joint_loss > grid_loss + objective_tolerance ) &
    &error stop 'synthetic refine3D joint refinement worsened the accepted objective'
if( joint_angle_error >= grid_angle_error - 1.d-4 ) &
    &error stop 'synthetic refine3D joint refinement did not improve angular recovery'
if( shift_rms > 0.25d0 ) &
    &error stop 'synthetic refine3D joint refinement did not recover the known shift'

call joint_search%kill
call seed_search%kill
call build%kill_strategy3D_tbox
call build%kill_general_tbox
deallocate(raw_losses, sigma2_noise)
write(*,'(a)') 'REFINE3D_INPL_CONT_RECOVERY FIXED_PROJECTION/ANGLE/SHIFT/OBJECTIVE: PASS'
end subroutine run_synthetic_recovery_contract

subroutine parse_recovery_arguments(vol_file, mskdiam, smpd, lp, truth_angle, truth_shift)
character(len=*), intent(out) :: vol_file
real, intent(out) :: mskdiam, smpd, lp, truth_angle, truth_shift(2)
character(len=4096) :: argument, key, value
integer :: iarg, separator

vol_file = ''
mskdiam = 120.
smpd = 1.3
lp = 8.
truth_angle = 37.
truth_shift = [2., -1.5]
do iarg = 1, command_argument_count()
    call get_command_argument(iarg, argument)
    separator = index(argument, '=')
    if( separator <= 1 ) cycle
    key = adjustl(argument(:separator-1))
    value = adjustl(argument(separator+1:))
    select case(trim(key))
    case('vol1')
        vol_file = trim(value)
    case('mskdiam')
        read(value,*) mskdiam
    case('smpd')
        read(value,*) smpd
    case('lp')
        read(value,*) lp
    case('angerr')
        read(value,*) truth_angle
    case('xsh')
        read(value,*) truth_shift(1)
    case('ysh')
        read(value,*) truth_shift(2)
    case default
        continue
    end select
enddo
if( len_trim(vol_file) == 0 )then
    write(*,'(a)') 'synthetic_recovery requires vol1=/absolute/path/reference_volume.mrc'
    stop 2
endif
if( mskdiam <= 0. .or. smpd <= 0. .or. lp <= 2.*smpd ) &
    &error stop 'invalid synthetic refine3D recovery sampling parameters'
end subroutine parse_recovery_arguments

pure function periodic_angle_error(angle, target) result(error)
use simple_core_module_api, only: dp
real(dp), intent(in) :: angle, target
real(dp) :: error
error = abs(modulo(angle-target+180.d0,360.d0)-180.d0)
end function periodic_angle_error

end module continuous_inplane_refine3D_recovery_test
