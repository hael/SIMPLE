program simple_test_continuous_inplane_rotation2D_stage1_validation
use simple_pftc_srch_api
use simple_string, only: string
use simple_builder, only: builder
use simple_matcher_smpl_and_lplims, only: set_bp_range3D
use simple_projector_pft, only: fproject_polar
use simple_polarft_calc, only: polarft_calc
use simple_pftc_shsrch_grad, only: pftc_shsrch_grad
use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
implicit none
#include "simple_local_flags.inc"

type :: recovery_result
    real(dp) :: angle_grid = 0.d0
    real(dp) :: angle_parabola = 0.d0
    real(dp) :: angle_error_grid = 0.d0
    real(dp) :: angle_error_parabola = 0.d0
    real(dp) :: angle_joint = 0.d0
    real(dp) :: angle_error_joint = 0.d0
    real(dp) :: shift_rms_grid = 0.d0
    real(dp) :: shift_rms_parabola = 0.d0
    real(dp) :: shift_rms_joint = 0.d0
    real(dp) :: shift_est_grid(2) = 0.d0
    real(dp) :: shift_est_parabola(2) = 0.d0
    real(dp) :: shift_est_joint(2) = 0.d0
    real(dp) :: shift_expected(2) = 0.d0
    logical :: shift_grid_accepted = .false.
    logical :: shift_parabola_accepted = .false.
    logical :: joint_accepted = .false.
    real(dp) :: loss_grid = 0.d0
    real(dp) :: loss_parabola = 0.d0
    real(dp) :: loss_joint = 0.d0
    real(dp) :: joint_gradient_max_error = 0.d0
    real(dp) :: joint_gradient_tol = 0.d0
    logical :: joint_gradient_finite = .false.
    logical :: joint_gradient_ok = .false.
    real(dp) :: parity_max_error = 0.d0
    real(dp) :: parity_tol = 0.d0
    logical :: parity_ok = .false.
    real(dp) :: stress_series_min = 0.d0
    real(dp) :: stress_parity_max_error = 0.d0
    real(dp) :: stress_score_max = 0.d0
    logical :: stress_ok = .false.
    logical :: finite = .false.
end type recovery_result

type(cmdline) :: args
integer, parameter :: ncases = 3
! FD noise floor: single-precision coefficient series differentiated with
! h=1e-3 leaves ~1e-3 absolute noise; a broken gradient errs at O(|grad|)
real(dp), parameter :: GRAD_FD_RTOL = 1.d-2
type(recovery_result) :: low_band(ncases), full_band(ncases), recovery(ncases)
character(len=512) :: vol_file
real :: smpd, mskdiam, lp_low, lp_full, truth_angle, shift_truth(2)
real :: truth_cases(ncases)
real(dp) :: low_alias_err, full_alias_err, low_alias_max, full_alias_max
integer :: icase

call parse_arguments(args, vol_file, mskdiam, smpd, lp_low, lp_full, truth_angle, shift_truth)
truth_cases = [0., truth_angle, 359.375]
do icase = 1, ncases
    call run_fixture(vol_file, mskdiam, smpd, lp_low, truth_cases(icase), [0.,0.], .false., low_band(icase))
    call run_fixture(vol_file, mskdiam, smpd, lp_full, truth_cases(icase), [0.,0.], .true., full_band(icase))
    call run_fixture(vol_file, mskdiam, smpd, lp_low, truth_cases(icase), shift_truth, .false., &
        &recovery(icase))
enddo

low_alias_err  = sqrt(sum([(low_band(icase)%angle_error_joint**2, icase=1,ncases)])/real(ncases,dp))
full_alias_err = sqrt(sum([(full_band(icase)%angle_error_joint**2, icase=1,ncases)])/real(ncases,dp))
low_alias_max  = maxval([(low_band(icase)%angle_error_joint, icase=1,ncases)])
full_alias_max = maxval([(full_band(icase)%angle_error_joint, icase=1,ncases)])

write(logfhandle,'(a)') 'STAGE1_ALIASING_EXPERIMENT_BEGIN'
write(logfhandle,'(a,f10.4)') 'STAGE1_ALIASING_LOW_PASS_A: ', lp_low
write(logfhandle,'(a,f10.4)') 'STAGE1_ALIASING_FULL_BAND_A: ', lp_full
do icase = 1, ncases
    write(logfhandle,'(a,f10.4,3es16.8)') 'STAGE1_ALIASING_CASE_TRUTH/RMS_LOW/RMS_FULL/DELTA: ', &
        &truth_cases(icase), low_band(icase)%angle_error_joint, &
        &full_band(icase)%angle_error_joint, &
        &full_band(icase)%angle_error_joint-low_band(icase)%angle_error_joint
enddo
write(logfhandle,'(a,es16.8)') 'STAGE1_ALIASING_LOW_PASS_RMS_ERROR_DEG: ', low_alias_err
write(logfhandle,'(a,es16.8)') 'STAGE1_ALIASING_FULL_BAND_RMS_ERROR_DEG: ', full_alias_err
write(logfhandle,'(a,es16.8)') 'STAGE1_ALIASING_LOW_PASS_MAX_ERROR_DEG: ', low_alias_max
write(logfhandle,'(a,es16.8)') 'STAGE1_ALIASING_FULL_BAND_MAX_ERROR_DEG: ', full_alias_max
write(logfhandle,'(a,es16.8)') 'STAGE1_ALIASING_ERROR_CHANGE_DEG: ', full_alias_err-low_alias_err
write(logfhandle,'(a)') 'STAGE1_ALIASING_EXPERIMENT_END'

write(logfhandle,'(a)') 'SYNTHETIC_RECOVERY_TABLE_BEGIN'
write(logfhandle,'(a)') 'CASE  STAGE                    ANGLE_ERROR_DEG       SHIFT_RMS_PIX         RAW_LOSS'
do icase = 1, ncases
    write(logfhandle,'(a,i4,f10.4)') 'SYNTHETIC_CASE_BEGIN: ', icase, truth_cases(icase)
    write(logfhandle,'(i4,1x,a,3es24.8)') icase, 'GRID_ONLY              ', recovery(icase)%angle_error_grid, &
        &recovery(icase)%shift_rms_grid, recovery(icase)%loss_grid
    write(logfhandle,'(i4,1x,a,3es24.8)') icase, 'PARABOLIC               ', recovery(icase)%angle_error_parabola, &
        &recovery(icase)%shift_rms_parabola, recovery(icase)%loss_parabola
    write(logfhandle,'(i4,1x,a,3es24.8)') icase, 'INPL_CONT_YES           ', &
        &recovery(icase)%angle_error_joint, recovery(icase)%shift_rms_joint, recovery(icase)%loss_joint
    write(logfhandle,'(a,2f12.5)') 'SYNTHETIC_EXPECTED_SHIFT_XY: ', &
        &recovery(icase)%shift_expected
    write(logfhandle,'(a,2f12.5)') 'SYNTHETIC_GRID_SHIFT_ESTIMATE: ', &
        &recovery(icase)%shift_est_grid
    write(logfhandle,'(a,2f12.5)') 'SYNTHETIC_PARABOLIC_SHIFT_ESTIMATE: ', &
        &recovery(icase)%shift_est_parabola
    write(logfhandle,'(a,2f12.5)') 'SYNTHETIC_INPL_CONT_YES_SHIFT_ESTIMATE: ', &
        &recovery(icase)%shift_est_joint
    write(logfhandle,'(a,2es16.8,2l2)') 'PHASE3_GRADIENT_MAX_ERROR/TOL/FINITE/OK: ', &
        &recovery(icase)%joint_gradient_max_error, recovery(icase)%joint_gradient_tol, &
        &recovery(icase)%joint_gradient_finite, recovery(icase)%joint_gradient_ok
    write(logfhandle,'(a,2es16.8,1x,l1)') 'PHASE3_PARITY_MAX_ERROR/TOL/OK: ', &
        &recovery(icase)%parity_max_error, recovery(icase)%parity_tol, &
        &recovery(icase)%parity_ok
    write(logfhandle,'(a,3es16.8,1x,l1)') 'PHASE3_STRESS_SERIES_MIN/PARITY/SCORE_MAX/OK: ', &
        &recovery(icase)%stress_series_min, recovery(icase)%stress_parity_max_error, &
        &recovery(icase)%stress_score_max, &
        &recovery(icase)%stress_ok
    write(logfhandle,'(a,3l2)') 'SYNTHETIC_SHIFT_ACCEPTED_GRID/PARAB/INPL_CONT_YES: ', &
        &recovery(icase)%shift_grid_accepted, recovery(icase)%shift_parabola_accepted, &
        &recovery(icase)%joint_accepted
enddo
write(logfhandle,'(a,3es24.8)') 'RMS_OVER_CASES GRID/PARAB/INPL_CONT_YES ANGLE: ', &
    &sqrt(sum([(recovery(icase)%angle_error_grid**2,icase=1,ncases)])/real(ncases,dp)), &
    &sqrt(sum([(recovery(icase)%angle_error_parabola**2,icase=1,ncases)])/real(ncases,dp)), &
    &sqrt(sum([(recovery(icase)%angle_error_joint**2,icase=1,ncases)])/real(ncases,dp))
write(logfhandle,'(a,3es24.8)') 'RMS_OVER_CASES GRID/PARAB/INPL_CONT_YES SHIFT: ', &
    &sqrt(sum([(recovery(icase)%shift_rms_grid**2,icase=1,ncases)])/real(ncases,dp)), &
    &sqrt(sum([(recovery(icase)%shift_rms_parabola**2,icase=1,ncases)])/real(ncases,dp)), &
    &sqrt(sum([(recovery(icase)%shift_rms_joint**2,icase=1,ncases)])/real(ncases,dp))
write(logfhandle,'(a)') 'SYNTHETIC_RECOVERY_TABLE_END'

if( any(.not. low_band%finite) .or. any(.not. full_band%finite) .or. any(.not. recovery%finite) .or. &
    &.not. ieee_is_finite(low_alias_err) .or. .not. ieee_is_finite(full_alias_err) )then
    error stop 'Stage 1 validation produced a non-finite result'
endif
if( any(.not. low_band%joint_gradient_ok) .or. any(.not. full_band%joint_gradient_ok) .or. &
    &any(.not. recovery%joint_gradient_ok) )then
    error stop 'Stage 1 validation: joint analytic gradient disagrees with central differences'
endif
if( any(.not. low_band%parity_ok) .or. any(.not. full_band%parity_ok) .or. &
    &any(.not. recovery%parity_ok) )then
    error stop 'Stage 1 validation: continuous evaluator disagrees with the discrete reference at integer angles'
endif
if( any(.not. low_band%stress_ok) .or. any(.not. full_band%stress_ok) .or. &
    &any(.not. recovery%stress_ok) )then
    error stop 'Stage 1 validation: joint objective series unphysical under near-noiseless sigma2 stress'
endif
write(logfhandle,'(a)') 'SIMPLE_TEST_CONT_INPL_ROT2D_STAGE1_VALIDATION NORMAL STOP'

contains

subroutine parse_arguments(args, vol_file, mskdiam, smpd, lp_low, lp_full, truth_angle, shift_truth)
    type(cmdline), intent(inout) :: args
    character(len=*), intent(out) :: vol_file
    real, intent(out) :: mskdiam, smpd, lp_low, lp_full, truth_angle, shift_truth(2)
    type(string) :: vol_arg
    vol_file = ''
    mskdiam = 120.
    smpd = 1.3
    lp_low = 8.
    lp_full = 2.7
    truth_angle = 37.
    shift_truth = [2., -1.5]
    call args%parse_oldschool
    if( .not. args%defined('vol1') )then
        write(logfhandle,'(a)') 'simple_test_continuous_inplane_rotation2D_stage1_validation vol1=xx mskdiam=120 smpd=1.3 '// &
            &'lp=8 lpstop=2.7 angerr=37 xsh=2 ysh=-1.5'
        error stop 'vol1 is required'
    endif
    vol_arg = args%get_carg('vol1')
    vol_file = vol_arg%to_char()
    if( args%defined('mskdiam') ) mskdiam = args%get_rarg('mskdiam')
    if( args%defined('smpd') )     smpd       = args%get_rarg('smpd')
    if( args%defined('lp') )       lp_low     = args%get_rarg('lp')
    if( args%defined('lpstop') )   lp_full    = args%get_rarg('lpstop')
    if( args%defined('angerr') )   truth_angle = args%get_rarg('angerr')
    if( args%defined('xsh') )      shift_truth(1) = args%get_rarg('xsh')
    if( args%defined('ysh') )      shift_truth(2) = args%get_rarg('ysh')
    if( mskdiam <= 0. .or. smpd <= 0. .or. lp_low <= 2.*smpd .or. lp_full < 2.*smpd )then
        error stop 'invalid mskdiam, smpd, or low-pass limits'
    endif
end subroutine parse_arguments

subroutine run_fixture(vol_file, mskdiam, smpd, lp, truth_angle, shift_truth, hard_edge, result)
    character(len=*), intent(in) :: vol_file
    real, intent(in) :: mskdiam, smpd, lp, truth_angle, shift_truth(2)
    logical, intent(in) :: hard_edge
    type(recovery_result), intent(out) :: result
    type(parameters) :: p
    type(cmdline) :: cline
    type(builder) :: b
    type(ori) :: o_ref, o_particle
    type(pftc_shsrch_grad) :: direct_search
    type(pftc_shsrch_grad) :: joint_search
    real(sp), allocatable, target :: sigma2_noise(:,:)
    real(sp), allocatable :: raw_losses(:), scores(:)
    real :: limits(2,2), joint_limits(3,2), cxy(3), joint_cxy(3)
    real(dp) :: rotind_grid, rotind_parab
    real(dp) :: shift_grid(2), shift_parabola(2), shift_joint(2)
    real(dp) :: expected_shift(2), objective_final
    real(dp) :: objective_final_grid, objective_final_parabola
    real(dp) :: theta_rad
    real(dp) :: rotind_joint, joint_loss
    real(dp) :: joint_grad(3)
    real(dp) :: fd_error, fd_scale, grad_scale, parity_scale
    real(dp) :: stale_parity, fresh_min
    integer :: nrots, igrid, irot_direct, irot_joint, k
    logical :: fd_ok, parity_finite
    logical :: shift_grid_accepted, shift_parabola_accepted, joint_accepted

    call cline%set('vol1', vol_file)
    call cline%set('mskdiam', mskdiam)
    call cline%set('smpd', smpd)
    call cline%set('lp', lp)
    call cline%set('nptcls', 1.)
    call cline%set('ctf', 'no')
    call cline%set('objfun', 'euclid')
    call cline%check
    call b%init_params_and_build_strategy3D_tbox(cline, p)
    call set_bp_range3D(p, b, cline)
    call b%pftc%new(p, 1, [1,1], p%kfromto)
    call b%vol%read(p%vols(1))
    if( hard_edge )then
        call b%vol%mask3D_hard(p%msk)
    else
        call b%vol%mask3D_soft(p%msk)
    endif
    call b%vol%fft()
    call b%vol%expand_cmat()

    call b%eulspace%get_ori(1, o_ref)
    call o_ref%e3set(0.0)
    o_particle = o_ref
    call o_particle%e3set(truth_angle)
    call fproject_polar(b%vol, 1, o_particle, b%pftc, iseven=.true.)
    call b%pftc%cp_even_ref2ptcl(1, 1)
    call fproject_polar(b%vol, 1, o_ref, b%pftc, iseven=.true.)
    call b%pftc%set_eo(1, .true.)
    ! The particle is first rotated by truth_angle, then its polar Fourier
    ! coefficients are multiplied by the production shift phase.  Consequently
    ! the recovered candidate is expressed in the rotated reference frame:
    ! expected_shift = R(truth_angle) * shift_truth.
    if( sum(abs(shift_truth)) > 0. ) call b%pftc%shift_ptcl(1, shift_truth)
    call b%pftc%memoize_refs
    call b%pftc%memoize_ptcls
    allocate(sigma2_noise(p%kfromto(1):p%kfromto(2),1), source=1.0_sp)
    call b%pftc%assign_sigma2_noise(sigma2_noise)
    call b%pftc%memoize_sqsum_ptcl(1)

    nrots = b%pftc%get_nrots()
    allocate(raw_losses(nrots), scores(nrots))
    call b%pftc%gen_raw_euclid_vals(1, 1, [0._sp,0._sp], raw_losses)
    igrid = minloc(raw_losses, dim=1)
    rotind_grid = real(igrid,dp)
    rotind_parab = rotind_grid + real(parabolic_peak_offset(-raw_losses, igrid),dp)
    limits(:,1) = -5.; limits(:,2) = 5.
    call direct_search%new_direct(b, limits)
    call direct_search%set_indices(1,1)
    shift_grid = 0.; shift_parabola = 0.
    irot_direct = igrid
    cxy = direct_search%minimize_direct(irot_direct, [0.,0.], .5, 16, sh_rot=.false., &
        &objective_final=objective_final, raw_euclid=.true.)
    shift_grid_accepted = irot_direct > 0
    if( shift_grid_accepted ) shift_grid = real(cxy(2:),dp)
    objective_final_grid = objective_final
    irot_direct = modulo(nint(rotind_parab)-1,nrots)+1
    cxy = direct_search%minimize_direct(irot_direct, [0.,0.], .5, 16, sh_rot=.false., &
        &objective_final=objective_final, raw_euclid=.true.)
    shift_parabola_accepted = irot_direct > 0
    if( shift_parabola_accepted ) shift_parabola = real(cxy(2:),dp)
    objective_final_parabola = objective_final
    call direct_search%kill

    ! The inpl_cont=yes route optimizes (sx,sy,rotind_frac) directly from
    ! the same selected discrete candidate and zero native shift.
    joint_limits(:,1) = [-5., -5., real(rotind_grid)-2.]
    joint_limits(:,2) = [ 5.,  5., real(rotind_grid)+2.]
    call joint_search%new_joint(b, joint_limits, 100)
    call joint_search%set_indices(1,1)
    irot_joint = igrid
    joint_cxy = joint_search%minimize_joint(irot_joint, [0.,0.], &
        &sh_rot=.false., rotind_frac=rotind_joint, irot_in=igrid)
    joint_accepted = irot_joint > 0
    if( joint_accepted )then
        shift_joint = real(joint_cxy(2:),dp)
    else
        shift_joint = 0.d0
        rotind_joint = rotind_grid
    endif
    call b%pftc%gen_raw_euclid_grad_at_angle(1, 1, shift_joint, rotind_joint, joint_loss, joint_grad)
    ! central-difference check at the optimized pose AND at a deliberately
    ! non-stationary fractional pose, so near-zero gradients at the optimum
    ! cannot mask a broken component
    call fd_gradient_check(b, shift_joint, rotind_joint, fd_error, fd_scale, fd_ok)
    result%joint_gradient_max_error = fd_error
    result%joint_gradient_finite    = fd_ok .and. &
        &ieee_is_finite(joint_loss) .and. all(ieee_is_finite(joint_grad))
    grad_scale = fd_scale
    call fd_gradient_check(b, shift_joint + [0.3d0,-0.2d0], rotind_joint + 0.37d0, &
        &fd_error, fd_scale, fd_ok)
    result%joint_gradient_max_error = max(result%joint_gradient_max_error, fd_error)
    result%joint_gradient_finite    = result%joint_gradient_finite .and. fd_ok
    grad_scale = max(grad_scale, fd_scale)
    result%joint_gradient_tol = GRAD_FD_RTOL * max(1.d0, grad_scale)
    result%joint_gradient_ok  = result%joint_gradient_finite .and. &
        &result%joint_gradient_max_error <= result%joint_gradient_tol
    ! integer-angle parity: continuous evaluator vs the discrete-rotation
    ! reference implementation at several nonzero shifts
    call parity_check(b, nrots, igrid, result%parity_max_error, parity_scale, parity_finite)
    result%parity_tol = GRAD_FD_RTOL * max(1.d0, parity_scale)
    result%parity_ok  = parity_finite .and. &
        &result%parity_max_error <= result%parity_tol
    ! stress: cavgs-like near-noiseless, shell-dependent sigma2 weighting.
    ! The raw loss is nonnegative by construction; scan the coefficient
    ! series densely over the joint (shift, fractional-angle) search domain
    ! and assert it never goes materially negative and stays in parity with
    ! the double-precision discrete evaluator at grid angles
    do k = p%kfromto(1), p%kfromto(2)
        sigma2_noise(k,1) = 1.e-4_sp / real(k,sp)
    enddo
    call b%pftc%assign_sigma2_noise(sigma2_noise)
    ! deliberately do NOT re-memoize the weighted square sums yet: production
    ! updates sigma2 per iteration without re-memoizing, and both the joint
    ! evaluator and the vector scoring path must normalize with the CURRENT
    ! sigma2 and stay physical under stale wsqsums_ptcls
    call series_floor_scan(b, igrid, result%stress_series_min, stale_parity)
    call b%pftc%gen_raw_euclid_vals(1, 1, [1.7_sp, -2.3_sp], raw_losses)
    call b%pftc%gen_objfun_vals(1, 1, [1.7_sp, -2.3_sp], scores)
    result%stress_series_min = min(result%stress_series_min, real(minval(raw_losses),dp))
    result%stress_score_max = real(maxval(scores),dp)
    ! re-memoize for the cross-route parity check (the discrete reference
    ! normalizes with the memoized wsqsums)
    call b%pftc%memoize_sqsum_ptcl(1)
    call series_floor_scan(b, igrid, fresh_min, result%stress_parity_max_error)
    result%stress_series_min = min(result%stress_series_min, fresh_min)
    result%stress_ok = ieee_is_finite(result%stress_series_min) .and. &
        &ieee_is_finite(result%stress_parity_max_error) .and. &
        &all(ieee_is_finite(scores)) .and. minval(scores) >= 0._sp .and. &
        &result%stress_score_max <= 1.d0 .and. &
        &result%stress_series_min > -GRAD_FD_RTOL .and. &
        &result%stress_parity_max_error <= GRAD_FD_RTOL

    call joint_search%kill

    theta_rad = real(truth_angle,dp) * acos(-1.d0) / 180.d0
    expected_shift = [cos(theta_rad)*real(shift_truth(1),dp) - sin(theta_rad)*real(shift_truth(2),dp), &
        &sin(theta_rad)*real(shift_truth(1),dp) + cos(theta_rad)*real(shift_truth(2),dp)]

    result%angle_grid = grid_index_to_angle(b%pftc, rotind_grid)
    result%angle_parabola = grid_index_to_angle(b%pftc, rotind_parab)
    result%angle_joint = grid_index_to_angle(b%pftc, rotind_joint)
    result%angle_error_grid = angular_error(result%angle_grid, -real(truth_angle,dp))
    result%angle_error_parabola = angular_error(result%angle_parabola, -real(truth_angle,dp))
    result%angle_error_joint = angular_error(result%angle_joint, -real(truth_angle,dp))
    result%shift_rms_grid = sqrt(sum((shift_grid-expected_shift)**2)/2.d0)
    result%shift_rms_parabola = sqrt(sum((shift_parabola-expected_shift)**2)/2.d0)
    result%shift_rms_joint = sqrt(sum((shift_joint-expected_shift)**2)/2.d0)
    result%shift_est_grid = shift_grid
    result%shift_est_parabola = shift_parabola
    result%shift_est_joint = shift_joint
    result%shift_expected = expected_shift
    result%shift_grid_accepted = shift_grid_accepted
    result%shift_parabola_accepted = shift_parabola_accepted
    result%joint_accepted = joint_accepted
    result%loss_grid = objective_final_grid
    result%loss_parabola = objective_final_parabola
    result%loss_joint = joint_loss
    result%finite = ieee_is_finite(result%angle_error_grid) .and. &
        &ieee_is_finite(result%angle_error_parabola) .and. &
        &ieee_is_finite(result%angle_error_joint) .and. &
        &ieee_is_finite(result%shift_rms_grid) .and. &
        &ieee_is_finite(result%shift_rms_parabola) .and. &
        &ieee_is_finite(result%shift_rms_joint) .and. &
        &ieee_is_finite(result%loss_grid) .and. &
        &ieee_is_finite(result%loss_parabola) .and. &
        &ieee_is_finite(result%loss_joint) .and. &
        &ieee_is_finite(result%joint_gradient_max_error) .and. &
        &result%joint_gradient_finite

    call b%kill_strategy3D_tbox
    call b%kill_general_tbox
    deallocate(raw_losses, scores, sigma2_noise)
end subroutine run_fixture

! Central-difference check of the joint analytic gradient at one pose.
! Returns the max component error, the gradient magnitude scale for
! tolerance construction, and whether all evaluations were finite.
subroutine fd_gradient_check(b, at_shift, at_rotind, max_error, grad_scale, ok)
    type(builder), intent(inout) :: b
    real(dp), intent(in)  :: at_shift(2), at_rotind
    real(dp), intent(out) :: max_error, grad_scale
    logical,  intent(out) :: ok
    real(dp), parameter :: FD_STEP = 1.d-3
    real(dp) :: f0, g0(3), gtmp(3), fd_plus, fd_minus
    real(dp) :: probe_shift(2), probe_rotind
    integer :: idim
    call b%pftc%gen_raw_euclid_grad_at_angle(1, 1, at_shift, at_rotind, f0, g0)
    ok         = ieee_is_finite(f0) .and. all(ieee_is_finite(g0))
    max_error  = 0.d0
    grad_scale = maxval(abs(g0))
    do idim = 1, 3
        probe_shift  = at_shift
        probe_rotind = at_rotind
        if( idim <= 2 )then
            probe_shift(idim) = probe_shift(idim) + FD_STEP
        else
            probe_rotind = probe_rotind + FD_STEP
        endif
        call b%pftc%gen_raw_euclid_grad_at_angle(1, 1, probe_shift, probe_rotind, fd_plus, gtmp)
        probe_shift  = at_shift
        probe_rotind = at_rotind
        if( idim <= 2 )then
            probe_shift(idim) = probe_shift(idim) - FD_STEP
        else
            probe_rotind = probe_rotind - FD_STEP
        endif
        call b%pftc%gen_raw_euclid_grad_at_angle(1, 1, probe_shift, probe_rotind, fd_minus, gtmp)
        max_error = max(max_error, abs((fd_plus-fd_minus)/(2.d0*FD_STEP) - g0(idim)))
    enddo
    ok = ok .and. ieee_is_finite(max_error)
end subroutine fd_gradient_check

! Integer-angle parity: at grid rotation indices the continuous evaluator's
! loss and x/y gradient must match the discrete-rotation reference
! implementation gen_raw_euclid_grad_for_rot_8 at nonzero shifts.
subroutine parity_check(b, nrots, igrid, max_error, scale, ok)
    type(builder), intent(inout) :: b
    integer,  intent(in)  :: nrots, igrid
    real(dp), intent(out) :: max_error, scale
    logical,  intent(out) :: ok
    real(dp), parameter :: probe_shifts(2,3) = reshape( &
        &[1.7d0,-2.3d0, -0.6d0,0.4d0, 3.1d0,1.2d0], [2,3])
    real(dp) :: f_cont, g_cont(3), f_disc, g_disc(2)
    integer :: irots(3), j, m
    irots = [igrid, modulo(igrid-1+nrots/3, nrots)+1, modulo(igrid-1+(2*nrots)/3, nrots)+1]
    max_error = 0.d0
    scale     = 0.d0
    ok        = .true.
    do j = 1, 3
        do m = 1, 3
            call b%pftc%gen_raw_euclid_grad_at_angle(1, 1, probe_shifts(:,m), &
                &real(irots(j),dp), f_cont, g_cont)
            call b%pftc%gen_raw_euclid_grad_for_rot_8(1, 1, probe_shifts(:,m), &
                &irots(j), f_disc, g_disc)
            ok = ok .and. ieee_is_finite(f_cont) .and. all(ieee_is_finite(g_cont)) .and. &
                &ieee_is_finite(f_disc) .and. all(ieee_is_finite(g_disc))
            max_error = max(max_error, abs(f_cont-f_disc), &
                &abs(g_cont(1)-g_disc(1)), abs(g_cont(2)-g_disc(2)))
            scale = max(scale, abs(f_disc), maxval(abs(g_disc)))
        enddo
    enddo
    ok = ok .and. ieee_is_finite(max_error)
end subroutine parity_check

! Dense scan of the continuous joint objective over the (shift, fractional
! angle) domain the optimizer searches.  Reports the minimum series value
! (a materially negative minimum means the truncated series is unphysical
! and L-BFGS-B will descend into the artifact) and the worst integer-angle
! disagreement with the double-precision discrete evaluator.
subroutine series_floor_scan( b, igrid, series_min, parity_max )
    type(builder), intent(inout) :: b
    integer,  intent(in)  :: igrid
    real(dp), intent(out) :: series_min, parity_max
    real(dp) :: f, g(3), f_disc, g_disc(2), shift(2), rotind
    integer :: ish, jsh, ith
    series_min = huge(0.d0)
    parity_max = 0.d0
    do ish = -2, 2
        do jsh = -2, 2
            shift = [real(ish,dp)*1.25d0, real(jsh,dp)*1.25d0]
            do ith = -40, 40
                rotind = real(igrid,dp) + real(ith,dp)*0.05d0
                call b%pftc%gen_raw_euclid_grad_at_angle(1, 1, shift, rotind, f, g)
                series_min = min(series_min, f)
            enddo
            call b%pftc%gen_raw_euclid_grad_at_angle(1, 1, shift, real(igrid,dp), f, g)
            call b%pftc%gen_raw_euclid_grad_for_rot_8(1, 1, shift, igrid, f_disc, g_disc)
            parity_max = max(parity_max, abs(f - f_disc))
        enddo
    enddo
end subroutine series_floor_scan

pure real function parabolic_peak_offset(vals, j) result(offset)
    real,    intent(in) :: vals(:)
    integer, intent(in) :: j
    integer :: jm, jp, n
    real    :: denom, scale
    n = size(vals)
    offset = 0.
    if( n < 3 .or. j < 1 .or. j > n ) return
    jm = modulo(j - 2, n) + 1
    jp = modulo(j,     n) + 1
    denom = vals(jm) - 2.*vals(j) + vals(jp)
    scale = max(1., abs(vals(jm)), abs(vals(j)), abs(vals(jp)))
    if( denom >= -10.*epsilon(1.)*scale ) return
    offset = 0.5 * (vals(jm) - vals(jp)) / denom
    if( .not. ieee_is_finite(offset) .or. abs(offset) > 0.5 ) offset = 0.
end function parabolic_peak_offset

real(dp) function grid_index_to_angle(pftc, rotind_frac) result(angle)
    class(*), intent(in) :: pftc
    real(dp), intent(in) :: rotind_frac
    select type(pftc)
        class is (polarft_calc)
            angle = (rotind_frac - 1.d0) * real(pftc%get_dang(),dp)
        class default
            error stop 'invalid polarft type in grid_index_to_angle'
    end select
end function grid_index_to_angle

real(dp) function angular_error(angle, target) result(err)
    real(dp), intent(in) :: angle, target
    err = abs(modulo(angle-target+180.d0,360.d0)-180.d0)
end function angular_error

end program simple_test_continuous_inplane_rotation2D_stage1_validation
