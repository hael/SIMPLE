program simple_test_euclid_stage1_validation
use simple_pftc_srch_api
use simple_string, only: string
use simple_builder, only: builder
use simple_matcher_smpl_and_lplims, only: set_bp_range3D
use simple_projector_pft, only: fproject_polar
use simple_polarft_calc, only: polarft_calc
use simple_pftc_shsrch_grad, only: parabolic_peak_offset, pftc_shsrch_grad
use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
implicit none
#include "simple_local_flags.inc"

type :: recovery_result
    real(dp) :: angle_grid = 0.d0
    real(dp) :: angle_parabola = 0.d0
    real(dp) :: angle_continuous = 0.d0
    real(dp) :: angle_error_grid = 0.d0
    real(dp) :: angle_error_parabola = 0.d0
    real(dp) :: angle_error_continuous = 0.d0
    real(dp) :: shift_rms_grid = 0.d0
    real(dp) :: shift_rms_parabola = 0.d0
    real(dp) :: shift_rms_continuous = 0.d0
    real(dp) :: shift_est_grid(2) = 0.d0
    real(dp) :: shift_est_parabola(2) = 0.d0
    real(dp) :: shift_est_continuous(2) = 0.d0
    real(dp) :: shift_expected(2) = 0.d0
    logical :: shift_grid_accepted = .false.
    logical :: shift_parabola_accepted = .false.
    logical :: shift_continuous_accepted = .false.
    real(dp) :: loss_grid = 0.d0
    real(dp) :: loss_parabola = 0.d0
    real(dp) :: loss_continuous = 0.d0
    logical :: finite = .false.
end type recovery_result

type(cmdline) :: args
integer, parameter :: ncases = 3
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

low_alias_err  = sqrt(sum([(low_band(icase)%angle_error_continuous**2, icase=1,ncases)])/real(ncases,dp))
full_alias_err = sqrt(sum([(full_band(icase)%angle_error_continuous**2, icase=1,ncases)])/real(ncases,dp))
low_alias_max  = maxval([(low_band(icase)%angle_error_continuous, icase=1,ncases)])
full_alias_max = maxval([(full_band(icase)%angle_error_continuous, icase=1,ncases)])

write(logfhandle,'(a)') 'STAGE1_ALIASING_EXPERIMENT_BEGIN'
write(logfhandle,'(a,f10.4)') 'STAGE1_ALIASING_LOW_PASS_A: ', lp_low
write(logfhandle,'(a,f10.4)') 'STAGE1_ALIASING_FULL_BAND_A: ', lp_full
do icase = 1, ncases
    write(logfhandle,'(a,f10.4,3es16.8)') 'STAGE1_ALIASING_CASE_TRUTH/RMS_LOW/RMS_FULL/DELTA: ', &
        &truth_cases(icase), low_band(icase)%angle_error_continuous, &
        &full_band(icase)%angle_error_continuous, &
        &full_band(icase)%angle_error_continuous-low_band(icase)%angle_error_continuous
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
    write(logfhandle,'(i4,1x,a,3es24.8)') icase, 'CONTINUOUS              ', recovery(icase)%angle_error_continuous, &
        &recovery(icase)%shift_rms_continuous, recovery(icase)%loss_continuous
    write(logfhandle,'(a,2f12.5)') 'SYNTHETIC_EXPECTED_SHIFT_XY: ', &
        &recovery(icase)%shift_expected
    write(logfhandle,'(a,2f12.5)') 'SYNTHETIC_GRID_SHIFT_ESTIMATE: ', &
        &recovery(icase)%shift_est_grid
    write(logfhandle,'(a,2f12.5)') 'SYNTHETIC_PARABOLIC_SHIFT_ESTIMATE: ', &
        &recovery(icase)%shift_est_parabola
    write(logfhandle,'(a,2f12.5)') 'SYNTHETIC_CONTINUOUS_SHIFT_ESTIMATE: ', &
        &recovery(icase)%shift_est_continuous
    write(logfhandle,'(a,3l2)') 'SYNTHETIC_SHIFT_ACCEPTED_GRID/PARAB/CONT: ', &
        &recovery(icase)%shift_grid_accepted, recovery(icase)%shift_parabola_accepted, &
        &recovery(icase)%shift_continuous_accepted
enddo
write(logfhandle,'(a,3es24.8)') 'RMS_OVER_CASES GRID/PARAB/CONT ANGLE: ', &
    &sqrt(sum([(recovery(icase)%angle_error_grid**2,icase=1,ncases)])/real(ncases,dp)), &
    &sqrt(sum([(recovery(icase)%angle_error_parabola**2,icase=1,ncases)])/real(ncases,dp)), &
    &sqrt(sum([(recovery(icase)%angle_error_continuous**2,icase=1,ncases)])/real(ncases,dp))
write(logfhandle,'(a,3es24.8)') 'RMS_OVER_CASES GRID/PARAB/CONT SHIFT: ', &
    &sqrt(sum([(recovery(icase)%shift_rms_grid**2,icase=1,ncases)])/real(ncases,dp)), &
    &sqrt(sum([(recovery(icase)%shift_rms_parabola**2,icase=1,ncases)])/real(ncases,dp)), &
    &sqrt(sum([(recovery(icase)%shift_rms_continuous**2,icase=1,ncases)])/real(ncases,dp))
write(logfhandle,'(a)') 'JOINT_PHASE3             NOT_IMPLEMENTED'
write(logfhandle,'(a)') 'SYNTHETIC_RECOVERY_TABLE_END'

if( any(.not. low_band%finite) .or. any(.not. full_band%finite) .or. any(.not. recovery%finite) .or. &
    &.not. ieee_is_finite(low_alias_err) .or. .not. ieee_is_finite(full_alias_err) )then
    error stop 'Stage 1 validation produced a non-finite result'
endif
write(logfhandle,'(a)') 'SIMPLE_TEST_EUCLID_STAGE1_VALIDATION NORMAL STOP'

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
        write(logfhandle,'(a)') 'simple_test_euclid_stage1_validation vol1=xx mskdiam=120 smpd=1.3 '// &
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
    real(sp), allocatable, target :: sigma2_noise(:,:)
    real(sp), allocatable :: raw_losses(:)
    complex(sp), allocatable :: coeffs(:)
    real :: limits(2,2), cxy(3)
    real(dp) :: theta_grid, theta_parab, theta_cont, residual, dtheta, ddtheta
    real(dp) :: trial_theta, trial_residual, trial_dtheta, trial_ddtheta, step
    real(dp) :: shift_grid(2), shift_parabola(2), shift_continuous(2)
    real(dp) :: expected_shift(2), objective_initial, objective_final
    real(dp) :: objective_final_grid, objective_final_parabola, objective_final_continuous
    real(dp) :: theta_rad
    integer :: nrots, igrid, irot_direct, iter, iback, pftsz
    logical :: accepted, shift_grid_accepted, shift_parabola_accepted, shift_continuous_accepted

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
    call b%vol%expand_cmat(p%box)

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
    pftsz = b%pftc%get_pftsz()
    allocate(raw_losses(nrots), coeffs(pftsz+1))
    call b%pftc%gen_raw_euclid_vals(1, 1, [0._sp,0._sp], raw_losses)
    igrid = minloc(raw_losses, dim=1)
    theta_grid = real(igrid,dp)
    theta_parab = theta_grid + real(parabolic_peak_offset(-raw_losses, igrid),dp)
    call b%pftc%gen_euclid_angular_coeffs(1, 1, [0._sp,0._sp], coeffs, irot_direct)
    theta_cont = real(irot_direct,dp)
    call b%pftc%eval_euclid_resid_at_angle(coeffs, theta_cont, residual, dtheta, ddtheta)
    if( ieee_is_finite(residual) .and. ieee_is_finite(dtheta) .and. ieee_is_finite(ddtheta) )then
        do iter = 1, 3
            if( ddtheta <= 0.d0 ) exit
            step = -dtheta/ddtheta
            if( .not. ieee_is_finite(step) ) exit
            step = max(-0.5d0,min(0.5d0,step))
            accepted = .false.
            do iback = 0, 3
                trial_theta = theta_cont + step
                call b%pftc%eval_euclid_resid_at_angle(coeffs, trial_theta, &
                    &trial_residual, trial_dtheta, trial_ddtheta)
                if( ieee_is_finite(trial_residual) .and. trial_residual <= residual )then
                    theta_cont = trial_theta
                    residual = trial_residual
                    dtheta = trial_dtheta
                    ddtheta = trial_ddtheta
                    accepted = .true.
                    exit
                endif
                step = 0.5d0*step
            enddo
            if( .not. accepted ) exit
        enddo
    endif

    limits(:,1) = -5.; limits(:,2) = 5.
    call direct_search%new(b, limits, opt_angle=.false., direct_only=.true.)
    call direct_search%set_indices(1,1)
    shift_grid = 0.; shift_parabola = 0.
    irot_direct = igrid
    cxy = direct_search%minimize_direct(irot_direct, [0.,0.], .5, 16, sh_rot=.false., &
        &objective_initial=objective_initial, objective_final=objective_final, raw_euclid=.true.)
    shift_grid_accepted = irot_direct > 0
    if( shift_grid_accepted ) shift_grid = real(cxy(2:),dp)
    objective_final_grid = objective_final
    irot_direct = modulo(nint(theta_parab)-1,nrots)+1
    cxy = direct_search%minimize_direct(irot_direct, [0.,0.], .5, 16, sh_rot=.false., &
        &objective_initial=objective_initial, objective_final=objective_final, raw_euclid=.true.)
    shift_parabola_accepted = irot_direct > 0
    if( shift_parabola_accepted ) shift_parabola = real(cxy(2:),dp)
    objective_final_parabola = objective_final

    shift_continuous = 0.
    irot_direct = modulo(nint(theta_cont)-1,nrots)+1
    cxy = direct_search%minimize_direct(irot_direct, [0.,0.], .5, 16, sh_rot=.false., &
        &objective_initial=objective_initial, objective_final=objective_final, raw_euclid=.true.)
    shift_continuous_accepted = irot_direct > 0
    if( shift_continuous_accepted ) shift_continuous = real(cxy(2:),dp)
    objective_final_continuous = objective_final
    call direct_search%kill

    theta_rad = real(truth_angle,dp) * acos(-1.d0) / 180.d0
    expected_shift = [cos(theta_rad)*real(shift_truth(1),dp) - sin(theta_rad)*real(shift_truth(2),dp), &
        &sin(theta_rad)*real(shift_truth(1),dp) + cos(theta_rad)*real(shift_truth(2),dp)]

    result%angle_grid = grid_index_to_angle(b%pftc, theta_grid)
    result%angle_parabola = grid_index_to_angle(b%pftc, theta_parab)
    result%angle_continuous = grid_index_to_angle(b%pftc, theta_cont)
    result%angle_error_grid = angular_error(result%angle_grid, -real(truth_angle,dp))
    result%angle_error_parabola = angular_error(result%angle_parabola, -real(truth_angle,dp))
    result%angle_error_continuous = angular_error(result%angle_continuous, -real(truth_angle,dp))
    result%shift_rms_grid = sqrt(sum((shift_grid-expected_shift)**2)/2.d0)
    result%shift_rms_parabola = sqrt(sum((shift_parabola-expected_shift)**2)/2.d0)
    result%shift_rms_continuous = sqrt(sum((shift_continuous-expected_shift)**2)/2.d0)
    result%shift_est_grid = shift_grid
    result%shift_est_parabola = shift_parabola
    result%shift_est_continuous = shift_continuous
    result%shift_expected = expected_shift
    result%shift_grid_accepted = shift_grid_accepted
    result%shift_parabola_accepted = shift_parabola_accepted
    result%shift_continuous_accepted = shift_continuous_accepted
    result%loss_grid = objective_final_grid
    result%loss_parabola = objective_final_parabola
    result%loss_continuous = objective_final_continuous
    result%finite = all(ieee_is_finite([result%angle_error_grid, result%angle_error_parabola, &
        &result%angle_error_continuous, result%shift_rms_grid, result%shift_rms_parabola, &
        &result%shift_rms_continuous, result%loss_grid, result%loss_parabola, result%loss_continuous]))

    call b%kill_strategy3D_tbox
    call b%kill_general_tbox
    deallocate(raw_losses, coeffs, sigma2_noise)
end subroutine run_fixture

real(dp) function grid_index_to_angle(pftc, theta) result(angle)
    class(*), intent(in) :: pftc
    real(dp), intent(in) :: theta
    integer :: iang, nrots
    select type(pftc)
        class is (polarft_calc)
            nrots = pftc%get_nrots()
            iang = modulo(nint(theta)-1,nrots)+1
            angle = real(pftc%get_rot(iang),dp) + (theta-real(iang,dp))*real(pftc%get_dang(),dp)
        class default
            error stop 'invalid polarft type in grid_index_to_angle'
    end select
end function grid_index_to_angle

real(dp) function residual_at(coeffs, pftc, theta) result(residual)
    complex(sp), intent(in) :: coeffs(:)
    class(*), intent(in) :: pftc
    real(dp), intent(in) :: theta
    real(dp) :: dtheta, ddtheta
    select type(pftc)
        class is (polarft_calc)
            call pftc%eval_euclid_resid_at_angle(coeffs, theta, residual, dtheta, ddtheta)
        class default
            error stop 'invalid polarft type in residual_at'
    end select
end function residual_at

real(dp) function angular_error(angle, target) result(err)
    real(dp), intent(in) :: angle, target
    err = abs(modulo(angle-target+180.d0,360.d0)-180.d0)
end function angular_error

end program simple_test_euclid_stage1_validation
