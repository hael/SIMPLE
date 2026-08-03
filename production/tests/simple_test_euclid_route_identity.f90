program simple_test_euclid_route_identity
use simple_pftc_srch_api
use simple_builder,         only: builder
use simple_matcher_smpl_and_lplims, only: set_bp_range3D
use simple_projector_pft,   only: fproject_polar
use simple_pftc_shsrch_grad, only: parabolic_peak_offset
use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
implicit none

type(parameters) :: p
type(cmdline)   :: cline
type(builder)   :: b
type(ori)       :: o
real(sp), allocatable :: legacy_scores(:), raw_losses(:)
real(dp), allocatable :: scalar_losses(:)
real(sp), allocatable :: sigma2_noise(:,:)
real(dp) :: max_legacy_error, max_route_error, tol, scalar_loss, grad(2)
complex(sp), allocatable :: angular_coeffs(:)
real(dp) :: max_angular_error, max_period_error, fd_dtheta_error, fd_ddtheta_error
real(dp) :: angular_loss, angular_dtheta, angular_ddtheta, loss_minus, loss_plus
real(dp) :: fd_dtheta, fd_ddtheta, theta_eval, h, angle_tol
real :: parabola_vals(5), parabola_offset
integer :: irot, nrots
integer :: angular_best_rot

parabola_vals = [ -1.5625, -0.0625, -0.5625, -3.0625, -5.0625 ]
parabola_offset = parabolic_peak_offset(parabola_vals, 2)
if( abs(parabola_offset - 0.25) > 10.*epsilon(1.) )then
    error stop 'Parabolic interpolation unit check failed'
endif

if( command_argument_count() < 4 )then
    write(logfhandle,'(a)',advance='no') &
        'simple_test_euclid_route_identity vol1=xx mskdiam=xx lp=xx'
    write(logfhandle,'(a)') ' smpd=xx>'
    stop 2
endif

call cline%parse_oldschool
call cline%checkvar('vol1', 1)
call cline%checkvar('mskdiam', 2)
call cline%checkvar('smpd', 3)
call cline%checkvar('lp', 4)
call cline%set('nptcls', 1.0)
call cline%set('ctf', 'no')
call cline%set('objfun', 'euclid')
call cline%check
call b%init_params_and_build_strategy3D_tbox(cline, p)
call set_bp_range3D(p, b, cline)

call b%pftc%new(p, p%nptcls, [1, p%nptcls], p%kfromto)
call b%vol%read(p%vols(1))
call b%vol%mask3D_soft(p%msk)
call b%vol%fft()
call b%vol%expand_cmat(p%box)
call b%eulspace%get_ori(irnd_uni(p%nspace), o)
call fproject_polar(b%vol, 1, o, b%pftc, iseven=.true.)
call b%pftc%cp_even_ref2ptcl(1, 1)
call b%pftc%set_eo(1, .true.)
call b%pftc%memoize_refs
call b%pftc%memoize_ptcls

allocate(sigma2_noise(p%kfromto(1):p%kfromto(2), 1), source=1.0_sp)
call b%pftc%assign_sigma2_noise(sigma2_noise)
call b%pftc%memoize_sqsum_ptcl(1)

nrots = b%pftc%get_nrots()
allocate(legacy_scores(nrots), raw_losses(nrots), scalar_losses(nrots))
call b%pftc%gen_objfun_vals(1, 1, [0.0_sp, 0.0_sp], legacy_scores)
call b%pftc%gen_raw_euclid_vals(1, 1, [0.0_sp, 0.0_sp], raw_losses)
allocate(angular_coeffs(b%pftc%get_pftsz()+1))
call b%pftc%gen_euclid_angular_coeffs(1, 1, [0.0_sp, 0.0_sp], &
    &angular_coeffs, angular_best_rot)

max_angular_error = 0.d0
do irot = 1, nrots
    call b%pftc%eval_euclid_resid_at_angle(angular_coeffs, real(irot,dp), &
        &angular_loss, angular_dtheta, angular_ddtheta)
    max_angular_error = max(max_angular_error, &
        &abs(angular_loss - real(raw_losses(irot),dp)))
enddo
if( angular_best_rot /= minloc(raw_losses, dim=1) )then
    error stop 'Angular coefficient best-rotation check failed'
endif

theta_eval = 0.37_dp
h = 1.e-2_dp
call b%pftc%eval_euclid_resid_at_angle(angular_coeffs, theta_eval, &
    &angular_loss, angular_dtheta, angular_ddtheta)
call b%pftc%eval_euclid_resid_at_angle(angular_coeffs, theta_eval-h, &
    &loss_minus, fd_dtheta, fd_ddtheta)
call b%pftc%eval_euclid_resid_at_angle(angular_coeffs, theta_eval+h, &
    &loss_plus, fd_dtheta, fd_ddtheta)
fd_dtheta  = (loss_plus - loss_minus) / (2.d0*h)
fd_ddtheta = (loss_plus - 2.d0*angular_loss + loss_minus) / (h*h)
fd_dtheta_error  = abs(fd_dtheta - angular_dtheta)
fd_ddtheta_error = abs(fd_ddtheta - angular_ddtheta)
call b%pftc%eval_euclid_resid_at_angle(angular_coeffs, theta_eval+real(nrots,dp), &
    &loss_plus, fd_dtheta, fd_ddtheta)
max_period_error = max(abs(loss_plus-angular_loss), abs(fd_dtheta-angular_dtheta), &
    &abs(fd_ddtheta-angular_ddtheta))
angle_tol = 5.e-4_dp * (1.d0 + maxval(abs(real(raw_losses,dp))))

max_legacy_error = 0.0_dp
max_route_error  = 0.0_dp
do irot = 1, nrots
    max_legacy_error = max(max_legacy_error, &
        abs(real(legacy_scores(irot), dp) - exp(-real(raw_losses(irot), dp))))
    call b%pftc%gen_raw_euclid_grad_for_rot_8(1, 1, [0.0_dp, 0.0_dp], &
        irot, scalar_loss, grad)
    scalar_losses(irot) = scalar_loss
    max_route_error = max(max_route_error, &
        abs(scalar_loss - real(raw_losses(irot), dp)))
enddo

tol = 5.0e-4_dp * (1.0_dp + maxval(abs(real(raw_losses, dp))))
write(logfhandle,'(a,i0)')    'ROUTE_IDENTITY_NROTS: ', nrots
write(logfhandle,'(a,es16.8)') 'ROUTE_IDENTITY_MAX_LEGACY_ERROR: ', max_legacy_error
write(logfhandle,'(a,es16.8)') 'ROUTE_IDENTITY_MAX_SCALAR_ERROR: ', max_route_error
write(logfhandle,'(a,es16.8)') 'ROUTE_IDENTITY_TOLERANCE: ', tol
write(logfhandle,'(a,i0)')    'ANGULAR_COEFF_BEST_ROT: ', angular_best_rot
write(logfhandle,'(a,es16.8)') 'ANGULAR_COEFF_MAX_GRID_ERROR: ', max_angular_error
write(logfhandle,'(a,es16.8)') 'ANGULAR_COEFF_MAX_PERIOD_ERROR: ', max_period_error
write(logfhandle,'(a,es16.8)') 'ANGULAR_COEFF_DTHETA_ERROR: ', fd_dtheta_error
write(logfhandle,'(a,es16.8)') 'ANGULAR_COEFF_DDTHETA_ERROR: ', fd_ddtheta_error
write(logfhandle,'(a,es16.8)') 'ANGULAR_COEFF_TOLERANCE: ', angle_tol

if( .not. ieee_is_finite(max_legacy_error) .or. &
    .not. ieee_is_finite(max_route_error) .or. &
    .not. ieee_is_finite(max_angular_error) .or. &
    .not. ieee_is_finite(max_period_error) .or. &
    .not. ieee_is_finite(fd_dtheta_error) .or. &
    .not. ieee_is_finite(fd_ddtheta_error) .or. &
    max_legacy_error > tol .or. max_route_error > tol .or. &
    max_angular_error > angle_tol .or. max_period_error > angle_tol .or. &
    fd_dtheta_error > 10.d0*angle_tol .or. fd_ddtheta_error > 100.d0*angle_tol )then
    call b%kill_strategy3D_tbox
    call b%kill_general_tbox
    deallocate(legacy_scores, raw_losses, scalar_losses, sigma2_noise, angular_coeffs)
    error stop 'Euclidean or angular coefficient route identity failed'
endif

call b%kill_strategy3D_tbox
call b%kill_general_tbox
deallocate(legacy_scores, raw_losses, scalar_losses, sigma2_noise, angular_coeffs)
write(logfhandle,'(a)') 'SIMPLE_TEST_EUCLID_ROUTE_IDENTITY NORMAL STOP'
end program simple_test_euclid_route_identity
