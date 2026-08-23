program simple_test_continuous_inplane_rotation2D_route_identity
use simple_pftc_srch_api
use simple_builder,         only: builder
use simple_matcher_smpl_and_lplims, only: set_bp_range3D
use simple_projector_pft,   only: fproject_polar
use simple_pftc_shsrch_grad, only: pftc_shsrch_grad
use simple_strategy2D_srch, only: strategy2D_srch, strategy2D_spec
use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
implicit none

type(parameters), target :: p
type(cmdline)   :: cline
type(builder)   :: b
type(ori)       :: o
type(pftc_shsrch_grad) :: legacy_shift_search
type(pftc_shsrch_grad) :: joint_search
type(strategy2D_srch) :: probabilistic_refine_search
type(strategy2D_spec) :: probabilistic_refine_spec
real(sp), allocatable :: legacy_scores(:), raw_losses(:)
real(dp), allocatable :: scalar_losses(:)
real(sp), allocatable :: sigma2_noise(:,:)
real(dp) :: max_legacy_error, max_route_error, tol, scalar_loss, grad(2)
real :: shift_limits(2,2), joint_limits(3,2), seed_shift(2), selected_corr
integer :: irot, nrots, selected_irot, expected_irot
integer :: env_status
character(len=16) :: rejection_mode
logical :: expect_hybrid_rejection

if( command_argument_count() < 4 )then
    write(logfhandle,'(a)',advance='no') &
        'simple_test_continuous_inplane_rotation2D_route_identity vol1=xx mskdiam=xx lp=xx'
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
expect_hybrid_rejection = .false.
rejection_mode = ''
call get_environment_variable('SIMPLE_TEST_EXPECT_HYBRID_INPL_REJECTION', &
    &rejection_mode, status=env_status)
if( env_status == 0 ) expect_hybrid_rejection = trim(rejection_mode) == 'yes'
call b%init_params_and_build_strategy3D_tbox(cline, p)
call set_bp_range3D(p, b, cline)

call b%pftc%new(p, p%nptcls, [1, p%nptcls], p%kfromto)
call b%vol%read(p%vols(1))
call b%vol%mask3D_soft(p%msk)
call b%vol%fft()
call b%vol%expand_cmat()
call b%eulspace%get_ori(irnd_uni(p%nspace), o)
call fproject_polar(b%vol, 1, o, b%pftc, iseven=.true.)
call b%pftc%cp_even_ref2ptcl(1, 1)
call b%pftc%set_eo(1, .true.)
call b%pftc%memoize_refs
call b%pftc%memoize_ptcls

allocate(sigma2_noise(p%kfromto(1):p%kfromto(2), 1), source=1.0_sp)
call b%pftc%assign_sigma2_noise(sigma2_noise)
call b%pftc%memoize_sqsum_ptcl(1)

shift_limits(:,1) = -1.
shift_limits(:,2) =  1.
call legacy_shift_search%new_legacy(b, shift_limits)
call legacy_shift_search%kill

probabilistic_refine_spec%iptcl       = 1
probabilistic_refine_spec%iptcl_batch = 1
probabilistic_refine_spec%iptcl_map   = 1
p%l_prob_align_mode = .false.
p%l_objfun_den      = .false.
p%inpl_cont         = 'no'
call probabilistic_refine_search%new(p, probabilistic_refine_spec, b)
if( probabilistic_refine_search%uses_continuous_refinement() )then
    error stop 'inpl_cont=no unexpectedly enables probabilistic continuous polish'
endif
if( probabilistic_refine_search%joint_inpl_optimizer%uses_joint_inplane() )then
    error stop 'inpl_cont=no unexpectedly constructed the joint optimizer'
endif
if( .not. probabilistic_refine_search%grad_shsrch_first_obj%does_opt_angle() )then
    error stop 'inpl_cont=no did not retain the legacy seed-search angle update'
endif
call probabilistic_refine_search%kill
p%inpl_cont = 'yes'
call probabilistic_refine_search%new(p, probabilistic_refine_spec, b)
if( .not. probabilistic_refine_search%uses_continuous_refinement() )then
    error stop 'inpl_cont=yes did not enable continuous refinement'
endif
if( .not. probabilistic_refine_search%joint_inpl_optimizer%uses_joint_inplane() )then
    error stop 'inpl_cont=yes did not construct the joint optimizer'
endif
! selection parity: shift seeding and candidate scoring stay on the legacy
! callback route regardless of inpl_cont; the joint optimizer only polishes
! the committed assignment
if( .not. probabilistic_refine_search%grad_shsrch_first_obj%does_opt_angle() )then
    error stop 'inpl_cont=yes must retain the legacy seed-search angle update (selection parity)'
endif
call probabilistic_refine_search%kill
p%l_prob_align_mode = .true.
p%inpl_cont = 'yes'
call probabilistic_refine_search%new(p, probabilistic_refine_spec, b)
if( .not. probabilistic_refine_search%uses_continuous_refinement() )then
    error stop 'Probabilistic continuous route did not enable selected-candidate refinement'
endif
if( .not. probabilistic_refine_search%joint_inpl_optimizer%uses_joint_inplane() )then
    error stop 'inpl_cont=yes probabilistic route did not construct the joint optimizer'
endif
call probabilistic_refine_search%kill
p%inpl_cont = 'yes'
p%l_objfun_den = .true.
if( b%pftc%is_raw_euclid_objfun() )then
    error stop 'Hybrid Euclidean objective unexpectedly reports raw continuous derivative support'
endif
call legacy_shift_search%new_legacy(b, shift_limits)
call legacy_shift_search%kill

! The strategy constructor deliberately rejects this policy with ERROR STOP.
! A child-process invocation reaches the call below and lets the mother test
! verify both the nonzero exit and the exact policy diagnostic.
if( expect_hybrid_rejection )then
    call probabilistic_refine_search%new(p, probabilistic_refine_spec, b)
    error stop 'Hybrid Euclidean objective unexpectedly accepted inpl_cont=yes'
endif
p%l_prob_align_mode = .false.
p%l_objfun_den      = .false.
p%inpl_cont         = 'no'

nrots = b%pftc%get_nrots()
allocate(legacy_scores(nrots), raw_losses(nrots), scalar_losses(nrots))
call b%pftc%gen_objfun_vals(1, 1, [0.0_sp, 0.0_sp], legacy_scores)
call b%pftc%gen_raw_euclid_vals(1, 1, [0.0_sp, 0.0_sp], raw_losses)
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

! Continuous refinement must initialize exactly like one invocation of the
! legacy discrete callback at the caller's x/y seed. It must not use the
! incoming particle angle or search a 5x5 grid of alternative shifts.
seed_shift = [0.5, -0.25]
call b%pftc%gen_objfun_vals(1, 1, seed_shift, legacy_scores)
expected_irot = maxloc(legacy_scores, dim=1)
joint_limits(1:2,1) = -1.
joint_limits(1:2,2) =  1.
joint_limits(3,:) = [1.-2., real(nrots)+2.]
call joint_search%new_joint(b, joint_limits, 100)
call joint_search%set_indices(1, 1)
call joint_search%select_best_discrete_angle(seed_shift, selected_irot, selected_corr)
if( joint_search%does_opt_angle() ) &
    &error stop 'joint seed selector unexpectedly attached the legacy callback'
if( selected_irot /= expected_irot ) &
    &error stop 'joint seed selector did not reproduce the callback index at x/y seed'
if( abs(selected_corr-legacy_scores(expected_irot)) > 5.e-6 ) &
    &error stop 'joint seed selector did not reproduce the callback score at x/y seed'
call joint_search%kill

tol = 5.0e-4_dp * (1.0_dp + maxval(abs(real(raw_losses, dp))))
write(logfhandle,'(a,i0)')    'ROUTE_IDENTITY_NROTS: ', nrots
write(logfhandle,'(a,es16.8)') 'ROUTE_IDENTITY_MAX_LEGACY_ERROR: ', max_legacy_error
write(logfhandle,'(a,es16.8)') 'ROUTE_IDENTITY_MAX_SCALAR_ERROR: ', max_route_error
write(logfhandle,'(a,es16.8)') 'ROUTE_IDENTITY_TOLERANCE: ', tol
if( .not. ieee_is_finite(max_legacy_error) .or. &
    .not. ieee_is_finite(max_route_error) .or. &
    max_legacy_error > tol .or. max_route_error > tol )then
    call b%kill_strategy3D_tbox
    call b%kill_general_tbox
    deallocate(legacy_scores, raw_losses, scalar_losses, sigma2_noise)
    error stop 'Euclidean route identity failed'
endif

call b%kill_strategy3D_tbox
call b%kill_general_tbox
deallocate(legacy_scores, raw_losses, scalar_losses, sigma2_noise)
write(logfhandle,'(a)') 'SIMPLE_TEST_CONT_INPL_ROT2D_ROUTE_IDENTITY NORMAL STOP'
end program simple_test_continuous_inplane_rotation2D_route_identity
