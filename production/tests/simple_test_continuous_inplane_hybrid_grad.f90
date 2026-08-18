! Validation harness for the continuous-angle hybrid evaluator used when
! objfun=euclid and objfun_den=yes.
!
! Checks integer-grid identity with gen_hybrid_scores through the public
! gen_objfun_vals dispatch, analytic versus central-difference derivatives in
! (sx,sy,rotind_frac), and joint-route capability advertisement.
program simple_test_continuous_inplane_hybrid_grad
use simple_pftc_srch_api
use simple_string,  only: string
use simple_builder, only: builder
use simple_pftc_shsrch_grad, only: pftc_shsrch_grad
use simple_matcher_smpl_and_lplims, only: set_bp_range3D
use simple_projector_pft, only: fproject_polar
use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
implicit none
#include "simple_local_flags.inc"

real(dp), parameter :: IDENTITY_ATOL  = 1.d-4
real(dp), parameter :: GRAD_FD_RTOL   = 1.d-2
real(dp), parameter :: GRAD_FD_AFLOOR = 3.d-3
real(dp), parameter :: PROBE_SHIFTS(2,3) = reshape( &
    &[0.d0,0.d0, 1.7d0,-2.3d0, -0.6d0,0.4d0], [2,3])
real(dp), parameter :: PROBE_ROTOFFS(4) = [0.37d0, -1.42d0, 0.5d0, 1.93d0]

type(cmdline)    :: args
type(parameters), target :: p
type(cmdline)    :: cline
type(builder), target :: b
type(pftc_shsrch_grad) :: joint_srch
type(ori)        :: o_ref, o_particle
type(string)     :: vol_arg
real(sp), allocatable, target :: sigma2_noise(:,:)
real(sp), allocatable :: scores(:)
complex(sp), allocatable :: den_pft(:,:)
character(len=512) :: vol_file
real     :: mskdiam, smpd, lp, truth_angle, shift_truth(2)
real(dp) :: f, grad(3), identity_err, identity_max, grad_err, grad_max, grad_tol
real(dp) :: probe_rotind
real     :: joint_lims(3,2)
integer  :: nrots, igrid, irot, ishift, ioff, failures
logical  :: ok

call parse_arguments(args, vol_file, mskdiam, smpd, lp, truth_angle, shift_truth)
failures = 0

call cline%set('vol1', trim(vol_file))
call cline%set('mskdiam', mskdiam)
call cline%set('smpd', smpd)
call cline%set('lp', lp)
call cline%set('nptcls', 1.)
call cline%set('ctf', 'no')
call cline%set('objfun', 'euclid')
call cline%set('objfun_den', 'yes')
call cline%set('objfun_den_w', 0.37)
call cline%check
call b%init_params_and_build_strategy3D_tbox(cline, p)
call set_bp_range3D(p, b, cline)
call b%pftc%new(p, 1, [1,1], p%kfromto)
call b%vol%read(p%vols(1))
call b%vol%mask3D_soft(p%msk)
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
if( sum(abs(shift_truth)) > 0. ) call b%pftc%shift_ptcl(1, shift_truth)

! Use the shifted raw particle as the denoised fixture. Production supplies a
! separately denoised PFT; equality here isolates the hybrid calculus while
! still exercising its distinct storage, memoization, and normalization.
den_pft = b%pftc%allocate_pft()
call b%pftc%get_ptcl_pft(1, den_pft(:,p%kfromto(1):p%kfromto(2)))
call b%pftc%set_ptcl_den_pft(1, den_pft)
call b%pftc%memoize_refs
call b%pftc%memoize_ptcls
call b%pftc%memoize_ptcls_den
allocate(sigma2_noise(p%kfromto(1):p%kfromto(2),1), source=1.0_sp)
call b%pftc%assign_sigma2_noise(sigma2_noise)
call b%pftc%memoize_sqsum_ptcl(1)
nrots = b%pftc%get_nrots()
allocate(scores(nrots))

ok = b%pftc%is_hybrid_objfun() .and. b%pftc%is_joint_grad_objfun()
write(logfhandle,'(a,l1)') 'HYBRID_JOINT_CAPABILITY_OK: ', ok
if( .not. ok ) failures = failures + 1
joint_lims(1:2,1) = -p%trs
joint_lims(1:2,2) =  p%trs
joint_lims(3,:)   = [1.-2., real(nrots)+2.]
call joint_srch%new_joint(b, joint_lims, p%maxits_sh)
call joint_srch%kill

identity_max = 0.d0
do ishift = 1, size(PROBE_SHIFTS, 2)
    call b%pftc%gen_objfun_vals(1, 1, real(PROBE_SHIFTS(:,ishift),sp), scores)
    do irot = 1, nrots
        call b%pftc%gen_hybrid_grad_at_angle(1, 1, PROBE_SHIFTS(:,ishift), &
            &real(irot,dp), f, grad)
        identity_err = abs(-f - real(scores(irot),dp))
        identity_max = max(identity_max, identity_err)
        if( .not. (ieee_is_finite(f) .and. all(ieee_is_finite(grad))) ) failures = failures + 1
    enddo
enddo
ok = identity_max <= IDENTITY_ATOL
write(logfhandle,'(a,2es16.8,1x,l1)') 'HYBRID_GRID_IDENTITY_MAX_ERROR/TOL/OK: ', &
    &identity_max, IDENTITY_ATOL, ok
if( .not. ok ) failures = failures + 1

call b%pftc%gen_objfun_vals(1, 1, [0._sp,0._sp], scores)
igrid = maxloc(scores, dim=1)
grad_max = 0.d0
do ishift = 1, size(PROBE_SHIFTS, 2)
    do ioff = 1, size(PROBE_ROTOFFS)
        probe_rotind = real(igrid,dp) + PROBE_ROTOFFS(ioff)
        call fd_gradient_check(b, PROBE_SHIFTS(:,ishift), probe_rotind, grad_err, grad_tol, ok)
        grad_max = max(grad_max, grad_err)
        if( .not. ok )then
            write(logfhandle,'(a,3f10.4,2es16.8)') 'HYBRID_GRADIENT_CHECK_FAILED_AT: ', &
                &PROBE_SHIFTS(:,ishift), probe_rotind, grad_err, grad_tol
            failures = failures + 1
        endif
    enddo
enddo
write(logfhandle,'(a,es16.8)') 'HYBRID_GRADIENT_MAX_ERROR: ', grad_max

call b%kill_strategy3D_tbox
call b%kill_general_tbox
deallocate(scores, sigma2_noise, den_pft)

if( failures > 0 )then
    write(logfhandle,'(a,i3)') 'SIMPLE_TEST_CONTINUOUS_INPLANE_HYBRID_GRAD FAILURES: ', failures
    error stop 1
endif
write(logfhandle,'(a)') 'SIMPLE_TEST_CONTINUOUS_INPLANE_HYBRID_GRAD: PASS'

contains

    subroutine parse_arguments(args, vol_file, mskdiam, smpd, lp, truth_angle, shift_truth)
        type(cmdline), intent(inout) :: args
        character(len=*), intent(out) :: vol_file
        real, intent(out) :: mskdiam, smpd, lp, truth_angle, shift_truth(2)
        vol_file    = ''
        mskdiam     = 120.
        smpd        = 1.3
        lp          = 8.
        truth_angle = 37.
        shift_truth = [2., -1.5]
        call args%parse_oldschool
        if( .not. args%defined('vol1') )then
            write(logfhandle,'(a)') 'simple_test_continuous_inplane_hybrid_grad vol1=xx '// &
                &'mskdiam=120 smpd=1.3 lp=8 angerr=37 xsh=2 ysh=-1.5'
            error stop 'vol1 is required'
        endif
        vol_arg  = args%get_carg('vol1')
        vol_file = vol_arg%to_char()
        if( args%defined('mskdiam') ) mskdiam        = args%get_rarg('mskdiam')
        if( args%defined('smpd')    ) smpd           = args%get_rarg('smpd')
        if( args%defined('lp')      ) lp             = args%get_rarg('lp')
        if( args%defined('angerr')  ) truth_angle    = args%get_rarg('angerr')
        if( args%defined('xsh')     ) shift_truth(1) = args%get_rarg('xsh')
        if( args%defined('ysh')     ) shift_truth(2) = args%get_rarg('ysh')
        if( mskdiam <= 0. .or. smpd <= 0. .or. lp <= 2.*smpd )then
            error stop 'invalid mskdiam, smpd, or low-pass limit'
        endif
    end subroutine parse_arguments

    subroutine fd_gradient_check(b, at_shift, at_rotind, max_error, tol, ok)
        type(builder), intent(inout) :: b
        real(dp), intent(in)  :: at_shift(2), at_rotind
        real(dp), intent(out) :: max_error, tol
        logical,  intent(out) :: ok
        real(dp), parameter :: FD_STEP = 1.d-3
        real(dp) :: f0, g0(3), gtmp(3), fd_plus, fd_minus
        real(dp) :: probe_shift(2), probe_rotind
        integer  :: idim
        call b%pftc%gen_hybrid_grad_at_angle(1, 1, at_shift, at_rotind, f0, g0)
        ok        = ieee_is_finite(f0) .and. all(ieee_is_finite(g0))
        max_error = 0.d0
        do idim = 1, 3
            probe_shift  = at_shift
            probe_rotind = at_rotind
            if( idim <= 2 )then
                probe_shift(idim) = probe_shift(idim) + FD_STEP
            else
                probe_rotind = probe_rotind + FD_STEP
            endif
            call b%pftc%gen_hybrid_grad_at_angle(1, 1, probe_shift, probe_rotind, fd_plus, gtmp)
            probe_shift  = at_shift
            probe_rotind = at_rotind
            if( idim <= 2 )then
                probe_shift(idim) = probe_shift(idim) - FD_STEP
            else
                probe_rotind = probe_rotind - FD_STEP
            endif
            call b%pftc%gen_hybrid_grad_at_angle(1, 1, probe_shift, probe_rotind, fd_minus, gtmp)
            max_error = max(max_error, abs((fd_plus-fd_minus)/(2.d0*FD_STEP) - g0(idim)))
        enddo
        tol = max(GRAD_FD_RTOL*maxval(abs(g0)), GRAD_FD_AFLOOR)
        ok  = ok .and. ieee_is_finite(max_error) .and. max_error <= tol
    end subroutine fd_gradient_check

end program simple_test_continuous_inplane_hybrid_grad
