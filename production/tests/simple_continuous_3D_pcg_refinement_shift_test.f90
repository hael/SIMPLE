module continuous_3D_pcg_refinement_shift_test
use continuous_3D_pcg_refinement_test_helpers, only: assert_true, build_truth_volume, TRUTH_VOLUME_BOX
use ieee_arithmetic, only: ieee_is_finite
use simple_defs, only: dp, sp, DPI
use simple_reconstructor_pcg, only: pcg_fourier_workspace, reconstructor_pcg, &
    &SHIFT_LM_ACCEPTED_IMPROVEMENT, SHIFT_LM_FINITE_NO_IMPROVEMENT, SHIFT_LM_NO_RELIABLE_UPDATE
use simple_type_defs, only: ctfparams, CTFFLAG_NO, CTFFLAG_YES
implicit none
private
public :: run_shift_gradient

integer, parameter :: N_FD_STEPS = 3
integer, parameter :: MAX_RECOVERY_ITERATIONS = 20
real(dp), parameter :: FD_STEPS(N_FD_STEPS) = [2.e-2_dp,1.e-2_dp,5.e-3_dp]
real(dp), parameter :: DIRECTIONAL_TOL = 3.e-3_dp
real(dp), parameter :: ADJOINT_TOL = 2.e-6_dp
real(dp), parameter :: OBJECTIVE_TOL = 3.e-3_dp
real(dp), parameter :: SIGN_TOL = 3.e-5_dp
real(dp), parameter :: RECOVERY_TOL = 2.e-3_dp

contains

subroutine run_shift_gradient()
    type(pcg_fourier_workspace) :: workspace
    type(reconstructor_pcg) :: pcgop
    real, allocatable :: phantom(:,:,:)
    complex, allocatable :: observed(:,:), residual(:,:), residual_minus(:,:), residual_plus(:,:)
    complex, allocatable :: zero_model(:,:), jv(:,:), z(:,:), transfer(:,:), ctf_transfer(:,:)
    complex, allocatable :: raw_ctf_observed(:,:), weighted_observed(:,:)
    real, allocatable :: sig2(:)
    type(ctfparams) :: weighted_ctf
    real(dp) :: rotmat(3,3), true_shift(2), estimate(2), direction(2), gradient(2), jhz(2)
    real(dp) :: objective, objective_minus, objective_plus, ignored_objective
    real(dp) :: directional_errors(N_FD_STEPS), objective_errors(N_FD_STEPS)
    real(dp) :: adjoint_error, sign_error, recovery_error, initial_error, inner_left, inner_right
    real(dp) :: max_trial_step, recovery_max_trial_step, weighted_gradient(2)
    real(dp) :: weighted_fd_error, weighted_recovery_error, whitening_error
    real(dp) :: accepted_objectives(0:MAX_RECOVERY_ITERATIONS)
    integer :: h, istep, k, lims2(2,2), naccepted, nattempted, status
    integer :: recovery_naccepted, recovery_nattempted

    call build_truth_volume(phantom)
    call pcgop%new(TRUTH_VOLUME_BOX,1._sp)
    call pcgop%set_volume(phantom)
    call pcgop%begin_fourier_workspace(workspace)
    lims2 = pcgop%get_lims2()
    allocate(observed(lims2(1,1):lims2(1,2),lims2(2,1):lims2(2,2)))
    allocate(residual, mold=observed)
    allocate(residual_minus, mold=observed)
    allocate(residual_plus, mold=observed)
    allocate(zero_model, mold=observed)
    allocate(jv, mold=observed)
    allocate(z, mold=observed)

    rotmat = 0._dp
    rotmat(1,1) = 1._dp
    rotmat(2,2) = 1._dp
    rotmat(3,3) = 1._dp
    true_shift = [0.35_dp,-0.27_dp]
    zero_model = cmplx(0.,0.)
    call workspace%shift_residual(rotmat,true_shift,zero_model,observed,ignored_objective)

    ! A separately coded real-space translation and direct DFT fixes the
    ! physical sign. It does not call the workspace gather or phase routine.
    sign_error = independent_shift_sign_error(pcgop,workspace,rotmat,true_shift,zero_model,observed,lims2)
    call assert_true(sign_error < SIGN_TOL, &
        &'positive build_transfer/shift-residual phase disagrees with the independent real-space DFT oracle')

    estimate = [0.12_dp,-0.08_dp]
    direction = [0.6_dp,-0.8_dp]
    call workspace%shift_residual(rotmat,estimate,observed,residual,objective)
    call workspace%shift_jvp(rotmat,estimate,direction,jv)
    do istep = 1, N_FD_STEPS
        call workspace%shift_residual(rotmat,estimate-FD_STEPS(istep)*direction,observed,&
            &residual_minus,objective_minus)
        call workspace%shift_residual(rotmat,estimate+FD_STEPS(istep)*direction,observed,&
            &residual_plus,objective_plus)
        directional_errors(istep) = plane_relative_error(&
            &(residual_plus-residual_minus)/real(2._dp*FD_STEPS(istep),sp),jv,lims2)
    enddo
    call assert_true(minval(directional_errors) < DIRECTIONAL_TOL, &
        &'shift residual Jv disagrees with centred directional differences')

    z = cmplx(0.,0.)
    do k = lims2(2,1), lims2(2,2)
        do h = lims2(1,1), lims2(1,2)
            if( h*h + k*k > (TRUTH_VOLUME_BOX/2)**2 ) cycle
            z(h,k) = cmplx(sin(0.17_dp*real(h,dp)+0.11_dp*real(k,dp)),&
                &cos(0.07_dp*real(h,dp)-0.13_dp*real(k,dp)),kind=sp)
        enddo
    enddo
    call workspace%shift_jhz(rotmat,estimate,z,jhz)
    inner_left = real_plane_inner_product(jv,z,lims2)
    inner_right = dot_product(direction,jhz)
    adjoint_error = abs(inner_left-inner_right) / max(1._dp,abs(inner_left),abs(inner_right))
    call assert_true(adjoint_error < ADJOINT_TOL, &
        &'real-parameter shift Jacobian adjoint identity failed')

    call workspace%shift_jhz(rotmat,estimate,residual,gradient)
    do istep = 1, N_FD_STEPS
        call workspace%shift_residual(rotmat,estimate-FD_STEPS(istep)*direction,observed,&
            &residual_minus,objective_minus)
        call workspace%shift_residual(rotmat,estimate+FD_STEPS(istep)*direction,observed,&
            &residual_plus,objective_plus)
        objective_errors(istep) = abs((objective_plus-objective_minus)/(2._dp*FD_STEPS(istep)) - &
            &dot_product(gradient,direction)) / max(1._dp,abs(dot_product(gradient,direction)))
    enddo
    call assert_true(minval(objective_errors) < OBJECTIVE_TOL, &
        &'shift objective gradient disagrees with centred directional differences')

    estimate = [-0.15_dp,0.12_dp]
    initial_error = sqrt(sum((estimate-true_shift)**2))
    call workspace%refine_shift_lm(rotmat,observed,estimate,MAX_RECOVERY_ITERATIONS,&
        &accepted_objectives,naccepted,status,nattempted,max_trial_step)
    recovery_error = sqrt(sum((estimate-true_shift)**2))
    call assert_true(naccepted > 0, 'shift recovery accepted no steps')
    call assert_true(nattempted >= naccepted, 'shift recovery reported fewer trials than accepted steps')
    call assert_true(status == SHIFT_LM_ACCEPTED_IMPROVEMENT, 'shift recovery returned the wrong terminal outcome')
    call assert_true(max_trial_step <= 1._dp+epsilon(1._dp), 'shift recovery exceeded the one-pixel step bound')
    call assert_true(all(accepted_objectives(1:naccepted) < accepted_objectives(0:naccepted-1)), &
        &'an accepted shift-refinement step did not lower the fully recomputed objective')
    call assert_true(recovery_error < initial_error .and. recovery_error < RECOVERY_TOL, &
        &'known injected shift was not recovered')
    call assert_true(all(ieee_is_finite(accepted_objectives(0:naccepted))), &
        &'shift recovery produced a non-finite accepted objective')
    recovery_naccepted = naccepted
    recovery_nattempted = nattempted
    recovery_max_trial_step = max_trial_step

    ! An exact solution is a finite no-op, not a fabricated update.
    estimate = true_shift
    call workspace%refine_shift_lm(rotmat,observed,estimate,MAX_RECOVERY_ITERATIONS,&
        &accepted_objectives,naccepted,status,nattempted,max_trial_step)
    call assert_true(status == SHIFT_LM_FINITE_NO_IMPROVEMENT .and. naccepted == 0, &
        &'an exact shift did not return finite_no_improvement')

    ! Exercise the production CTF/sigma transfer in the executed shift objective.
    allocate(sig2(0:lims2(1,2)))
    sig2 = [(1.0+0.03*real(h),h=0,lims2(1,2))]
    allocate(transfer,mold=observed)
    allocate(ctf_transfer,mold=observed)
    allocate(raw_ctf_observed,mold=observed)
    allocate(weighted_observed,mold=observed)
    weighted_ctf%smpd = 1.0
    weighted_ctf%kv = 300.0
    weighted_ctf%cs = 2.7
    weighted_ctf%fraca = 0.1
    weighted_ctf%dfx = 1.4
    weighted_ctf%dfy = 1.55
    weighted_ctf%angast = 23.0
    weighted_ctf%phshift = 0.0
    weighted_ctf%ctfflag = CTFFLAG_YES
    ctf_transfer = pcgop%build_transfer(weighted_ctf,[0.,0.])
    transfer = pcgop%build_transfer(weighted_ctf,[0.,0.],sig2)
    raw_ctf_observed = ctf_transfer*observed
    weighted_observed = pcgop%whiten_observation(raw_ctf_observed,sig2)
    whitening_error = plane_relative_error(weighted_observed,transfer*observed,lims2)
    call assert_true(whitening_error < ADJOINT_TOL, &
        &'production observation whitening disagrees with the weighted transfer')
    estimate = [0.09_dp,-0.04_dp]
    direction = [0.6_dp,-0.8_dp]
    call workspace%shift_objective_gradient(rotmat,estimate,weighted_observed, &
        &objective,weighted_gradient,transfer)
    call workspace%shift_objective_gradient(rotmat,estimate-FD_STEPS(2)*direction, &
        &weighted_observed,objective_minus,gradient,transfer)
    call workspace%shift_objective_gradient(rotmat,estimate+FD_STEPS(2)*direction, &
        &weighted_observed,objective_plus,gradient,transfer)
    weighted_fd_error = abs((objective_plus-objective_minus)/(2._dp*FD_STEPS(2)) - &
        &dot_product(weighted_gradient,direction)) / max(1._dp,abs(dot_product(weighted_gradient,direction)))
    call assert_true(weighted_fd_error < OBJECTIVE_TOL, &
        &'CTF/sigma-weighted shift gradient disagrees with centred differences')
    call workspace%refine_shift_lm(rotmat,weighted_observed,estimate,MAX_RECOVERY_ITERATIONS, &
        &accepted_objectives,naccepted,status,nattempted,max_trial_step,transfer)
    weighted_recovery_error = sqrt(sum((estimate-true_shift)**2))
    call assert_true(status == SHIFT_LM_ACCEPTED_IMPROVEMENT .and. weighted_recovery_error < RECOVERY_TOL, &
        &'CTF/sigma-weighted shift recovery failed')

    ! A zero Fourier volume has no observable shift direction.
    phantom = 0.
    call pcgop%set_volume(phantom)
    call pcgop%begin_fourier_workspace(workspace)
    observed = cmplx(0.,0.)
    estimate = [0.2_dp,-0.1_dp]
    call workspace%refine_shift_lm(rotmat,observed,estimate,MAX_RECOVERY_ITERATIONS,&
        &accepted_objectives,naccepted,status,nattempted,max_trial_step)
    call assert_true(status == SHIFT_LM_NO_RELIABLE_UPDATE .and. naccepted == 0, &
        &'an unobservable shift did not return no_reliable_update')

    write(*,'(a,3(es14.6,1x))') 'CONTINUOUS_3D_PCG_SHIFT min-Jv/adjoint/min-objective error: ', &
        &minval(directional_errors), adjoint_error, minval(objective_errors)
    write(*,'(a,3(es14.6,1x))') 'CONTINUOUS_3D_PCG_SHIFT Jv centred-FD curve: ', directional_errors
    write(*,'(a,3(es14.6,1x))') 'CONTINUOUS_3D_PCG_SHIFT objective centred-FD curve: ', objective_errors
    write(*,'(a,3(es14.6,1x),i0)') 'CONTINUOUS_3D_PCG_SHIFT sign/initial/final/accepted: ', &
        &sign_error, initial_error, recovery_error, recovery_naccepted
    write(*,'(a,es14.6,1x,i0)') 'CONTINUOUS_3D_PCG_SHIFT max-step/attempted: ', &
        &recovery_max_trial_step, recovery_nattempted
    write(*,'(a,3(es14.6,1x))') 'CONTINUOUS_3D_PCG_SHIFT weighted-gradient/recovery/whitening error: ', &
        &weighted_fd_error, weighted_recovery_error, whitening_error
    write(*,'(a)') 'CONTINUOUS_3D_PCG_SHIFT_GRADIENT: PASS'
    deallocate(observed,residual,residual_minus,residual_plus,zero_model,jv,z,phantom)
    deallocate(transfer,ctf_transfer,raw_ctf_observed,weighted_observed,sig2)
    call workspace%kill
    call pcgop%kill
end subroutine run_shift_gradient

function independent_shift_sign_error(pcgop,workspace,rotmat,shift,zero_model,shifted_model,lims2) result(max_error)
    type(reconstructor_pcg), intent(in) :: pcgop
    type(pcg_fourier_workspace), intent(in) :: workspace
    real(dp), intent(in) :: rotmat(3,3), shift(2)
    integer, intent(in) :: lims2(2,2)
    complex, intent(inout) :: zero_model(lims2(1,1):lims2(1,2),lims2(2,1):lims2(2,2))
    complex, intent(in) :: shifted_model(lims2(1,1):lims2(1,2),lims2(2,1):lims2(2,2))
    real(dp) :: max_error
    integer, parameter :: NMODES = 3
    integer, parameter :: modes(2,NMODES) = reshape([1,0, 0,1, 2,-1],[2,NMODES])
    type(ctfparams) :: ctfparms
    complex, allocatable :: transfer(:,:), scratch(:,:)
    complex(dp) :: oracle_ratio, workspace_ratio
    real(dp) :: ignored_objective
    integer :: h, imode, k

    ctfparms%ctfflag = CTFFLAG_NO
    allocate(transfer(lims2(1,1):lims2(1,2),lims2(2,1):lims2(2,2)))
    transfer = pcgop%build_transfer(ctfparms,real(shift,sp))
    allocate(scratch, mold=zero_model)
    scratch = cmplx(0.,0.)
    call workspace%shift_residual(rotmat,[0._dp,0._dp],scratch,zero_model,ignored_objective)
    max_error = 0._dp
    do imode = 1, NMODES
        h = modes(1,imode)
        k = modes(2,imode)
        oracle_ratio = translated_direct_dft_ratio(TRUTH_VOLUME_BOX,h,k,shift)
        max_error = max(max_error,relative_complex_error(cmplx(transfer(h,k),kind=dp),oracle_ratio))
        call assert_true(abs(zero_model(h,k)) > 1.e-6_sp, &
            &'shift sign oracle selected a near-zero workspace Fourier sample')
        workspace_ratio = cmplx(shifted_model(h,k),kind=dp) / cmplx(zero_model(h,k),kind=dp)
        max_error = max(max_error,relative_complex_error(workspace_ratio,oracle_ratio))
    enddo
    deallocate(transfer,scratch)
end function independent_shift_sign_error

function translated_direct_dft_ratio(box,h,k,shift) result(ratio)
    integer, intent(in) :: box, h, k
    real(dp), intent(in) :: shift(2)
    complex(dp) :: ratio, base_dft, shifted_dft, exponential
    real(dp) :: x, y, base_value, shifted_value, arg
    integer :: ix, iy
    base_dft = cmplx(0._dp,0._dp,kind=dp)
    shifted_dft = cmplx(0._dp,0._dp,kind=dp)
    do iy = 0, box-1
        y = real(iy,dp)
        do ix = 0, box-1
            x = real(ix,dp)
            base_value = independent_periodic_field(x,y,box)
            ! Positive SIMPLE shift is tested as the physical periodic
            ! resampling y(x)=x(x+t), which must have a positive Fourier phase.
            shifted_value = independent_periodic_field(x+shift(1),y+shift(2),box)
            arg = -2._dp*DPI*(real(h,dp)*x+real(k,dp)*y)/real(box,dp)
            exponential = cmplx(cos(arg),sin(arg),kind=dp)
            base_dft = base_dft + base_value*exponential
            shifted_dft = shifted_dft + shifted_value*exponential
        enddo
    enddo
    if( abs(base_dft) <= 100._dp*epsilon(1._dp) )then
        error stop 'independent direct-DFT oracle selected an absent Fourier mode'
    endif
    ratio = shifted_dft/base_dft
end function translated_direct_dft_ratio

pure real(dp) function independent_periodic_field(x,y,box) result(value)
    real(dp), intent(in) :: x, y
    integer, intent(in) :: box
    real(dp) :: scale
    scale = 2._dp*DPI/real(box,dp)
    value = 1.1_dp*cos(scale*x+0.2_dp) + 0.8_dp*sin(scale*y-0.4_dp) + &
        &0.6_dp*cos(scale*(2._dp*x-y)+0.3_dp)
end function independent_periodic_field

function plane_relative_error(actual,expected,lims2) result(error)
    integer, intent(in) :: lims2(2,2)
    complex, intent(in) :: actual(lims2(1,1):lims2(1,2),lims2(2,1):lims2(2,2))
    complex, intent(in) :: expected(lims2(1,1):lims2(1,2),lims2(2,1):lims2(2,2))
    real(dp) :: error, numerator, denominator
    integer :: h, k
    numerator = 0._dp
    denominator = 0._dp
    do k = lims2(2,1), lims2(2,2)
        do h = lims2(1,1), lims2(1,2)
            if( h*h + k*k > (TRUTH_VOLUME_BOX/2)**2 ) cycle
            numerator = numerator + real(abs(actual(h,k)-expected(h,k)),dp)**2
            denominator = denominator + real(abs(expected(h,k)),dp)**2
        enddo
    enddo
    error = sqrt(numerator/max(denominator,epsilon(1._dp)))
end function plane_relative_error

function real_plane_inner_product(left,right,lims2) result(product)
    integer, intent(in) :: lims2(2,2)
    complex, intent(in) :: left(lims2(1,1):lims2(1,2),lims2(2,1):lims2(2,2))
    complex, intent(in) :: right(lims2(1,1):lims2(1,2),lims2(2,1):lims2(2,2))
    real(dp) :: product
    integer :: h, k
    product = 0._dp
    do k = lims2(2,1), lims2(2,2)
        do h = lims2(1,1), lims2(1,2)
            if( h*h + k*k > (TRUTH_VOLUME_BOX/2)**2 ) cycle
            product = product + real(conjg(cmplx(left(h,k),kind=dp))*cmplx(right(h,k),kind=dp),dp)
        enddo
    enddo
end function real_plane_inner_product

pure real(dp) function relative_complex_error(actual,expected) result(error)
    complex(dp), intent(in) :: actual, expected
    error = abs(actual-expected)/max(1._dp,abs(actual),abs(expected))
end function relative_complex_error

end module continuous_3D_pcg_refinement_shift_test
