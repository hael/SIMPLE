module simple_pcg_pose_polisher
use ieee_arithmetic, only: ieee_is_finite
use simple_core_module_api
use simple_builder, only: builder
use simple_cmdline, only: cmdline
use simple_parameters, only: parameters
use simple_image, only: image
use simple_matcher_ptcl_io, only: prepimgbatch, discrete_read_imgbatch, &
    &discrete_read_imgbatch_source, killimgbatch
use simple_ptcl_cache, only: ptcl_cache_in_use, ptcl_cache_assert_ready, ptcl_cache_read_batch
use simple_sigma2_files, only: load_sigma2_groups
use simple_math_ft, only: resample_sigma2
use simple_refine3D_fnames, only: refine3D_state_halfvol_fname
use simple_reconstructor_pcg, only: reconstructor_pcg, pcg_fourier_workspace, &
    &SHIFT_LM_ACCEPTED_IMPROVEMENT, SHIFT_LM_FINITE_NO_IMPROVEMENT, &
    &SHIFT_LM_NO_RELIABLE_UPDATE, SHIFT_LM_STEP_BOUND_REJECTED, &
    &SHIFT_LM_INVALID_NUMERICS, SHIFT_LM_ITERATION_LIMIT, &
    &POSE_LM_ACCEPTED_IMPROVEMENT, POSE_LM_FINITE_NO_IMPROVEMENT, &
    &POSE_LM_NO_RELIABLE_UPDATE, POSE_LM_STEP_BOUND_REJECTED, &
    &POSE_LM_INVALID_NUMERICS, POSE_LM_ITERATION_LIMIT
implicit none
private
#include "simple_local_flags.inc"
public :: pcg_pose_polish_summary, polish_fixed_volume_shifts, polish_fixed_volume_poses
public :: execute_final_pcg_pose_polish

integer, parameter :: POLISH_BATCH_SIZE = min(32,MAXIMGBATCHSZ)
integer, parameter :: POLISH_MAX_ITERATIONS = 8

type :: pcg_pose_polish_summary
    integer :: nparticles = 0
    integer :: nimproved = 0
    integer :: nunchanged = 0
    integer :: nunreliable = 0
    integer :: nstep_bound = 0
    integer :: ninvalid = 0
    integer :: niteration_limit = 0
    integer :: naccepted_steps = 0
    integer :: nattempted_steps = 0
    real(dp) :: objective_before = 0._dp
    real(dp) :: objective_after = 0._dp
    real(dp) :: max_trial_step = 0._dp
    real(dp) :: max_rotation_step = 0._dp
    real(dp) :: max_shift_step = 0._dp
    integer :: nstencil_switches = 0
end type pcg_pose_polish_summary

contains

subroutine polish_fixed_volume_shifts(workspace, rotmats, observed, shifts, &
    &max_iterations, statuses, summary, transfers)
    type(pcg_fourier_workspace), intent(in) :: workspace
    real(dp), intent(in) :: rotmats(:,:,:)
    complex, intent(in) :: observed(:,:,:)
    real(dp), intent(inout) :: shifts(:,:)
    integer, intent(in) :: max_iterations
    integer, intent(out) :: statuses(:)
    type(pcg_pose_polish_summary), intent(out) :: summary
    complex, optional, intent(in) :: transfers(:,:,:)
    real(dp), allocatable :: accepted_objectives(:)
    real(dp) :: input_shift(2), max_trial_step
    integer :: iparticle, naccepted, nattempted, status

    call validate_batch_shapes(workspace, rotmats, observed, shifts, statuses, max_iterations, transfers)
    summary = pcg_pose_polish_summary()
    summary%nparticles = size(shifts,2)
    allocate(accepted_objectives(0:max_iterations))

    do iparticle = 1, summary%nparticles
        ! Preserve the incoming pose unless LM returns an accepted improvement.
        input_shift = shifts(:,iparticle)
        if( present(transfers) )then
            call workspace%refine_shift_lm(rotmats(:,:,iparticle), observed(:,:,iparticle), &
                &shifts(:,iparticle), max_iterations, accepted_objectives, naccepted, status, &
                &nattempted, max_trial_step, transfers(:,:,iparticle))
        else
            call workspace%refine_shift_lm(rotmats(:,:,iparticle), observed(:,:,iparticle), &
                &shifts(:,iparticle), max_iterations, accepted_objectives, naccepted, status, &
                &nattempted, max_trial_step)
        endif
        statuses(iparticle) = status
        summary%naccepted_steps = summary%naccepted_steps + naccepted
        summary%nattempted_steps = summary%nattempted_steps + nattempted
        summary%max_trial_step = max(summary%max_trial_step,max_trial_step)

        if( ieee_is_finite(accepted_objectives(0)) )then
            summary%objective_before = summary%objective_before + accepted_objectives(0)
        endif
        select case(status)
        case(SHIFT_LM_ACCEPTED_IMPROVEMENT)
            if( naccepted < 1 ) error stop 'accepted shift polish has no accepted LM step'
            if( .not. ieee_is_finite(accepted_objectives(naccepted)) ) &
                &error stop 'accepted shift polish has a non-finite objective'
            if( accepted_objectives(naccepted) >= accepted_objectives(0) ) &
                &error stop 'accepted shift polish did not reduce the objective'
            summary%nimproved = summary%nimproved + 1
            summary%objective_after = summary%objective_after + accepted_objectives(naccepted)
        case(SHIFT_LM_FINITE_NO_IMPROVEMENT)
            ! All nonaccepted terminal outcomes restore the caller's shift.
            shifts(:,iparticle) = input_shift
            summary%nunchanged = summary%nunchanged + 1
            summary%objective_after = summary%objective_after + accepted_objectives(0)
        case(SHIFT_LM_NO_RELIABLE_UPDATE)
            shifts(:,iparticle) = input_shift
            summary%nunreliable = summary%nunreliable + 1
            summary%objective_after = summary%objective_after + accepted_objectives(0)
        case(SHIFT_LM_STEP_BOUND_REJECTED)
            shifts(:,iparticle) = input_shift
            summary%nstep_bound = summary%nstep_bound + 1
            summary%objective_after = summary%objective_after + accepted_objectives(0)
        case(SHIFT_LM_INVALID_NUMERICS)
            shifts(:,iparticle) = input_shift
            summary%ninvalid = summary%ninvalid + 1
            if( ieee_is_finite(accepted_objectives(0)) ) &
                &summary%objective_after = summary%objective_after + accepted_objectives(0)
        case(SHIFT_LM_ITERATION_LIMIT)
            shifts(:,iparticle) = input_shift
            summary%niteration_limit = summary%niteration_limit + 1
            summary%objective_after = summary%objective_after + accepted_objectives(0)
        case default
            error stop 'shift polish returned an unknown LM outcome'
        end select
    enddo
    deallocate(accepted_objectives)
end subroutine polish_fixed_volume_shifts

subroutine polish_fixed_volume_poses(workspace, rotmats, observed, shifts, rotation_scale, &
    &max_iterations, statuses, summary, transfers)
    type(pcg_fourier_workspace), intent(in) :: workspace
    real(dp), intent(inout) :: rotmats(:,:,:), shifts(:,:)
    complex, intent(in) :: observed(:,:,:)
    real(dp), intent(in) :: rotation_scale
    integer, intent(in) :: max_iterations
    integer, intent(out) :: statuses(:)
    type(pcg_pose_polish_summary), intent(out) :: summary
    complex, optional, intent(in) :: transfers(:,:,:)
    real(dp), allocatable :: accepted_objectives(:)
    real(dp) :: input_rotmat(3,3), input_shift(2), max_rotation_step, max_shift_step
    integer :: iparticle, naccepted, nattempted, status, nstencil_switches

    call validate_batch_shapes(workspace,rotmats,observed,shifts,statuses,max_iterations,transfers)
    if( rotation_scale <= 0._dp .or. .not. ieee_is_finite(rotation_scale) ) &
        &error stop 'pose polish requires a positive finite rotation scale'
    summary = pcg_pose_polish_summary()
    summary%nparticles = size(shifts,2)
    allocate(accepted_objectives(0:max_iterations))

    do iparticle = 1, summary%nparticles
        ! A particle pose is one transaction: retain both R and t or restore both.
        input_rotmat = rotmats(:,:,iparticle)
        input_shift = shifts(:,iparticle)
        if( present(transfers) )then
            call workspace%refine_pose_lm(rotmats(:,:,iparticle),observed(:,:,iparticle), &
                &shifts(:,iparticle),rotation_scale,max_iterations,accepted_objectives,naccepted, &
                &status,nattempted,max_rotation_step,max_shift_step,nstencil_switches, &
                &transfers(:,:,iparticle))
        else
            call workspace%refine_pose_lm(rotmats(:,:,iparticle),observed(:,:,iparticle), &
                &shifts(:,iparticle),rotation_scale,max_iterations,accepted_objectives,naccepted, &
                &status,nattempted,max_rotation_step,max_shift_step,nstencil_switches)
        endif
        statuses(iparticle) = status
        summary%naccepted_steps = summary%naccepted_steps+naccepted
        summary%nattempted_steps = summary%nattempted_steps+nattempted
        summary%max_rotation_step = max(summary%max_rotation_step,max_rotation_step)
        summary%max_shift_step = max(summary%max_shift_step,max_shift_step)
        summary%nstencil_switches = summary%nstencil_switches+nstencil_switches
        if( ieee_is_finite(accepted_objectives(0)) ) &
            &summary%objective_before = summary%objective_before+accepted_objectives(0)

        select case(status)
        case(POSE_LM_ACCEPTED_IMPROVEMENT)
            if( naccepted < 1 ) error stop 'accepted pose polish has no accepted LM step'
            if( .not. ieee_is_finite(accepted_objectives(naccepted)) ) &
                &error stop 'accepted pose polish has a non-finite objective'
            if( accepted_objectives(naccepted) >= accepted_objectives(0) ) &
                &error stop 'accepted pose polish did not reduce the objective'
            summary%nimproved = summary%nimproved+1
            summary%objective_after = summary%objective_after+accepted_objectives(naccepted)
        case(POSE_LM_FINITE_NO_IMPROVEMENT)
            summary%nunchanged = summary%nunchanged+1
            call restore_input_pose
        case(POSE_LM_NO_RELIABLE_UPDATE)
            summary%nunreliable = summary%nunreliable+1
            call restore_input_pose
        case(POSE_LM_STEP_BOUND_REJECTED)
            summary%nstep_bound = summary%nstep_bound+1
            call restore_input_pose
        case(POSE_LM_INVALID_NUMERICS)
            summary%ninvalid = summary%ninvalid+1
            rotmats(:,:,iparticle) = input_rotmat
            shifts(:,iparticle) = input_shift
            if( ieee_is_finite(accepted_objectives(0)) ) &
                &summary%objective_after = summary%objective_after+accepted_objectives(0)
        case(POSE_LM_ITERATION_LIMIT)
            summary%niteration_limit = summary%niteration_limit+1
            call restore_input_pose
        case default
            error stop 'pose polish returned an unknown LM outcome'
        end select
    enddo
    deallocate(accepted_objectives)

contains

    subroutine restore_input_pose
        rotmats(:,:,iparticle) = input_rotmat
        shifts(:,iparticle) = input_shift
        summary%objective_after = summary%objective_after+accepted_objectives(0)
    end subroutine restore_input_pose

end subroutine polish_fixed_volume_poses

subroutine execute_final_pcg_pose_polish(params, build, cline, total_summary)
    type(parameters), intent(inout) :: params
    type(builder), intent(inout) :: build
    class(cmdline), intent(inout) :: cline
    type(pcg_pose_polish_summary), intent(out) :: total_summary
    type(pcg_pose_polish_summary) :: batch_summary
    type(reconstructor_pcg) :: pcgop
    type(pcg_fourier_workspace) :: workspace
    type(image) :: halfmap, mskimg
    type(ctfparams) :: ctfparms
    type(ori) :: orientation
    type(string) :: halfmap_name
    integer, allocatable :: selected_pinds(:), half_pinds(:), statuses(:)
    complex, allocatable :: observed(:,:,:), transfers(:,:,:)
    real, allocatable :: sig2(:,:), volume(:,:,:)
    logical, allocatable :: normalization_mask(:,:,:)
    real(dp), allocatable :: rotmats(:,:,:), shifts(:,:)
    integer :: state, eo, ibatch, batchlims(2), batchsz, i, iptcl, nselected
    integer :: lims2(2,2), r, kfromto(2)
    real :: crop_factor, sdev_noise, edge_mean
    real(dp) :: rotation_scale
    logical :: l_cached, l_sigma_loaded

    if( trim(params%pcg_pose_polish) /= 'yes' )then
        total_summary = pcg_pose_polish_summary()
        return
    endif
    if( trim(params%rec_backend) /= 'pcg' )then
        THROW_HARD('final PCG pose polish requires an executed PCG reconstruction route')
    endif
    if( params%box < 1 .or. params%box_crop < 1 )then
        THROW_HARD('final PCG pose polish requires positive native and cropped boxes')
    endif

    total_summary = pcg_pose_polish_summary()
    nselected = 0
    call build%spproj_field%sample4rec([params%fromp,params%top],nselected,selected_pinds)
    if( nselected < 1 ) THROW_HARD('final PCG pose polish found no active reconstruction particles')
    l_sigma_loaded = .false.
    call ptcl_cache_assert_ready(params,build)
    l_cached = ptcl_cache_in_use(params,build)
    if( params%cc_objfun == OBJFUN_EUCLID )then
        l_sigma_loaded = allocated(build%esig%sigma2_noise)
        if( .not. l_sigma_loaded )then
            call load_sigma2_groups(params,build%pftc,build%esig,build%spproj_field,cline,l_sigma_loaded)
        endif
        if( .not. l_sigma_loaded ) THROW_HARD('final PCG pose polish requires sigma2 files for objfun=euclid')
    endif
    if( l_cached )then
        call prepimgbatch(params,build,POLISH_BATCH_SIZE,box=params%box_crop,smpd=params%smpd_crop)
    else
        call prepimgbatch(params,build,POLISH_BATCH_SIZE)
        ! Distributed masters own the project but need not build the general-toolbox mask.
        if( allocated(build%lmsk) )then
            allocate(normalization_mask,source=build%lmsk)
        else
            call mskimg%disc([params%box,params%box,1],params%smpd,params%msk,normalization_mask)
        endif
        if( .not. allocated(normalization_mask) )then
            THROW_HARD('final PCG pose polish could not build its normalization mask')
        endif
        if( any(shape(normalization_mask) /= [params%box,params%box,1]) )then
            THROW_HARD('final PCG pose polish normalization mask has the wrong shape')
        endif
    endif
    call orientation%new(.false.)
    crop_factor = real(params%box_crop) / real(params%box)
    ! One rotation step moves the mask edge by at most about one pixel.
    rotation_scale = 1._dp/max(real(params%msk_crop,dp),1._dp)

    do state = 1, params%nstates
        do eo = 0, 1
            call collect_state_half_particles(state,eo,half_pinds)
            if( size(half_pinds) == 0 )then
                deallocate(half_pinds)
                cycle
            endif
            ! Use the final FSC-regularized half-map, not its unfiltered diagnostic copy.
            halfmap_name = refine3D_state_halfvol_fname(state,merge('odd ','even',eo==1))
            if( .not. file_exists(halfmap_name) ) THROW_HARD('final PCG pose polish half-map is missing')
            call halfmap%new([params%box_crop,params%box_crop,params%box_crop],params%smpd_crop)
            call halfmap%read(halfmap_name)
            if( halfmap%get_box() /= params%box_crop )then
                THROW_HARD('final PCG pose polish half-map box does not match box_crop')
            endif
            if( abs(halfmap%get_smpd()-params%smpd_crop) > 1.e-4*max(1.,params%smpd_crop) )then
                THROW_HARD('final PCG pose polish half-map sampling does not match smpd_crop')
            endif
            volume = halfmap%get_rmat()
            call pcgop%new(params%box_crop,params%smpd_crop,1.e-3)
            call pcgop%set_volume(volume)
            call pcgop%begin_fourier_workspace(workspace)
            call workspace%set_shell_range(params%kfromto)
            lims2 = pcgop%get_lims2()
            r = lims2(1,2)

            do ibatch = 1, size(half_pinds), POLISH_BATCH_SIZE
                batchlims = [ibatch,min(size(half_pinds),ibatch+POLISH_BATCH_SIZE-1)]
                batchsz = batchlims(2)-batchlims(1)+1
                allocate(observed(lims2(1,1):lims2(1,2),lims2(2,1):lims2(2,2),batchsz))
                allocate(transfers(lims2(1,1):lims2(1,2),lims2(2,1):lims2(2,2),batchsz))
                allocate(sig2(0:r,batchsz),source=1.0)
                allocate(rotmats(3,3,batchsz),shifts(2,batchsz),statuses(batchsz))
                if( params%cc_objfun == OBJFUN_EUCLID )then
                    kfromto = build%esig%get_kfromto()
                    do i = 1, batchsz
                        iptcl = half_pinds(batchlims(1)+i-1)
                        call resample_sigma2(kfromto(1),kfromto(2), &
                            &build%esig%sigma2_noise(kfromto(1):kfromto(2),iptcl),r,1.0,sig2(:,i))
                    enddo
                endif
                if( l_cached )then
                    call ptcl_cache_read_batch(params,build,size(half_pinds),half_pinds,batchlims)
                else if( params%l_ptcl_src_den )then
                    call discrete_read_imgbatch_source(params,build,'den',size(half_pinds),half_pinds, &
                        &batchlims,build%imgbatch(:batchsz))
                else
                    call discrete_read_imgbatch(params,build,size(half_pinds),half_pinds,batchlims)
                endif
                sdev_noise = 0.0
                edge_mean = 0.0
                do i = 1, batchsz
                    iptcl = half_pinds(batchlims(1)+i-1)
                    if( .not. l_cached ) call build%imgbatch(i)%norm_noise(normalization_mask,sdev_noise)
                    call build%imgbatch(i)%taper_edges_particle(nint(COSMSKHALFWIDTH),edge_mean)
                    call build%imgbatch(i)%fft
                    observed(:,:,i) = pcgop%whiten_observation( &
                        &pcgop%extract_native_plane(build%imgbatch(i)),sig2(:,i))
                    call build%spproj_field%get_ori(iptcl,orientation)
                    rotmats(:,:,i) = real(orientation%get_mat(),dp)
                    shifts(:,i) = real(build%spproj_field%get_2Dshift(iptcl)*crop_factor,dp)
                    ctfparms = build%spproj%get_ctfparams(params%oritype,iptcl)
                    ctfparms%smpd = params%smpd_crop
                    transfers(:,:,i) = pcgop%build_transfer(ctfparms,[0.,0.],sig2(:,i))
                enddo
                call polish_fixed_volume_poses(workspace,rotmats,observed,shifts,rotation_scale, &
                    &POLISH_MAX_ITERATIONS,statuses,batch_summary,transfers)
                call add_summary(total_summary,batch_summary)
                do i = 1, batchsz
                    iptcl = half_pinds(batchlims(1)+i-1)
                    ! Persist the complete pose only after an accepted objective reduction.
                    if( statuses(i) == POSE_LM_ACCEPTED_IMPROVEMENT )then
                        call build%spproj_field%set_euler(iptcl,real(dm2euler(rotmats(:,:,i)),kind=sp))
                        call build%spproj_field%set_shift(iptcl,real(shifts(:,i)/real(crop_factor,dp)))
                    endif
                enddo
                deallocate(observed,transfers,sig2,rotmats,shifts,statuses)
            enddo
            call workspace%kill
            call pcgop%kill
            call halfmap%kill
            call halfmap_name%kill
            deallocate(volume,half_pinds)
        enddo
    enddo
    call orientation%kill
    call killimgbatch(build)
    call mskimg%kill
    if( allocated(normalization_mask) ) deallocate(normalization_mask)
    deallocate(selected_pinds)
    write(logfhandle,'(A,I0,A,I0,A,I0,A,I0)') '>>> PCG POSE POLISH: PARTICLES ', &
        &total_summary%nparticles, ' IMPROVED ', total_summary%nimproved, &
        &' UNCHANGED ', total_summary%nunchanged, ' UNRELIABLE ', total_summary%nunreliable
    write(logfhandle,'(A,I0,A,I0,A,I0)') '>>> PCG POSE POLISH: STEP-BOUND ', &
        &total_summary%nstep_bound, ' INVALID ', total_summary%ninvalid, &
        &' ITERATION-LIMIT ', total_summary%niteration_limit
    write(logfhandle,'(A,I0,A,I0,A,I0)') '>>> PCG POSE POLISH: ATTEMPTED STEPS ', &
        &total_summary%nattempted_steps, ' ACCEPTED STEPS ', total_summary%naccepted_steps, &
        &' STENCIL SWITCHES ', total_summary%nstencil_switches
    write(logfhandle,'(A,ES12.5,A,ES12.5)') '>>> PCG POSE POLISH: MAX ROTATION/SHIFT STEP ', &
        &total_summary%max_rotation_step, ' / ', total_summary%max_shift_step
    write(logfhandle,'(A,ES12.5,A,ES12.5)') '>>> PCG POSE POLISH: OBJECTIVE ', &
        &total_summary%objective_before, ' -> ', total_summary%objective_after

contains

    subroutine collect_state_half_particles(state_here,eo_here,pinds)
        integer, intent(in) :: state_here, eo_here
        integer, allocatable, intent(out) :: pinds(:)
        integer :: ip, n
        n = 0
        do ip = 1, nselected
            if( build%spproj_field%get_state(selected_pinds(ip)) /= state_here ) cycle
            if( build%spproj_field%get_eo(selected_pinds(ip)) /= eo_here ) cycle
            n = n + 1
        enddo
        allocate(pinds(n))
        n = 0
        do ip = 1, nselected
            if( build%spproj_field%get_state(selected_pinds(ip)) /= state_here ) cycle
            if( build%spproj_field%get_eo(selected_pinds(ip)) /= eo_here ) cycle
            n = n + 1
            pinds(n) = selected_pinds(ip)
        enddo
    end subroutine collect_state_half_particles

    subroutine add_summary(total,part)
        type(pcg_pose_polish_summary), intent(inout) :: total
        type(pcg_pose_polish_summary), intent(in) :: part
        total%nparticles = total%nparticles + part%nparticles
        total%nimproved = total%nimproved + part%nimproved
        total%nunchanged = total%nunchanged + part%nunchanged
        total%nunreliable = total%nunreliable + part%nunreliable
        total%nstep_bound = total%nstep_bound + part%nstep_bound
        total%ninvalid = total%ninvalid + part%ninvalid
        total%niteration_limit = total%niteration_limit + part%niteration_limit
        total%naccepted_steps = total%naccepted_steps + part%naccepted_steps
        total%nattempted_steps = total%nattempted_steps + part%nattempted_steps
        total%objective_before = total%objective_before + part%objective_before
        total%objective_after = total%objective_after + part%objective_after
        total%max_trial_step = max(total%max_trial_step,part%max_trial_step)
        total%max_rotation_step = max(total%max_rotation_step,part%max_rotation_step)
        total%max_shift_step = max(total%max_shift_step,part%max_shift_step)
        total%nstencil_switches = total%nstencil_switches+part%nstencil_switches
    end subroutine add_summary

end subroutine execute_final_pcg_pose_polish

subroutine validate_batch_shapes(workspace, rotmats, observed, shifts, statuses, max_iterations, transfers)
    type(pcg_fourier_workspace), intent(in) :: workspace
    real(dp), intent(in) :: rotmats(:,:,:), shifts(:,:)
    complex, intent(in) :: observed(:,:,:)
    integer, intent(in) :: statuses(:), max_iterations
    complex, optional, intent(in) :: transfers(:,:,:)
    integer :: lims2(2,2), nparticles

    nparticles = size(shifts,2)
    lims2 = workspace%get_lims2()
    if( max_iterations < 1 ) error stop 'pose-polish batch requires at least one LM iteration'
    if( size(shifts,1) /= 2 ) error stop 'pose-polish batch requires two shift coordinates'
    if( size(rotmats,1) /= 3 .or. size(rotmats,2) /= 3 ) &
        &error stop 'pose-polish batch requires 3-by-3 rotation matrices'
    if( size(rotmats,3) /= nparticles ) error stop 'pose-polish batch rotation count mismatch'
    if( size(observed,1) /= lims2(1,2)-lims2(1,1)+1 .or. &
        &size(observed,2) /= lims2(2,2)-lims2(2,1)+1 ) &
        &error stop 'pose-polish batch Fourier-plane shape mismatch'
    if( size(observed,3) /= nparticles ) error stop 'pose-polish batch observation count mismatch'
    if( present(transfers) )then
        if( any(shape(transfers) /= shape(observed)) ) &
            &error stop 'pose-polish batch transfer-plane shape mismatch'
    endif
    if( size(statuses) /= nparticles ) error stop 'pose-polish batch status count mismatch'
end subroutine validate_batch_shapes

end module simple_pcg_pose_polisher
