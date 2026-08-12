!@descr: shared-memory production strategy body for kernel PCG reconstruct3D
module simple_rec3D_pcg_strategy
use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
use simple_core_module_api
use simple_builder,           only: builder
use simple_cmdline,           only: cmdline
use simple_parameters,        only: parameters
use simple_reconstructor_pcg, only: reconstructor_pcg, pcg_solver_outcome, PCG_OP_KERNEL
use simple_matcher_ptcl_io,   only: prepimgbatch, discrete_read_imgbatch, &
    &discrete_read_imgbatch_source, killimgbatch
use simple_sigma2_files,      only: load_sigma2_groups
use simple_math_ft,           only: upsample_sigma2
use simple_fsc,               only: phase_rand_fsc, fsc_area_score_result
use simple_image,             only: image
use simple_image_msk,         only: image_msk
use simple_refine3D_fnames,   only: refine3D_state_halfvol_fname, refine3D_state_vol_fname, &
    &refine3D_fsc_fname, refine3D_resolution_txt_fbody, refine3D_pcg_raw_accum_fname
implicit none

public :: execute_rec3D_pcg_shared, execute_rec3D_pcg_worker, execute_rec3D_pcg_distributed_master
private
#include "simple_local_flags.inc"

real,    parameter :: PCG_LAMBDA = 1.0e-3

contains

    subroutine execute_rec3D_pcg_shared( params, build, cline )
        type(parameters), intent(inout) :: params
        type(builder),    intent(inout) :: build
        class(cmdline),   intent(inout) :: cline
        type(image) :: half_even, half_odd, merged
        type(string) :: fname_even, fname_odd, fname_vol, fname_fsc
        integer, allocatable :: selected_pinds(:), half_pinds(:)
        real, allocatable :: fsc(:), res0143s(:)
        logical, allocatable :: state_written(:)
        integer :: nselected, state, n_state, n_even, n_odd, iptcl, istate
        real :: res05, cfar
        integer(timer_int_kind) :: t_state_phase
        real(dp) :: time_map_output, time_fsc_output
        logical :: l_sigma_loaded

        call validate_supported_mode()
        nselected = 0
        call build%spproj_field%sample4rec([params%fromp,params%top], nselected, selected_pinds)
        if( nselected < 1 ) THROW_HARD('no active particles selected for PCG reconstruct3D')

        if( params%cc_objfun == OBJFUN_EUCLID )then
            call load_sigma2_groups(params, build%pftc, build%esig, build%spproj_field, cline, l_sigma_loaded)
            if( .not. l_sigma_loaded ) THROW_HARD('PCG objfun=euclid requires sigma2 files')
        endif

        allocate(res0143s(params%nstates), source=0.0)
        allocate(state_written(params%nstates), source=.false.)
        call prepimgbatch(params, build, MAXIMGBATCHSZ)

        do state = 1, params%nstates
            n_state = count_state(state)
            n_even = count_state_half(state, 0)
            n_odd  = count_state_half(state, 1)
            if( n_state == 0 )then
                write(logfhandle,'(A,I0,A)') '>>> PCG RECONSTRUCT3D: STATE ', state, ' HAS NO SELECTED PARTICLES; SKIPPING'
                cycle
            endif
            if( n_even + n_odd /= n_state ) THROW_HARD('PCG reconstruct3D found invalid halfset labels')
            if( n_even < 1 .or. n_odd < 1 ) THROW_HARD('PCG reconstruct3D requires particles in both halfsets')

            call collect_state_half(state, 0, n_even, half_pinds)
            call solve_state_half(state, 'even', half_pinds, half_even)
            deallocate(half_pinds)
            call collect_state_half(state, 1, n_odd, half_pinds)
            call solve_state_half(state, 'odd', half_pinds, half_odd)
            deallocate(half_pinds)

            fname_even = refine3D_state_halfvol_fname(state, 'even')
            fname_odd  = refine3D_state_halfvol_fname(state, 'odd')
            fname_vol  = refine3D_state_vol_fname(state)
            fname_fsc  = refine3D_fsc_fname(state)
            t_state_phase = tic()
            call half_even%write(fname_even, del_if_exists=.true.)
            call half_odd%write(fname_odd, del_if_exists=.true.)
            call merged%copy(half_even)
            call merged%add(half_odd)
            call merged%mul(0.5)
            call merged%write(fname_vol, del_if_exists=.true.)
            time_map_output = real(toc(t_state_phase),dp)
            t_state_phase = tic()
            call calculate_state_fsc(state, half_even, half_odd, merged, fsc, res05, res0143s(state), cfar)
            call arr2file(fsc, fname_fsc)
            call write_fsc_summary(state, merged, fsc, res05, res0143s(state), cfar)
            time_fsc_output = real(toc(t_state_phase),dp)
            write(logfhandle,'(A,I0)') '>>> PCG SHARED OUTPUT PHASES: STATE ', state
            write(logfhandle,'(A,F9.3)') '    halfmap + merged-map output : ', time_map_output
            write(logfhandle,'(A,F9.3)') '    FSC + cFAR + summary        : ', time_fsc_output
            call flush(logfhandle)

            params%vols(state)      = fname_vol
            params%vols_even(state) = fname_even
            params%vols_odd(state)  = fname_odd
            call cline%set('vol'//int2str(state), fname_vol)
            state_written(state) = .true.
            call half_even%kill
            call half_odd%kill
            call merged%kill
            call fname_even%kill
            call fname_odd%kill
            call fname_vol%kill
            call fname_fsc%kill
            if( allocated(fsc) ) deallocate(fsc)
        enddo

        call killimgbatch(build)
        if( .not. any(state_written) ) THROW_HARD('PCG reconstruct3D produced no populated states')
        if( params%nstates == 1 )then
            call build%spproj_field%set_all2single('res', res0143s(1))
        else
            do iptcl = 1, build%spproj_field%get_noris()
                istate = build%spproj_field%get_state(iptcl)
                if( istate > 0 .and. istate <= params%nstates )then
                    if( state_written(istate) ) call build%spproj_field%set(iptcl, 'res', res0143s(istate))
                endif
            enddo
        endif
        call build%spproj%write_segment_inside(params%oritype, params%projfile)
        call register_project_outputs()

        deallocate(selected_pinds, res0143s, state_written)

    contains

        subroutine validate_supported_mode()
            if( params%nparts > 1 .or. cline%defined('part') )then
                THROW_HARD('shared PCG entry received distributed parameters')
            endif
            if( trim(params%pcgop) /= 'kernel' ) THROW_HARD('production rec_backend=pcg requires pcgop=kernel')
            if( params%maxits > 8 ) THROW_HARD('production rec_backend=pcg requires maxits<=8')
            if( trim(params%projrec) /= 'no' ) THROW_HARD('rec_backend=pcg does not yet support projrec=yes')
            if( params%box_crop /= params%box )then
                THROW_HARD('rec_backend=pcg box cropping awaits a scale-equivalence acceptance test')
            endif
            if( params%l_update_frac .or. params%l_trail_rec )then
                THROW_HARD('PCG fractional and trailing reconstruction require accumulator-domain integration')
            endif
            if( trim(params%conical_fsc) == 'yes' ) THROW_HARD('PCG conical FSC integration is not implemented')
            if( params%msk <= 0.5 .or. params%msk_crop <= 0.5 )then
                THROW_HARD('rec_backend=pcg requires mskdiam for normalization and solve support')
            endif
            if( params%maxits < 1 ) THROW_HARD('PCG maxits must be at least 1')
            if( .not. ieee_is_finite(params%rtol) ) THROW_HARD('PCG rtol must be finite')
        end subroutine validate_supported_mode

        integer function count_state( state_here ) result(n)
            integer, intent(in) :: state_here
            integer :: i
            n = 0
            do i = 1, nselected
                if( build%spproj_field%get_state(selected_pinds(i)) == state_here ) n = n + 1
            enddo
        end function count_state

        integer function count_state_half( state_here, eo_here ) result(n)
            integer, intent(in) :: state_here, eo_here
            integer :: i, p
            n = 0
            do i = 1, nselected
                p = selected_pinds(i)
                if( build%spproj_field%get_state(p) /= state_here ) cycle
                if( build%spproj_field%get_eo(p) == eo_here ) n = n + 1
            enddo
        end function count_state_half

        subroutine collect_state_half( state_here, eo_here, n, pinds )
            integer,              intent(in)  :: state_here, eo_here, n
            integer, allocatable, intent(out) :: pinds(:)
            integer :: i, p, cnt
            allocate(pinds(n))
            cnt = 0
            do i = 1, nselected
                p = selected_pinds(i)
                if( build%spproj_field%get_state(p) /= state_here ) cycle
                if( build%spproj_field%get_eo(p) /= eo_here ) cycle
                cnt = cnt + 1
                pinds(cnt) = p
            enddo
            if( cnt /= n ) THROW_HARD('inconsistent PCG state-half particle count')
        end subroutine collect_state_half

        subroutine solve_state_half( state_here, half, pinds, volume, outcome )
            integer,          intent(in)    :: state_here
            character(len=*), intent(in)    :: half
            integer,          intent(in)    :: pinds(:)
            type(image),      intent(inout) :: volume
            type(pcg_solver_outcome), optional, intent(out) :: outcome
            type(reconstructor_pcg) :: pcgop
            type(pcg_solver_outcome) :: result
            type(oris)      :: selection
            type(ori)       :: orientation
            type(ctfparams) :: ctfparms
            complex, allocatable :: y_batch(:,:,:)
            real,    allocatable :: sig2(:,:), x(:,:,:), rel_res_hist(:)
            integer :: lims2(2,2), R, kfromto(2), batchlims(2), batchsz
            integer :: i, ii, iptcl, ibatch, niters
            real    :: shift(2), crop_factor, sdev_noise, edge_mean
            integer(timer_int_kind) :: t_half, t_phase
            real(dp) :: time_metadata, time_particles, time_accum_init, time_accum
            real(dp) :: time_finalize, time_solve, time_total

            t_half = tic()
            time_metadata  = 0.0_dp
            time_particles = 0.0_dp
            time_accum_init = 0.0_dp
            time_accum     = 0.0_dp
            time_finalize  = 0.0_dp
            time_solve     = 0.0_dp
            write(logfhandle,'(A,I0,3A,I0)') '>>> PCG RECONSTRUCT3D: STATE ', state_here, &
                &' ', trim(half), ' PARTICLES = ', size(pinds)
            t_phase = tic()
            crop_factor = real(params%box_crop) / real(params%box)
            call pcgop%new(params%box_crop, params%smpd_crop, PCG_LAMBDA)
            call pcgop%set_sym(build%pgrpsyms)
            call pcgop%set_mask(params%msk_crop)
            lims2 = pcgop%get_lims2()
            R     = lims2(1,2)
            allocate(sig2(0:R,size(pinds)), source=1.0)
            if( params%cc_objfun == OBJFUN_EUCLID )then
                kfromto = build%esig%get_kfromto()
                do i = 1, size(pinds)
                    call upsample_sigma2(kfromto(1), kfromto(2), &
                        &build%esig%sigma2_noise(kfromto(1):kfromto(2),pinds(i)), R, sig2(0:R,i))
                enddo
            endif

            call selection%new(size(pinds), .true.)
            call orientation%new(.false.)
            do i = 1, size(pinds)
                iptcl = pinds(i)
                call build%spproj_field%get_ori(iptcl, orientation)
                ctfparms      = build%spproj%get_ctfparams(params%oritype, iptcl)
                ctfparms%smpd = params%smpd_crop
                shift         = build%spproj_field%get_2Dshift(iptcl) * crop_factor
                call orientation%set_ctfvars(ctfparms)
                call orientation%set_shift(shift)
                call selection%set_ori(i, orientation)
            enddo
            call pcgop%prep_particles(selection, use_ctf=.true., sig2=sig2)
            time_metadata = real(toc(t_phase),dp)
            t_phase = tic()
            allocate(y_batch(lims2(1,1):lims2(1,2), lims2(2,1):lims2(2,2), MAXIMGBATCHSZ))
            call pcgop%begin_accum
            time_accum_init = real(toc(t_phase),dp)
            sdev_noise = 0.0
            edge_mean  = 0.0
            do ibatch = 1, size(pinds), MAXIMGBATCHSZ
                batchlims = [ibatch, min(size(pinds),ibatch+MAXIMGBATCHSZ-1)]
                batchsz   = batchlims(2) - batchlims(1) + 1
                t_phase = tic()
                if( params%l_ptcl_src_den )then
                    call discrete_read_imgbatch_source(params, build, 'den', size(pinds), pinds, &
                        &batchlims, build%imgbatch(:batchsz))
                else
                    call discrete_read_imgbatch(params, build, size(pinds), pinds, batchlims)
                endif
                do ii = 1, batchsz
                    call build%imgbatch(ii)%norm_noise(build%lmsk, sdev_noise)
                    call build%imgbatch(ii)%taper_edges_particle(nint(COSMSKHALFWIDTH), edge_mean)
                    call build%imgbatch(ii)%fft
                    y_batch(:,:,ii) = pcgop%extract_native_plane(build%imgbatch(ii))
                enddo
                time_particles = time_particles + real(toc(t_phase),dp)
                t_phase = tic()
                call pcgop%accumulate_batch(y_batch, batchsz, batchlims(1))
                time_accum = time_accum + real(toc(t_phase),dp)
            enddo
            t_phase = tic()
            call pcgop%end_accum(.true.)
            call pcgop%set_op_mode(PCG_OP_KERNEL)
            time_finalize = real(toc(t_phase),dp)
            allocate(x(params%box_crop,params%box_crop,params%box_crop), source=0.0)
            t_phase = tic()
            call pcgop%solve_accum(x, maxits=params%maxits, rtol=params%rtol, &
                &rel_res_hist=rel_res_hist, niters=niters, outcome=result)
            time_solve = real(toc(t_phase),dp)
            call volume%new([params%box_crop,params%box_crop,params%box_crop], params%smpd_crop)
            call volume%set_rmat(x, .false.)
            time_total = real(toc(t_half),dp)
            call report_half_timings(state_here, half, time_metadata, time_particles, time_accum_init, &
                &time_accum, time_finalize, time_solve, time_total)
            call write_half_diagnostics(state_here, half, size(pinds), result, rel_res_hist, &
                &time_metadata, time_particles, time_accum_init, time_accum, time_finalize, time_solve, time_total)
            write(logfhandle,'(A,I0,3A,I0,2A)') '>>> PCG RECONSTRUCT3D: STATE ', state_here, ' ', trim(half), &
                &' FINISHED AFTER ', niters, ' ITERATIONS, STOP=', trim(result%stop_reason)
            if( present(outcome) ) outcome = result

            call pcgop%kill
            call selection%kill
            call orientation%kill
            deallocate(y_batch, sig2, x, rel_res_hist)
        end subroutine solve_state_half

        subroutine calculate_state_fsc( state_here, even, odd, avg, fsc_out, res05, res0143, cfar )
            integer,     intent(in)    :: state_here
            type(image), intent(in)    :: even, odd, avg
            real, allocatable, intent(out) :: fsc_out(:)
            real,                    intent(out) :: res05, res0143, cfar
            type(image)     :: work_even, work_odd
            type(image_msk) :: envmask
            type(fsc_area_score_result) :: cones
            real, allocatable :: fsc_t(:), fsc_n(:), res(:)
            integer :: nyq
            nyq = even%get_filtsz()
            if( params%l_envfsc )then
                call envmask%automask3D(params, avg, .false., lp_override=params%envmsklp)
                call envmask%write(string(AUTOMASK_FBODY//int2str_pad(state_here,2)//MRC_EXT))
                call phase_rand_fsc(even, odd, envmask, state_here, nyq, fsc_out, fsc_t, fsc_n)
                call work_even%copy(even)
                call work_odd%copy(odd)
                call work_even%zero_env_background(envmask)
                call work_odd%zero_env_background(envmask)
                call work_even%mul(envmask)
                call work_odd%mul(envmask)
                deallocate(fsc_t, fsc_n)
                call envmask%kill_bimg
            else
                call work_even%copy(even)
                call work_odd%copy(odd)
                call work_even%mask3D_soft(params%msk_crop, backgr=0.0)
                call work_odd%mask3D_soft(params%msk_crop, backgr=0.0)
                allocate(fsc_out(nyq), source=0.0)
            endif
            call cones%new(work_even, 256, 20.0, 0.143, 1)
            call cones%calc_fsc_area_score(work_even, work_odd, state=state_here)
            cfar = cones%cfar
            if( .not. params%l_envfsc ) call work_even%fsc(work_odd, fsc_out)
            res = avg%get_res()
            call get_resolution(fsc_out, res, res05, res0143)
            res05   = max(res05,   2.0 * params%smpd_crop)
            res0143 = max(res0143, 2.0 * params%smpd_crop)
            write(logfhandle,'(A,I0,A,F8.3)') '>>> PCG RECONSTRUCT3D: STATE ', state_here, &
                &' FSC=0.500 RESOLUTION = ', res05
            write(logfhandle,'(A,I0,A,F8.3)') '>>> PCG RECONSTRUCT3D: STATE ', state_here, &
                &' FSC=0.143 RESOLUTION = ', res0143
            call cones%kill
            call work_even%kill
            call work_odd%kill
            deallocate(res)
        end subroutine calculate_state_fsc

        subroutine write_fsc_summary( state_here, avg, fsc_in, res05, res0143, cfar )
            integer,     intent(in) :: state_here
            type(image), intent(in) :: avg
            real,        intent(in) :: fsc_in(:), res05, res0143, cfar
            type(string) :: fname
            real, allocatable :: res(:)
            integer :: funit, k
            fname = refine3D_resolution_txt_fbody(state_here)
            res   = avg%get_res()
            call fopen(funit, file=fname, status='replace', action='write')
            do k = 1, min(size(res),size(fsc_in))
                write(funit,'(A,1X,F6.2,1X,A,1X,F7.3)') &
                    &'>>> RESOLUTION:', res(k), '>>> CORRELATION:', fsc_in(k)
            enddo
            write(funit,'(A,1X,F6.2)') '>>> RESOLUTION AT FSC=0.500 DETERMINED TO:', res05
            write(funit,'(A,1X,F6.2)') '>>> RESOLUTION AT FSC=0.143 DETERMINED TO:', res0143
            write(funit,'(A,1X,F6.2)') '>>> CONICAL FSC AREA RATIO (cFAR) SCORE  :', cfar
            call fclose(funit)
            call fname%kill
            deallocate(res)
        end subroutine write_fsc_summary

        subroutine report_half_timings( state_here, half, metadata, particles, accum_init, accum, &
                &finalize, solve_time, total )
            integer,          intent(in) :: state_here
            character(len=*), intent(in) :: half
            real(dp),         intent(in) :: metadata, particles, accum_init, accum, finalize, solve_time, total
            real(dp) :: other
            other = max(0.0_dp, total - metadata - particles - accum_init - accum - finalize - solve_time)
            write(logfhandle,'(A,I0,2A)') '>>> PCG SHARED PHASES: STATE ', state_here, ' ', trim(half)
            write(logfhandle,'(A,F9.3)') '    metadata/operator preparation : ', metadata
            write(logfhandle,'(A,F9.3)') '    particle I/O + preprocessing  : ', particles
            write(logfhandle,'(A,F9.3)') '    accumulator initialization    : ', accum_init
            write(logfhandle,'(A,F9.3)') '    fused (B,D) accumulation      : ', accum
            write(logfhandle,'(A,F9.3)') '    accumulator finalization      : ', finalize
            write(logfhandle,'(A,F9.3)') '    PCG solve                     : ', solve_time
            write(logfhandle,'(A,F9.3)') '    other                         : ', other
            write(logfhandle,'(A,F9.3)') '    total half                    : ', total
            call flush(logfhandle)
        end subroutine report_half_timings

        subroutine write_half_diagnostics( state_here, half, nptcls, result, history, &
                &metadata, particles, accum_init, accum, finalize, solve_time, total )
            integer,                  intent(in) :: state_here, nptcls
            character(len=*),         intent(in) :: half
            type(pcg_solver_outcome), intent(in) :: result
            real,                     intent(in) :: history(:)
            real(dp),                 intent(in) :: metadata, particles, accum_init, accum
            real(dp),                 intent(in) :: finalize, solve_time, total
            type(string) :: fname
            integer :: funit, i
            fname = 'reconstruct3D_pcg_state'//int2str_pad(state_here,2)//'_'//trim(half)//'_diagnostics.txt'
            call fopen(funit, file=fname, status='replace', action='write')
            write(funit,'(A,I0)')     'nptcls=',               nptcls
            write(funit,'(A,I0)')     'requested_maxits=',     result%requested_maxits
            write(funit,'(A,I0)')     'iteration_count=',      result%iteration_count
            write(funit,'(A,A)')      'stop_reason=',          trim(result%stop_reason)
            write(funit,'(A,L1)')     'converged=',            result%converged
            write(funit,'(A,ES14.6)') 'initial_rel_resid_l2=', result%initial_rel_residual
            write(funit,'(A,ES14.6)') 'final_rel_resid_l2=',   result%final_rel_residual
            write(funit,'(A,ES14.6)') 'final_rel_update=',     result%final_rel_update
            write(funit,'(A,F12.6)')  'metadata_seconds=',     metadata
            write(funit,'(A,F12.6)')  'particle_io_prep_seconds=', particles
            write(funit,'(A,F12.6)')  'accum_init_seconds=',   accum_init
            write(funit,'(A,F12.6)')  'fused_B_D_accum_seconds=', accum
            write(funit,'(A,F12.6)')  'accum_finalize_seconds=', finalize
            write(funit,'(A,F12.6)')  'solve_seconds=',        solve_time
            write(funit,'(A,F12.6)')  'total_half_seconds=',   total
            do i = 1, size(history)
                write(funit,'(A,I0,A,ES14.6)') 'iter', i, '_rel_resid_l2=', history(i)
            enddo
            call fclose(funit)
            call fname%kill
        end subroutine write_half_diagnostics

        subroutine register_project_outputs()
            character(len=16) :: imgkind
            integer :: s
            if( trim(params%mkdir) /= 'yes' ) return
            if( trim(params%oritype) == 'cls3D' )then
                imgkind = 'vol_cavg'
            else
                imgkind = 'vol'
            endif
            do s = 1, params%nstates
                if( .not. state_written(s) ) cycle
                fname_vol = refine3D_state_vol_fname(s)
                fname_fsc = refine3D_fsc_fname(s)
                call build%spproj%add_vol2os_out(fname_vol, params%smpd_crop, s, trim(imgkind), box=params%box_crop)
                call build%spproj%add_fsc2os_out(fname_fsc, s, params%box_crop)
                call fname_vol%kill
                call fname_fsc%kill
            enddo
            call build%spproj%write_segment_inside('out', params%projfile)
        end subroutine register_project_outputs

    end subroutine execute_rec3D_pcg_shared

    !> Distributed worker: accumulate and atomically publish raw full-range B
    !! and real D for every local (state,half). Workers never call end_accum;
    !! folding and every nonlinear finalization step belong to the master.
    subroutine execute_rec3D_pcg_worker( params, build, cline, selected_pinds )
        type(parameters), intent(inout) :: params
        type(builder),    intent(inout) :: build
        class(cmdline),   intent(inout) :: cline
        integer,          intent(in)    :: selected_pinds(:)
        integer, allocatable :: half_pinds(:)
        character(len=256) :: provenance
        integer :: state, eo, n_half
        logical :: l_sigma_loaded

        call validate_pcg_common(params)
        provenance = pcg_raw_provenance(params)
        if( params%cc_objfun == OBJFUN_EUCLID )then
            call load_sigma2_groups(params, build%pftc, build%esig, build%spproj_field, cline, l_sigma_loaded)
            if( .not. l_sigma_loaded ) THROW_HARD('PCG objfun=euclid requires sigma2 files')
        endif
        if( size(selected_pinds) > 0 ) call prepimgbatch(params, build, MAXIMGBATCHSZ)
        do state = 1, params%nstates
            do eo = 0, 1
                call collect_worker_state_half(state, eo, selected_pinds, half_pinds)
                n_half = size(half_pinds)
                call accumulate_worker_state_half(state, eo, half_pinds, provenance)
                deallocate(half_pinds)
                write(logfhandle,'(A,I0,A,I0,A,I0,A,I0)') '>>> PCG RAW WORKER: PART ', params%part, &
                    &' STATE ', state, ' HALF ', eo, ' PARTICLES ', n_half
            enddo
        enddo
        if( size(selected_pinds) > 0 ) call killimgbatch(build)

    contains

        subroutine collect_worker_state_half( state_here, eo_here, pinds, selected )
            integer,              intent(in)  :: state_here, eo_here, pinds(:)
            integer, allocatable, intent(out) :: selected(:)
            integer :: i, n, p
            n = 0
            do i = 1, size(pinds)
                p = pinds(i)
                if( build%spproj_field%get_state(p) /= state_here ) cycle
                if( build%spproj_field%get_eo(p) /= eo_here ) cycle
                n = n + 1
            enddo
            allocate(selected(n))
            n = 0
            do i = 1, size(pinds)
                p = pinds(i)
                if( build%spproj_field%get_state(p) /= state_here ) cycle
                if( build%spproj_field%get_eo(p) /= eo_here ) cycle
                n = n + 1
                selected(n) = p
            enddo
        end subroutine collect_worker_state_half

        subroutine accumulate_worker_state_half( state_here, eo_here, pinds, provenance_here )
            integer,          intent(in) :: state_here, eo_here, pinds(:)
            character(len=*), intent(in) :: provenance_here
            type(reconstructor_pcg) :: pcgop
            type(oris)      :: selection
            type(ori)       :: orientation
            type(ctfparams) :: ctfparms
            type(string)    :: fname
            complex, allocatable :: y_batch(:,:,:)
            real,    allocatable :: sig2(:,:)
            integer :: lims2(2,2), R, kfromto(2), batchlims(2), batchsz
            integer :: i, ii, iptcl, ibatch
            real    :: shift(2), crop_factor, sdev_noise, edge_mean

            call pcgop%new(params%box_crop, params%smpd_crop, PCG_LAMBDA)
            fname = refine3D_pcg_raw_accum_fname(state_here, params%part, params%numlen, &
                &merge('odd ', 'even', eo_here == 1))
            if( size(pinds) == 0 )then
                call pcgop%write_raw_accum(fname, state_here, eo_here, params%part, &
                    &params%nparts, 0, provenance_here)
                call pcgop%kill
                call fname%kill
                return
            endif
            call pcgop%set_sym(build%pgrpsyms)
            lims2 = pcgop%get_lims2()
            R     = lims2(1,2)
            allocate(sig2(0:R,size(pinds)), source=1.0)
            if( params%cc_objfun == OBJFUN_EUCLID )then
                kfromto = build%esig%get_kfromto()
                do i = 1, size(pinds)
                    call upsample_sigma2(kfromto(1), kfromto(2), &
                        &build%esig%sigma2_noise(kfromto(1):kfromto(2),pinds(i)), R, sig2(0:R,i))
                enddo
            endif
            call selection%new(size(pinds), .true.)
            call orientation%new(.false.)
            crop_factor = real(params%box_crop) / real(params%box)
            do i = 1, size(pinds)
                iptcl = pinds(i)
                call build%spproj_field%get_ori(iptcl, orientation)
                ctfparms      = build%spproj%get_ctfparams(params%oritype, iptcl)
                ctfparms%smpd = params%smpd_crop
                shift         = build%spproj_field%get_2Dshift(iptcl) * crop_factor
                call orientation%set_ctfvars(ctfparms)
                call orientation%set_shift(shift)
                call selection%set_ori(i, orientation)
            enddo
            call pcgop%prep_particles(selection, use_ctf=.true., sig2=sig2)
            allocate(y_batch(lims2(1,1):lims2(1,2), lims2(2,1):lims2(2,2), MAXIMGBATCHSZ))
            call pcgop%begin_accum
            sdev_noise = 0.0
            edge_mean  = 0.0
            do ibatch = 1, size(pinds), MAXIMGBATCHSZ
                batchlims = [ibatch, min(size(pinds),ibatch+MAXIMGBATCHSZ-1)]
                batchsz   = batchlims(2) - batchlims(1) + 1
                if( params%l_ptcl_src_den )then
                    call discrete_read_imgbatch_source(params, build, 'den', size(pinds), pinds, &
                        &batchlims, build%imgbatch(:batchsz))
                else
                    call discrete_read_imgbatch(params, build, size(pinds), pinds, batchlims)
                endif
                do ii = 1, batchsz
                    call build%imgbatch(ii)%norm_noise(build%lmsk, sdev_noise)
                    call build%imgbatch(ii)%taper_edges_particle(nint(COSMSKHALFWIDTH), edge_mean)
                    call build%imgbatch(ii)%fft
                    y_batch(:,:,ii) = pcgop%extract_native_plane(build%imgbatch(ii))
                enddo
                call pcgop%accumulate_batch(y_batch, batchsz, batchlims(1))
            enddo
            call pcgop%write_raw_accum(fname, state_here, eo_here, params%part, &
                &params%nparts, size(pinds), provenance_here)
            call pcgop%kill
            call selection%kill
            call orientation%kill
            call fname%kill
            deallocate(y_batch, sig2)
        end subroutine accumulate_worker_state_half

    end subroutine execute_rec3D_pcg_worker

    !> Distributed master: reduce raw worker B,D artifacts in ascending part
    !! order, then perform all folding, finalization and PCG locally. Independent
    !! state/half reductions are completed and released one at a time.
    subroutine execute_rec3D_pcg_distributed_master( params, build, cline )
        type(parameters), intent(inout) :: params
        type(builder),    intent(inout) :: build
        class(cmdline),   intent(inout) :: cline
        type(image) :: half_even, half_odd, merged
        type(string) :: fname_even, fname_odd, fname_vol, fname_fsc, raw_fname
        real, allocatable :: fsc(:), res0143s(:)
        logical, allocatable :: state_written(:)
        character(len=256) :: provenance
        integer :: state, part, eo, n_even, n_odd, iptcl, istate
        real :: res05, cfar
        integer(timer_int_kind) :: t_state_phase
        real(dp) :: time_map_output, time_fsc_output

        call validate_pcg_common(params)
        if( params%nparts < 2 ) THROW_HARD('distributed PCG requires nparts>1')
        provenance = pcg_raw_provenance(params)
        allocate(res0143s(params%nstates), source=0.0)
        allocate(state_written(params%nstates), source=.false.)
        do state = 1, params%nstates
            call reduce_solve_state_half(state, 0, 'even', half_even, n_even)
            call reduce_solve_state_half(state, 1, 'odd',  half_odd,  n_odd)
            if( n_even == 0 .and. n_odd == 0 )then
                write(logfhandle,'(A,I0,A)') '>>> PCG DISTRIBUTED: STATE ', state, &
                    &' HAS NO SELECTED PARTICLES; SKIPPING'
                cycle
            endif
            if( n_even < 1 .or. n_odd < 1 ) THROW_HARD('distributed PCG requires both halfsets')
            fname_even = refine3D_state_halfvol_fname(state, 'even')
            fname_odd  = refine3D_state_halfvol_fname(state, 'odd')
            fname_vol  = refine3D_state_vol_fname(state)
            fname_fsc  = refine3D_fsc_fname(state)
            t_state_phase = tic()
            call half_even%write(fname_even, del_if_exists=.true.)
            call half_odd%write(fname_odd, del_if_exists=.true.)
            call merged%copy(half_even)
            call merged%add(half_odd)
            call merged%mul(0.5)
            call merged%write(fname_vol, del_if_exists=.true.)
            time_map_output = real(toc(t_state_phase),dp)
            t_state_phase = tic()
            call calculate_distributed_fsc(state, half_even, half_odd, merged, &
                &fsc, res05, res0143s(state), cfar)
            call arr2file(fsc, fname_fsc)
            call write_distributed_fsc_summary(state, merged, fsc, res05, res0143s(state), cfar)
            time_fsc_output = real(toc(t_state_phase),dp)
            write(logfhandle,'(A,I0)') '>>> PCG DISTRIBUTED OUTPUT PHASES: STATE ', state
            write(logfhandle,'(A,F9.3)') '    halfmap + merged-map output : ', time_map_output
            write(logfhandle,'(A,F9.3)') '    FSC + cFAR + summary        : ', time_fsc_output
            params%vols(state)      = fname_vol
            params%vols_even(state) = fname_even
            params%vols_odd(state)  = fname_odd
            call cline%set('vol'//int2str(state), fname_vol)
            state_written(state) = .true.
            call half_even%kill
            call half_odd%kill
            call merged%kill
            call fname_even%kill
            call fname_odd%kill
            call fname_vol%kill
            call fname_fsc%kill
            if( allocated(fsc) ) deallocate(fsc)
        enddo
        if( .not. any(state_written) ) THROW_HARD('distributed PCG produced no populated states')
        if( params%nstates == 1 )then
            call build%spproj_field%set_all2single('res', res0143s(1))
        else
            do iptcl = 1, build%spproj_field%get_noris()
                istate = build%spproj_field%get_state(iptcl)
                if( istate > 0 .and. istate <= params%nstates )then
                    if( state_written(istate) ) call build%spproj_field%set(iptcl, 'res', res0143s(istate))
                endif
            enddo
        endif
        call build%spproj%write_segment_inside(params%oritype, params%projfile)
        ! Delete raw artifacts only after every state has completed. Until this
        ! point they remain a restart/debug boundary for a failed master solve.
        do state = 1, params%nstates
            do eo = 0, 1
                do part = 1, params%nparts
                    raw_fname = refine3D_pcg_raw_accum_fname(state, part, params%numlen, &
                        &merge('odd ', 'even', eo == 1))
                    call del_file(raw_fname)
                enddo
            enddo
        enddo
        call raw_fname%kill
        deallocate(res0143s, state_written)

    contains

        subroutine reduce_solve_state_half( state_here, eo_here, half, volume, nptcls )
            integer,          intent(in)    :: state_here, eo_here
            character(len=*), intent(in)    :: half
            type(image),      intent(inout) :: volume
            integer,          intent(out)   :: nptcls
            type(reconstructor_pcg) :: pcgop
            type(pcg_solver_outcome) :: result
            type(string) :: fname
            real, allocatable :: x(:,:,:), rel_res_hist(:)
            integer :: part_here, n_part, niters
            integer(timer_int_kind) :: t_phase
            real(dp) :: time_reduce, time_finalize, time_solve

            call pcgop%new(params%box_crop, params%smpd_crop, PCG_LAMBDA)
            call pcgop%set_mask(params%msk_crop)
            call pcgop%begin_reduction
            nptcls = 0
            t_phase = tic()
            do part_here = 1, params%nparts
                fname = refine3D_pcg_raw_accum_fname(state_here, part_here, params%numlen, half)
                call pcgop%add_raw_accum(fname, state_here, eo_here, part_here, params%nparts, &
                    &provenance, n_part)
                nptcls = nptcls + n_part
                call fname%kill
            enddo
            time_reduce = real(toc(t_phase),dp)
            write(logfhandle,'(A,I0,3A,I0)') '>>> PCG DISTRIBUTED: STATE ', state_here, &
                &' ', trim(half), ' PARTICLES = ', nptcls
            if( nptcls == 0 )then
                call pcgop%kill
                return
            endif
            t_phase = tic()
            call pcgop%end_accum(.true.)
            call pcgop%set_op_mode(PCG_OP_KERNEL)
            time_finalize = real(toc(t_phase),dp)
            allocate(x(params%box_crop,params%box_crop,params%box_crop), source=0.0)
            t_phase = tic()
            call pcgop%solve_accum(x, maxits=params%maxits, rtol=params%rtol, &
                &rel_res_hist=rel_res_hist, niters=niters, outcome=result)
            time_solve = real(toc(t_phase),dp)
            call volume%new([params%box_crop,params%box_crop,params%box_crop], params%smpd_crop)
            call volume%set_rmat(x, .false.)
            write(logfhandle,'(A,F9.3)') '    fixed-order raw reduction    : ', time_reduce
            write(logfhandle,'(A,F9.3)') '    master finalization          : ', time_finalize
            write(logfhandle,'(A,F9.3)') '    master PCG solve             : ', time_solve
            write(logfhandle,'(A,I0,3A,I0,2A)') '>>> PCG DISTRIBUTED: STATE ', state_here, ' ', trim(half), &
                &' FINISHED AFTER ', niters, ' ITERATIONS, STOP=', trim(result%stop_reason)
            call write_distributed_diagnostics(state_here, half, nptcls, result, rel_res_hist, &
                &time_reduce, time_finalize, time_solve)
            call pcgop%kill
            deallocate(x, rel_res_hist)
        end subroutine reduce_solve_state_half

        subroutine calculate_distributed_fsc( state_here, even, odd, avg, fsc_out, res05, res0143, cfar )
            integer,     intent(in)    :: state_here
            type(image), intent(in)    :: even, odd, avg
            real, allocatable, intent(out) :: fsc_out(:)
            real,                    intent(out) :: res05, res0143, cfar
            type(image)     :: work_even, work_odd
            type(image_msk) :: envmask
            type(fsc_area_score_result) :: cones
            real, allocatable :: fsc_t(:), fsc_n(:), res(:)
            integer :: nyq
            nyq = even%get_filtsz()
            if( params%l_envfsc )then
                call envmask%automask3D(params, avg, .false., lp_override=params%envmsklp)
                call envmask%write(string(AUTOMASK_FBODY//int2str_pad(state_here,2)//MRC_EXT))
                call phase_rand_fsc(even, odd, envmask, state_here, nyq, fsc_out, fsc_t, fsc_n)
                call work_even%copy(even)
                call work_odd%copy(odd)
                call work_even%zero_env_background(envmask)
                call work_odd%zero_env_background(envmask)
                call work_even%mul(envmask)
                call work_odd%mul(envmask)
                deallocate(fsc_t, fsc_n)
                call envmask%kill_bimg
            else
                call work_even%copy(even)
                call work_odd%copy(odd)
                call work_even%mask3D_soft(params%msk_crop, backgr=0.0)
                call work_odd%mask3D_soft(params%msk_crop, backgr=0.0)
                allocate(fsc_out(nyq), source=0.0)
            endif
            call cones%new(work_even, 256, 20.0, 0.143, 1)
            call cones%calc_fsc_area_score(work_even, work_odd, state=state_here)
            cfar = cones%cfar
            if( .not. params%l_envfsc ) call work_even%fsc(work_odd, fsc_out)
            res = avg%get_res()
            call get_resolution(fsc_out, res, res05, res0143)
            res05   = max(res05,   2.0 * params%smpd_crop)
            res0143 = max(res0143, 2.0 * params%smpd_crop)
            write(logfhandle,'(A,I0,A,F8.3)') '>>> PCG DISTRIBUTED: STATE ', state_here, &
                &' FSC=0.500 RESOLUTION = ', res05
            write(logfhandle,'(A,I0,A,F8.3)') '>>> PCG DISTRIBUTED: STATE ', state_here, &
                &' FSC=0.143 RESOLUTION = ', res0143
            call cones%kill
            call work_even%kill
            call work_odd%kill
            deallocate(res)
        end subroutine calculate_distributed_fsc

        subroutine write_distributed_fsc_summary( state_here, avg, fsc_in, res05, res0143, cfar )
            integer,     intent(in) :: state_here
            type(image), intent(in) :: avg
            real,        intent(in) :: fsc_in(:), res05, res0143, cfar
            type(string) :: fname
            real, allocatable :: res(:)
            integer :: funit, k
            fname = refine3D_resolution_txt_fbody(state_here)
            res   = avg%get_res()
            call fopen(funit, file=fname, status='replace', action='write')
            do k = 1, min(size(res),size(fsc_in))
                write(funit,'(A,1X,F6.2,1X,A,1X,F7.3)') &
                    &'>>> RESOLUTION:', res(k), '>>> CORRELATION:', fsc_in(k)
            enddo
            write(funit,'(A,1X,F6.2)') '>>> RESOLUTION AT FSC=0.500 DETERMINED TO:', res05
            write(funit,'(A,1X,F6.2)') '>>> RESOLUTION AT FSC=0.143 DETERMINED TO:', res0143
            write(funit,'(A,1X,F6.2)') '>>> CONICAL FSC AREA RATIO (cFAR) SCORE  :', cfar
            call fclose(funit)
            call fname%kill
            deallocate(res)
        end subroutine write_distributed_fsc_summary

        subroutine write_distributed_diagnostics( state_here, half, nptcls, result, history, &
                &reduce_time, finalize_time, solve_time )
            integer,                  intent(in) :: state_here, nptcls
            character(len=*),         intent(in) :: half
            type(pcg_solver_outcome), intent(in) :: result
            real,                     intent(in) :: history(:)
            real(dp),                 intent(in) :: reduce_time, finalize_time, solve_time
            type(string) :: fname
            integer :: funit, i
            fname = 'reconstruct3D_pcg_state'//int2str_pad(state_here,2)//'_'//trim(half)//'_diagnostics.txt'
            call fopen(funit, file=fname, status='replace', action='write')
            write(funit,'(A,A)')      'execution_mode=',        'distributed'
            write(funit,'(A,I0)')     'nparts=',                params%nparts
            write(funit,'(A,I0)')     'nptcls=',                nptcls
            write(funit,'(A,I0)')     'requested_maxits=',      result%requested_maxits
            write(funit,'(A,I0)')     'iteration_count=',       result%iteration_count
            write(funit,'(A,A)')      'stop_reason=',           trim(result%stop_reason)
            write(funit,'(A,L1)')     'converged=',             result%converged
            write(funit,'(A,ES14.6)') 'initial_rel_resid_l2=',  result%initial_rel_residual
            write(funit,'(A,ES14.6)') 'final_rel_resid_l2=',    result%final_rel_residual
            write(funit,'(A,ES14.6)') 'final_rel_update=',      result%final_rel_update
            write(funit,'(A,F12.6)')  'raw_reduce_seconds=',    reduce_time
            write(funit,'(A,F12.6)')  'master_finalize_seconds=', finalize_time
            write(funit,'(A,F12.6)')  'solve_seconds=',         solve_time
            do i = 1, size(history)
                write(funit,'(A,I0,A,ES14.6)') 'iter', i, '_rel_resid_l2=', history(i)
            enddo
            call fclose(funit)
            call fname%kill
        end subroutine write_distributed_diagnostics

    end subroutine execute_rec3D_pcg_distributed_master

    subroutine validate_pcg_common( params )
        type(parameters), intent(in) :: params
        if( trim(params%pcgop) /= 'kernel' ) THROW_HARD('production rec_backend=pcg requires pcgop=kernel')
        if( params%maxits < 1 .or. params%maxits > 8 ) THROW_HARD('production PCG requires 1<=maxits<=8')
        if( trim(params%projrec) /= 'no' ) THROW_HARD('rec_backend=pcg does not yet support projrec=yes')
        if( params%box_crop /= params%box ) THROW_HARD('rec_backend=pcg does not yet support box cropping')
        if( params%l_update_frac .or. params%l_trail_rec )then
            THROW_HARD('PCG fractional and trailing reconstruction are not implemented')
        endif
        if( trim(params%conical_fsc) == 'yes' ) THROW_HARD('PCG conical FSC integration is not implemented')
        if( params%msk <= 0.5 .or. params%msk_crop <= 0.5 ) THROW_HARD('rec_backend=pcg requires mskdiam')
        if( .not. ieee_is_finite(params%rtol) ) THROW_HARD('PCG rtol must be finite')
    end subroutine validate_pcg_common

    function pcg_raw_provenance( params ) result(provenance)
        type(parameters), intent(in) :: params
        character(len=256) :: provenance
        provenance = 'pcgraw-v1|pgrp='//trim(params%pgrp)//'|objfun='//trim(params%objfun)// &
            &'|ptcl_src='//trim(params%ptcl_src)//'|iter='//trim(int2str(params%which_iter))
    end function pcg_raw_provenance

end module simple_rec3D_pcg_strategy
