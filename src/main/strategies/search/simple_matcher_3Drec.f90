!@descr: Cartesian online/offline 3D reconstruction module
module simple_matcher_3Drec
use simple_core_module_api
use simple_timer
use simple_builder,         only: builder
use simple_classaverager,  only: fourier_2d_accumulator
use simple_cmdline,         only: cmdline
use simple_matcher_ptcl_io, only: discrete_read_imgbatch, discrete_read_imgbatch_source, prepimgbatch, killimgbatch
use simple_memoize_ft_maps, only: memoize_ft_maps, forget_ft_maps
use simple_parameters,      only: parameters
use simple_reconstructor,   only: reconstructor
use simple_refine3D_fnames, only: refine3D_partial_rec_fbody, refine3D_state_vol_fname
implicit none

public :: init_rec, prep_imgs4rec, cleanup_rec_buffers, write_state_partial, calc_3Drec, calc_projdir3Drec
private
#include "simple_local_flags.inc"

contains

    !> volumetric 3d reconstruction
    subroutine calc_3Drec( params, build, cline, nptcls, pinds )
        class(parameters), intent(inout) :: params
        class(builder),    intent(inout) :: build
        class(cmdline),    intent(inout) :: cline
        integer,           intent(in)    :: nptcls
        integer,           intent(in)    :: pinds(nptcls)
        type(fplane_type), allocatable   :: fpls(:)
        type(reconstructor) :: recvol
        integer, allocatable :: grouped_pinds(:), state_eo_offsets(:)
        integer :: batchlims(2), ibatch, batchsz, state, eo, group
        logical :: l_den_src
        logical :: DEBUG = .false.
        integer(timer_int_kind) :: t, t0
        real(timer_int_kind)    :: t_init, t_read, t_prep, t_grid, t_tot
        if( nptcls < 1 ) return
        if( DEBUG ) t0 = tic()
        call group_pinds_by_state_eo(params, build, nptcls, pinds, grouped_pinds, state_eo_offsets)
        ! Initialize state-independent reconstruction buffers only after
        ! registration and assignment are complete.
        if( DEBUG ) t = tic()
        call init_rec(params, build, MAXIMGBATCHSZ, fpls)
        call prepimgbatch(params, build, MAXIMGBATCHSZ)
        l_den_src = params%l_ptcl_src_den
        if( DEBUG ) t_init = toc(t)
        if( DEBUG ) then
            t_read = 0.d0
            t_prep = 0.d0
            t_grid = 0.d0
        endif
        do state = 1,params%nstates
            if( state_eo_offsets(2*state+1) <= state_eo_offsets(2*state-1) )then
                call mark_empty_state(build, state)
                cycle
            endif
            do eo = 0,1
                group = state_eo_group(state, eo)
                call init_state_half_rec(params, build, recvol)
                do ibatch = state_eo_offsets(group),state_eo_offsets(group+1)-1,MAXIMGBATCHSZ
                    batchlims = [ibatch, min(state_eo_offsets(group+1)-1, ibatch+MAXIMGBATCHSZ-1)]
                    batchsz   = batchlims(2) - batchlims(1) + 1
                    if( DEBUG ) t = tic()
                    if( l_den_src )then
                        call discrete_read_imgbatch_source(params, build, 'den', batchsz, &
                            grouped_pinds(batchlims(1):batchlims(2)), [1,batchsz], build%imgbatch(:batchsz))
                    else
                        call discrete_read_imgbatch(params, build, size(grouped_pinds), grouped_pinds, batchlims)
                    endif
                    if( DEBUG ) t_read = t_read + toc(t)
                    if( DEBUG ) t = tic()
                    call prep_imgs4rec(params, build, batchsz, build%imgbatch(:batchsz), &
                        grouped_pinds(batchlims(1):batchlims(2)), fpls(:batchsz))
                    if( DEBUG ) t_prep = t_prep + toc(t)
                    if( DEBUG ) t = tic()
                    call update_state_half_rec(state, eo, build, batchsz, &
                        &grouped_pinds(batchlims(1):batchlims(2)), fpls(:batchsz), recvol)
                    if( DEBUG ) t_grid = t_grid + toc(t)
                enddo
                ! Preserve the paired Cartesian partial contract even when a
                ! populated state has no particles in one half.
                call write_state_half_partial(params, recvol, state, eo)
                call kill_state_half_rec(recvol)
            enddo
            call set_state_vol_output(params, cline, state)
        enddo
        call cleanup_rec_buffers(build, fpls)
        deallocate(grouped_pinds, state_eo_offsets)
        if( DEBUG .and. (params%part==1) )then
            t_tot = toc(t0)
            print *,'Init          : ', t_init
            print *,'Read          : ', t_read
            print *,'Prep          : ', t_prep
            print *,'Grid          : ', t_grid
            print *,'Total rec time: ', t_tot
        endif
    end subroutine calc_3Drec

    !> Volumetric 3D reconstruction from compact projection-direction sums.
    !! Particles are first accumulated on the native 2D Fourier grid with the
    !! exact KB numerator/CTF^2 machinery used for class averages.  Those raw
    !! sums are then inserted directly into the 3D reconstruction; they are
    !! never CTF-density corrected and never transformed through real space.
    subroutine calc_projdir3Drec( params, build, cline, nptcls, pinds )
        class(parameters), intent(inout) :: params
        class(builder),    intent(inout) :: build
        class(cmdline),    intent(inout) :: cline
        integer,           intent(in)    :: nptcls
        integer,           intent(in)    :: pinds(nptcls)
        type(fourier_2d_accumulator) :: projdir_sums
        type(fplane_type), allocatable :: fpls(:)
        type(fplane_type) :: compact_fpl
        type(reconstructor) :: recvol
        type(ori) :: orientation
        integer, allocatable :: eopops(:), grouped_pinds(:), state_eo_offsets(:), proj2slice(:)
        integer :: batchlims(2), batchsz, ibatch, i, iptcl, iproj, eo, state, nproj, group
        logical :: l_den_src
        if( nptcls < 1 ) return
        if( params%nspace /= build%eulspace%get_noris() )then
            THROW_HARD('nspace/eulspace mismatch; calc_projdir3Drec')
        endif
        ! Standalone reconstruct3D callers do not have the matcher phase that
        ! finalizes these labels before writing orientation metadata.
        call build%spproj_field%set_projs(build%eulspace)
        call group_pinds_by_state_eo(params, build, nptcls, pinds, grouped_pinds, state_eo_offsets)
        call init_rec(params, build, MAXIMGBATCHSZ, fpls)
        call prepimgbatch(params, build, MAXIMGBATCHSZ)
        allocate(eopops(params%nspace), proj2slice(params%nspace), source=0)
        l_den_src = params%l_ptcl_src_den
        do state = 1,params%nstates
            if( state_eo_offsets(2*state+1) <= state_eo_offsets(2*state-1) )then
                call mark_empty_state(build, state)
                cycle
            endif
            do eo = 0,1
                group = state_eo_group(state, eo)
                eopops = 0
                !$omp parallel do default(shared) private(i,iptcl,iproj) &
                !$omp schedule(static) proc_bind(close) reduction(+:eopops)
                do i = state_eo_offsets(group),state_eo_offsets(group+1)-1
                    iptcl = grouped_pinds(i)
                    iproj = build%spproj_field%get_int(iptcl, 'proj')
                    if( iproj < 1 .or. iproj > params%nspace ) cycle
                    eopops(iproj) = eopops(iproj) + 1
                enddo
                !$omp end parallel do
                call init_state_half_rec(params, build, recvol)
                ! Allocate only populated projection directions for this
                ! state/half, avoiding a dense nspace*box^2 allocation.
                proj2slice = 0
                nproj = 0
                do iproj = 1,params%nspace
                    if( eopops(iproj) == 0 ) cycle
                    nproj = nproj + 1
                    proj2slice(iproj) = nproj
                enddo
                call projdir_sums%new([params%box_crop,params%box_crop], max(1,nproj))
                do ibatch = state_eo_offsets(group),state_eo_offsets(group+1)-1,MAXIMGBATCHSZ
                    batchlims = [ibatch, min(state_eo_offsets(group+1)-1,ibatch+MAXIMGBATCHSZ-1)]
                    batchsz   = batchlims(2) - batchlims(1) + 1
                    if( l_den_src )then
                        call discrete_read_imgbatch_source(params, build, 'den', batchsz, &
                            &grouped_pinds(batchlims(1):batchlims(2)), [1,batchsz], build%imgbatch(:batchsz))
                    else
                        call discrete_read_imgbatch(params, build, size(grouped_pinds), grouped_pinds, batchlims)
                    endif
                    call prep_imgs4rec(params, build, batchsz, build%imgbatch(:batchsz), &
                        &grouped_pinds(batchlims(1):batchlims(2)), fpls(:batchsz))
                    ! Each OpenMP iteration owns one projection-direction
                    ! accumulator, matching the race-free class-average policy.
                    !$omp parallel do default(shared) private(iproj,i,iptcl) &
                    !$omp schedule(dynamic) proc_bind(close)
                    do iproj = 1,params%nspace
                        if( eopops(iproj) == 0 ) cycle
                        do i = batchlims(1),batchlims(2)
                            iptcl = grouped_pinds(i)
                            if( build%spproj_field%get_int(iptcl,'proj') /= iproj ) cycle
                            call projdir_sums%add_fplane(build%spproj_field%e3get(iptcl), &
                                &fpls(i-batchlims(1)+1), proj2slice(iproj))
                        enddo
                    enddo
                    !$omp end parallel do
                enddo
                do iproj = 1,params%nspace
                    if( eopops(iproj) == 0 ) cycle
                    call build%eulspace%get_ori(iproj, orientation)
                    call orientation%set_state(state)
                    call orientation%set('eo',0)
                    if( eo == 1 ) call orientation%set('eo',1)
                    call projdir_sums%export_fplane(proj2slice(iproj), compact_fpl)
                    call recvol%insert_plane_oversamp(build%pgrpsyms, orientation, compact_fpl, compact_source=.true.)
                enddo
                call projdir_sums%kill
                call write_state_half_partial(params, recvol, state, eo)
                call kill_state_half_rec(recvol)
            enddo
            call set_state_vol_output(params, cline, state)
        enddo
        if( allocated(compact_fpl%cmplx_plane) ) deallocate(compact_fpl%cmplx_plane)
        if( allocated(compact_fpl%ctfsq_plane) ) deallocate(compact_fpl%ctfsq_plane)
        deallocate(eopops, grouped_pinds, state_eo_offsets, proj2slice)
        call orientation%kill
        call cleanup_rec_buffers(build, fpls)
    end subroutine calc_projdir3Drec

    pure integer function state_eo_group( state, eo )
        integer, intent(in) :: state, eo
        state_eo_group = 2 * (state - 1) + eo + 1
    end function state_eo_group

    pure integer function normalized_eo( eo )
        integer, intent(in) :: eo
        select case(eo)
            case(-1,0)
                normalized_eo = 0
            case(1)
                normalized_eo = 1
            case default
                normalized_eo = -1
        end select
    end function normalized_eo

    !> Group the selected reconstruction particles by final hard state and
    !! even/odd half.  This lets a worker construct only one half-volume at a
    !! time while preserving the existing paired partial-file contract.
    subroutine group_pinds_by_state_eo( params, build, nptcls, pinds, grouped_pinds, state_eo_offsets )
        class(parameters),              intent(in)  :: params
        class(builder),                 intent(in)  :: build
        integer,                        intent(in)  :: nptcls, pinds(nptcls)
        integer, allocatable,           intent(out) :: grouped_pinds(:), state_eo_offsets(:)
        integer, allocatable :: group_counts(:), next_pos(:)
        integer :: i, iptcl, state, eo, group, nvalid, ninvalid
        allocate(group_counts(2*params%nstates), source=0)
        ninvalid = 0
        do i = 1,nptcls
            iptcl = pinds(i)
            state = build%spproj_field%get_state(iptcl)
            if( state < 1 .or. state > params%nstates )then
                ninvalid = ninvalid + 1
                cycle
            endif
            eo = normalized_eo(build%spproj_field%get_eo(iptcl))
            if( eo < 0 )then
                ninvalid = ninvalid + 1
                cycle
            endif
            group = state_eo_group(state, eo)
            group_counts(group) = group_counts(group) + 1
        enddo
        allocate(state_eo_offsets(2*params%nstates+1))
        state_eo_offsets(1) = 1
        do group = 1,2*params%nstates
            state_eo_offsets(group+1) = state_eo_offsets(group) + group_counts(group)
        enddo
        nvalid = state_eo_offsets(2*params%nstates+1) - 1
        allocate(grouped_pinds(nvalid), next_pos(2*params%nstates))
        next_pos = state_eo_offsets(1:2*params%nstates)
        do i = 1,nptcls
            iptcl = pinds(i)
            state = build%spproj_field%get_state(iptcl)
            if( state < 1 .or. state > params%nstates ) cycle
            eo = normalized_eo(build%spproj_field%get_eo(iptcl))
            if( eo < 0 ) cycle
            group = state_eo_group(state, eo)
            grouped_pinds(next_pos(group)) = iptcl
            next_pos(group) = next_pos(group) + 1
        enddo
        if( nvalid + ninvalid /= nptcls ) THROW_HARD('invalid state/even-odd grouping count; group_pinds_by_state_eo')
        if( ninvalid > 0 )then
            write(logfhandle,'(A,I0)') '>>> RECONSTRUCTION: SKIPPING PARTICLES WITH INVALID STATE OR EVEN/ODD LABELS: ', ninvalid
        endif
        do state = 1,params%nstates
            write(logfhandle,'(A,I0,A,I0,A,I0,A,I0)') '>>> RECONSTRUCTION STATE ', state, ' PARTICLE COUNT: ', &
                &group_counts(state_eo_group(state,0)) + group_counts(state_eo_group(state,1)), &
                &' (EVEN/ODD): ', group_counts(state_eo_group(state,0)), '/', group_counts(state_eo_group(state,1))
        enddo
        deallocate(group_counts, next_pos)
    end subroutine group_pinds_by_state_eo

    !> Initialize one half-map worker reconstructor.  The EO composite remains
    !! owned by volassemble, where both halves are required for FSC and merge.
    subroutine init_state_half_rec( params, build, recvol )
        class(parameters),     intent(inout) :: params
        class(builder),        intent(inout) :: build
        class(reconstructor),  intent(inout) :: recvol
        call recvol%dealloc_rho
        call recvol%kill
        call recvol%new([params%box_crop,params%box_crop,params%box_crop], params%smpd_crop)
        call recvol%alloc_rho(params, build%spproj, expand=.true.)
        call recvol%reset_exp
    end subroutine init_state_half_rec

    subroutine kill_state_half_rec( recvol )
        class(reconstructor), intent(inout) :: recvol
        call recvol%dealloc_rho
        call recvol%kill
    end subroutine kill_state_half_rec

    !> Insert a state/even-odd homogeneous particle batch into one half-map.
    subroutine update_state_half_rec( state, eo, build, nptcls, pinds, fplanes, recvol )
        integer,              intent(in)    :: state, eo, nptcls, pinds(nptcls)
        class(builder),       intent(inout) :: build
        type(fplane_type),    intent(inout) :: fplanes(nptcls)
        class(reconstructor), intent(inout) :: recvol
        type(ori) :: orientation
        integer :: iptcl, i
        do i = 1,nptcls
            iptcl = pinds(i)
            call build%spproj_field%get_ori(iptcl, orientation)
            if( orientation%get_state() /= state )then
                THROW_HARD('non-homogeneous state reconstruction batch; update_state_half_rec')
            endif
            if( normalized_eo(orientation%get_eo()) /= eo )then
                THROW_HARD('non-homogeneous even-odd reconstruction batch; update_state_half_rec')
            endif
            call recvol%insert_plane_oversamp(build%pgrpsyms, orientation, fplanes(i))
        enddo
        call orientation%kill
    end subroutine update_state_half_rec

    !> Write one half-map partial in the existing Cartesian volume/rho format.
    subroutine write_state_half_partial( params, recvol, state, eo )
        class(parameters),     intent(in)    :: params
        class(reconstructor),  intent(inout) :: recvol
        integer,               intent(in)    :: state, eo
        type(string) :: fbody
        integer :: numlen_part
        numlen_part = max(1, params%numlen)
        fbody = refine3D_partial_rec_fbody(state, params%part, numlen_part)
        call recvol%compress_exp
        select case(eo)
            case(0)
                call recvol%write(fbody//'_even'//MRC_EXT, del_if_exists=.true.)
                call recvol%write_rho(string('rho_')//fbody//'_even'//MRC_EXT)
            case(1)
                call recvol%write(fbody//'_odd'//MRC_EXT, del_if_exists=.true.)
                call recvol%write_rho(string('rho_')//fbody//'_odd'//MRC_EXT)
            case default
                THROW_HARD('unsupported even-odd half; write_state_half_partial')
        end select
    end subroutine write_state_half_partial

    subroutine set_state_vol_output( params, cline, state )
        class(parameters), intent(inout) :: params
        class(cmdline),    intent(inout) :: cline
        integer,           intent(in)    :: state
        if( .not. cline%defined('force_volassemble') )then
            params%vols(state) = refine3D_state_vol_fname(state)
            call cline%set('vol'//int2str(state), params%vols(state))
        endif
    end subroutine set_state_vol_output

    !> Write one current-state Cartesian partial using the EO composite.  This
    !! remains the OpenMP-offload writer; CPU workers write one half at a time.
    subroutine write_state_partial( params, build, cline, state )
        class(parameters), intent(inout) :: params
        class(builder),    intent(inout) :: build
        class(cmdline),    intent(inout) :: cline
        integer,           intent(in)    :: state
        integer :: numlen_part
        numlen_part = max(1, params%numlen)
        call build%eorecvol%compress_exp
        call build%eorecvol%write_eos(refine3D_partial_rec_fbody(state, params%part, numlen_part))
        if( .not. cline%defined('force_volassemble') )then
            params%vols(state) = refine3D_state_vol_fname(state)
            call cline%set('vol'//int2str(state), params%vols(state))
        endif
    end subroutine write_state_partial

    subroutine mark_empty_state( build, state )
        class(builder),    intent(inout) :: build
        integer,           intent(in)    :: state
        if( allocated(build%fsc) ) build%fsc(state,:) = 0.
    end subroutine mark_empty_state

    !> Release reconstruction-only buffers after the final state is written.
    subroutine cleanup_rec_buffers( build, fplanes )
        use simple_imgarr_utils, only: dealloc_imgarr
        class(builder),                 intent(inout) :: build
        type(fplane_type), allocatable, intent(inout) :: fplanes(:)
        integer :: i
        if( allocated(fplanes) )then
            do i = 1,size(fplanes)
                if( allocated(fplanes(i)%cmplx_plane) )    deallocate(fplanes(i)%cmplx_plane)
                if( allocated(fplanes(i)%ctfsq_plane) )    deallocate(fplanes(i)%ctfsq_plane)
                if( allocated(fplanes(i)%transfer_plane) ) deallocate(fplanes(i)%transfer_plane)
            enddo
            deallocate(fplanes)
        endif
        call dealloc_imgarr(build%img_pad_heap)
        call forget_ft_maps
        call killimgbatch(build)
    end subroutine cleanup_rec_buffers

    !>  Initiates objects required for online volumetric 3d reconstruction
    !>  Does not read images
    subroutine init_rec( params, build, maxbatchsz, fplanes )
        use simple_imgarr_utils, only: alloc_imgarr
        class(parameters),              intent(in)    :: params
        class(builder),                 intent(inout) :: build
        integer,                        intent(in)    :: maxbatchsz
        type(fplane_type), allocatable, intent(inout) :: fplanes(:)
        ! sanity check for ml_reg
        if( params%l_ml_reg )then
            if( .not. allocated(build%esig%sigma2_noise) )then
                THROW_HARD('build%esig%sigma2_noise is not allocated while ml_reg is enabled; calc_3Drec')
            endif
        endif
        ! allocate convenience CTF & memory aligned objects
        if( allocated(fplanes) )  deallocate(fplanes)
        allocate(fplanes(maxbatchsz))
        ! heap of padded images
        call alloc_imgarr(nthr_glob, [params%boxpd, params%boxpd, 1], params%smpd, build%img_pad_heap)
    end subroutine init_rec

    !> Preprocess particle images for online volumetric 3d reconstruction
    subroutine prep_imgs4rec( params, build, nptcls, ptcl_imgs, pinds, fplanes )
        use simple_image, only: image
        class(parameters), intent(in)    :: params
        class(builder),    intent(inout) :: build
        integer,           intent(in)    :: nptcls
        class(image),      intent(inout) :: ptcl_imgs(nptcls)
        integer,           intent(in)    :: pinds(nptcls)
        type(fplane_type), intent(inout) :: fplanes(nptcls)
        type(ctfparams) :: ctfparms(nthr_glob)
        real      :: shift(2)
        integer   :: iptcl, i, ithr, kfromto(2)
        ! logical/physical address mapping for padded Fourier planes
        call memoize_ft_maps([params%boxpd, params%boxpd, 1], params%smpd)
        ! gridding batch loop
        kfromto = build%esig%get_kfromto()
        !$omp parallel do default(shared) private(i,ithr,iptcl,shift) schedule(static) proc_bind(close)
        do i = 1,nptcls
            ithr   = omp_get_thread_num() + 1
            iptcl  = pinds(i)
            call ptcl_imgs(i)%norm_noise_taper_edge_pad_fft(build%lmsk, build%img_pad_heap(ithr))
            ctfparms(ithr) = build%spproj%get_ctfparams(params%oritype, iptcl)
            shift = build%spproj_field%get_2Dshift(iptcl)
            if( params%l_ml_reg )then
                call build%img_pad_heap(ithr)%gen_fplane4rec(kfromto, params%smpd_crop, ctfparms(ithr),&
                &shift, fplanes(i), build%esig%sigma2_noise(kfromto(1):kfromto(2),iptcl))
            else
                call build%img_pad_heap(ithr)%gen_fplane4rec(kfromto, params%smpd_crop, ctfparms(ithr),&
                &shift, fplanes(i))
            endif
        end do
        !$omp end parallel do
    end subroutine prep_imgs4rec

end module simple_matcher_3Drec
