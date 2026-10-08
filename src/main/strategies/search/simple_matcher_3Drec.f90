!@descr: Cartesian online/offline 3D reconstruction module
module simple_matcher_3Drec
use simple_core_module_api
use simple_timer
use simple_builder,         only: builder
use simple_matcher_ptcl_io, only: discrete_read_imgbatch, prepimgbatch, prep_rec_observation, killimgbatch
use simple_memoize_ft_maps, only: memoize_ft_maps, forget_ft_maps
use simple_parameters,      only: parameters
use simple_reconstructor,   only: reconstructor
use simple_refine3D_fnames, only: refine3D_partial_rec_fbody
use simple_state_weight_set, only: state_weight_set
implicit none

public :: init_rec, prep_imgs4rec, gen_rec_plane, cleanup_rec_buffers, write_state_half_partial, calc_3Drec
private
#include "simple_local_flags.inc"

contains

    !> volumetric 3d reconstruction
    !> Writes the even/odd partial accumulators of every state; the caller
    !! assembles them (volassemble) and names the output volumes. With a state weight set, or a
    !! caller's weight table wtab(nptcls,nstates) aligned with pinds, a particle enters every state it
    !! weighs into (weight above zero), with that weight.
    subroutine calc_3Drec( params, build, nptcls, pinds, wset, wtab )
        use simple_image,        only: image
        use simple_imgarr_utils, only: alloc_imgarr, dealloc_imgarr
        class(parameters), target,         intent(inout) :: params
        class(builder),                    intent(inout) :: build
        integer,                           intent(in)    :: nptcls
        integer,                           intent(in)    :: pinds(nptcls)
        class(state_weight_set), optional, intent(inout) :: wset
        real,                    optional, intent(in)    :: wtab(:,:)
        type(fplane_type), allocatable   :: fpls(:)
        type(image),       allocatable   :: crop_imgs(:)
        type(reconstructor) :: recvol
        integer, allocatable :: grouped_pinds(:), state_eo_offsets(:)
        real,    allocatable :: grouped_w(:)
        integer :: batchlims(2), ibatch, batchsz, state, eo, group
        logical :: DEBUG = .false.
        integer(timer_int_kind) :: t, t0
        real(timer_int_kind)    :: t_init, t_read, t_prep, t_grid, t_tot
        if( nptcls < 1 ) return
        if( DEBUG ) t0 = tic()
        if( present(wset) .or. present(wtab) )then
            call group_pinds_by_weights(params, build, nptcls, pinds, grouped_pinds, grouped_w, state_eo_offsets, &
                &wset=wset, wtab=wtab)
        else
            call group_pinds_by_state_eo(params, build, nptcls, pinds, grouped_pinds, state_eo_offsets)
        endif
        ! Initialize state-independent reconstruction buffers only after
        ! registration and assignment are complete.
        if( DEBUG ) t = tic()
        call init_rec(params, build, MAXIMGBATCHSZ, fpls, &
            &cropped=params%box_crop < params%box)
        if( params%box_crop < params%box )then
            call alloc_imgarr(nthr_glob, [params%box_crop,params%box_crop,1], &
                &params%smpd_crop, crop_imgs)
        endif
        call prepimgbatch(params, build, MAXIMGBATCHSZ)
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
                    call discrete_read_imgbatch(params, build, size(grouped_pinds), grouped_pinds, batchlims)
                    if( DEBUG ) t_read = t_read + toc(t)
                    if( DEBUG ) t = tic()
                    call prep_imgs4rec(params, build, batchsz, build%imgbatch(:batchsz), &
                        &grouped_pinds(batchlims(1):batchlims(2)), fpls(:batchsz), &
                        &crop_imgs=crop_imgs)
                    if( DEBUG ) t_prep = t_prep + toc(t)
                    if( DEBUG ) t = tic()
                    if( allocated(grouped_w) )then
                        call update_state_half_rec(state, eo, build, batchsz, &
                            &grouped_pinds(batchlims(1):batchlims(2)), fpls(:batchsz), recvol, &
                            &grouped_w(batchlims(1):batchlims(2)))
                    else
                        call update_state_half_rec(state, eo, build, batchsz, &
                            &grouped_pinds(batchlims(1):batchlims(2)), fpls(:batchsz), recvol)
                    endif
                    if( DEBUG ) t_grid = t_grid + toc(t)
                enddo
                ! Preserve the paired Cartesian partial contract even when a
                ! populated state has no particles in one half.
                call write_state_half_partial(params, recvol, state, eo)
                call kill_state_half_rec(recvol)
            enddo
        enddo
        call cleanup_rec_buffers(build, fpls)
        if( allocated(crop_imgs) ) call dealloc_imgarr(crop_imgs)
        deallocate(grouped_pinds, state_eo_offsets)
        if( allocated(grouped_w) ) deallocate(grouped_w)
        if( DEBUG .and. (params%part==1) )then
            t_tot = toc(t0)
            print *,'Init          : ', t_init
            print *,'Read          : ', t_read
            print *,'Prep          : ', t_prep
            print *,'Grid          : ', t_grid
            print *,'Total rec time: ', t_tot
        endif
    end subroutine calc_3Drec

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
        deallocate(group_counts, next_pos)
    end subroutine group_pinds_by_state_eo

    !> Group the selected particles by (state, even/odd) from the weight set: a particle is a member
    !! of every state whose weight for it is above zero (only params%state when state= was given),
    !! in selection order within a group, with its weight beside it
    subroutine group_pinds_by_weights( params, build, nptcls, pinds, grouped_pinds, grouped_w, state_eo_offsets, wset, wtab )
        class(parameters),                 intent(in)    :: params
        class(builder),                    intent(in)    :: build
        integer,                           intent(in)    :: nptcls, pinds(nptcls)
        integer, allocatable,              intent(out)   :: grouped_pinds(:), state_eo_offsets(:)
        real,    allocatable,              intent(out)   :: grouped_w(:)
        class(state_weight_set), optional, intent(inout) :: wset
        real,                    optional, intent(in)    :: wtab(:,:)
        type :: state_members
            integer, allocatable :: rows(:)
            real,    allocatable :: w(:)
        end type state_members
        type(state_members), allocatable :: mem(:)
        integer, allocatable :: next_pos(:)
        integer :: i, state, eo, group, ninvalid, ntot
        if( present(wset) )then
            if( wset%get_nstates() /= params%nstates ) THROW_HARD('nstates differs from the state weight set; group_pinds_by_weights')
        else
            if( size(wtab,1) /= nptcls .or. size(wtab,2) /= params%nstates ) THROW_HARD('weight table shape; group_pinds_by_weights')
        endif
        allocate(mem(params%nstates), state_eo_offsets(2*params%nstates+1), next_pos(2*params%nstates))
        ! counts per (state, half)
        state_eo_offsets = 0
        ninvalid = 0
        do state = 1,params%nstates
            if( params%l_state_defined .and. state /= params%state )then
                allocate(mem(state)%rows(0), mem(state)%w(0))
                cycle
            endif
            if( present(wset) )then
                call wset%get_members(state, 0., pinds, mem(state)%rows, mem(state)%w)
            else
                mem(state)%rows = pack(pinds,        wtab(:,state) > 0.)
                mem(state)%w    = pack(wtab(:,state), wtab(:,state) > 0.)
            endif
            do i = 1,size(mem(state)%rows)
                eo = normalized_eo(build%spproj_field%get_eo(mem(state)%rows(i)))
                if( eo < 0 )then
                    ninvalid = ninvalid + 1
                    cycle
                endif
                group = state_eo_group(state, eo)
                state_eo_offsets(group+1) = state_eo_offsets(group+1) + 1
            enddo
        enddo
        state_eo_offsets(1) = 1
        do group = 1,2*params%nstates
            state_eo_offsets(group+1) = state_eo_offsets(group) + state_eo_offsets(group+1)
        enddo
        ntot = state_eo_offsets(2*params%nstates+1) - 1
        allocate(grouped_pinds(ntot), grouped_w(ntot))
        next_pos = state_eo_offsets(1:2*params%nstates)
        do state = 1,params%nstates
            do i = 1,size(mem(state)%rows)
                eo = normalized_eo(build%spproj_field%get_eo(mem(state)%rows(i)))
                if( eo < 0 ) cycle
                group = state_eo_group(state, eo)
                grouped_pinds(next_pos(group)) = mem(state)%rows(i)
                grouped_w(next_pos(group))     = mem(state)%w(i)
                next_pos(group) = next_pos(group) + 1
            enddo
        enddo
        if( ninvalid > 0 )then
            write(logfhandle,'(A,I0)') '>>> RECONSTRUCTION: SKIPPING WEIGHTED MEMBERS WITH INVALID EVEN/ODD LABELS: ', ninvalid
        endif
        deallocate(mem, next_pos)
    end subroutine group_pinds_by_weights

    !> Initialize one half-map worker reconstructor. Pair assembly owns explicit
    !! even and odd operands only where both halves are required.
    subroutine init_state_half_rec( params, build, recvol )
        class(parameters), target, intent(inout) :: params
        class(builder),        intent(inout) :: build
        class(reconstructor),  intent(inout) :: recvol
        call recvol%new_accumulator(params, build%spproj, expand=.true.)
    end subroutine init_state_half_rec

    subroutine kill_state_half_rec( recvol )
        class(reconstructor), intent(inout) :: recvol
        call recvol%kill
    end subroutine kill_state_half_rec

    !> Insert a state/even-odd homogeneous particle batch into one half-map. With weights the batch
    !! is the state's weighted membership (labels need not match) and every plane enters with its weight.
    subroutine update_state_half_rec( state, eo, build, nptcls, pinds, fplanes, recvol, weights )
        integer,              intent(in)    :: state, eo, nptcls, pinds(nptcls)
        class(builder),       intent(inout) :: build
        type(fplane_type),    intent(inout) :: fplanes(nptcls)
        class(reconstructor), intent(inout) :: recvol
        real, optional,       intent(in)    :: weights(nptcls)
        type(ori) :: orientation
        integer :: iptcl, i
        do i = 1,nptcls
            iptcl = pinds(i)
            call build%spproj_field%get_ori(iptcl, orientation)
            if( .not. present(weights) )then
                if( orientation%get_state() /= state )then
                    THROW_HARD('non-homogeneous state reconstruction batch; update_state_half_rec')
                endif
            endif
            if( normalized_eo(orientation%get_eo()) /= eo )then
                THROW_HARD('non-homogeneous even-odd reconstruction batch; update_state_half_rec')
            endif
            if( present(weights) )then
                call recvol%insert_plane_oversamp(build%pgrpsyms, orientation, fplanes(i), w=weights(i))
            else
                call recvol%insert_plane_oversamp(build%pgrpsyms, orientation, fplanes(i))
            endif
        enddo
        call orientation%kill
    end subroutine update_state_half_rec

    !> Write one half-map partial in the existing Cartesian volume/rho format.
    subroutine write_state_half_partial( params, recvol, state, eo )
        class(parameters),     intent(in)    :: params
        class(reconstructor),  intent(inout) :: recvol
        integer,               intent(in)    :: state, eo
        type(string) :: fbody, numerator_fname, density_fname
        integer :: numlen_part
        numlen_part = max(1, params%numlen)
        fbody = refine3D_partial_rec_fbody(state, params%part, numlen_part)
        call recvol%compress_exp
        select case(eo)
            case(0)
                numerator_fname = fbody//'_even'//MRC_EXT
                density_fname   = string('rho_')//fbody//'_even'//MRC_EXT
            case(1)
                numerator_fname = fbody//'_odd'//MRC_EXT
                density_fname   = string('rho_')//fbody//'_odd'//MRC_EXT
            case default
                THROW_HARD('unsupported even-odd half; write_state_half_partial')
        end select
        call recvol%write_raw_accum(numerator_fname, density_fname)
        call fbody%kill
        call numerator_fname%kill
        call density_fname%kill
    end subroutine write_state_half_partial

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
    !>  cropped explicitly selects the box_crop pad-heap grid for the Fourier-crop
    !!  path. It is explicit rather than inferred from params, so external init_rec
    !!  callers (flex, offload) keep full-size buffers unless they opt in.
    subroutine init_rec( params, build, maxbatchsz, fplanes, cropped )
        use simple_imgarr_utils, only: alloc_imgarr
        class(parameters),              intent(in)    :: params
        class(builder),                 intent(inout) :: build
        integer,                        intent(in)    :: maxbatchsz
        type(fplane_type), allocatable, intent(inout) :: fplanes(:)
        logical, optional,              intent(in)    :: cropped
        logical :: l_cropped
        l_cropped = .false.
        if( present(cropped) ) l_cropped = cropped
        ! Sigma weighting is part of the Euclidean data objective, independent
        ! of whether volassemble subsequently adds the ML prior.
        if( params%cc_objfun == OBJFUN_EUCLID )then
            if( .not. allocated(build%esig%sigma2_noise) )then
                THROW_HARD('sigma2_noise is not allocated for Euclidean reconstruction')
            endif
        endif
        ! allocate convenience CTF & memory aligned objects
        if( allocated(fplanes) )  deallocate(fplanes)
        allocate(fplanes(maxbatchsz))
        ! Heap of padded images. boxpd and box_croppd cover the same physical
        ! extent (box*smpd == box_crop*smpd_crop), so a given Fourier index means
        ! the same spatial frequency on either grid; the cropped one simply stops
        ! at the crop Nyquist, which is all the box_crop reconstructor can hold.
        if( l_cropped )then
            call alloc_imgarr(nthr_glob, [params%box_croppd, params%box_croppd, 1], params%smpd_crop, build%img_pad_heap)
        else
            call alloc_imgarr(nthr_glob, [params%boxpd, params%boxpd, 1], params%smpd, build%img_pad_heap)
        endif
    end subroutine init_rec

    !> Preprocess particle images for online 3D reconstruction. With crop_imgs, as the cropped 2D restoration
    !! (cavger_update_sums): normalize at the native box, Fourier crop, taper and pad at box_crop without renorm,
    !! shift and CTF in cropped pixels, and the sigma2 range capped at the crop Nyquist the cropped grid spans.
    !! For model comparisons (FLEX's projected model): kfromto, the band (default the sigma2 range of build%esig;
    !! the planes' Nyquist is then capped at it); observation_model, planes holding the whitened observation and
    !! the forward transfer; whiten, whitening by the sigma2 spectra (default with the Euclidean objective).
    subroutine prep_imgs4rec( params, build, nptcls, ptcl_imgs, pinds, fplanes, crop_imgs, kfromto, &
            &observation_model, whiten )
        use simple_image, only: image
        class(parameters), intent(in)    :: params
        class(builder),    intent(inout) :: build
        integer,           intent(in)    :: nptcls
        class(image),      intent(inout) :: ptcl_imgs(nptcls)
        integer,           intent(in)    :: pinds(nptcls)
        type(fplane_type), intent(inout) :: fplanes(nptcls)
        type(image), allocatable, optional, intent(inout) :: crop_imgs(:)
        integer,                  optional, intent(in)    :: kfromto(2)
        logical,                  optional, intent(in)    :: observation_model, whiten
        integer :: i, ithr, kfromto_here(2)
        logical :: l_crop, l_band, l_obs, l_whiten
        l_crop = present(crop_imgs)
        if( l_crop ) l_crop = allocated(crop_imgs)
        l_band = present(kfromto)
        l_obs  = .false.
        if( present(observation_model) ) l_obs = observation_model
        l_whiten = params%cc_objfun == OBJFUN_EUCLID
        if( present(whiten) ) l_whiten = whiten
        if( l_whiten .and. .not. allocated(build%esig%sigma2_noise) )then
            THROW_HARD('whitened reconstruction planes need the sigma2 spectra; prep_imgs4rec')
        endif
        ! logical/physical address mapping for padded Fourier planes
        if( l_crop )then
            call memoize_ft_maps([params%box_croppd, params%box_croppd, 1], params%smpd_crop)
        else
            call memoize_ft_maps([params%boxpd, params%boxpd, 1], params%smpd)
        endif
        ! gridding batch loop
        if( l_band )then
            kfromto_here = kfromto
        else
            kfromto_here = build%esig%get_kfromto()
        endif
        if( l_crop ) kfromto_here(2) = min(kfromto_here(2), params%box_crop/2)
        !$omp parallel do default(shared) private(i,ithr) schedule(static) proc_bind(close)
        do i = 1,nptcls
            ithr = omp_get_thread_num() + 1
            if( l_crop )then
                ! backend-neutral crop step (PCG prepares through the same call); the
                ! fused routine below tapers at the cropped box, pads and transforms
                call prep_rec_observation(ptcl_imgs(i), build%lmsk, crop_imgs(ithr), .false.)
                call crop_imgs(ithr)%norm_noise_taper_edge_pad_fft(build%lmsk_crop, &
                    &build%img_pad_heap(ithr), renorm=.false.)
            else
                call ptcl_imgs(i)%norm_noise_taper_edge_pad_fft(build%lmsk, build%img_pad_heap(ithr))
            endif
            call gen_rec_plane(params, build, pinds(i), build%img_pad_heap(ithr), l_crop, kfromto_here, &
                &l_whiten, l_obs, l_band, fplanes(i))
        end do
        !$omp end parallel do
    end subroutine prep_imgs4rec

    !> The reconstruction plane of particle iptcl from its padded transform pad_img, prep_imgs4rec's
    !! per-particle step: CTF and shift scaled to the cropped grid when cropped; whitened by the particle's
    !! sigma2 spectrum when whiten; with observation_model the plane holds the whitened observation and the
    !! forward transfer; with cap its Nyquist is capped at the band kfromto(2).
    subroutine gen_rec_plane( params, build, iptcl, pad_img, cropped, kfromto, whiten, observation_model, cap, fplane )
        use simple_image, only: image
        class(parameters), intent(in)    :: params
        class(builder),    intent(inout) :: build
        integer,           intent(in)    :: iptcl, kfromto(2)
        class(image),      intent(inout) :: pad_img
        logical,           intent(in)    :: cropped, whiten, observation_model, cap
        type(fplane_type), intent(inout) :: fplane
        type(ctfparams) :: ctfparms
        real :: shift(2), crop_factor
        ctfparms = build%spproj%get_ctfparams(params%oritype, iptcl)
        shift    = build%spproj_field%get_2Dshift(iptcl)
        if( cropped )then
            ! shconst is in pixel units of the padded box the image actually has,
            ! and the CTF kernel reads cycles/pixel of the current grid
            crop_factor   = real(params%box_crop) / real(params%box)
            ctfparms%smpd = ctfparms%smpd / crop_factor ! = smpd_crop
            shift         = shift * crop_factor
        endif
        if( whiten )then
            if( iptcl < lbound(build%esig%sigma2_noise,2) .or. iptcl > ubound(build%esig%sigma2_noise,2) )then
                THROW_HARD('particle index outside the sigma2 table; gen_rec_plane')
            endif
            call pad_img%gen_fplane4rec(kfromto, params%smpd_crop, ctfparms, shift, fplane, &
                &build%esig%sigma2_noise(kfromto(1):kfromto(2),iptcl), store_transfer=observation_model, &
                &observation_model=observation_model)
        else
            call pad_img%gen_fplane4rec(kfromto, params%smpd_crop, ctfparms, shift, fplane, &
                &store_transfer=observation_model, observation_model=observation_model)
        endif
        ! the plane claims no frequency beyond the band (padded units)
        if( cap .and. fplane%nyq > 0 ) fplane%nyq = min(fplane%nyq, max(OSMPL_PAD_FAC, OSMPL_PAD_FAC*kfromto(2)))
    end subroutine gen_rec_plane

end module simple_matcher_3Drec
