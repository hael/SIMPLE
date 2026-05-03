!@descr: Cartesian online/offline 3D reconstruction module
module simple_matcher_3Drec
use simple_core_module_api
use simple_timer
use simple_builder,         only: builder
use simple_cmdline,         only: cmdline
use simple_matcher_ptcl_io, only: discrete_read_imgbatch, prepimgbatch
use simple_memoize_ft_maps, only: memoize_ft_maps, forget_ft_maps
use simple_parameters,      only: parameters
use simple_refine3D_fnames, only: refine3D_partial_rec_fbody, refine3D_state_vol_fname
implicit none

public :: init_rec, prep_imgs4rec, update_rec, write_partial_recs, finalize_rec_objs, calc_3Drec
private
#include "simple_local_flags.inc"

! Experimental Cartesian reconstruction mode. When enabled, particle Fourier
! planes are inserted with weighted nearest-cell gridding instead of the
! standard full KB splat in reconstructor%insert_plane_oversamp.
logical, parameter :: RECON_USE_GRID_PLANE_NN = .false.

contains

    !>  \brief  initializes all volumes for reconstruction
    subroutine preprecvols( params, build )
        class(parameters), intent(in)    :: params
        class(builder),    intent(inout) :: build
        integer, allocatable :: pops(:)
        integer :: istate
        call build%spproj_field%get_pops(pops, 'state', maxn=params%nstates)
        do istate = 1, params%nstates
            if( pops(istate) > 0)then
                call build%eorecvols(istate)%new(params, build%spproj)
                call build%eorecvols(istate)%reset_all
            endif
        end do
        deallocate(pops)
    end subroutine preprecvols

    !>  \brief  destructs all volumes for reconstruction
    subroutine killrecvols( params, build )
        class(parameters), intent(in) :: params
        class(builder),    intent(inout) :: build
        integer :: istate
        do istate = 1, params%nstates
            call build%eorecvols(istate)%kill
        end do
    end subroutine killrecvols

    !>  \brief  grids one particle image to the volume
    subroutine grid_ptcl( build, fpl, se, o )
        class(builder),     intent(inout) :: build
        class(fplane_type), intent(in)    :: fpl
        class(sym),         intent(inout) :: se
        class(ori),         intent(inout) :: o
        real    :: pw
        integer :: s, eo
        ! state flag
        s = o%get_state()
        if( s == 0 ) return
        ! eo flag
        eo = o%get_eo()
        ! particle-weight
        pw = 1.0
        if( o%isthere('w') ) pw = o%get('w')
        if( pw > TINY )then
            if( RECON_USE_GRID_PLANE_NN )then
                call build%eorecvols(s)%grid_plane_nn(se, o, fpl, eo, pw)
            else
                call build%eorecvols(s)%grid_plane(se, o, fpl, eo, pw)
            endif
        endif
    end subroutine grid_ptcl

    !> volumetric 3d reconstruction
    subroutine calc_3Drec( params, build, cline, nptcls, pinds )
        use simple_imgarr_utils, only: alloc_imgarr, dealloc_imgarr
        class(parameters), intent(inout) :: params
        class(builder),    intent(inout) :: build
        class(cmdline),    intent(inout) :: cline
        integer,           intent(in)    :: nptcls
        integer,           intent(in)    :: pinds(nptcls)
        type(fplane_type), allocatable   :: fpls(:)
        integer :: batchlims(2), ibatch, batchsz
        logical :: DEBUG = .false.
        integer(timer_int_kind) :: t, t0
        real(timer_int_kind)    :: t_init, t_read, t_prep, t_grid, t_tot
        if( DEBUG ) t0 = tic()
        ! Initialize objects for recontruction
        if( DEBUG ) t = tic()
        call init_rec(params, build, MAXIMGBATCHSZ, fpls)
        ! Prep batch image objects
        call prepimgbatch(params, build, MAXIMGBATCHSZ)
        if( DEBUG ) t_init = toc(t)
        ! gridding batch loop
        if( DEBUG ) then
            t_read = 0.d0
            t_prep = 0.d0
            t_grid = 0.d0
        endif
        do ibatch = 1,nptcls,MAXIMGBATCHSZ
            batchlims = [ibatch, min(nptcls, ibatch+MAXIMGBATCHSZ-1)]
            batchsz   = batchlims(2) - batchlims(1) + 1
            ! read images
            if( DEBUG ) t = tic()
            call discrete_read_imgbatch(params, build, nptcls, pinds, batchlims)
            if( DEBUG ) t_read = t_read + toc(t)
            ! preprocess images into padded objects
            if( DEBUG ) t = tic()
            call prep_imgs4rec(params, build, batchsz, build%imgbatch(:batchsz),&
                                &pinds(batchlims(1):batchlims(2)), fpls(:batchsz))
            if( DEBUG ) t_prep = t_prep + toc(t)
            ! insert padded slices into lattice
            if( DEBUG ) t = tic()
            call update_rec(params, build, batchsz, pinds(batchlims(1):batchlims(2)), fpls(:batchsz))
            if( DEBUG ) t_grid = t_grid + toc(t)
        end do
        ! Write partial reconstructions and clean up reconstruction objects
        call write_partial_recs(params, build, cline, fpls)
        call finalize_rec_objs(params, build)
        if( DEBUG .and. (params%part==1) )then
            t_tot = toc(t0)
            print *,'Init          : ', t_init
            print *,'Read          : ', t_read
            print *,'Prep          : ', t_prep
            print *,'Grid          : ', t_grid
            print *,'Total rec time: ', t_tot
        endif
    end subroutine calc_3Drec

    !>  Initiates objects required for online volumetric 3d reconstruction
    !>  Does not read images
    subroutine init_rec( params, build, maxbatchsz, fplanes, init_volumes )
        use simple_imgarr_utils, only: alloc_imgarr
        class(parameters),              intent(in)    :: params
        class(builder),                 intent(inout) :: build
        integer,                        intent(in)    :: maxbatchsz
        type(fplane_type), allocatable, intent(inout) :: fplanes(:)
        logical, optional,              intent(in)    :: init_volumes
        logical :: l_init_volumes
        l_init_volumes = .true.
        if( present(init_volumes) ) l_init_volumes = init_volumes
        ! sanity check for ml_reg
        if( params%l_ml_reg )then
            if( .not. allocated(build%esig%sigma2_noise) )then
                THROW_HARD('build%esig%sigma2_noise is not allocated while ml_reg is enabled; calc_3Drec')
            endif
        endif
        ! init volumes
        if( l_init_volumes ) call preprecvols(params, build)
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

    !> volumetric 3d reconstruction
    subroutine update_rec( params, build, nptcls, pinds, fplanes )
        class(parameters), intent(in)    :: params
        class(builder),    intent(inout) :: build
        integer,           intent(in)    :: nptcls
        integer,           intent(in)    :: pinds(nptcls)
        type(fplane_type), intent(inout) :: fplanes(nptcls)
        type(ori) :: orientation
        integer   :: iptcl, i
        call memoize_ft_maps([params%boxpd, params%boxpd, 1], params%smpd)
        ! gridding
        do i = 1,nptcls
            iptcl  = pinds(i)
            call build%spproj_field%get_ori(iptcl, orientation)
            if( orientation%isstatezero() ) cycle
            call grid_ptcl(build, fplanes(i), build%pgrpsyms, orientation)
        end do
        call orientation%kill
    end subroutine update_rec

    !> volumetric 3d reconstruction
    subroutine write_partial_recs( params, build, cline, fplanes )
        use simple_imgarr_utils, only: alloc_imgarr, dealloc_imgarr
        class(parameters),              intent(inout) :: params
        class(builder),                 intent(inout) :: build
        class(cmdline),                 intent(inout) :: cline
        type(fplane_type), allocatable, intent(inout) :: fplanes(:)
        integer :: i, s, numlen_part
        ! deallocate convenience objects
        do i = 1,size(fplanes)
            if( allocated(fplanes(i)%cmplx_plane) ) deallocate(fplanes(i)%cmplx_plane)
            if( allocated(fplanes(i)%ctfsq_plane) ) deallocate(fplanes(i)%ctfsq_plane)
        end do
        deallocate(fplanes)
        call dealloc_imgarr(build%img_pad_heap)
        ! write partial reconstructions for downstream volassemble
        numlen_part = max(1, params%numlen)
        do s=1,params%nstates
            if( build%spproj_field%get_pop(s, 'state') == 0 )then
                build%fsc(s,:) = 0.
                cycle
            endif
            call build%eorecvols(s)%compress_exp
            call build%eorecvols(s)%write_eos(refine3D_partial_rec_fbody(s, params%part, numlen_part))
            if( .not. cline%defined('force_volassemble') )then
                params%vols(s) = refine3D_state_vol_fname(s)
                call cline%set('vol'//int2str(s), params%vols(s))
            endif
        end do
    end subroutine write_partial_recs

    subroutine finalize_rec_objs( params, build )
        use simple_imgarr_utils, only: dealloc_imgarr
        class(parameters), intent(in)    :: params
        class(builder),    intent(inout) :: build
        call dealloc_imgarr(build%img_pad_heap)
        call forget_ft_maps
        call killrecvols(params, build)
    end subroutine finalize_rec_objs

end module simple_matcher_3Drec
