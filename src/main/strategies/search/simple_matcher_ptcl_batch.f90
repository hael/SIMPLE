!@descr: particle batch preparation helpers for matcher workflows
module simple_matcher_ptcl_batch
use simple_pftc_srch_api
use simple_builder,         only: builder
use simple_matcher_ptcl_io, only: prepimgbatch, discrete_read_imgbatch, killimgbatch
!$ use omp_lib, only: omp_get_max_threads
use simple_matcher_2Dprep,  only: prepimg4align, prepimg4align_cached, prepimg4align_cart
use simple_cartft_calc,     only: cartft_calc
use simple_ptcl_cache,      only: ptcl_cache_in_use, ptcl_cache_read_batch
implicit none

public :: prep_sigmas_objfun, alloc_ptcl_imgs
public :: build_batch_particles3D, build_batch_particles2D, prep_cart_batch
public :: clean_batch_particles2D, clean_batch_particles3D
private
#include "simple_local_flags.inc"

contains

    !> The sigma2 state of a pass that needs it (Euclidean scoring, or the CC residual update of
    !! the polar CC pose initialization); a Cartesian cc pass reads none (C5). A polish starts from
    !! the discrete pass's rows, not the committed generation's.
    subroutine prep_sigmas_objfun( params, build )
        class(parameters), intent(inout) :: params
        class(builder),    intent(inout) :: build
        type(string)      :: fname
        logical           :: found
        ! cc_emit_sigma is a CC-only update path. CC does not consume sigma,
        ! while Euclidean scoring requires populated grouped sigma values.
        if( trim(params%cc_emit_sigma) == 'yes' .and. params%cc_objfun == OBJFUN_EUCLID )then
            THROW_HARD('cc_emit_sigma=yes requires objfun=cc; euclid scoring would use unpopulated sigmas')
        endif
        if( params%cc_objfun == OBJFUN_EUCLID .or. trim(params%cc_emit_sigma) == 'yes' )then
            call build%spproj%get_sigma2_state_path(fname, found)
            if( .not. found ) THROW_HARD('particle project has no canonical sigma2 state path')
            if( trim(params%cc_emit_sigma) == 'yes' .and. .not. file_exists(fname) )then
                THROW_HARD('CC residual sigma update requires image-bootstrap sigma2')
            endif
            if( params%l_cart_refine )then
                call build%esig%new(params, fname, params%box)
            else
                call build%esig%new(params, build%pftc, fname, params%box)
            endif
            call build%esig%read_part(  build%spproj_field)
            ! a polish overlays its Cartesian residuals on those of the discrete pass it follows
            if( params%l_cont_polish ) call build%esig%read_pending_range
            call build%esig%read_groups(build%spproj_field)
            call fname%kill
        end if
    end subroutine prep_sigmas_objfun

    !>  imgbatch_box/imgbatch_smpd size the raw read buffer; pass params%box_crop and
    !!  params%smpd_crop when the batch will be filled from the downscaled particle
    !!  cache rather than the originals. The sampling distance has to travel with the
    !!  box: prepimg4align_cached derives img_out's smpd from the input image, so a
    !!  buffer left at params%smpd would stamp the wrong sampling onto ptcl_match_imgs.
    subroutine alloc_ptcl_imgs( params, build, ptcl_imgs, ptcl_imgs_pad, batchsz, imgbatch_box, imgbatch_smpd )
        use simple_imgarr_utils, only: alloc_imgarr, dealloc_imgarr
        class(parameters),        intent(inout) :: params
        class(builder),           intent(inout) :: build
        type(image), allocatable, intent(inout) :: ptcl_imgs(:)
        type(image), allocatable, intent(inout) :: ptcl_imgs_pad(:)
        integer,                  intent(in)    :: batchsz
        integer, optional,        intent(in)    :: imgbatch_box
        real,    optional,        intent(in)    :: imgbatch_smpd
        integer           :: ithr
        if( present(imgbatch_box) )then
            if( present(imgbatch_smpd) )then
                call prepimgbatch(params, build, batchsz, box=imgbatch_box, smpd=imgbatch_smpd)
            else
                call prepimgbatch(params, build, batchsz, box=imgbatch_box)
            endif
        else
            call prepimgbatch(params, build, batchsz)
        endif
        allocate(ptcl_imgs(nthr_glob), ptcl_imgs_pad(nthr_glob))
        !$omp parallel do default(shared) private(ithr) schedule(static) proc_bind(close)
        do ithr = 1,nthr_glob
            call ptcl_imgs(ithr)%new(    [params%box_crop,  params%box_crop,  1], params%smpd_crop, wthreads=.false.)
            call ptcl_imgs_pad(ithr)%new([params%box_croppd,params%box_croppd,1], params%smpd_crop, wthreads=.false.)
        enddo
        !$omp end parallel do
    end subroutine alloc_ptcl_imgs

    !> Prepare one 3-D particle batch from a single read for the representation(s) of the pass:
    !! the particle slots of the Cartesian calculator when the pass has built it (build%cftc
    !! holds the references a Cartesian pass or the in-matcher polish read; other 3-D batch
    !! users, the probability-table commanders, build none) and the polar particles unless the
    !! pass is Cartesian (l_cart_refine). The Cartesian slots are prepared first, from copies,
    !! because the polar preparation normalizes the images in place.
    subroutine build_batch_particles3D( params, build, nptcls_here, pinds_here, tmp_imgs, tmp_imgs_pad )
        class(parameters),      intent(in)    :: params
        class(builder),         intent(inout) :: build
        integer,                intent(in)    :: nptcls_here
        integer,                intent(in)    :: pinds_here(nptcls_here)
        class(image),           intent(inout) :: tmp_imgs(params%nthr), tmp_imgs_pad(params%nthr)
        logical :: l_polar, l_cart
        l_polar   = .not. params%l_cart_refine
        l_cart    = build%cftc%has_refs()
        if( l_polar ) call build%pftc%reallocate_ptcls(nptcls_here, pinds_here)
        call discrete_read_imgbatch(params, build, nptcls_here, pinds_here, [1,nptcls_here])
        if( l_cart ) call prep_cart_batch_of_build(params, build, nptcls_here, pinds_here)
        if( params%l_cart_refine .and. .not. l_cart ) THROW_HARD('a Cartesian pass requires its references in build%cftc')
        if( .not. l_polar ) return
        call polarize_batch_particles3D(params, build, nptcls_here, pinds_here, build%imgbatch(:nptcls_here), &
            tmp_imgs, tmp_imgs_pad)
    end subroutine build_batch_particles3D

    !> The inputs of prep_cart_batch from the project: the stored shifts on the cropped grid, the
    !! CTF parameters, the noise mask and, under euclid only, the particles' sigma2 (C5).
    subroutine prep_cart_batch_of_build( params, build, nptcls_here, pinds_here )
        class(parameters), intent(in)    :: params
        class(builder),    intent(inout) :: build
        integer,           intent(in)    :: nptcls_here
        integer,           intent(in)    :: pinds_here(nptcls_here)
        type(ctfparams), allocatable :: ctfparms(:)
        real,            allocatable :: crop_shifts(:,:), sigma2(:,:)
        real    :: crop_factor
        integer :: i, iptcl
        crop_factor = real(params%box_crop) / real(params%box)
        allocate(ctfparms(nptcls_here), crop_shifts(2,nptcls_here))
        do i = 1, nptcls_here
            iptcl            = pinds_here(i)
            crop_shifts(:,i) = build%spproj_field%get_2Dshift(iptcl) * crop_factor
            ctfparms(i)      = build%spproj%get_ctfparams(params%oritype, iptcl)
        end do
        if( params%cc_objfun == OBJFUN_EUCLID )then
            if( .not. allocated(build%esig%sigma2_noise) ) THROW_HARD('a Euclidean Cartesian batch requires sigma2')
            allocate(sigma2(0:params%kfromto(2),nptcls_here), source=1.)
            do i = 1, nptcls_here
                sigma2(params%kfromto(1):params%kfromto(2),i) = &
                    &build%esig%sigma2_noise(params%kfromto(1):params%kfromto(2),pinds_here(i))
            end do
            call prep_cart_batch(build%cftc, nptcls_here, build%imgbatch(:nptcls_here), build%lmsk, params%box_crop, &
                &params%smpd_crop, params%msk_crop, crop_shifts, ctfparms, params%kfromto, sigma2)
        else
            call prep_cart_batch(build%cftc, nptcls_here, build%imgbatch(:nptcls_here), build%lmsk, params%box_crop, &
                &params%smpd_crop, params%msk_crop, crop_shifts, ctfparms, params%kfromto)
        endif
    end subroutine prep_cart_batch_of_build

    !> Fill particle slots 1..nptcls of the Cartesian calculator, slot i from raw_imgs(i), in
    !! one parallel loop, one particle per iteration (plan section 6.2): prepimg4align_cart on a
    !! thread copy of the raw image (the input is left as read) with the stored shift on the
    !! cropped grid and the calculator's taper, then set_ptcl for euclid when sigma2 (0:kto,
    !! nptcls) is present, for cc otherwise (C5). The soft-mask coordinates and Fourier maps of
    !! the cropped box are memoized here, serially. THROW_HARD when the calculator has fewer
    !! slots than nptcls.
    subroutine prep_cart_batch( cftc, nptcls, raw_imgs, noise_mask, box_crop, smpd_crop, mskrad, crop_shifts, &
        &ctfparms, kfromto, sigma2 )
        class(cartft_calc), intent(inout) :: cftc
        integer,            intent(in)    :: nptcls, box_crop
        class(image),       intent(in)    :: raw_imgs(nptcls)
        logical,            intent(in)    :: noise_mask(:,:,:)
        real,               intent(in)    :: smpd_crop, mskrad, crop_shifts(2,nptcls)
        type(ctfparams),    intent(in)    :: ctfparms(nptcls)
        integer,            intent(in)    :: kfromto(2)
        real, optional,     intent(in)    :: sigma2(0:,:)
        type(image), allocatable :: raw_work(:), out_work(:)
        type(ctfparams)          :: ctf_crop
        complex,     allocatable :: observed(:,:)
        real,        allocatable :: taper(:)
        integer :: i, ithr, nworkers, ldim_raw(3)
        if( nptcls < 1 ) return
        if( cftc%get_nptcls() < nptcls ) THROW_HARD('Cartesian calculator has fewer particle slots than the batch')
        if( present(sigma2) )then
            if( size(sigma2,2) < nptcls ) THROW_HARD('Cartesian batch sigma2 does not cover the batch')
        endif
        nworkers = 1
        !$ nworkers = omp_get_max_threads()
        ldim_raw = raw_imgs(1)%get_ldim()
        taper    = cftc%get_ptcl_taper()
        allocate(raw_work(nworkers), out_work(nworkers))
        do ithr = 1, nworkers
            call raw_work(ithr)%new(ldim_raw, raw_imgs(1)%get_smpd(), wthreads=.false.)
            call out_work(ithr)%new([box_crop, box_crop, 1], smpd_crop, wthreads=.false.)
        end do
        call out_work(1)%memoize_mask_coords
        call memoize_ft_maps([box_crop, box_crop, 1], smpd_crop)
        !$omp parallel do default(shared) private(i,ithr,observed,ctf_crop) schedule(static) proc_bind(close)
        do i = 1, nptcls
            ithr = omp_get_thread_num() + 1
            call raw_work(ithr)%copy_fast(raw_imgs(i))
            call prepimg4align_cart(raw_work(ithr), noise_mask, out_work(ithr), mskrad, smpd_crop, crop_shifts(:,i), &
                &taper, ctfparms(i), observed, ctf_crop)
            if( present(sigma2) )then
                call cftc%set_ptcl(i, observed, ctf_crop, sigma2(:,i), kfromto)
            else
                call cftc%set_ptcl(i, observed, ctf_crop, kfromto)
            endif
        end do
        !$omp end parallel do
        call forget_ft_maps
        do ithr = 1, nworkers
            call raw_work(ithr)%kill
            call out_work(ithr)%kill
        end do
        deallocate(raw_work, out_work)
    end subroutine prep_cart_batch

    subroutine polarize_batch_particles3D( params, build, nptcls_here, pinds_here, src_imgs, tmp_imgs, &
        &tmp_imgs_pad )
        class(parameters),      intent(in)    :: params
        class(builder),         intent(inout) :: build
        integer,                intent(in)    :: nptcls_here
        integer,                intent(in)    :: pinds_here(nptcls_here)
        class(image),           intent(inout) :: src_imgs(nptcls_here)
        class(image),           intent(inout) :: tmp_imgs(params%nthr), tmp_imgs_pad(params%nthr)
        integer :: iptcl_batch, iptcl, ithr, pdim_interp(3)
        call tmp_imgs(1)%memoize_mask_coords
        call memoize_ft_maps(tmp_imgs(1)%get_ldim(), tmp_imgs(1)%get_smpd())
        pdim_interp = build%pftc%get_pdim_interp()
        call tmp_imgs_pad(1)%memoize4polarize_oversamp(pdim_interp)
        !$omp parallel do default(shared) private(iptcl,iptcl_batch,ithr) schedule(static) proc_bind(close)
        do iptcl_batch = 1,nptcls_here
            ithr  = omp_get_thread_num() + 1
            iptcl = pinds_here(iptcl_batch)
            call prepimg4align(params, build, iptcl, src_imgs(iptcl_batch), tmp_imgs(ithr), tmp_imgs_pad(ithr))
            call build%pftc%polarize_ptcl_pft(tmp_imgs_pad(ithr), iptcl, pdim=pdim_interp, oversamp=.true.)
            call build%pftc%set_eo(iptcl, nint(build%spproj_field%get(iptcl,'eo'))<=0 )
        end do
        !$omp end parallel do
        call forget_ft_maps
        call build%pftc%create_polar_absctfmats(build%spproj, 'ptcl3D')
        call build%pftc%memoize_ptcls
    end subroutine polarize_batch_particles3D

    !>  ptcl_imgs receives the raw full-size images, which only callers that restore
    !!  class averages from them need; omit it to skip both the buffer and the copy.
    subroutine build_batch_particles2D( params, build, nptcls_here, pinds, ptcl_imgs, ptcl_match_imgs, ptcl_match_imgs_pad )
        class(parameters),      intent(in)    :: params
        class(builder),         intent(inout) :: build
        integer,                intent(in)    :: nptcls_here
        integer,                intent(in)    :: pinds(nptcls_here)
        class(image), optional, intent(inout) :: ptcl_imgs(nptcls_here)
        class(image),           intent(inout) :: ptcl_match_imgs(params%nthr)
        class(image),           intent(inout) :: ptcl_match_imgs_pad(params%nthr)
        integer :: iptcl_batch, iptcl, ithr, pdim_interp(3)
        logical :: l_keep_raw, l_cached
        l_keep_raw = present(ptcl_imgs)
        ! When cached, ptcl_imgs receives the Fourier-cropped particle rather than the
        ! full-size original, so callers that pass it must have sized it at box_crop and
        ! told cavger_init_online to expect cropped particles.
        l_cached   = ptcl_cache_in_use(params, build)
        if( l_cached )then
            call ptcl_cache_read_batch(params, build, nptcls_here, pinds, [1,nptcls_here])
        else
            call discrete_read_imgbatch(params, build, nptcls_here, pinds, [1,nptcls_here])
        endif
        call build%pftc%reallocate_ptcls(nptcls_here, pinds)
        pdim_interp = build%pftc%get_pdim_interp()
        call ptcl_match_imgs_pad(1)%memoize4polarize_oversamp(pdim_interp)
        call ptcl_match_imgs(1)%memoize_mask_coords
        call memoize_ft_maps(ptcl_match_imgs(1)%get_ldim(), ptcl_match_imgs(1)%get_smpd())
        !$omp parallel do default(shared) private(iptcl,iptcl_batch,ithr) schedule(static) proc_bind(close)
        do iptcl_batch = 1,nptcls_here
            ithr  = omp_get_thread_num() + 1
            iptcl = pinds(iptcl_batch)
            if( l_keep_raw ) call ptcl_imgs(iptcl_batch)%copy_fast(build%imgbatch(iptcl_batch))
            if( l_cached )then
                call prepimg4align_cached(params, build, iptcl, build%imgbatch(iptcl_batch), &
                    &ptcl_match_imgs(ithr), ptcl_match_imgs_pad(ithr))
            else
                call prepimg4align(params, build, iptcl, build%imgbatch(iptcl_batch), &
                    &ptcl_match_imgs(ithr), ptcl_match_imgs_pad(ithr))
            endif
            call build%pftc%polarize_ptcl_pft(ptcl_match_imgs_pad(ithr), iptcl, pdim=pdim_interp, oversamp=.true.)
            call build%pftc%set_eo(iptcl, nint(build%spproj_field%get(iptcl,'eo'))<=0 )
        end do
        !$omp end parallel do
        call build%pftc%create_polar_absctfmats(build%spproj, 'ptcl2D')
        call build%pftc%memoize_ptcls
        call forget_ft_maps
    end subroutine build_batch_particles2D

    subroutine clean_batch_particles2D( build, ptcl_imgs, ptcl_match_imgs, ptcl_match_imgs_pad )
        use simple_imgarr_utils, only: dealloc_imgarr
        class(builder),           intent(inout) :: build
        type(image), allocatable, intent(inout) :: ptcl_imgs(:), ptcl_match_imgs(:), ptcl_match_imgs_pad(:)
        call killimgbatch(build)
        call dealloc_imgarr(ptcl_imgs)
        call dealloc_imgarr(ptcl_match_imgs)
        call dealloc_imgarr(ptcl_match_imgs_pad)
    end subroutine clean_batch_particles2D

    subroutine clean_batch_particles3D( build, ptcl_imgs, ptcl_imgs_pad )
        use simple_imgarr_utils, only: dealloc_imgarr
        class(builder),           intent(inout) :: build
        type(image), allocatable, intent(inout) :: ptcl_imgs(:)
        type(image), allocatable, intent(inout) :: ptcl_imgs_pad(:)
        call killimgbatch(build)
        call dealloc_imgarr(ptcl_imgs)
        call dealloc_imgarr(ptcl_imgs_pad)
    end subroutine clean_batch_particles3D

end module simple_matcher_ptcl_batch
