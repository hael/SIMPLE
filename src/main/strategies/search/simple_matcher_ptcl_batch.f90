!@descr: particle batch preparation helpers for matcher workflows
module simple_matcher_ptcl_batch
use simple_pftc_srch_api
use simple_builder,         only: builder
use simple_euclid_sigma2,   only: sigma2_star_from_iter
use simple_matcher_ptcl_io, only: prepimgbatch, discrete_read_imgbatch, discrete_read_imgbatch_source, killimgbatch
use simple_matcher_2Dprep,  only: prepimg4align
implicit none

public :: prep_sigmas_objfun, alloc_ptcl_imgs
public :: build_batch_particles3D, build_batch_particles2D
public :: clean_batch_particles2D, clean_batch_particles3D
private
#include "simple_local_flags.inc"

contains

    subroutine prep_sigmas_objfun( params, build, l_stream )
        class(parameters), intent(inout) :: params
        class(builder),    intent(inout) :: build
        logical,           intent(in)    :: l_stream
        type(string)      :: fname
        logical           :: l_group_only_init
        if( params%cc_objfun == OBJFUN_EUCLID )then
            fname = SIGMA2_FBODY//int2str_pad(params%part,params%numlen)//'.dat'
            call build%esig%new(params, build%pftc, fname, params%box)
            l_group_only_init = (.not. file_exists(fname)) .and. file_exists(sigma2_star_from_iter(params%which_iter))
            if( l_stream .or. l_group_only_init )then
                call build%esig%read_groups(build%spproj_field)
                call build%esig%allocate_ptcls
            else
                call build%esig%read_part(  build%spproj_field)
                call build%esig%read_groups(build%spproj_field)
            endif
            call fname%kill
        end if
    end subroutine prep_sigmas_objfun

    subroutine alloc_ptcl_imgs( params, build, ptcl_imgs, ptcl_imgs_pad, batchsz )
        class(parameters),        intent(inout) :: params
        class(builder),           intent(inout) :: build
        type(image), allocatable, intent(inout) :: ptcl_imgs(:)
        type(image), allocatable, intent(inout) :: ptcl_imgs_pad(:)
        integer,                  intent(in)    :: batchsz
        integer           :: ithr
        call prepimgbatch(params, build, batchsz)
        allocate(ptcl_imgs(nthr_glob), ptcl_imgs_pad(nthr_glob))
        !$omp parallel do default(shared) private(ithr) schedule(static) proc_bind(close)
        do ithr = 1,nthr_glob
            call ptcl_imgs(ithr)%new(    [params%box_crop,  params%box_crop,  1], params%smpd_crop, wthreads=.false.)
            call ptcl_imgs_pad(ithr)%new([params%box_croppd,params%box_croppd,1], params%smpd_crop, wthreads=.false.)
        enddo
        !$omp end parallel do
    end subroutine alloc_ptcl_imgs

    subroutine build_batch_particles3D( params, build, nptcls_here, pinds_here, tmp_imgs, tmp_imgs_pad, imgs4rec )
        class(parameters),      intent(in)    :: params
        class(builder),         intent(inout) :: build
        integer,                intent(in)    :: nptcls_here
        integer,                intent(in)    :: pinds_here(nptcls_here)
        class(image),           intent(inout) :: tmp_imgs(params%nthr), tmp_imgs_pad(params%nthr)
        class(image), optional, intent(inout) :: imgs4rec(nptcls_here)
        integer :: iptcl_batch, iptcl, ithr, pdim_interp(3)
        logical :: l_backup_imgs, l_copy_rec_from_match
        l_backup_imgs = present(imgs4rec)
        l_copy_rec_from_match = l_backup_imgs .and. trim(params%recimg_source) == trim(params%matchimg_source)
        call build%pftc%reallocate_ptcls(nptcls_here, pinds_here)
        if( trim(params%matchimg_source) == 'raw' )then
            call discrete_read_imgbatch(params, build, nptcls_here, pinds_here, [1,nptcls_here])
        else
            call discrete_read_imgbatch_source(params, build, trim(params%matchimg_source), &
                nptcls_here, pinds_here, [1,nptcls_here], build%imgbatch(:nptcls_here))
        endif
        if( l_backup_imgs .and. .not. l_copy_rec_from_match )then
            call discrete_read_imgbatch_source(params, build, trim(params%recimg_source), &
                nptcls_here, pinds_here, [1,nptcls_here], imgs4rec)
        endif
        call tmp_imgs(1)%memoize_mask_coords
        call memoize_ft_maps(tmp_imgs(1)%get_ldim(), tmp_imgs(1)%get_smpd())
        pdim_interp = build%pftc%get_pdim_interp()
        call tmp_imgs_pad(1)%memoize4polarize_oversamp(pdim_interp)
        !$omp parallel do default(shared) private(iptcl,iptcl_batch,ithr) schedule(static) proc_bind(close)
        do iptcl_batch = 1,nptcls_here
            ithr  = omp_get_thread_num() + 1
            iptcl = pinds_here(iptcl_batch)
            if( l_copy_rec_from_match ) call imgs4rec(iptcl_batch)%copy_fast(build%imgbatch(iptcl_batch))
            call prepimg4align(params, build, iptcl, build%imgbatch(iptcl_batch), tmp_imgs(ithr), tmp_imgs_pad(ithr))
            call build%pftc%polarize_ptcl_pft(tmp_imgs_pad(ithr), iptcl, pdim=pdim_interp, oversamp=.true.)
            call build%pftc%set_eo(iptcl, nint(build%spproj_field%get(iptcl,'eo'))<=0 )
        end do
        !$omp end parallel do
        call forget_ft_maps
        call build%pftc%create_polar_absctfmats(build%spproj, 'ptcl3D')
        call build%pftc%memoize_ptcls
    end subroutine build_batch_particles3D

    subroutine build_batch_particles2D( params, build, nptcls_here, pinds, ptcl_imgs, ptcl_match_imgs, ptcl_match_imgs_pad )
        class(parameters), intent(in)    :: params
        class(builder),    intent(inout) :: build
        integer,           intent(in)    :: nptcls_here
        integer,           intent(in)    :: pinds(nptcls_here)
        class(image),      intent(inout) :: ptcl_imgs(nptcls_here)
        class(image),      intent(inout) :: ptcl_match_imgs(params%nthr)
        class(image),      intent(inout) :: ptcl_match_imgs_pad(params%nthr)
        integer :: iptcl_batch, iptcl, ithr, pdim_interp(3)
        call discrete_read_imgbatch(params, build, nptcls_here, pinds, [1,nptcls_here])
        call build%pftc%reallocate_ptcls(nptcls_here, pinds)
        pdim_interp = build%pftc%get_pdim_interp()
        call ptcl_match_imgs_pad(1)%memoize4polarize_oversamp(pdim_interp)
        call ptcl_match_imgs(1)%memoize_mask_coords
        call memoize_ft_maps(ptcl_match_imgs(1)%get_ldim(), ptcl_match_imgs(1)%get_smpd())
        !$omp parallel do default(shared) private(iptcl,iptcl_batch,ithr) schedule(static) proc_bind(close)
        do iptcl_batch = 1,nptcls_here
            ithr  = omp_get_thread_num() + 1
            iptcl = pinds(iptcl_batch)
            call ptcl_imgs(iptcl_batch)%copy_fast(build%imgbatch(iptcl_batch))
            call prepimg4align(params, build, iptcl, build%imgbatch(iptcl_batch), ptcl_match_imgs(ithr), ptcl_match_imgs_pad(ithr))
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

    subroutine clean_batch_particles3D( build, ptcl_imgs, ptcl_imgs_pad, imgs4rec )
        use simple_imgarr_utils, only: dealloc_imgarr
        class(builder),                     intent(inout) :: build
        type(image), allocatable,           intent(inout) :: ptcl_imgs(:)
        type(image), allocatable,           intent(inout) :: ptcl_imgs_pad(:)
        type(image), allocatable, optional, intent(inout) :: imgs4rec(:)
        call killimgbatch(build)
        if( present(imgs4rec) ) call dealloc_imgarr(imgs4rec)
        call dealloc_imgarr(ptcl_imgs)
        call dealloc_imgarr(ptcl_imgs_pad)
    end subroutine clean_batch_particles3D

end module simple_matcher_ptcl_batch
