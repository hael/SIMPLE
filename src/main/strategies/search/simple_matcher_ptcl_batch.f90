!@descr: particle batch preparation helpers for matcher workflows
module simple_matcher_ptcl_batch
use simple_pftc_srch_api
use simple_builder,         only: builder
use simple_matcher_ptcl_io, only: prepimgbatch, discrete_read_imgbatch, killimgbatch
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
        if( params%cc_objfun == OBJFUN_EUCLID )then
            fname = SIGMA2_FBODY//int2str_pad(params%part,params%numlen)//'.dat'
            call build%esig%new(params, build%pftc, fname, params%box)
            if( l_stream )then
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
        complex, allocatable :: pft(:,:)
        integer :: iptcl_batch, iptcl, ithr
        logical :: l_backup_imgs
        l_backup_imgs = present(imgs4rec)
        call build%pftc%reallocate_ptcls(nptcls_here, pinds_here)
        call discrete_read_imgbatch(params, build, nptcls_here, pinds_here, [1,nptcls_here])
        call tmp_imgs(1)%memoize_mask_coords
        call memoize_ft_maps(tmp_imgs(1)%get_ldim(), tmp_imgs(1)%get_smpd())
        call tmp_imgs_pad(1)%memoize4polarize_oversamp(build%pftc%get_pdim_interp())
        !$omp parallel do default(shared) private(iptcl,iptcl_batch,ithr,pft) schedule(static) proc_bind(close)
        do iptcl_batch = 1,nptcls_here
            ithr  = omp_get_thread_num() + 1
            iptcl = pinds_here(iptcl_batch)
            if( l_backup_imgs ) call imgs4rec(iptcl_batch)%copy_fast(build%imgbatch(iptcl_batch))
            call prepimg4align(params, build, iptcl, build%imgbatch(iptcl_batch), tmp_imgs(ithr), tmp_imgs_pad(ithr))
            pft = build%pftc%allocate_ptcl_pft()
            call tmp_imgs_pad(ithr)%polarize_oversamp(pft)
            call build%pftc%set_ptcl_pft(iptcl, pft)
            deallocate(pft)
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
        complex, allocatable :: pft(:,:)
        integer :: iptcl_batch, iptcl, ithr
        call discrete_read_imgbatch(params, build, nptcls_here, pinds, [1,nptcls_here])
        call build%pftc%reallocate_ptcls(nptcls_here, pinds)
        call ptcl_match_imgs_pad(1)%memoize4polarize_oversamp(build%pftc%get_pdim_interp())
        call ptcl_match_imgs(1)%memoize_mask_coords
        call memoize_ft_maps(ptcl_match_imgs(1)%get_ldim(), ptcl_match_imgs(1)%get_smpd())
        !$omp parallel do default(shared) private(iptcl,iptcl_batch,ithr,pft) schedule(static) proc_bind(close)
        do iptcl_batch = 1,nptcls_here
            ithr  = omp_get_thread_num() + 1
            iptcl = pinds(iptcl_batch)
            call ptcl_imgs(iptcl_batch)%copy_fast(build%imgbatch(iptcl_batch))
            call prepimg4align(params, build, iptcl, build%imgbatch(iptcl_batch), ptcl_match_imgs(ithr), ptcl_match_imgs_pad(ithr))
            pft = build%pftc%allocate_pft()
            call ptcl_match_imgs_pad(ithr)%polarize_oversamp(pft)
            call build%pftc%set_ptcl_pft(iptcl, pft)
            deallocate(pft)
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
