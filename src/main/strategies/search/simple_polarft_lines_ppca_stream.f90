!@descr: streaming helpers for feeding prepared particle polar Fourier lines into complex PPCA
module simple_polarft_lines_ppca_stream
use simple_pftc_srch_api
use simple_builder,             only: builder
use simple_complex_ppca,        only: complex_ppca
use simple_matcher_ptcl_batch,  only: alloc_ptcl_imgs, build_batch_particles2D, clean_batch_particles2D
use simple_imgarr_utils,        only: alloc_imgarr
implicit none

public :: stream_pft_lines_ppca
public :: denoise_write_pft_lines_ppca
private
#include "simple_local_flags.inc"

contains

    subroutine stream_pft_lines_ppca( params, build, q, model, nlines )
        class(parameters),        intent(inout) :: params
        class(builder),           intent(inout) :: build
        integer,                  intent(in)    :: q
        class(complex_ppca),      intent(inout) :: model
        integer(kind=8), optional,intent(out)   :: nlines
        type(image), allocatable :: ptcl_imgs(:), ptcl_match_imgs(:), ptcl_match_imgs_pad(:)
        integer, allocatable :: pinds(:), batches(:,:)
        integer :: batchsz_max, nbatches, batch_start, batch_end, batchsz
        integer :: kfromto(2), nk, ibatch
        integer(kind=8) :: nlines_here
        ! The caller owns application-specific pftc preparation. This routine only
        ! prepares the particle batch scratch needed to populate the already-configured
        ! polar Fourier container and stream its lines into the PPCA model.
        call get_stream_batch_layout(params, build, pinds, batches, batchsz_max)
        kfromto = params%kfromto
        nk = kfromto(2) - kfromto(1) + 1
        call model%new(0, nk, q)
        call model%stream_reset()
        call alloc_ptcl_imgs(params, build, ptcl_match_imgs, ptcl_match_imgs_pad, batchsz_max)
        call alloc_imgarr(batchsz_max, [params%box, params%box, 1], params%smpd, ptcl_imgs)
        nlines_here = 0_8
        nbatches = size(batches,1)
        do ibatch = 1,nbatches
            batch_start = batches(ibatch,1)
            batch_end   = batches(ibatch,2)
            batchsz     = batch_end - batch_start + 1
            call build_batch_particles2D(params, build, batchsz, pinds(batch_start:batch_end), &
                &ptcl_imgs(1:batchsz), ptcl_match_imgs, ptcl_match_imgs_pad)
            call update_model_from_batch_lines(build%pftc, pinds(batch_start:batch_end), model, nlines_here)
        enddo
        call model%stream_finalize()
        call clean_batch_particles2D(build, ptcl_imgs, ptcl_match_imgs, ptcl_match_imgs_pad)
        deallocate(pinds, batches)
        if( present(nlines) ) nlines = nlines_here
    end subroutine stream_pft_lines_ppca

    subroutine denoise_write_pft_lines_ppca( params, build, model, fname, nlines )
        class(parameters),        intent(inout) :: params
        class(builder),           intent(inout) :: build
        class(complex_ppca),      intent(in)    :: model
        class(string),            intent(in)    :: fname
        integer(kind=8), optional,intent(out)   :: nlines
        type(image), allocatable :: ptcl_imgs(:), ptcl_match_imgs(:), ptcl_match_imgs_pad(:)
        integer, allocatable :: pinds(:), batches(:,:)
        integer :: batchsz_max, nbatches, batch_start, batch_end, batchsz
        integer :: ibatch
        integer(kind=8) :: nlines_here
        call get_stream_batch_layout(params, build, pinds, batches, batchsz_max)
        call alloc_ptcl_imgs(params, build, ptcl_match_imgs, ptcl_match_imgs_pad, batchsz_max)
        call alloc_imgarr(batchsz_max, [params%box, params%box, 1], params%smpd, ptcl_imgs)
        nlines_here = 0_8
        nbatches = size(batches,1)
        do ibatch = 1,nbatches
            batch_start = batches(ibatch,1)
            batch_end   = batches(ibatch,2)
            batchsz     = batch_end - batch_start + 1
            call build_batch_particles2D(params, build, batchsz, pinds(batch_start:batch_end), &
                &ptcl_imgs(1:batchsz), ptcl_match_imgs, ptcl_match_imgs_pad)
            call denoise_batch_lines(build%pftc, pinds(batch_start:batch_end), model, nlines_here)
            call build%pftc%write_ptcl_pft_range(fname, size(pinds), batch_start, batch_end)
        enddo
        call clean_batch_particles2D(build, ptcl_imgs, ptcl_match_imgs, ptcl_match_imgs_pad)
        deallocate(pinds, batches)
        if( present(nlines) ) nlines = nlines_here
    end subroutine denoise_write_pft_lines_ppca

    subroutine get_stream_batch_layout( params, build, pinds, batches, batchsz_max )
        class(parameters),        intent(in)    :: params
        class(builder),           intent(inout) :: build
        integer, allocatable,     intent(out)   :: pinds(:), batches(:,:)
        integer,                  intent(out)   :: batchsz_max
        integer :: nptcls2update, nbatches
        if( allocated(pinds) ) deallocate(pinds)
        call build%spproj_field%sample4update_all([params%fromp,params%top], nptcls2update, pinds, .true.)
        batchsz_max = min(nptcls2update, params%nthr * BATCHTHRSZ)
        nbatches    = ceiling(real(nptcls2update) / real(batchsz_max))
        batches     = split_nobjs_even(nptcls2update, nbatches)
        batchsz_max = maxval(batches(:,2) - batches(:,1) + 1)
    end subroutine get_stream_batch_layout

    subroutine update_model_from_batch_lines( pftc, pinds_here, model, nlines )
        class(polarft_calc),      intent(inout) :: pftc
        integer,                  intent(in)    :: pinds_here(:)
        class(complex_ppca),      intent(inout) :: model
        integer(kind=8),          intent(inout) :: nlines
        integer :: iptcl_map, iptcl, irot, kfromto(2), nk
        complex(sp), allocatable :: line_sp(:)
        complex(dp), allocatable :: line_dp(:)
        kfromto = pftc%get_kfromto()
        nk = kfromto(2) - kfromto(1) + 1
        allocate(line_sp(kfromto(1):kfromto(2)), line_dp(nk))
        do iptcl_map = 1,size(pinds_here)
            iptcl = pinds_here(iptcl_map)
            do irot = 1,pftc%get_pftsz()
                call pftc%get_ptcl_line(iptcl, irot, line_sp)
                line_dp = cmplx(real(line_sp,dp), real(aimag(line_sp),dp), kind=dp)
                call model%stream_update(line_dp)
                nlines = nlines + 1_8
            enddo
        enddo
        deallocate(line_sp, line_dp)
    end subroutine update_model_from_batch_lines

    subroutine denoise_batch_lines( pftc, pinds_here, model, nlines )
        class(polarft_calc),      intent(inout) :: pftc
        integer,                  intent(in)    :: pinds_here(:)
        class(complex_ppca),      intent(in)    :: model
        integer(kind=8),          intent(inout) :: nlines
        integer :: iptcl_map, iptcl, irot, kfromto(2), nk
        complex(sp), allocatable :: pft_sp(:,:), denoised_sp(:,:)
        complex(dp), allocatable :: line_dp(:), denoised_dp(:)
        kfromto = pftc%get_kfromto()
        nk = kfromto(2) - kfromto(1) + 1
        allocate(line_dp(nk), denoised_dp(nk))
        do iptcl_map = 1,size(pinds_here)
            iptcl = pinds_here(iptcl_map)
            pft_sp = pftc%allocate_ptcl_pft()
            denoised_sp = pftc%allocate_ptcl_pft()
            call pftc%get_ptcl_pft(iptcl, pft_sp)
            ! Preserve the exact polar Fourier geometry and any shells outside the
            ! denoising band. We only replace the modeled kfromto segment.
            denoised_sp = pft_sp
            do irot = 1,pftc%get_pftsz()
                line_dp = cmplx(real(pft_sp(irot,kfromto(1):kfromto(2)),dp), &
                    &real(aimag(pft_sp(irot,kfromto(1):kfromto(2))),dp), kind=dp)
                call model%denoise(line_dp, denoised_dp)
                denoised_sp(irot,kfromto(1):kfromto(2)) = cmplx(real(denoised_dp,sp), real(aimag(denoised_dp),sp), kind=sp)
                nlines = nlines + 1_8
            enddo
            call pftc%set_ptcl_pft(iptcl, denoised_sp)
            deallocate(pft_sp, denoised_sp)
        enddo
        deallocate(line_dp, denoised_dp)
    end subroutine denoise_batch_lines

end module simple_polarft_lines_ppca_stream
