module simple_gui_utils
include 'simple_lib.f08'
use simple_image,      only: image
use simple_stack_io,   only: stack_io
#include "simple_local_flags.inc"

contains

    subroutine mic2thumb( mic, jpgname, l_neg )
        class(image),      intent(inout) :: mic
        class(string),     intent(in)    :: jpgname
        logical, optional, intent(in)    :: l_neg
        type(image) :: thumb, padded
        real, parameter :: BACKGR = 128. ! taken as centre of [0.255] for jpegs
        integer :: ldim(3), ldim_thumb(3), ldim_pad(3)
        real    :: scale, smpd_thumb
        logical :: ll_neg
        ll_neg = .false.
        if( present(l_neg) ) ll_neg = l_neg
        ldim            = mic%get_ldim()
        scale           = real(GUI_PSPECSZ)/maxval(ldim(1:2))
        ldim_thumb(1:2) = round2even(real(ldim(1:2))*scale)
        ldim_pad(1:2)   = GUI_PSPECSZ
        ldim_thumb(3)   = 1
        ldim_pad(3)     = 1
        smpd_thumb      = mic%get_smpd() / scale
        call thumb%new(ldim_thumb, smpd_thumb)
        call padded%new(ldim_pad, smpd_thumb)
        padded = BACKGR
        call mic%fft()
        call mic%clip(thumb)
        if( ll_neg ) call thumb%neg
        call thumb%ifft()
        call mic%ifft()
        call thumb%norm4viz
        call thumb%pad(padded, backgr=BACKGR)
        call padded%write_jpg(jpgname, norm=.true., quality=92)
        call thumb%kill
        call padded%kill
    end subroutine mic2thumb

    ! write tiled jpeg of mrc file
    subroutine mrc2jpeg_tiled(mrcfile, outfile, scale, ntiles, msk, n_xtiles, n_ytiles, mskdiam_px)
        class(string),                  intent(in)  :: mrcfile, outfile
        real,    optional,              intent(out) :: scale
        integer, optional,              intent(out) :: ntiles, n_xtiles, n_ytiles
        logical, optional, allocatable, intent(in)  :: msk(:)
        integer, optional,              intent(in)  :: mskdiam_px
        type(image)    :: img, img_pad, img_jpeg
        type(stack_io) :: stkio_r
        integer        :: ldim_stk(3)
        integer        :: ncls_here, xtiles, ytiles, icls, ix, iy, l_ntiles
        real           :: smpd
        smpd = 1.0
        if(.not. file_exists(mrcfile)) return
        call find_ldim_nptcls(mrcfile, ldim_stk, ncls_here)
        xtiles = floor(sqrt(real(ncls_here)))
        ytiles = ceiling(real(ncls_here) / real(xtiles))
        call img%new([ldim_stk(1), ldim_stk(1), 1], smpd)
        call img_pad%new([JPEG_DIM, JPEG_DIM, 1], smpd)
        call img_jpeg%new([xtiles * JPEG_DIM, ytiles * JPEG_DIM, 1], smpd)
        call stkio_r%open(mrcfile, smpd, 'read', bufsz=ncls_here)
        call stkio_r%read_whole
        ix = 1
        iy = 1
        l_ntiles = 0
        do icls=1, ncls_here
            if(present(msk)) then
                if( .not. msk(icls)) cycle
            end if
            call img%zero_and_unflag_ft
            call stkio_r%get_image(icls, img)
            if(present(mskdiam_px)) call img%mask(mskdiam_px / 2.0, 'softavg')
            call img%fft
            if(ldim_stk(1) > JPEG_DIM) then
                call img%clip(img_pad)
            else
                call img%pad(img_pad, backgr=0., antialiasing=.false.)
            end if
            call img_pad%ifft
            call img_jpeg%tile(img_pad, ix, iy)
            l_ntiles = l_ntiles + 1
            ix = ix + 1
            if(ix > xtiles) then
                ix = 1
                iy = iy + 1
            end if
        enddo
        call stkio_r%close()
        call img_jpeg%write_jpg(string(outfile%to_char() // ".tmp"))
        call simple_rename(outfile%to_char()// ".tmp", outfile, overwrite=.true.)
        if(present(scale))  scale  = real(JPEG_DIM) / ldim_stk(1)
        if(present(ntiles)) ntiles = l_ntiles
        if(present(n_xtiles)) n_xtiles = xtiles
        if(present(n_ytiles)) n_ytiles = ytiles
        call img%kill()
        call img_pad%kill()
        call img_jpeg%kill()
    end subroutine mrc2jpeg_tiled

end module simple_gui_utils
