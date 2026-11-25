!@descr: various GUI utility functions
module simple_gui_utils
use simple_core_module_api
use simple_image,      only: image
use simple_stack_io,   only: stack_io
#include "simple_local_flags.inc"

contains

    subroutine mic2thumb( mic, jpgname, l_neg )
        class(image),      intent(inout)        :: mic
        class(string),     intent(in)           :: jpgname
        logical,           optional, intent(in) :: l_neg
        type(image)     :: thumb, padded
        real, parameter :: BACKGR = 128. ! taken as centre of [0.255] for jpegs
        integer         :: ldim(3), ldim_thumb(3), ldim_pad(3)
        real            :: scale, smpd_thumb
        logical         :: ll_neg
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
        if( any(ldim_thumb > ldim) )then
            call mic%pad(thumb, antialiasing=.false.)
        else
            call mic%clip(thumb)
        endif
        if( ll_neg ) call thumb%neg
        call thumb%ifft()
        call mic%ifft()
        call thumb%norm4viz
        call thumb%pad(padded, backgr=BACKGR)
        call padded%write_jpg(jpgname, norm=.true., quality=92)
        call thumb%kill
        call padded%kill
    end subroutine mic2thumb

    !> write tiled jpeg of mrc file; frames (one code per image: 0 none, 1 red, 2 orange) writes it in colour with a
    !! frame around each tile of a coded image
    subroutine mrc2jpeg_tiled(mrcfile, outfile, scale, ntiles, msk, n_xtiles, n_ytiles, mskdiam_px, frames)
        class(string),                  intent(in)                         :: mrcfile, outfile
        real,                           optional,              intent(out) :: scale
        integer,                        optional,              intent(out) :: ntiles, n_xtiles, n_ytiles
        logical,                        optional, allocatable, intent(in)  :: msk(:)
        integer,                        optional,              intent(in)  :: mskdiam_px
        integer,                        optional,              intent(in)  :: frames(:)
        type(image)          :: img, img_pad, img_jpeg
        type(stack_io)       :: stkio_r
        integer, allocatable :: tile_xy(:,:)
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
        ! mask memoization
        call img%memoize_mask_coords
        allocate(tile_xy(2,ncls_here), source=0)
        do icls=1, ncls_here
            if(present(msk)) then
                if( .not. msk(icls)) cycle
            end if
            call img%zero_and_unflag_ft
            call stkio_r%get_image(icls, img)
            if(present(mskdiam_px)) call img%mask2D_softavg(mskdiam_px / 2.0)
            call img%fft
            if(ldim_stk(1) > JPEG_DIM) then
                call img%clip(img_pad)
            else
                call img%pad(img_pad, backgr=0., antialiasing=.false.)
            end if
            call img_pad%ifft
            call img_jpeg%tile(img_pad, ix, iy)
            tile_xy(:,icls) = [ix, iy]
            l_ntiles = l_ntiles + 1
            ix = ix + 1
            if(ix > xtiles) then
                ix = 1
                iy = iy + 1
            end if
        enddo
        call stkio_r%close()
        if( present(frames) )then
            call write_framed(string(outfile%to_char() // ".tmp"))
        else
            call img_jpeg%write_jpg(string(outfile%to_char() // ".tmp"))
        endif
        call simple_rename(outfile%to_char()// ".tmp", outfile)
        if(present(scale))  scale  = real(JPEG_DIM) / ldim_stk(1)
        if(present(ntiles)) ntiles = l_ntiles
        if(present(n_xtiles)) n_xtiles = xtiles
        if(present(n_ytiles)) n_ytiles = ytiles
        call img%kill()
        call img_pad%kill()
        call img_jpeg%kill()

      contains

        !> the montage in colour, scaled between its extremes as the grey writer does, a frame of FRAME_PX around
        !! each placed tile whose image has a code
        subroutine write_framed( fname )
            use simple_jpg, only: write_rgb_jpeg
            class(string), intent(in) :: fname
            integer, parameter :: FRAME_PX = 3
            real,    parameter :: RGB_RED(3) = [230., 40., 40.] / 255., RGB_ORANGE(3) = [255., 160., 0.] / 255.
            real, allocatable  :: rmat(:,:,:), rgb(:,:,:)
            real    :: lo, hi, col(3)
            integer :: w, h, i, j, k, x0, y0, d
            if( size(frames) /= ncls_here ) THROW_HARD('one frame code per image; mrc2jpeg_tiled')
            rmat = img_jpeg%get_rmat()
            w    = xtiles * JPEG_DIM
            h    = ytiles * JPEG_DIM
            lo   = minval(rmat(:w,:h,1))
            hi   = maxval(rmat(:w,:h,1))
            allocate(rgb(3,w,h))
            do j = 1,h
                do i = 1,w
                    rgb(:,i,j) = 0.5
                    if( hi > lo ) rgb(:,i,j) = (rmat(i,j,1) - lo) / (hi - lo)
                enddo
            enddo
            do k = 1,ncls_here
                if( frames(k) < 1 .or. tile_xy(1,k) < 1 ) cycle
                col = RGB_RED
                if( frames(k) == 2 ) col = RGB_ORANGE
                x0 = (tile_xy(1,k) - 1) * JPEG_DIM
                y0 = (tile_xy(2,k) - 1) * JPEG_DIM
                do d = 1,FRAME_PX
                    do i = x0 + 1,x0 + JPEG_DIM
                        rgb(:,i,y0 + d)                = col
                        rgb(:,i,y0 + JPEG_DIM + 1 - d) = col
                    enddo
                    do j = y0 + 1,y0 + JPEG_DIM
                        rgb(:,x0 + d,j)                = col
                        rgb(:,x0 + JPEG_DIM + 1 - d,j) = col
                    enddo
                enddo
            enddo
            call write_rgb_jpeg(fname%to_char(), rgb, quality=90)
        end subroutine write_framed

    end subroutine mrc2jpeg_tiled

end module simple_gui_utils
