module simple_gui_utils
include 'simple_lib.f08'
use simple_image,      only: image
use simple_stack_io,   only: stack_io
use simple_oris,       only: oris
use simple_sp_project, only: sp_project
implicit none
#include "simple_local_flags.inc"

interface cavgs2jpg
    module procedure cavgs2jpg_1
    module procedure cavgs2jpg_2
end interface cavgs2jpg

contains

    subroutine mic2thumb( mic, jpgname )
        class(image),     intent(inout) :: mic
        character(len=*), intent(in)    :: jpgname
        type(image) :: thumb
        integer :: ldim(3), ldim_thumb(3)
        real    :: scale, smpd_thumb
        ldim            = mic%get_ldim()
        scale           = real(GUI_PSPECSZ)/maxval(ldim(1:2))
        ldim_thumb(1:2) = round2even(real(ldim(1:2))*scale)
        ldim_thumb(3)   = 1
        smpd_thumb      = mic%get_smpd() / scale
        call thumb%new(ldim_thumb, smpd_thumb)
        call mic%fft()
        call mic%clip(thumb)
        call thumb%ifft()
        call mic%ifft()
        call thumb%write_jpg(jpgname, norm=.true., quality=92)
        call thumb%kill
    end subroutine mic2thumb

    subroutine cavgs2jpg_1( spproj, cavg_inds, jpgname, xtiles, ytiles )
        class(sp_project),    intent(inout) :: spproj
        integer, allocatable, intent(inout) :: cavg_inds(:)
        character(len=*),     intent(in)    :: jpgname
        integer,              intent(out)   :: xtiles, ytiles
        logical,          allocatable :: cls_mask(:)
        character(len=:), allocatable :: cavgsstk, stkpath
        type(image)    :: img, jpegimg
        type(stack_io) :: stkio_r
        integer        :: i, icls, isel, ncls, ncls_stk, ldim_read(3), ncls_sel, ix, iy, ntiles
        real           :: smpd
        ncls = spproj%os_cls2D%get_noris()
        if( ncls == 0 ) return
        if( allocated(cavg_inds) ) deallocate(cavg_inds)
        allocate(cls_mask(ncls), cavg_inds(ncls))
        cls_mask = .true.
        do icls = 1, ncls
            if( spproj%os_cls2D%get_int(icls,'pop') == 0 ) cls_mask(icls) = .false.
            if( spproj%os_cls2D%get_state(icls)     == 0 ) cls_mask(icls) = .false.
            cavg_inds(icls) = icls
        enddo
        ncls_sel = count(cls_mask)
        if( ncls_sel == 0 ) return
        call spproj%get_cavgs_stk(cavgsstk, ncls_stk, smpd, imgkind='cavg', stkpath=stkpath)
        if(.not. file_exists(cavgsstk)) cavgsstk = trim(stkpath) // '/' // trim(cavgsstk)
        if(.not. file_exists(cavgsstk)) THROW_HARD('cavgs stk does not exist')
        if( ncls /= ncls_stk ) THROW_HARD('Inconsistent # cavgs in spproj and stack file')
        call stkio_r%open(trim(cavgsstk), smpd, 'read', bufsz=ncls)
        ldim_read    = stkio_r%get_ldim()
        ldim_read(3) = 1
        call stkio_r%read_whole
        xtiles = floor(sqrt(real(ncls_sel)))
        ytiles = ceiling(real(ncls_sel) / real(xtiles))
        call jpegimg%new([xtiles * JPEG_DIM, ytiles * JPEG_DIM, 1], smpd)
        isel   = 0
        ix     = 1
        iy     = 1
        ntiles = 0
        do icls = 1,ncls
            if( cls_mask(icls) ) then
                isel = isel + 1
                call img%new(ldim_read, smpd)
                call stkio_r%get_image(icls, img)
                call img%fft
                if(ldim_read(1) > JPEG_DIM) then
                    call img%clip_inplace([JPEG_DIM,JPEG_DIM,1])
                else
                    call img%pad_inplace([JPEG_DIM,JPEG_DIM,1], backgr=0., antialiasing=.false.)
                end if
                call img%ifft
                call jpegimg%tile(img, ix, iy)
                ntiles = ntiles + 1
                ix = ix + 1
                if(ix > xtiles) then
                    ix = 1
                    iy = iy + 1
                end if
            endif
        enddo
        call jpegimg%write_jpg(jpgname)
        call stkio_r%close
        call img%kill
        call jpegimg%kill
    end subroutine cavgs2jpg_1

    subroutine cavgs2jpg_2( cavgsstk, jpgname, xtiles, ytiles )
        character(len=*), intent(in)  :: cavgsstk
        character(len=*), intent(in)  :: jpgname
        integer,          intent(out) :: xtiles, ytiles
        type(image)    :: img, jpegimg
        type(stack_io) :: stkio_r
        integer        :: i, icls, ncls, ldim_read(3), ix, iy, ntiles
        real           :: smpd
        if(.not. file_exists(cavgsstk)) THROW_HARD('cavgs stk does not exist')
        call find_ldim_nptcls(cavgsstk, ldim_read, ncls, smpd)
        call stkio_r%open(trim(cavgsstk), smpd, 'read', bufsz=ncls)
        ldim_read    = stkio_r%get_ldim()
        ldim_read(3) = 1
        call stkio_r%read_whole
        xtiles = floor(sqrt(real(ncls)))
        ytiles = ceiling(real(ncls) / real(xtiles))
        call jpegimg%new([xtiles * JPEG_DIM, ytiles * JPEG_DIM, 1], smpd)
        ix     = 1
        iy     = 1
        ntiles = 0
        do icls = 1,ncls
            call img%new(ldim_read, smpd)
            call stkio_r%get_image(icls, img)
            call img%fft
            if(ldim_read(1) > JPEG_DIM) then
                call img%clip_inplace([JPEG_DIM,JPEG_DIM,1])
            else
                call img%pad_inplace([JPEG_DIM,JPEG_DIM,1], backgr=0., antialiasing=.false.)
            end if
            call img%ifft
            call jpegimg%tile(img, ix, iy)
            ntiles = ntiles + 1
            ix = ix + 1
            if(ix > xtiles) then
                ix = 1
                iy = iy + 1
            end if
        enddo
        call jpegimg%write_jpg(jpgname)
        call stkio_r%close
        call img%kill
        call jpegimg%kill
    end subroutine cavgs2jpg_2

    subroutine shape_ranked_cavgs2jpg( spproj, cavg_inds, jpgname, xtiles, ytiles )
        class(sp_project),    intent(inout) :: spproj
        integer, allocatable, intent(inout) :: cavg_inds(:)
        character(len=*),     intent(in)    :: jpgname
        integer,              intent(out)   :: xtiles, ytiles
        integer,          allocatable :: shape_ranks(:)
        character(len=:), allocatable :: cavgsstk, stkpath
        type(image)    :: img, jpegimg
        type(stack_io) :: stkio_r
        integer        :: i, icls, ncls, ncls_stk, ldim_read(3), ncls_sel, ix, iy, ntiles
        real           :: smpd
        ncls = spproj%os_cls2D%get_noris()
        if( ncls == 0 ) return
        if( .not. spproj%os_cls2D%isthere('shape_rank') ) THROW_HARD('No shape ranking available!')
        if( allocated(cavg_inds) ) deallocate(cavg_inds)
        allocate(shape_ranks(ncls), cavg_inds(ncls))
        do icls = 1, ncls
            shape_ranks(icls) = spproj%os_cls2D%get_int(icls,'shape_rank')
            cavg_inds(icls)   = icls
        end do
        cavg_inds   = pack(cavg_inds,   mask=shape_ranks > 0)
        shape_ranks = pack(shape_ranks, mask=shape_ranks > 0)
        call hpsort(shape_ranks, cavg_inds)
        call spproj%get_cavgs_stk(cavgsstk, ncls_stk, smpd, imgkind='cavg', stkpath=stkpath)
        if(.not. file_exists(cavgsstk)) cavgsstk = trim(stkpath) // '/' // trim(cavgsstk)
        if(.not. file_exists(cavgsstk)) THROW_HARD('cavgs stk does not exist')
        if( ncls /= ncls_stk ) THROW_HARD('Inconsistent # cavgs in spproj and stack file')
        call stkio_r%open(trim(cavgsstk), smpd, 'read', bufsz=ncls)
        ldim_read    = stkio_r%get_ldim()
        ldim_read(3) = 1
        call stkio_r%read_whole
        ncls_sel = size(cavg_inds)
        xtiles   = floor(sqrt(real(ncls_sel)))
        ytiles   = ceiling(real(ncls_sel) / real(xtiles))
        call jpegimg%new([xtiles * JPEG_DIM, ytiles * JPEG_DIM, 1], smpd)
        ix     = 1
        iy     = 1
        ntiles = 0
        do icls = 1,ncls_sel
            call img%new(ldim_read, smpd)
            call stkio_r%get_image(cavg_inds(icls), img)
            call img%fft
            if(ldim_read(1) > JPEG_DIM) then
                call img%clip_inplace([JPEG_DIM,JPEG_DIM,1])
            else
                call img%pad_inplace([JPEG_DIM,JPEG_DIM,1], backgr=0., antialiasing=.false.)
            end if
            call img%ifft
            call jpegimg%tile(img, ix, iy)
            ntiles = ntiles + 1
            ix = ix + 1
            if(ix > xtiles) then
                ix = 1
                iy = iy + 1
            end if
        enddo
        call jpegimg%write_jpg(jpgname)
        call stkio_r%close
        call img%kill
        call jpegimg%kill
    end subroutine shape_ranked_cavgs2jpg

end module simple_gui_utils
