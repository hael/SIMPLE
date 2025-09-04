! stack image processing routines for SPIDER/MRC files and other useful img jiffys
module simple_imgproc
include 'simple_lib.f08'
use simple_image,    only: image
use simple_stack_io, only: stack_io
implicit none

public :: make_pcavecs, make_pcavol, mrc2jpeg_tiled

private

#include "simple_local_flags.inc"

contains

    subroutine make_pcavecs( imgs, D, avg, pcavecs, l_mask, transp, avg_pix, do_fft )
        class(image),                intent(inout) :: imgs(:)       !< images to serialize
        integer,                     intent(out)   :: D             !< vector dimension
        real, allocatable,           intent(inout) :: avg(:)        !< average of pcavecs
        real, allocatable,           intent(inout) :: pcavecs(:,:)  !< PCA vectors
        logical,           optional, intent(in)    :: l_mask(:,:,:) !< true for pixels to extract
        logical,           optional, intent(in)    :: transp        !< pixel-wise learning
        real, allocatable, optional, intent(inout) :: avg_pix(:)    !< pixel average
        logical,           optional, intent(in)    :: do_fft
        logical, allocatable :: ll_mask(:,:,:)
        logical              :: ttransp, l_fft
        integer              :: n, i, ldim(3), ldim_mask(3)
        n       = size(imgs)
        l_fft   = .false.
        ttransp = .false.
        if( present(transp) ) ttransp = transp
        if( present(do_fft) ) l_fft = do_fft
        if( l_fft )then
            do i = 1, n
                call imgs(i)%fft
            enddo
            ldim = imgs(1)%get_array_shape()
            D    = product(ldim)
        else
            ldim = imgs(1)%get_ldim()
            if( present(l_mask) )then
                ldim_mask = [size(l_mask, dim=1),size(l_mask, dim=2),size(l_mask, dim=3)]
                if( .not. all(ldim_mask .eq. ldim) )then
                    write(logfhandle,*) 'ERROR! nonconforming matrix sizes'
                    write(logfhandle,*) 'ldim of image: ', ldim
                    write(logfhandle,*) 'ldim of mask : ', ldim_mask
                    THROW_HARD('make_pcavec_stack')
                endif
                allocate(ll_mask(ldim_mask(1),ldim_mask(2),ldim_mask(3)), source=l_mask)
            else
                ldim_mask = ldim
                allocate(ll_mask(ldim_mask(1),ldim_mask(2),ldim_mask(3)), source=.true.)
            endif
            D = count(ll_mask)
        endif
        if( allocated(pcavecs) ) deallocate(pcavecs)
        if( allocated(avg)     ) deallocate(avg)
        allocate(pcavecs(D,n), source=0.)
        do i = 1, n
            pcavecs(:,i) = imgs(i)%serialize(ll_mask)
        end do
        if( ttransp )then
            if( present(avg_pix) ) avg_pix = sum(pcavecs, dim=2) / real(n)
            pcavecs = transpose(pcavecs)
            avg = sum(pcavecs, dim=2) / real(D)
            do i = 1, D
                pcavecs(:,i) = pcavecs(:,i) - avg
            end do
        else
            avg = sum(pcavecs, dim=2) / real(n)
            if( present(avg_pix) ) avg_pix = avg
            do i = 1, n
                pcavecs(:,i) = pcavecs(:,i) - avg
            end do
        endif
        if( l_fft )then
            do i = 1, n
                call imgs(i)%ifft
            enddo
        endif
    end subroutine make_pcavecs

    subroutine make_pcavol( vol, D, avg, pcavec )
        class(image),      intent(in)    :: vol           !< vol to serialize
        integer,           intent(out)   :: D             !< vector dimension
        real, allocatable, intent(inout) :: avg(:)        !< average of pcavecs
        real, allocatable, intent(inout) :: pcavec(:,:)   !< PCA vector
        logical, allocatable :: ll_mask(:,:,:)
        real,    allocatable :: volvec(:)
        integer              :: i, j, ldim(3)
        ldim = vol%get_ldim()
        allocate(ll_mask(ldim(1),ldim(2),ldim(3)), source=.true.)
        D = count(ll_mask)
        if( allocated(pcavec) ) deallocate(pcavec)
        allocate(pcavec(D,D), volvec(D), avg(D), source=0.)
        volvec = vol%serialize(ll_mask)
        do i = 1, D
            do j = 1, D
                pcavec(i,j) = abs(volvec(i) - volvec(j))
            enddo
        enddo
        avg = sum(pcavec, dim=2) / real(D)
        do i = 1, D
            pcavec(:,i) = pcavec(:,i) - avg
        enddo
    end subroutine make_pcavol

    ! write tiled jpeg of mrc file
    subroutine mrc2jpeg_tiled(mrcfile, outfile, scale, ntiles, msk, n_xtiles, n_ytiles)
        character(len=*),               intent(in)  :: mrcfile, outfile
        real,    optional,              intent(out) :: scale
        integer, optional,              intent(out) :: ntiles, n_xtiles, n_ytiles
        logical, optional, allocatable, intent(in)  :: msk(:)
        type(image)    :: img, img_pad, img_jpeg
        type(stack_io) :: stkio_r
        integer        :: ldim_stk(3)
        integer        :: ncls_here, xtiles, ytiles, icls, ix, iy, l_ntiles
        real           :: smpd
        smpd = 1.0
        if(.not. file_exists(trim(mrcfile))) return
        call find_ldim_nptcls(trim(mrcfile), ldim_stk, ncls_here)
        xtiles = floor(sqrt(real(ncls_here)))
        ytiles = ceiling(real(ncls_here) / real(xtiles))
        call img%new([ldim_stk(1), ldim_stk(1), 1], smpd)
        call img_pad%new([JPEG_DIM, JPEG_DIM, 1], smpd)
        call img_jpeg%new([xtiles * JPEG_DIM, ytiles * JPEG_DIM, 1], smpd)
        call stkio_r%open(trim(mrcfile), smpd, 'read', bufsz=ncls_here)
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
        call img_jpeg%write_jpg(trim(outfile) // ".tmp")
        call simple_rename(trim(outfile) // ".tmp", trim(outfile), overwrite=.true.)
        if(present(scale))  scale  = real(JPEG_DIM) / ldim_stk(1)
        if(present(ntiles)) ntiles = l_ntiles
        if(present(n_xtiles)) n_xtiles = xtiles
        if(present(n_ytiles)) n_ytiles = ytiles
        call img%kill()
        call img_pad%kill()
        call img_jpeg%kill()
    end subroutine mrc2jpeg_tiled

end module simple_imgproc