! stack image processing routines for SPIDER/MRC files
module simple_imgproc
include 'simple_lib.f08'
use simple_image, only: image
implicit none

public :: make_pcavecs
private

#include "simple_local_flags.inc"

contains

    subroutine make_pcavecs( imgs, D, avg, pcavecs, l_mask )
        class(image),      intent(in)    :: imgs(:)       !< images to serialize
        integer,           intent(out)   :: D             !< vector dimension
        real, allocatable, intent(inout) :: avg(:)        !< average of pcavecs 
        real, allocatable, intent(inout) :: pcavecs(:,:)  !< PCA vectors
        logical, optional, intent(in)    :: l_mask(:,:,:) !< true for pixels to extract
        type(image)          :: img
        logical, allocatable :: ll_mask(:,:,:)
        integer              :: n, i, ldim(3), ldim_mask(3)
        real                 :: smpd
        ldim = imgs(1)%get_ldim()
        n    = size(imgs)
        smpd = imgs(1)%get_smpd()
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
        call img%new(ldim,smpd)
        D = count(ll_mask)
        if( allocated(pcavecs) ) deallocate(pcavecs)
        allocate(pcavecs(n,D), source=0.)
        if( allocated(avg) ) deallocate(avg)
        allocate(avg(D), source=0.)
        do i = 1, n
            pcavecs(i,:) = imgs(i)%serialize(ll_mask)
            avg(:) = avg(:) + pcavecs(i,:)
        end do
        avg = avg / real(n)
        do i=1,n
            pcavecs(i,:) = pcavecs(i,:) - avg(:)
        end do
    end subroutine make_pcavecs

end module simple_imgproc