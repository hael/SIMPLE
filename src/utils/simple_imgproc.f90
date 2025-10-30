! stack image processing routines for SPIDER/MRC files and other useful img jiffys
module simple_imgproc
include 'simple_lib.f08'
use simple_cmdline,        only: cmdline
use simple_default_clines, only: set_automask2D_defaults
use simple_image,          only: image
use simple_masker,         only: automask2D
use simple_parameters,     only: parameters
use simple_stack_io,       only: stack_io
implicit none
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

    subroutine mrc2mskdiam( mrcfile, smpd, mskdiam, box_for_pick, box_for_extract )
        character(len=*), intent(in)  :: mrcfile
        real,             intent(in)  :: smpd
        real,             intent(out) :: mskdiam
        integer,          intent(out) :: box_for_pick, box_for_extract
        type(parameters)         :: params
        type(stack_io)           :: stkio_r
        type(image)              :: ref2D, ref2D_clip
        type(cmdline)            :: cline
        type(image), allocatable :: masks(:)
        real,        allocatable :: diams(:), shifts(:,:)
        real,    parameter :: MSKDIAM2LP = 0.15, lP_LB = 30., LP_UB = 15.
        integer, parameter :: NREFS=100
        real    :: diam_max, maxdiam
        integer :: ldim(3), ncavgs, icavg
        call cline%set('pickrefs', trim(mrcfile))
        call cline%set('smpd',     smpd)
        ! set defaults
        call set_automask2D_defaults(cline)
        ! parse parameters
        call params%new(cline)
        ! read selected cavgs
        call find_ldim_nptcls(params%pickrefs, ldim, ncavgs)
        ldim(3) = 1
        params%msk = real(ldim(1)/2) - COSMSKHALFWIDTH ! for automasking
        ! read
        allocate( masks(ncavgs) )
        call stkio_r%open(params%pickrefs, params%smpd, 'read', bufsz=ncavgs)
        do icavg = 1, ncavgs
            call stkio_r%read(icavg, masks(icavg))
        end do
        call stkio_r%close
        call automask2D(masks, params%ngrow, nint(params%winsz), params%edge, diams, shifts)
        box_for_pick    = min(round2even(diam_max / params%smpd + 2. * COSMSKHALFWIDTH), ldim(1))
        mskdiam         = params%smpd * box_for_pick
        maxdiam         = mskdiam + mskdiam * BOX_EXP_FAC
        box_for_extract = find_larger_magic_box(round2even(maxdiam / params%smpd))
    end subroutine mrc2mskdiam

end module simple_imgproc