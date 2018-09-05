! operations on micrographs
module simple_micops
include 'simple_lib.f08'
use simple_image, only: image
implicit none

public :: read_micrograph, shrink_micrograph, set_box, extract_boxes2file
public :: extract_windows2file, generate_sampling_coordinates
private
#include "simple_local_flags.inc"

logical, parameter :: DEBUG_HERE = .true.
type(image) :: micrograph, micrograph_shrunken, micrograph_denoised
integer     :: ldim(3), ldim_shrunken(3), box, box_shrunken, nx, ny, border_offset
real        :: shrink_factor

contains
    subroutine read_micrograph( micfname, smpd, save_copy )
        character(len=*),  intent(in) :: micfname
        real,              intent(in) :: smpd
        logical, optional, intent(in) :: save_copy
        logical :: ssave_copy
        integer :: nframes
        ssave_copy = .false.
        if(present(save_copy)) ssave_copy = save_copy
        call find_ldim_nptcls(trim(micfname), ldim, nframes)
        if( nframes > 1 )then
            write(*,*) '# frames: ', nframes
            THROW_HARD('simple_micops only intended for nframes=1 micrographs; read_micrograph')
        endif
        call micrograph%new(ldim, smpd)
        call micrograph%read(trim(micfname))
        if(ssave_copy) call micrograph%write('Original_micrograph.mrc')
        if( DEBUG ) print *, 'debug(micops); read micrograph'
    end subroutine read_micrograph

    subroutine shrink_micrograph( shrink_fact, ldim_out, smpd_out )
        real,              intent(in)  :: shrink_fact
        integer, optional, intent(out) :: ldim_out(3)
        real,    optional, intent(out) :: smpd_out
        shrink_factor    = shrink_fact
        ldim_shrunken(1) = round2even(real(ldim(1))/shrink_factor)
        ldim_shrunken(2) = round2even(real(ldim(2))/shrink_factor)
        ldim_shrunken(3) = 1
        call micrograph_shrunken%new(ldim_shrunken, micrograph%get_smpd() * shrink_factor)
        call micrograph%fft()
        call micrograph_shrunken%set_ft(.true.)
        call micrograph%clip(micrograph_shrunken)
        if( present(ldim_out) ) ldim_out = ldim_shrunken
        if( present(smpd_out) ) smpd_out = micrograph%get_smpd() * shrink_factor
        if( DEBUG_HERE ) print *, 'DEBUG_HERE(micops); did shrink micrograph to logical dimension: ', ldim_shrunken
    end subroutine shrink_micrograph

    subroutine set_box( box_in, box_out, snr )
        integer,           intent(in)  :: box_in
        integer, optional, intent(out) :: box_out
        real,    optional, intent(in)  :: snr
        box = box_in
        box_shrunken = round2even(real(box)/shrink_factor)
        !high-pass filter shrunken micrograph according to box_shrunken
        call micrograph_shrunken%bp(box_shrunken * micrograph_shrunken%get_smpd(), 0.)!, width=real(box_shrunken/2.)) Chiara
        ! return filtered micrograph in real-space
        call micrograph_shrunken%ifft()
        if(present(snr)) call micrograph_shrunken%add_gauran(snr)
        if( DEBUG_HERE ) call micrograph_shrunken%write('shrunken_hpassfiltered.mrc')
        ! loop dimensions for target extraction will be 0:nx and 0:ny
        nx = ldim_shrunken(1) - box_shrunken
        ny = ldim_shrunken(2) - box_shrunken
        if( present(box_out) ) box_out = box_shrunken
        if( DEBUG_HERE ) print *, 'DEBUG_HERE(micops); did set box_shrunken to: ', box_shrunken
    end subroutine set_box

    subroutine extract_boxes2file( offset, outstk, n_images, boffset, coord )
        integer,           intent(in)               :: offset
        character(len=*),  intent(in)               :: outstk
        integer,           intent(out)              :: n_images
        integer, optional, intent(in)               :: boffset
        integer, optional, allocatable, intent(out) :: coord(:,:)
        type(image) :: imgwin
        integer     :: xind, yind, cnt, toc(2)
        logical     :: outside
        if( present(boffset) )then
          border_offset = boffset
        else
          border_offset = 0
        endif
        call imgwin%new([box_shrunken,box_shrunken,1], micrograph_shrunken%get_smpd())
        cnt = 0
        do xind=border_offset,nx-border_offset,offset
            do yind=border_offset,ny-border_offset,offset
                call micrograph_shrunken%window_slim([xind,yind], box_shrunken, imgwin, outside)
                if( .not. outside )then
                    cnt = cnt + 1
                    call imgwin%write(outstk, cnt)
                endif
            end do
        end do
        n_images = cnt
        if(present(coord)) then
              allocate(coord(cnt,2))
              cnt = 0
              do xind=border_offset,nx-border_offset,offset
                  do yind=border_offset,ny-border_offset,offset
                    outside = .false.
                    toc = [xind,yind] + box_shrunken
                    if(toc(1) > ldim_shrunken(1) .or. toc(2) > ldim_shrunken(2))  outside = .true.
                      if( .not. outside )then
                          cnt = cnt + 1
                          coord(cnt,:) = [xind,yind]
                      endif
                  end do
              end do

        endif
        if( DEBUG_HERE ) print *, 'DEBUG_HERE(micops); wrote # images to stack: ', n_images
    end subroutine extract_boxes2file

    ! This subroutine does exactly what 'extract_boxes2file' in simple_micops does, but it works on
    ! whatever input image (not necessary a mic)
    subroutine extract_windows2file( micrograph_shrunken, offset, outstk, n_images, boffset, coord )
        type(image), intent(inout) :: micrograph_shrunken
        integer,           intent(in)               :: offset
        character(len=*),  intent(in)               :: outstk
        integer,           intent(out)              :: n_images
        integer, optional, intent(in)               :: boffset
        integer, optional, allocatable, intent(out) :: coord(:,:)
        type(image) :: imgwin
        integer     :: xind, yind, cnt, toc(2)
        logical     :: outside
        if( present(boffset) )then
          border_offset = boffset
        else
          border_offset = 0
        endif
        call imgwin%new([box_shrunken,box_shrunken,1], micrograph_shrunken%get_smpd())
        cnt = 0
        do xind=border_offset,nx-border_offset,offset
            do yind=border_offset,ny-border_offset,offset
                call micrograph_shrunken%window_slim([xind,yind], box_shrunken, imgwin, outside)
                if( .not. outside )then
                    cnt = cnt + 1
                    call imgwin%write(outstk, cnt)
                endif
            end do
        end do
        n_images = cnt
        if(present(coord)) then
              allocate(coord(cnt,2))
              cnt = 0
              do xind=border_offset,nx-border_offset,offset
                  do yind=border_offset,ny-border_offset,offset
                    outside = .false.
                    toc = [xind,yind] + box_shrunken
                    if(toc(1) > ldim_shrunken(1) .or. toc(2) > ldim_shrunken(2))  outside = .true.
                      if( .not. outside )then
                          cnt = cnt + 1
                          coord(cnt,:) = [xind,yind]
                      endif
                  end do
              end do

        endif
        if( DEBUG_HERE ) print *, 'DEBUG_HERE(micops); wrote # images to stack: ', n_images
    end subroutine extract_windows2file

    function generate_sampling_coordinates( offset ) result( coords )
        integer, intent(in)  :: offset
        integer, allocatable :: coords(:,:)
        integer :: xind, yind, cnt
        ! generate sampling coordinates
        cnt = 0
        do xind=border_offset,nx-border_offset,offset
            do yind=border_offset,ny-border_offset,offset
                cnt = cnt + 1
            end do
        end do
        allocate(coords(cnt,2))
        cnt = 0
        do xind=border_offset,nx-border_offset,offset
            do yind=border_offset,ny-border_offset,offset
                cnt = cnt + 1
                coords(cnt,:) = [xind,yind]
            end do
        end do
    end function generate_sampling_coordinates

end module simple_micops
