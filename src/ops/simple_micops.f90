! operations on micrographs
module simple_micops
include 'simple_lib.f08'
use simple_image, only: image
implicit none

public :: read_micrograph, shrink_micrograph
private

type(image) :: micrograph, micrograph_shrunken
integer     :: ldim(3), ldim_shrunken(3), box, box_shrunken
real        :: shrink_factor


contains

    subroutine read_micrograph( micfname, smpd )
        character(len=*), intent(in) :: micfname
        real,             intent(in) :: smpd
        integer :: nframes
        call find_ldim_nptcls(trim(micfname), ldim, nframes)
        if( nframes > 1 )then
            write(*,*) '# frames: ', nframes
            write(*,*) 'ERROR! simple_micops only intended for nframes=1 micrographs'
            stop 'simple_micops :: read_micrograph'
        endif
        call micrograph%new(ldim, smpd)
        call micrograph%read(trim(micfname))
    end subroutine read_micrograph

    subroutine shrink_micrograph( shrink_fact, ldim_out )
        real,              intent(in)  :: shrink_fact
        integer, optional, intent(out) :: ldim_out(3)
        shrink_factor    = shrink_fact
        ldim_shrunken(1) = round2even(real(ldim(1))/shrink_factor)
        ldim_shrunken(2) = round2even(real(ldim(2))/shrink_factor)
        ldim_shrunken(3) = 1
        call micrograph_shrunken%new(ldim_shrunken, micrograph%get_smpd() * shrink_factor)
        call micrograph%fft()
        call micrograph_shrunken%set_ft(.true.)
        call micrograph%clip(micrograph_shrunken)
        if( present(ldim_out) ) ldim_out = ldim_shrunken
    end subroutine shrink_micrograph


    subroutine set_box( box_in )
        integer, intent(in) :: box_in
        box = box_in
        box_shrunken = round2even(real(box)/shrink_factor)
        ! high-pass filter shrunken micrograph according to box_shrunken
        call micrograph_shrunken%bp(box_shrunken * micrograph_shrunken%get_smpd(), 0.)
    end subroutine set_box


    ! subroutine extract_boxes2file

    ! end subroutine extract_boxes2file

    ! subroutine extract_box_features( box, offset )



    ! end subroutine extract_box_features







end module simple_micops

