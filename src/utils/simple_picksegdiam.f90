module simple_picksegdiam
include 'simple_lib.f08'
use simple_image,    only: image
use simple_binimage, only: binimage
use simple_micproc
use simple_linked_list
use simple_syslib
implicit none

public :: picksegdiam
private
#include "simple_local_flags.inc"

real, parameter :: SMPD_SHRINK1 = 4.0
real, parameter :: FRAC_FG      = 0.17

! instance
type picksegdiam
    private
    real              :: moldiam_max
    type(linked_list) :: diameters, xpos, ypos
contains
    generic            :: pick => pick_1, pick_2 
    procedure, private :: pick_1, pick_2
    procedure          :: get_nboxes
    procedure          :: get_diameters
    procedure          :: write_pos_and_diams
    procedure          :: kill
end type picksegdiam

contains

    subroutine pick_1( self, micname, smpd, moldiam_max, pcontrast, denfname, topofname, binfname )
        class(picksegdiam),         intent(inout) :: self
        character(len=*),           intent(in)    :: micname !< micrograph file name
        real,                       intent(in)    :: smpd, moldiam_max
        character(len=*),           intent(in)    :: pcontrast
        character(len=*), optional, intent(in)    :: denfname, topofname, binfname
        character(len=:), allocatable :: fname, output_dir
        integer,          allocatable :: cc_imat(:,:,:), cc_imat_copy(:,:,:)
        logical,          allocatable :: picking_mask(:,:)
        type(image)    :: mic_raw, mic_den
        type(binimage) :: mic_shrink, mic_bin, img_cc
        real           :: rpos(2), diam, scale
        integer        :: ldim_raw(3), ldim(3), pos(2), icc, nccs, nmasked
        logical        :: l_empty
        call self%kill
        self%moldiam_max = moldiam_max
        scale = smpd / SMPD_SHRINK1
        call read_mic_subtr_backgr_shrink(micname, smpd, scale, pcontrast, mic_raw, mic_shrink, l_empty)
        if( l_empty ) return
        ldim_raw = mic_raw%get_ldim()
        ldim     = mic_shrink%get_ldim()
        ! Prep segmentation
        call flag_amorphous_carbon(mic_shrink, picking_mask)
        nmasked = count(.not.picking_mask)
        if( real(nmasked) > 0.98 * real(product(ldim)) ) return
        call cascade_filter_biomol( mic_shrink, mic_den )
        if( present(denfname) ) call mic_den%write(denfname)
        if( present(topofname) ) call mic_shrink%write(topofname)
        call binarize_mic_den(mic_shrink, FRAC_FG, mic_bin)
        if( nmasked > 0 ) call mic_bin%apply_mask(picking_mask)
        ! identify connected components
        call mic_bin%find_ccs(img_cc)
        call img_cc%get_nccs(nccs)
        ! gather size info
        call img_cc%get_imat(cc_imat)
        call img_cc%get_imat(cc_imat_copy)
        do icc = 1, nccs
            call img_cc%diameter_cc(icc, diam)
            if( diam + 2. * SMPD_SHRINK1 > self%moldiam_max .or. diam < SMPD_SHRINK1 * 3. ) then
                ! remove connected component
                where ( cc_imat == icc ) cc_imat_copy = 0
            else
                ! stash diameter & box coordinate
                call self%diameters%push_back(diam)
                call img_cc%masscen_cc(icc, rpos)
                pos = nint(rpos/scale) + ldim_raw(1:2)/2 ! base 0
                call self%xpos%push_back(pos(1))
                call self%ypos%push_back(pos(2))
            endif
        end do
        ! binarize back
        cc_imat = cc_imat_copy
        where( cc_imat_copy > 0 )
            cc_imat = 1
        elsewhere
            cc_imat = 0
        endwhere
        call mic_bin%set_imat(cc_imat)
        if( present(binfname) ) call mic_bin%write(binfname)
        ! destruct
        call mic_shrink%kill
        call mic_den%kill
        call mic_bin%kill_bimg
        call img_cc%kill_bimg
        deallocate(cc_imat, cc_imat_copy)
    end subroutine pick_1

    subroutine pick_2( self, ldim_raw, smpd_raw, binmicname, diam_fromto )
        class(picksegdiam), intent(inout) :: self
        integer,            intent(in)    :: ldim_raw(3)
        real,               intent(in)    :: smpd_raw
        character(len=*),   intent(in)    :: binmicname !< binarized micrograph file name, smpd is SMPD_SHRINK1
        real,               intent(in)    :: diam_fromto(2)
        integer,  allocatable :: cc_imat(:,:,:), cc_imat_copy(:,:,:)
        real,     allocatable :: masscens(:,:)
        type(binimage) :: mic_bin, img_cc
        real           :: rpos(2), diam, diam_adj, scale
        integer        :: pos(2), icc, nccs, nboxes, ldim(3), nframes
        logical        :: l_empty
        call self%kill
        scale = smpd_raw / SMPD_SHRINK1
        ! read binary
        call find_ldim_nptcls(trim(binmicname), ldim, nframes)
        if( ldim(3) /= 1 .or. nframes /= 1 ) THROW_HARD('Only for 2D images')
        call mic_bin%new_bimg(ldim, SMPD_SHRINK1)
        call mic_bin%read_bimg(trim(binmicname), 1)
        call img_cc%new_bimg(ldim, SMPD_SHRINK1)
        ! identify connected components
        call mic_bin%find_ccs(img_cc, update_imat=.true.)
        call img_cc%get_nccs(nccs)
        ! gather size info
        call img_cc%get_imat(cc_imat)
        call img_cc%get_imat(cc_imat_copy)
        do icc = 1, nccs
            call img_cc%diameter_cc(icc, diam)
            diam_adj = diam + 2. * SMPD_SHRINK1 ! becasue of the 2X erosion
            if( diam_adj >= diam_fromto(1) .and. diam_adj <= diam_fromto(2) ) then
                ! stash diameter & box coordinate
                call self%diameters%push_back(diam)
                call img_cc%masscen_cc(icc, rpos)
                pos = nint(rpos/scale) + ldim_raw(1:2)/2 ! base 0
                call self%xpos%push_back(pos(1))
                call self%ypos%push_back(pos(2))
            else
                ! remove connected component
                where ( cc_imat == icc ) cc_imat_copy = 0
            endif
        end do
        ! destruct
        call mic_bin%kill_bimg
        call img_cc%kill_bimg
        deallocate(cc_imat, cc_imat_copy)
    end subroutine pick_2

    ! Getters

    integer function get_nboxes( self )
        class(picksegdiam), intent(in) :: self
        get_nboxes = self%diameters%size()
    end function get_nboxes

    subroutine get_diameters( self, arr )
        class(picksegdiam), intent(in) :: self
        real, allocatable, intent(inout) :: arr(:)
        type(list_iterator)   :: diams_iter
        class(*), allocatable :: any
        real    :: areal
        integer :: i, n
        if( allocated(arr) ) deallocate(arr)
        n = self%diameters%size()
        if( n == 0 ) return
        allocate(arr(n))
        i = 0
        diams_iter = self%diameters%begin()
        do while (diams_iter%has_next())
            i = i + 1
            call diams_iter%next(any)
            select type(any)
                type is (real(kind(areal)))
                    arr(i) = any
            end select
        end do
    end subroutine get_diameters

    ! for writing coordinates & diameters
    subroutine write_pos_and_diams( self, fname, nptcls, box )
        class(picksegdiam), intent(in)    :: self
        character(len=*),   intent(in)    :: fname
        integer,            intent(inout) :: nptcls
        integer, optional,  intent(in)    :: box
        type(list_iterator)   :: xpos_iter, ypos_iter
        real,     allocatable :: diams(:)
        class(*), allocatable :: any
        real    :: areal, diam
        integer :: i, funit, iostat, x,y, box_here
        nptcls = self%diameters%size()
        if( nptcls == 0 ) return
        if( present(box) )then
            box_here = box
        else
            box_here = round2even(BOX_EXP_FACTOR_DEFAULT*self%moldiam_max)
        endif
        allocate(diams(nptcls),source=0.)
        call self%get_diameters(diams)
        call fopen(funit, status='REPLACE', action='WRITE', file=trim(adjustl(fname)), iostat=iostat)
        xpos_iter = self%xpos%begin()
        ypos_iter = self%ypos%begin()
        i = 0
        do while (xpos_iter%has_next())
            i = i + 1
            call xpos_iter%next(any)
            select type(any)
                type is (integer(kind(i)))
                    x = any
            end select
            call ypos_iter%next(any)
            select type(any)
                type is (integer(kind(i)))
                    y = any
            end select
            x = x - box_here/2
            y = y - box_here/2
            write(funit,'(4I7,2F8.1)') x, y, box_here, box_here, diams(i)
        end do
        call fclose(funit)
        deallocate(diams,any)
    end subroutine write_pos_and_diams

    subroutine kill( self )
        class(picksegdiam), intent(inout) :: self
        self%moldiam_max = 0.
        call self%diameters%kill
        call self%xpos%kill
        call self%ypos%kill
    end subroutine kill

end module simple_picksegdiam
