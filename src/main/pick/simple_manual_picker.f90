!@descr: manual box selection with lower-origin image coordinates, a pan/zoom viewport and EMAN box output
module simple_manual_picker
use, intrinsic :: ieee_arithmetic, only: ieee_is_finite, ieee_rint
use simple_defs,  only: dp
use simple_error, only: simple_exception
implicit none
private
public :: manual_picker
#include "simple_local_flags.inc"

type :: manual_picker
    private
    real(dp), allocatable :: centres(:,:)
    integer :: dims(2) = 0, box = 0, selected = 0
    real(dp) :: scale = 1._dp, origin(2) = 0._dp
  contains
    procedure :: new
    procedure :: kill
    procedure :: count_boxes
    procedure :: get_selected
    procedure :: get_box
    procedure :: set_box
    procedure :: get_centres
    procedure :: click
    procedure :: move_selected
    procedure :: delete_selected
    procedure :: fit
    procedure :: pan
    procedure :: zoom
    procedure :: image_to_screen
    procedure :: screen_to_image
    procedure :: save
end type manual_picker

contains

    subroutine new( self, dims, box )
        class(manual_picker), intent(inout) :: self
        integer, intent(in) :: dims(2), box
        call self%kill
        if( box < 1 .or. any(dims < box) ) THROW_HARD('picking box must fit inside both image dimensions')
        self%dims = dims
        self%box = box
        allocate(self%centres(2,0))
    end subroutine new

    subroutine kill( self )
        class(manual_picker), intent(inout) :: self
        if( allocated(self%centres) ) deallocate(self%centres)
        self%dims = 0
        self%box = 0
        self%selected = 0
        self%scale = 1._dp
        self%origin = 0._dp
    end subroutine kill

    integer function count_boxes( self )
        class(manual_picker), intent(in) :: self
        count_boxes = 0
        if( allocated(self%centres) ) count_boxes = size(self%centres,2)
    end function count_boxes

    integer function get_selected( self )
        class(manual_picker), intent(in) :: self
        get_selected = self%selected
    end function get_selected

    integer function get_box( self )
        class(manual_picker), intent(in) :: self
        get_box = self%box
    end function get_box

    subroutine set_box( self, box )
        class(manual_picker), intent(inout) :: self
        integer, intent(in) :: box
        integer :: i
        real(dp) :: half
        self%box = max(1, min(minval(self%dims), box))
        half = real(self%box,dp)/2._dp
        do i = 1,self%count_boxes()
            self%centres(:,i) = max(half, min(real(self%dims,dp)-half, self%centres(:,i)))
        end do
    end subroutine set_box

    function get_centres( self ) result( xy )
        class(manual_picker), intent(in) :: self
        real(dp), allocatable :: xy(:,:)
        xy = self%centres
    end function get_centres

    ! topmost box wins; a click on empty image adds one in insertion order
    subroutine click( self, xy )
        class(manual_picker), intent(inout) :: self
        real(dp), intent(in) :: xy(2)
        integer :: i, n
        real(dp), allocatable :: grown(:,:)
        if( .not. all(ieee_is_finite(xy)) ) return
        if( any(xy < 0._dp) .or. any(xy > real(self%dims,dp)) ) return
        do i = self%count_boxes(),1,-1
            if( all(abs(xy - self%centres(:,i)) <= real(self%box,dp)/2._dp) )then
                self%selected = i
                return
            endif
        end do
        n = self%count_boxes()
        allocate(grown(2,n+1))
        grown(:,:n) = self%centres
        grown(:,n+1) = xy
        call move_alloc(grown, self%centres)
        self%selected = self%count_boxes()
        call self%move_selected(xy)
    end subroutine click

    subroutine move_selected( self, xy )
        class(manual_picker), intent(inout) :: self
        real(dp), intent(in) :: xy(2)
        real(dp) :: half
        if( self%selected == 0 ) return
        if( .not. all(ieee_is_finite(xy)) ) return
        half = real(self%box,dp)/2._dp
        self%centres(:,self%selected) = max(half, min(real(self%dims,dp)-half, xy))
    end subroutine move_selected

    subroutine delete_selected( self )
        class(manual_picker), intent(inout) :: self
        integer :: i, n
        i = self%selected
        if( i == 0 ) return
        n = self%count_boxes()
        self%centres = reshape([reshape(self%centres(:,:i-1),[2*(i-1)]), &
            &reshape(self%centres(:,i+1:n),[2*(n-i)])], [2,n-1])
        self%selected = 0
    end subroutine delete_selected

    subroutine fit( self, width, height )
        class(manual_picker), intent(inout) :: self
        integer, intent(in) :: width, height
        self%scale = min(real(max(1,width),dp)/self%dims(1), real(max(1,height),dp)/self%dims(2))
        self%origin = (real([width,height],dp)-self%scale*real(self%dims,dp))/2._dp
    end subroutine fit

    subroutine pan( self, delta )
        class(manual_picker), intent(inout) :: self
        real(dp), intent(in) :: delta(2)
        self%origin = self%origin + delta
    end subroutine pan

    subroutine zoom( self, screen, factor )
        class(manual_picker), intent(inout) :: self
        real(dp), intent(in) :: screen(2), factor
        real(dp) :: ratio, next_scale
        if( .not. ieee_is_finite(factor) ) return
        if( factor <= 0._dp ) return
        next_scale = min(256._dp, max(0.001_dp, self%scale*factor))
        ratio = next_scale/self%scale
        self%origin = screen + ratio*(self%origin-screen)
        self%scale = next_scale
    end subroutine zoom

    function image_to_screen( self, xy ) result( screen )
        class(manual_picker), intent(in) :: self
        real(dp), intent(in) :: xy(2)
        real(dp) :: screen(2)
        screen = self%origin + self%scale*[xy(1), real(self%dims(2),dp)-xy(2)]
    end function image_to_screen

    function screen_to_image( self, screen ) result( xy )
        class(manual_picker), intent(in) :: self
        real(dp), intent(in) :: screen(2)
        real(dp) :: xy(2)
        xy = (screen-self%origin)/self%scale
        xy(2) = real(self%dims(2),dp)-xy(2)
    end function screen_to_image

    ! status is returned to the window so a failed save never closes it
    subroutine save( self, filename, status, message )
        class(manual_picker), intent(in) :: self
        character(len=*), intent(in) :: filename
        integer, intent(out) :: status
        character(len=*), intent(out) :: message
        integer :: unit, i, corner(2), close_status
        message = ''
        open(newunit=unit, file=filename, status='replace', action='write', iostat=status, iomsg=message)
        if( status /= 0 ) return
        do i = 1,self%count_boxes()
            corner = int(ieee_rint(self%centres(:,i)-real(self%box,dp)/2._dp))
            write(unit,'(i0,3(a,i0))', iostat=status, iomsg=message) &
                &corner(1), achar(9), corner(2), achar(9), self%box, achar(9), self%box
            if( status /= 0 ) exit
        end do
        close(unit, iostat=close_status)
        if( status == 0 )then
            status = close_status
            if( status /= 0 ) message = 'Could not close the output box file'
        endif
    end subroutine save

end module simple_manual_picker
