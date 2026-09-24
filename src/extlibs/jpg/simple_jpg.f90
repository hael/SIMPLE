!! Michael Eager Feb 2018

!! Grayscale Bit representation in JPEG is limited to 24 bits
!! (better than 8-bit PNG)

module simple_jpg
use simple_core_module_api
implicit none

public :: jpg_img, write_rgb_jpeg
private
#include "simple_local_flags.inc"

integer,         parameter :: max_colors = 256
integer(kind=4), parameter :: boz_00ff = INT(z'000000ff',kind=4)
integer(kind=4), parameter :: boz_ffff = INT(z'00ffffff',kind=4)

type jpg_img
    private
    integer :: width      =  0
    integer :: height     =  0
    integer :: quality    =  90
    integer :: colorspace =  1
contains
    procedure, private :: save_jpeg_r4
    procedure, private :: save_jpeg_r4_3D
    generic            :: writeJpg => save_jpeg_r4, save_jpeg_r4_3D
end type jpg_img

interface

    integer function stbi_write_jpg (file_name, w, h, comp, data, quality ) bind ( c, name="stbi_write_jpg" )
        use,intrinsic :: iso_c_binding
        implicit none
        character(c_char),dimension(*),intent(in)    :: file_name
        integer(c_int), intent(in), VALUE            :: w
        integer(c_int), intent(in), VALUE            :: h
        integer(c_int), intent(in), VALUE            :: comp     ! Each pixel contains 'comp' channels of data stored interleaved with 8-bits
        !   per channel, in the following order: 1=Y, 2=YA, 3=RGB, 4=RGBA. (Y is  monochrome color.)
        type (C_PTR), VALUE             :: data ! (const void *)
        integer(c_int), intent(in), value            :: quality  ! 1 to 100. Higher quality looks better but results in a bigger image.
    end function stbi_write_jpg

end interface

contains

    !> True-colour JPEG from an interleaved RGB raster rgb(3,w,h) with channel values in [0,1]
    !! (row 1 = top of the image). The jpg_img colorspec=3 path packs ONE scalar into 24 bits, which
    !! is a false-colour ramp, not RGB; this writes the three channels as the stb writer expects them.
    subroutine write_rgb_jpeg( fname, rgb, quality )
        use, intrinsic :: iso_c_binding, only: c_ptr, c_loc, c_null_char
        character(len=*),  intent(in) :: fname
        real,              intent(in) :: rgb(:,:,:)
        integer, optional, intent(in) :: quality
        integer(1), pointer           :: buf(:) => NULL()
        character(len=:), allocatable :: fstr
        type(c_ptr) :: img
        integer     :: w, h, i, j, c, idx, q, status, v
        w = size(rgb,2); h = size(rgb,3)
        if( size(rgb,1) /= 3 .or. w < 1 .or. h < 1 ) THROW_HARD('write_rgb_jpeg: rgb must be (3,w,h)')
        q = 90; if( present(quality) ) q = quality
        allocate(buf(3*w*h)); allocate(fstr, source=trim(fname)//c_null_char)
        do j = 1, h
            do i = 1, w
                idx = 3*((j-1)*w + (i-1))
                do c = 1, 3
                    v = nint(255.0 * min(1.0, max(0.0, rgb(c,i,j))))
                    if( v > 127 ) v = v - 256          ! two's complement into integer(1)
                    buf(idx+c) = int(v, kind=1)
                end do
            end do
        end do
        img    = c_loc(buf)
        status = stbi_write_jpg(fstr, w, h, 3, img, q)
        if( status == 0 ) THROW_HARD('write_rgb_jpeg: stbi_write_jpg failed for '//trim(fname))
        deallocate(buf, fstr)
    end subroutine write_rgb_jpeg

    function save_jpeg_r4_3D (self, fname, in_buffer, quality, colorspec) result(status)
        class(jpg_img),   intent(inout)        :: self
        character(len=*), intent(in)           :: fname
        real,             intent(in)           :: in_buffer(:,:,:)
        integer,          intent(in), optional :: quality
        integer,          intent(in), optional :: colorspec
        character(len=:), allocatable :: fname_here
        integer               ::  w,h,c,slice
        integer               ::  status
        character(len=STDLEN) ::  fstr
        status = 1
        fname_here = trim(adjustl(fname))
        c=4
        w = size(in_buffer,1)
        h= size(in_buffer,2)
        if(w == 0 .or. h == 0) return
        do slice=1, size(in_buffer,3)
            write(fstr, '(a,i4.4,a)') fname_here(1:(len_trim(fname_here)-c))//'_',slice,'.jpg'//c_null_char
            status = self%save_jpeg_r4 (fstr, in_buffer(:,:,slice), quality, colorspec)
        end do
    end function save_jpeg_r4_3D

    function save_jpeg_r4(self, fname, in_buffer, quality, colorspec) result(status)
        class(jpg_img),    intent(inout) :: self
        character(len=*),  intent(in)    :: fname
        real,              intent(in)    :: in_buffer(:,:)
        integer, optional, intent(in)    :: quality, colorspec
        type(c_ptr)                      :: img
        real                          :: lo, hi
        integer                       :: status, w, h,c, i, j, idx
        integer(c_int)                :: pixel
        character(len=:), allocatable :: fstr
        integer(1),       pointer     :: img_buffer(:) => NULL()
        status = 1
        c      = 1
        self%width  = 0
        self%height = 0
        w = size(in_buffer,1)
        h = size(in_buffer,2)
        if(w == 0 .or. h == 0) return
        if(present(quality)) self%quality = quality
        if(present(colorspec)) self%colorspace = colorspec
        allocate(fstr, source=trim(fname)//c_null_char)
        lo = minval(in_buffer)
        hi = maxval(in_buffer)
        allocate(img_buffer(w*h*3))
        do j=0,h-1
            do i=0,w-1
                if  (self%colorspace == 3) then
                    c=3
                    pixel = NINT( real(boz_ffff) * (in_buffer(i+1,j+1)-lo)/(hi-lo),kind=4)
                    idx = i*c + (j*w*c) + 1
                    img_buffer(idx)     = INT( ISHFT( pixel , -16) ,kind=c_char)
                    img_buffer(idx + 1) = INT( IAND( ISHFT( pixel , -8_c_int) , boz_00ff) ,kind=c_char)
                    img_buffer(idx + 2) = INT( IAND( pixel , boz_00ff) ,kind=c_char)
                else
                    c=1
                    pixel =  INT( REAL( max_colors - 1)*REAL( (in_buffer(i+1,j+1)-lo)/REAL(hi - lo) ) ,kind=c_int)
                    pixel =  IAND( pixel , boz_ffff)
                    idx = i*c + (j*w*c) + 1
                    img_buffer(idx) = INT(pixel,kind=1)
                end if

            end do
        end do
        img = c_loc(img_buffer)
        self%width = w
        self%height = h
        status = stbi_write_jpg (fstr, self%width, self%height, self%colorspace, img, self%quality )
        if(status == 0 ) THROW_HARD('call to write_jpeg failed')
        status = 0
        deallocate(fstr,img_buffer)
    end function save_jpeg_r4

end module simple_jpg
