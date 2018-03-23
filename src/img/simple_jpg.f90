!! Michael Eager Feb 2018

!! Grayscale Bit representation in JPEG is limited to 24 bits
!! (better than 8-bit PNG)

module simple_jpg
include 'simple_lib.f08'
use,intrinsic                                        :: iso_c_binding
implicit none

public                                               :: jpg_img, test_jpg_export
private
#include "simple_local_flags.inc"
    !*** cptr should be have the same size as a c pointer
    !*** It doesn't matter whether it is an integer or a real
integer, parameter                                   :: cptr          = kind(5)
integer, parameter                                   :: numComponents = 3
!  integer, parameter                                   :: pixelsize     = selected_int_kind(3)
!   integer, parameter                                   :: longint       = selected_int_kind(9)
integer, public, parameter                           :: max_colors =  256
!buf_range(1)buf_range(1)im2 = gdImageScale(im, 1, 65535);
integer, parameter                                   :: GREYSCALE  =  1


type jpg_img
    private
    integer                                          :: width      =  0
    integer                                          :: height     =  0
    integer                                          :: quality    =  90
    integer                                          :: colorspace =  1
    integer(cptr)                                    :: ptr        =  0
    logical                                          :: fmode      = .false.
    real                                             :: gamma      = 0.707
    logical                                          :: normalised = .false.
    integer(1), dimension(:), pointer                :: img_buffer => NULL()
contains
    procedure :: destructor
    procedure :: constructor
    procedure                                        :: getWidth
    procedure                                        :: getHeight
    procedure                                        :: getQuality
    procedure                                        :: getColorspace
    procedure                                        :: setQuality
    procedure                                        :: setColorspace
    procedure                                        :: getGamma
    procedure                                        :: getNormalised
    procedure                                        :: setGamma
    procedure                                        :: setNormalised
    procedure :: montage
    procedure, private                               :: save_jpeg_r4_3D
    procedure, private                               :: save_jpeg_r4
    procedure, private                               :: save_jpeg_i4
    generic                                          :: writeJpg  => save_jpeg_r4_3D, save_jpeg_r4, save_jpeg_i4
    procedure, private                                        :: load_jpeg_r4
    procedure, private                                        :: load_jpeg_i4
    generic                                          :: loadJpg  => load_jpeg_r4, load_jpeg_i4

end type jpg_img

! interface jpg_img
!     ! Constructor declaration
!     module procedure constructor
! end interface jpg_img

interface
#ifdef HAVE_LIBJPEG
    subroutine setup_jpeg ( width, height, in, out ) bind ( c, name="setup_jpeg" )
        import
        integer(C_INT), intent(in), value            :: width
        integer(C_INT), intent(in), value            :: height
        integer(C_INT), dimension(*), intent(out)    :: in
        integer(C_INT), dimension(*), intent(out)    :: out
    end subroutine setup_jpeg
    subroutine read_jpeg (file_name, img, width, height, colorspec, status ) bind ( c, name="read_jpeg" )
        import
        character(c_char), dimension(*), intent(in)  :: file_name
        type(C_PTR), value            :: img
        integer(c_int), intent(inout)                :: width
        integer(c_int), intent(inout)                :: height
        integer(c_int), intent(inout)                :: colorspec
        integer(c_int), intent(inout)                :: status
    end subroutine read_jpeg

    subroutine write_jpeg (img, file_name, width, height, quality, colorspec, status ) bind ( c, name="write_jpeg" )
        import
        type (C_PTR), VALUE             :: img
        character(c_char),dimension(*),intent(in)    :: file_name
        integer(c_int), intent(in), VALUE            :: width
        integer(c_int), intent(in), VALUE            :: height
        integer(c_int), intent(in)                   :: quality
        integer(c_int), intent(in), VALUE            :: colorspec
        integer(c_int), intent(inout)                :: status
    end subroutine write_jpeg
#endif

    ! function fgltLoadTGA(FileName,width,height,components,eform,image)
    !     use, intrinsic :: iso_c_binding
    !     type(c_ptr), target :: fgltloadTGA
    !     character(len=*), intent(in) :: FileName
    !     integer(c_int), intent(out) :: width, height
    !     integer(c_int), intent(out) :: components, eform
    !     integer(c_char), dimension(:), allocatable, &
    !         intent(out), target :: image
    ! end function fgltLoadTGA
    !        STBIWDEF int stbi_write_jpg(char const *filename, int x, int y, int comp, const void  *data, int quality)
    integer function stbi_write_jpg (file_name, w, h, comp, data, quality ) bind ( c, name="stbi_write_jpg" )
        import
        character(c_char),dimension(*),intent(in)    :: file_name
        integer(c_int), intent(in), VALUE            :: w
        integer(c_int), intent(in), VALUE            :: h
        integer(c_int), intent(in), VALUE            :: comp     ! Each pixel contains 'comp' channels of data stored interleaved with 8-bits
        !   per channel, in the following order: 1=Y, 2=YA, 3=RGB, 4=RGBA. (Y is  monochrome color.)
        type (C_PTR), VALUE             :: data ! (const void *)
        integer(c_int), intent(in), value            :: quality  ! 1 to 100. Higher quality looks better but results in a bigger image.
    end function stbi_write_jpg
    integer function stbi_write_png (file_name, w, h, comp, data ) bind ( c, name="stbi_write_jpg" )
        import
        character(c_char),dimension(*),intent(in)    :: file_name
        integer(c_int), intent(in), VALUE            :: w
        integer(c_int), intent(in), VALUE            :: h
        integer(c_int), intent(in), VALUE            :: comp     ! Each pixel contains 'comp' channels of data stored interleaved with 8-bits
        !   per channel, in the following order: 1=Y, 2=YA, 3=RGB, 4=RGBA. (Y is  monochrome color.)
        type (C_PTR), VALUE             :: data ! (const void *)

    end function stbi_write_png
    integer function stbi_read_jpg (file_name, data, w, h, comp ) bind ( c, name="stbi_read_jpg" )
        use,intrinsic                                        :: iso_c_binding
        implicit none
        character(c_char),dimension(*),intent(in)    :: file_name
        type (C_PTR)             :: data ! (const void *)
        integer(c_int), intent(inout)            :: w
        integer(c_int), intent(inout)            :: h
        integer(c_int), intent(inout)            :: comp     ! Each pixel contains 'comp' channels of data stored interleaved with 8-bits
        !   per channel, in the following order: 1=Y, 2=YA, 3=RGB, 4=RGBA. (Y is  monochrome color.)


    end function stbi_read_jpg

end interface

contains

    subroutine constructor (self)
        class(jpg_img) :: self
#if !defined (PGI)
        DebugPrint 'In constructor address is: ', loc(self)
#endif
        if(associated(self%img_buffer)) nullify(self%img_buffer)
    end subroutine constructor

    integer function getWidth(self)
        class(jpg_img), intent(in)                       :: self
        getWidth =self%width
    end function getWidth
    integer function getHeight(self)
        class(jpg_img), intent(in)                       :: self
        getHeight = self%height
    end function getHeight
    integer function getQuality(self)
        class(jpg_img), intent(in)                       :: self
        getQuality =self%quality
    end function getQuality
    integer function getColorspace(self)
        class(jpg_img), intent(in)                       :: self
        getColorspace = self%colorspace
    end function getColorspace
    real function getGamma(self)
        class(jpg_img), intent(in)                       :: self
        getGamma =self%gamma
    end function getGamma
    logical function getNormalised(self)
        class(jpg_img), intent(in)                       :: self
        getNormalised = self%normalised
    end function getNormalised

    subroutine setWidth(self, width)
        class(jpg_img), intent(inout)                    :: self
        integer, intent(in)                              :: width
        self%width=width
    end subroutine setWidth
    subroutine setHeight(self,height)
        class(jpg_img), intent(inout)                    :: self
        integer, intent(in) :: height
        self%height=height
    end subroutine setHeight
    subroutine setQuality(self,quality)
        class(jpg_img), intent(inout)                    :: self
        integer, intent(in)                              :: quality
        self%quality=quality
    end subroutine setQuality
    subroutine setColorspace(self,colorspace)
        class(jpg_img), intent(inout)                    :: self
        integer, intent(in)                              :: colorspace
        self%colorspace=colorspace
    end subroutine setColorspace
    subroutine setGamma(self,gamma)
        class(jpg_img), intent(inout)                    :: self
        real,intent(in)                                  :: gamma
        self%gamma=gamma
    end subroutine setGamma
    subroutine setNormalised(self,normalised)
        class(jpg_img), intent(inout)                    :: self
        logical, intent(in)                              :: normalised
        self%normalised=normalised
    end subroutine setNormalised

    function montage (self, fname, in_buffer, quality, colorspec) result(status)
        class(jpg_img), intent(inout)    :: self
        character(len=*), intent(inout)  :: fname
        real,    intent(in)              :: in_buffer(:,:,:)
        integer, intent(in), optional    :: quality
        integer, intent(in), optional    :: colorspec
        type(c_ptr)                      :: img
        integer                          :: w,h,d,c, slice, dx, dy
        integer                          :: status, newshape(2)
        real, allocatable                :: new_buffer(:,:)
        character(len=STDLEN)            :: fstr
        status = 1
        VerbosePrint ' In save_jpeg_r4_3D '
        fname = trim(adjustl(fname))
        c=4
        !if(fname(len_trim(fname)) == c_null_char) c=5
        w = size(in_buffer,1)
        h= size(in_buffer,2)
        d = size(in_buffer,3)
        if(w == 0 .or. h == 0) return
        if (d == 1) then
            newshape = (/ w, h /)
            allocate(new_buffer(w,h), source=reshape(in_buffer, newshape))
        else
            dx = max( 1, ceiling( sqrt( real(d) ) ) )
            dy = min( dx*dx, floor( sqrt( real(d) ) ) )
            newshape = (/ w*dx, h*dy /)
            allocate(new_buffer(w*dx, h*dy), source=0.0)
            new_buffer = reshape(in_buffer, newshape)
        end if

        write(fstr, '(a)') fname//c_null_char
        status = self%save_jpeg_r4 (fstr, new_buffer, quality, colorspec)
        deallocate(new_buffer)
    end function montage


    function save_jpeg_r4_3D (self, fname, in_buffer, quality, colorspec) result(status)
        class(jpg_img),   intent(inout)        :: self
        character(len=*), intent(inout)        :: fname
        real,             intent(in)           :: in_buffer(:,:,:)
        integer,          intent(in), optional :: quality
        integer,          intent(in), optional :: colorspec
        type(c_ptr)           ::  img
        integer               ::  w,h,c,slice
        integer               ::  status
        character(len=STDLEN) ::  fstr
         status = 1
        VerbosePrint ' In save_jpeg_r4_3D '
        fname = trim(adjustl(fname))
        c=4
        !if(fname(len_trim(fname)) == c_null_char) c=5
        w = size(in_buffer,1)
        h= size(in_buffer,2)
        if(w == 0 .or. h == 0) return
        do slice=1, size(in_buffer,3)
            write(fstr, '(a,i4.4,a)') fname(1:(len_trim(fname)-c))//'_',slice,'.jpg'//c_null_char
            !            fstr = trim(adjustl(fstr))
            status = self%save_jpeg_r4 (fstr, in_buffer(:,:,slice), quality, colorspec)
            !    if(status /= 0 ) call simple_stop('save_jpeg_r4_3D call to save_jpeg_r4  failed' )
        end do
        !call write_jpeg(img, fname, w, h, quality_here, c, status)
        !if(status /= 0 ) call simple_stop('save_jpeg_r4 call to write_jpeg returned status==1' )
        !self%width = w
        !self%height = h
        !        deallocate(fstr)
    end function save_jpeg_r4_3D

    function save_jpeg_r4 (self, fname, in_buffer, quality, colorspec) result(status)
        class(jpg_img), intent(inout)   :: self
        character(len=*), intent(inout) :: fname
        real, intent(in)                :: in_buffer(:,:)
        integer, intent(in), optional   :: quality
        integer, intent(in), optional   :: colorspec
        type(c_ptr)                     :: img
        integer                         :: w, h,c, i, j
        real lo, hi
        integer                         :: status
        integer(c_int)                  :: pixel
        character(len=:), allocatable   :: fstr
        status = 1
        c=1
        ! verbose=.true.
        VerbosePrint '>>> In save_jpeg_r4 '
        self%width =0
        self%height = 0
        w = size(in_buffer,1)
        h = size(in_buffer,2)
        if(w == 0 .or. h == 0) return
        if(present(quality)) self%quality = quality
        if(present(colorspec)) self%colorspace = colorspec
        allocate(fstr, source=trim(fname)//c_null_char)
        VerbosePrint '>>> In save_jpeg_r4 input w x h ', w, h
        !        self%img_buffer = transfer(in_buffer,1_c_int8_t)
        lo = minval(in_buffer)
        hi = maxval(in_buffer)
        VerbosePrint '>>> In save_jpeg_r4 input lo hi ', lo, hi
        print *, in_buffer
        allocate(self%img_buffer(w*h*3))!, source=0_c_char)
        !  self%img_buffer = reshape(in_buffer, (/ w*h /) )
        do j=0,h-1
            do i=0,w-1
                if  (self%colorspace == 3) then
                    c=3
                    pixel = NINT( (2**24) * (in_buffer(i+1,j+1)-lo)/(hi-lo),kind=4)

                    self%img_buffer((i-1)*c + (j-1) * w * c+ 1) = INT( ISHFT( pixel , -16) ,kind=c_char)
                    self%img_buffer((i-1)*c + (j-1) * w * c + 2) =  INT( IAND( ISHFT( pixel , -8_c_int) , z'000000ff') ,kind=c_char)
                    self%img_buffer((i-1)*c + (j-1) * w * c + 3) =  INT( IAND( pixel , z'000000ff') ,kind=c_char)
                    DebugPrint (i-1)*c + (j-1) * w * c+ 1, self%img_buffer((i)*c + (j) * w * c+ 1)
                else
                    c=1
                    pixel =  INT( REAL( max_colors - 1)*REAL( (in_buffer(i+1,j+1)-lo)/REAL(hi - lo) ) ,kind=c_int)
                    pixel =  IAND( pixel , z'00ffffff')
                    self%img_buffer(i*c + (j*w*c) + 1) = INT(pixel,kind=1)
                    print *, (i*c)+(j*w*c) + 1, in_buffer(i+1,j+1), pixel, self%img_buffer(i*c + j*w*c + 1)
                end if

            end do
        end do
        do i=1,w
            do j=1,h
                print *,i + j*w,  self%img_buffer(i + (j-1)*w)
            end do
        end do

        img = c_loc(self%img_buffer)
#if !defined(PGI)
        VerbosePrint ">>>  save_jpeg_r4 img address: ", img
#endif
        !   if(c_associated(img, c_loc(self%img_buffer))) &
        !       call simple_stop( 'save_jpeg img buffer and C pointer do not point to same target')
        !   print *, "shape img ", shape(img)
        self%width = w
        self%height = h
        VerbosePrint '>>> Calling stbi_write_jpg fname ', fname
        VerbosePrint '>>> Calling stbi_write_jpg width ',self%width
        VerbosePrint '>>> Calling stbi_write_jpg height ',  self%height
        VerbosePrint '>>> Calling stbi_write_jpg comp ', self%colorspace
        VerbosePrint '>>> Calling stbi_write_jpg quality ', self%quality
        status = stbi_write_jpg (fstr, self%width, self%height, self%colorspace, img, self%quality )
        !call write_jpeg(img, fname, w, h, self%quality, self%colorspace, status)
        VerbosePrint '>>> stbi_write_jpg returned status ', status
        !STBI returns 0 on failure and non-0 on success.
        if(status == 0 ) call simple_stop('save_jpeg_r4 call to write_jpeg failed' )
        status = 0
        deallocate(fstr)
        VerbosePrint '>>> Done  save_jpeg_r4 '
    end function save_jpeg_r4

    function save_jpeg_i4 (self, fname, in_buffer, quality, colorspec) result(status)
        class(jpg_img), intent(inout)   :: self
        character(len=*), intent(inout) :: fname
        integer, intent(in)             :: in_buffer(:,:)
        integer, intent(in), optional   :: quality, colorspec
        type(c_ptr)                     :: img
        integer                         :: w,h,j,i,c, lo,hi
        integer                         :: status
        integer(c_int)                  :: pixel
        character(len=:), allocatable    :: fstr
        VerbosePrint '>>>  In save_jpeg_i4 '
        status       = 0
        w            = size(in_buffer,1)
        h            = size(in_buffer,2)
        if(w == 0 .or. h == 0) return
        if(present(quality)) self%quality = quality
        if(present(quality)) self%colorspace = colorspec
        self%width = w
        self%height = h
        lo = minval(in_buffer)
        hi = maxval(in_buffer)
        VerbosePrint '>>> In save_jpeg_i4 input lo hi ', lo, hi
        !        self%img_buffer = transfer(in_buffer,1_c_int8_t)
        !self%img_buffer = reshape(in_buffer, (/ w * h /))
        allocate(self%img_buffer(w*h*3))
        do j=0,h-1
            do i=0,w-1
                if  (self%colorspace == 3) then
                    c=3
                    pixel = NINT( REAL(2**24) * REAL(in_buffer(i+1,j+1)-lo)/REAL(hi-lo),kind=4)

                    self%img_buffer((i-1)*c + (j-1) * w * c+ 1) = INT( ISHFT( pixel , -16) ,kind=c_char)
                    self%img_buffer((i-1)*c + (j-1) * w * c + 2) =  INT( IAND( ISHFT( pixel , -8) , z'000000ff') ,kind=c_char)
                    self%img_buffer((i-1)*c + (j-1) * w * c + 3) =  INT( IAND( pixel , z'000000ff') ,kind=c_char)
                    print *, (i-1)*c + (j-1) * w * c+ 1, self%img_buffer((i-1)*c + (j-1) * w * c+ 1)
                else
                    c=1
                    pixel =  INT( REAL( max_colors - 1)*( REAL(in_buffer(i+1,j+1)-lo) / REAL(hi - lo) ) ,kind=c_int)
                    pixel =  IAND( pixel , z'00ffffff')
                    self%img_buffer(i*c + (j*w*c) + 1) = INT(pixel,kind=1)
                    print *, (i*c)+(j*w*c) + 1, in_buffer(i+1,j+1), pixel, self%img_buffer(i*c + j*w*c + 1)
                end if
            end do
        end do
        img=c_loc(self%img_buffer)
        ! if(c_associated(img, c_loc(self%img_buffer))) &
        !     stop 'save_jpeg_i4 img buffer and C pointer do not point to same target'
        allocate(fstr, source=trim(fname)//c_null_char)
        !call write_jpeg(img, fname, w, h, self%quality, self%colorspace, status)
        status = stbi_write_jpg (fstr, self%width, self%height, self%colorspace, img, self%quality )
        !STBI returns 0 on failure and non-0 on success.
        if( status == 0 ) call simple_stop('save_jpeg_i4 call to write_jpeg returned status ==1' )
        status = 0
        if(allocated(fstr)) deallocate(fstr)
    end function save_jpeg_i4

    function load_jpeg_r4(self, fname, out_buffer) result(status)
        class(jpg_img), intent(inout)    :: self
        character(len=*), intent(inout)  :: fname
        real, allocatable, intent(inout) :: out_buffer(:,:)
        type(c_ptr)                      :: img
        integer                          :: w,h,c, i,j
        integer                          :: status, bshape(1), pixel
        character(len=:), allocatable            :: fstr
         integer(1), dimension(:),pointer :: imgbuffer
        status = 1
        verbose =.true.
        allocate(fstr, source=trim(fname)//c_null_char)
        !        call read_jpeg(fstr, img, w, h, c, status)
        status = stbi_read_jpg(fstr, img, w, h, c )
        if(status ==0) call simple_stop ("simple_jpg::load_jpeg_r4 stbi_read_jpeg failed ")
        status=0
        print *, " read_jpeg returned img width x height ", w, h
        print *, "shape img ", shape(img)
        if(associated(imgbuffer)) nullify(imgbuffer)

        bshape = (/ w * h /)
        call c_f_pointer(img,imgbuffer, bshape)
        VerbosePrint 'load_jpeg_r4 img size / shape', size(imgbuffer), shape(imgbuffer)
        VerbosePrint 'load_jpeg_r4 img_buffer l/u bounds ', lbound(imgbuffer), ubound(imgbuffer)
        VerbosePrint 'load_jpeg_r4 img_buffer ', associated(imgbuffer)
        VerbosePrint 'load_jpeg_r4 img_buffer '
        !        out_buffer = transfer(self%img_buffer,1.0)
        allocate(out_buffer(w,h))
        do i=0,w-1
            do j=0,h-1
                pixel = 0
                if  (c > 1 ) then
                    pixel = ISHFT(INT(imgbuffer((i*c) + (j * w * c)+ 3), kind=4), 16)
                    pixel = pixel + ISHFT(INT(imgbuffer((i*c) + (j * w * c) + 2),kind=4), 8)
                endif
                pixel = pixel + INT(imgbuffer((i*c) + (j * w * c) + 1),kind=4)
                out_buffer(i+1,j+1) = REAL( pixel, kind=4 )
            end do
        end do
        self%width = w
        self%height = h
        self%colorspace = c
        status=0
        !  if(allocated(self%img_buffer)) deallocate(self%img_buffer)
        if(associated(imgbuffer)) nullify(imgbuffer)
        VerbosePrint 'load_jpeg_r4 done'
        if(allocated(fstr))deallocate(fstr)
    end function load_jpeg_r4

    function load_jpeg_i4(self, fname, out_buffer) result(status)
        use,intrinsic   :: iso_c_binding
        implicit none
        class(jpg_img), intent(inout)       :: self
        character(len=*), intent(inout)     :: fname
        integer, allocatable, intent(inout) :: out_buffer(:,:)
        type(c_ptr)                         :: img
        integer                             :: bshape(1)
        integer                             :: i,j,w,h,c
        integer                             :: status,pixel
        integer(1), dimension(:),pointer :: imgbuffer
        character(len=:),allocatable  :: fstr
        status = 1
        verbose =.true.
        allocate(fstr, source=trim(fname)//c_null_char)
        VerbosePrint 'load_jpeg_i4 ', fstr
        !call read_jpeg(fstr, img, w, h, c, status)
        status = stbi_read_jpg (fstr, img, w, h, c )
        if(status==0) call simple_stop ("simple_jpg::load_jpeg_i4 read_jpeg failed ")
        status=0
        !   w = size(self%img_buffer,1)
        !   h = size(self%img_buffer,2)
        VerbosePrint 'load_jpeg_i4 return values width, height, components', w, h, c
 !       VerbosePrint 'load_jpeg_i4 img_buffer l/u bounds ', lbound(img), ubound(img)
        if(associated(imgbuffer)) nullify(imgbuffer)
        bshape = (/ w * h /)
        call c_f_pointer(img,imgbuffer, bshape)
        VerbosePrint 'load_jpeg_i4 img size / shape', size(imgbuffer), shape(imgbuffer)
        VerbosePrint 'load_jpeg_i4 img_buffer l/u bounds ', lbound(imgbuffer), ubound(imgbuffer)
        VerbosePrint 'load_jpeg_i4 img_buffer ', associated(imgbuffer)
        VerbosePrint 'load_jpeg_i4 img_buffer '
        print *, imgbuffer
        allocate(out_buffer(w,h))
        do j=0,h-1
            do i=0,w-1
                VerbosePrint 'load_jpeg_i4 ', (i*c) + (w*j * c) + 1, imgbuffer((i*c) + (w*j*c) + 1)
                out_buffer(i+1,j+1) = imgbuffer((i * c) + (w * j * c) + 1)
            end do
        end do
        self%width = w
        self%height = h
        self%colorspace = c

        if(allocated(fstr)) deallocate(fstr)
        if(associated(imgbuffer)) nullify(imgbuffer)
        VerbosePrint 'load_jpeg_i4 done'
        status=0
    end function load_jpeg_i4

    subroutine destructor(self)
        class(jpg_img) :: self
#if !defined(PGI)
        DebugPrint 'Destructor of jpg_img object with address: ', loc(self)
#endif
        if(associated(self%img_buffer)) nullify(self%img_buffer)
    end subroutine destructor

    subroutine test_jpg_export()
        use,intrinsic                                        :: iso_c_binding
        implicit none
        character                          :: buf(1024)
        integer, allocatable               :: int32_buffer(:,:)
        real, allocatable                  :: real32_buffer(:,:)
        character(len=:), allocatable :: simple_path_str
        character(len=:), allocatable :: testimg
        integer status, bmp_size , width , height , pixel_size
        character(len=:), allocatable              :: cmd, fstr
        type(jpg_img)                      :: jpg
        logical passed
        debug           = .true.
        passed = .true.
        pixel_size      = 3
        width           = 227
        height          = 149
        status          = 1
        bmp_size        = width * height * pixel_size
        !        allocate(int32_buffer(bmp_size))
        status = simple_getenv('SIMPLE_PATH', simple_path_str)
        if(status/=0) call simple_stop(" SIMPLE_PATH not found in environment ")
        print *,"test_jpg_export: Starting"
        allocate(testimg, source=trim(simple_path_str)//"/bin/gui/ext/src/jpeg/testimg.jpg")
        allocate(cmd, source="cp "//trim(adjustl(testimg))//" ./; convert testimg.jpg -colorspace Gray  test.jpg")
        if(.not. file_exists("test.jpg")) then
            if(.not. file_exists(testimg)) call simple_stop("test_jpg_export needs the jpeglib testimg to run test")
            print *," Executing ", cmd
            call exec_cmdline(cmd)
        end if
        deallocate(cmd)
        allocate(cmd, source="identify  test.jpg")
        call exec_cmdline(cmd)

        print *,"test_jpg_export: Testing load_jpeg_i4 "
        allocate(fstr, source="test.jpg")
        status = jpg%loadJpg(fstr, int32_buffer)
        if(status == 0) then
            print *, " test_jpg_export: load_jpeg_i4 success "
            if(allocated(int32_buffer)) then
                width=size(int32_buffer, 1)
                height = size(int32_buffer, 2)
                print *, " test_jpg_export: load_jpeg int32 width x height ", width, height
                if(width /= jpg%getWidth() .or. height /= jpg%getHeight())then
                    print *,"test_jpg_export: FAILED: load_jpeg int32 width height not consistent "
                    status=1
                endif
                print *, "test_jpg_export:  sum | max ", sum(int32_buffer), maxval(int32_buffer)
            else
                print *,"test_jpg_export: FAILED: int32 image buffer not allocated"
                status=1
            end if
        end if
        if(status==1) call simple_stop("test_jpg_export: FAILED: load_jpeg int32 failed")

        print *,"test_jpg_export: Testing load_jpeg_r4 "
        status = jpg%loadJpg(fstr, real32_buffer)
        if(status == 0) then
            print *, "test_jpg_export: load_jpeg_r4 success "
            if(allocated(real32_buffer))then
                width=size(real32_buffer, 1)
                height = size(real32_buffer, 2)
                print *, " width x height ", width, height
                if(width /= jpg%getWidth() .or. height /= jpg%getHeight())then
                    print *,"test_jpg_export: FAILED: load_jpeg real32 width height not consistent "
                    status=1
                endif
                print *,"test_jpg_export:  sum | max ", sum(real32_buffer), maxval(real32_buffer)
            else
                print *,"test_jpg_export: FAILED: real32 image buffer not allocated"
                status=1
            end if
        end if
        if(status==1) call simple_stop("test_jpg_export: FAILED: load_jpeg real32 failed")


        deallocate(fstr)
        allocate(fstr, source="test-i4.jpg")
        call del_file(fstr)
        status = jpg%writeJpg(fstr, int32_buffer)
        if(status == 0) then
            print *, " test_jpg_export: save_jpeg_i4 success "
            if(file_exists("test-i4.jpg"))then
                call exec_cmdline("display test-i4.jpg")
            else
                status=1
            end if
        end if
        if(status==1) call simple_stop("test_jpg_export: FAILED: save_jpeg int32 failed")

        deallocate(fstr)
        allocate(fstr, source="test-r4.jpg")
        call del_file(fstr)
        print *, " test_jpg_export: save_jpeg_r4  "
        status = jpg%writeJpg(fstr, real32_buffer)
        if(status == 0) then
            print *, " test_jpg_export: save_jpeg_r4 success "
            if(file_exists("test-r4.jpg"))then
                call exec_cmdline("display test-r4.jpg")
            else
                status=1
            end if
        end if
        if(status==1) call simple_stop("test_jpg_export: FAILED: save_jpeg real32 failed")

        ! Write the decompressed bitmap out to a ppm file, just to make sure
        ! it worked.
        ! call fopen(fd, "output.ppm", status='old')
        ! write(fd, '("P6",1x,i0,1x,i0,1x,"255")') width, height
        ! write(fd, *) bmp_buffer
        ! call fclose(fd)
        if (allocated(int32_buffer)) deallocate(int32_buffer)
        if (allocated(real32_buffer)) deallocate(real32_buffer)
        if (allocated(simple_path_str)) deallocate( simple_path_str)
        if (allocated(testimg)) deallocate( testimg)
        if (allocated(cmd)) deallocate(cmd)
        if (allocated(fstr))  deallocate(fstr)
    end subroutine test_jpg_export


end module simple_jpg
