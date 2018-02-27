!! Michael Eager Feb 2018

!! Grayscale Bit representation in JPEG is limited to 24 bits
!! (better than 8-bit PNG)

module simple_jpg
    include 'simple_lib.f08'
    use,intrinsic ::  iso_c_binding
    implicit none

    public:: jpg_img, test_jpg_export
private
    !*** cptr should be have the same size as a c pointer
    !*** It doesn't matter whether it is an integer or a real
    integer, parameter :: cptr = kind(5)
    integer, parameter :: numComponents = 3;
    integer, parameter :: pixelsize = selected_int_kind(3)
    integer, public, parameter :: max_colors = 256
    !buf_range(1)buf_range(1)im2 = gdImageScale(im, 1, 65535);
    integer, parameter :: GREYSCALE = 1
#include "simple_local_flags.inc"
    type jpg_img
        private
        integer        :: width = 0
        integer        :: height = 0
        integer(cptr)  :: ptr = 0
        logical        :: fmode = .false.
    contains
        procedure          :: getWidth
        procedure          :: getHeight
        procedure, private :: save_jpeg_r4
        procedure, private :: save_jpeg_i4
        generic :: save_jpeg => save_jpeg_r4, save_jpeg_i4
        procedure :: load_jpeg_r4
        procedure :: load_jpeg_i4
        generic :: load_jpeg => load_jpeg_r4, load_jpeg_i4
    end type jpg_img


    interface
        subroutine read_jpeg (file_name, img, width, height, colorspec, status ) bind ( c, name="read_jpeg" )
            use iso_c_binding
            character ( c_char ), dimension(*), intent(in) :: file_name
            type ( c_ptr ) :: img
            integer ( c_int ) :: width
            integer ( c_int ) :: height
            integer ( c_int ) :: colorspec
            logical ( c_bool ) :: status
        end subroutine read_jpeg
        subroutine write_jpeg (img, file_name, width, height, quality, colorspec, status ) bind ( c, name="write_jpeg" )
            use iso_c_binding
            type ( c_ptr ), VALUE :: img
            character ( c_char), dimension(*), intent(in) :: file_name
            integer ( c_int ),VALUE :: width
            integer ( c_int ),VALUE :: height
            integer ( c_int ),VALUE :: quality
            integer ( c_int ),VALUE :: colorspec
            logical ( c_bool ) :: status
        end subroutine write_jpeg
    end interface

contains

    integer function getWidth(self)
        class(jpg_img), intent(in) :: self
        getWidth =self%width
    end function getWidth

    integer function getHeight(self)
        class(jpg_img), intent(in) :: self
        getHeight =self%height
    end function getHeight


    function save_jpeg_r4 (self, fname, in_buffer, quality) result(status)
        class(jpg_img), intent(inout) :: self
        character(len=*), intent(inout) :: fname
        real, intent(in) :: in_buffer(:,:)
        integer, intent(in), optional :: quality
        integer , allocatable, target :: img_buffer(:,:)
        type(c_ptr) :: img
        integer :: w,h,c, i,j, quality_here
        logical(1) :: status
        status = .false.
        c = -1
        quality_here = 100
        w = size(in_buffer,1)
        h= size(in_buffer,2)
        if(w == 0 .or. h == 0) return
        if(present(quality)) quality_here = quality
        !        img_buffer = transfer(in_buffer,1_c_int8_t)
        allocate(img_buffer(w,h))
        do i=1,w
            do j=1,h
                img_buffer(i,j) = INT(in_buffer(i,j),kind=4)
            end do
        end do
        img = c_loc(img_buffer)
        if(c_associated(img, c_loc(img_buffer))) &
                 stop 'save_jpeg img buffer and C pointer do not point to same target'
        call write_jpeg(img, trim(fname)//c_null_char, w, h, quality_here, c, status)

    end function save_jpeg_r4

    function save_jpeg_i4 (self, fname, in_buffer, quality) result(status)
        class(jpg_img), intent(inout) :: self
        character(len=*), intent(inout) :: fname
        integer, intent(in) :: in_buffer(:,:)
        integer, intent(in), optional :: quality
        integer , pointer :: img_buffer(:,:)
        type(c_ptr) :: img
        integer :: w,h,c, quality_here
        logical(1) :: status
        status = .false.
        c = -1
        quality_here = 100
        w = size(in_buffer,1)
        h = size(in_buffer,2)
        if(w == 0 .or. h == 0) return
        if(present(quality)) quality_here = quality
        !        img_buffer = transfer(in_buffer,1_c_int8_t)
        img_buffer = in_buffer
        img=c_loc(img_buffer)
        if(c_associated(img, c_loc(img_buffer))) &
            stop 'save_jpeg_i4 img buffer and C pointer do not point to same target'
        call write_jpeg(img, trim(fname)//c_null_char, w, h, quality_here, -1, status)
        self%width = w
        self%height = h
    end function save_jpeg_i4

    function load_jpeg_r4(self, fname, out_buffer) result(status)
        class(jpg_img), intent(inout) :: self
        character(len=*), intent(inout) :: fname
        real, allocatable, intent(inout) :: out_buffer(:,:)
        integer , pointer :: img_buffer(:,:)
        type(c_ptr) :: img
        integer :: w,h,c, i,j
        logical(1) :: status
        character(len=STDLEN) :: fstr
        status = .false.
        fstr = trim(fname)//c_null_char
        call read_jpeg(fstr, img, w, h,c, status)
        if(.not. status) call simple_stop ("simple_jpg::load_jpeg_r4 read_jpeg failed ")
        w = size(img_buffer,1)
        h=  size(img_buffer,2)
        allocate(out_buffer(w,h))
        call c_f_pointer(img,img_buffer, [w,h])
        !        out_buffer = transfer(img_buffer,1.0)
        do i=1,w
            do j=1,h
                out_buffer(i,j) = REAL(img_buffer(i,j))
            end do
        end do
         self%width = w
        self%height = h
      !  if(allocated(img_buffer)) deallocate(img_buffer)
    end function load_jpeg_r4

    function load_jpeg_i4(self, fname, out_buffer) result(status)
        class(jpg_img), intent(inout) :: self
        character(len=*), intent(inout) :: fname
        integer, allocatable, intent(inout) :: out_buffer(:,:)
        integer, pointer :: img_buffer(:,:)
        type(c_ptr) :: img
        integer :: shape(2)
        integer :: i,j,w,h,c
        logical(1) :: status
        character(len=STDLEN) :: fstr
        status = .false.
        fstr = trim(fname)//c_null_char
        DebugPrint 'load_jpeg_i4 ', fstr
        call read_jpeg(fstr, img, w, h, c, status)
        !if(.not. status) call simple_stop ("simple_jpg::load_jpeg_i4 read_jpeg failed ")
        w = size(img_buffer,1)
        h = size(img_buffer,2)
        DebugPrint 'load_jpeg_i4 returns img size ', w, h
        shape = (/ w, h /)
        call c_f_pointer(img,img_buffer, shape)
        DebugPrint 'load_jpeg_i4 c_f_pointer  ', w, h

        allocate(out_buffer(w,h))
        do i=1,w
            do j=1,h
                out_buffer(i,j) = INT(img_buffer(i, j))
            end do
        end do
        self%width = w
        self%height = h
    !    if(allocated(img_buffer)) deallocate(img_buffer)
    end function load_jpeg_i4


    subroutine test_jpg_export()
        character :: buf(1024)
        integer, allocatable :: int32_buffer(:,:)
        real, allocatable :: real32_buffer(:,:)
        character(len=STDLEN), allocatable :: simple_path_str
        character(len=STDLEN), allocatable :: testimg
        integer status, bmp_size , width , height , pixel_size
        character(len=STDLEN) :: cmd, fstr
        type(jpg_img) :: jpg
        debug=.true.
        pixel_size= 3
        width =  227
        height = 149
        status=1
        bmp_size = width * height * pixel_size
!        allocate(int32_buffer(bmp_size))
        simple_path_str = simple_getenv('SIMPLE_PATH')
        print *,"test_jpg_export: Starting"
        allocate(testimg, source=trim(simple_path_str)//"/bin/gui/ext/src/jpeg/testimg.jpg")
        cmd = "convert " //trim(adjustl(testimg))//" -colorspace Gray  test.jpg"
        if(.not. file_exists("test.jpg")) then
            if(.not. file_exists(testimg)) call simple_stop("test_jpg_export needs the jpeglib testimg to run test")
            print *," Executing ", cmd
            call exec_cmdline(cmd)
        end if
        print *,"test_jpg_export: Testing load_jpeg_i4 "
        fstr = "test.jpg"
        status = jpg%load_jpeg_i4(fstr, int32_buffer)
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


        status = jpg%load_jpeg_r4(fstr, real32_buffer)
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

        fstr="test-i4.jpg"
        call del_file(fstr)
        status = jpg%save_jpeg_i4(fstr, int32_buffer)
        if(status == 0) then
            print *, " test_jpg_export: save_jpeg_i4 success "
            if(file_exists("test-i4.jpg"))then
                call exec_cmdline("identify test-i4.jpg")
            else
                status=1
            end if
        end if
        if(status==1) call simple_stop("test_jpg_export: FAILED: save_jpeg int32 failed")

        fstr="test-r4.jpg"
        call del_file(fstr)
        status = jpg%save_jpeg_r4(fstr, real32_buffer)
        if(status == 0) then
            print *, " test_jpg_export: save_jpeg_i4 success "
            if(file_exists("test-i4.jpg"))then
                call exec_cmdline("identify test-i4.jpg")
            else
                status=1
            end if
        end if
        if(status==1) call simple_stop("test_jpg_export: FAILED: save_jpeg int32 failed")


        ! Write the decompressed bitmap out to a ppm file, just to make sure
        ! it worked.
        ! call fopen(fd, "output.ppm", status='old')
        ! write(fd, '("P6",1x,i0,1x,i0,1x,"255")') width, height
        ! write(fd, *) bmp_buffer
        ! call fclose(fd)
        if (allocated(int32_buffer)) deallocate(int32_buffer)
        if (allocated(real32_buffer)) deallocate(real32_buffer)

    end subroutine test_jpg_export


end module simple_jpg
