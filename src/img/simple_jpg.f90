!! Michael Eager Feb 2018

!! Grayscale Bit representation in JPEG is limited to 24 bits
!! (better than 8-bit PNG)

module simple_jpg
    include 'simple_lib.f08'
    use,intrinsic                                        :: iso_c_binding
    implicit none

    public                                               :: jpg_img, test_jpg_export
    private
    !*** cptr should be have the same size as a c pointer
    !*** It doesn't matter whether it is an integer or a real
    integer, parameter                                   :: cptr          = kind(5)
    integer, parameter                                   :: numComponents = 3;
    integer, parameter                                   :: pixelsize     = selected_int_kind(3)
    integer, parameter                                   :: longint       = selected_int_kind(9)

    integer, public, parameter                           :: max_colors =  256
    !buf_range(1)buf_range(1)im2 = gdImageScale(im, 1, 65535);
    integer, parameter                                   :: GREYSCALE  =  1
#include "simple_local_flags.inc"
    type jpg_img
        private
        integer                                          :: width      =  0
        integer                                          :: height     =  0
        integer(cptr)                                    :: ptr        =  0
        logical                                          :: fmode      =  .false.
    contains
        procedure                                        :: getWidth
        procedure                                        :: getHeight
        procedure, private                               :: save_jpeg_r4
        procedure, private                               :: save_jpeg_i4
        generic                                          :: save_jpeg  => save_jpeg_r4, save_jpeg_i4
        procedure                                        :: load_jpeg_r4
        procedure                                        :: load_jpeg_i4
        generic                                          :: load_jpeg  => load_jpeg_r4, load_jpeg_i4
    end type jpg_img


    interface
        subroutine setup_jpeg ( width, height, in, out ) bind ( c, name="setup_jpeg" )
            import
            integer(C_INT), value                        :: width
            integer(C_INT), value                        :: height
            integer(C_INT), dimension(*), intent(out)    :: in
            integer(C_INT), dimension(*), intent(out)    :: out
        end subroutine setup_jpeg
        subroutine read_jpeg (file_name, img, width, height, colorspec, status ) bind ( c, name="cread_jpeg" )
            import
            character(c_char ), dimension(*), intent(in) :: file_name
            type(C_PTR), value                           :: img
            integer(c_int), intent(inout)                :: width
            integer(c_int), intent(inout)                :: height
            integer(c_int), intent(inout)                :: colorspec
            integer(c_int), intent(inout)                :: status
        end subroutine read_jpeg
        subroutine write_jpeg (img, file_name, width, height, quality, colorspec, status ) bind ( c, name="write_jpeg" )
            import
            type (C_PTR), value                          :: img
            character(c_char),dimension(*),intent(in)    :: file_name
            integer(c_int), intent(in), VALUE            :: width
            integer(c_int), intent(in), VALUE            :: height
            integer(c_int), intent(in)                   :: quality
            integer(c_int), intent(in), VALUE            :: colorspec
            integer(c_int), intent(inout)                :: status
        end subroutine write_jpeg


            function fgltLoadTGA(FileName,width,height,components,eform,image)
                use, intrinsic :: iso_c_binding
                type(c_ptr), target :: fgltloadTGA
                character(len=*), intent(in) :: FileName
                integer(c_int), intent(out) :: width, height
                integer(c_int), intent(out) :: components, eform
                integer(c_char), dimension(:), allocatable, &
                    intent(out), target :: image
            end function fgltLoadTGA


    end interface

contains

    integer function getWidth(self)
        class(jpg_img), intent(in) :: self
        getWidth =self%width
    end function getWidth

    integer function getHeight(self)
        class(jpg_img), intent(in) :: self
        getHeight = self%height
    end function getHeight


    function save_jpeg_r4_3D (self, fname, in_buffer, quality) result(status)
        class(jpg_img), intent(inout)    :: self
        character(len=*), intent(inout)  :: fname
        real,    intent(in)              :: in_buffer(:,:,:)
        integer, intent(in), optional    :: quality
        integer(c_int), allocatable, target    :: img_buffer(:,:)
        type(c_ptr)                      :: img
        integer                          :: w,h,c, i,j, quality_here
        integer                          :: status

        status = 1
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
                img_buffer(i,j) = INT(in_buffer(i,j,1),kind=4)
            end do
        end do
         img = c_loc(img_buffer)
        ! if(c_associated(img, c_loc(img_buffer))) &
        !          stop 'save_jpeg img buffer and C pointer do not point to same target'
        call write_jpeg(img, fname, w, h, quality_here, c, status)
        if(status /= 0 ) call simple_stop('save_jpeg_r4 call to write_jpeg returned status==1' )
        self%width = w
        self%height = h
    end function save_jpeg_r4_3D

    function save_jpeg_r4 (self, fname, in_buffer, quality) result(status)
        class(jpg_img), intent(inout)   :: self
        character(len=*), intent(inout) :: fname
        real, intent(in)                :: in_buffer(:,:)
        integer, intent(in), optional   :: quality
        integer(c_int) , allocatable, target   :: img_buffer(:,:)
        type(c_ptr)                     :: img
        integer                         :: w,h,c, i,j, quality_here
        integer                         :: status

        status = 1
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
        ! if(c_associated(img, c_loc(img_buffer))) &
        !          stop 'save_jpeg img buffer and C pointer do not point to same target'
        call write_jpeg(img, fname, w, h, quality_here, c, status)
        if(status /= 0 ) call simple_stop('save_jpeg_r4 call to write_jpeg returned status==1' )
        self%width = w
        self%height = h
    end function save_jpeg_r4

    function save_jpeg_i4 (self, fname, in_buffer, quality) result(status)
        class(jpg_img), intent(inout)   :: self
        character(len=*), intent(inout) :: fname
        integer, intent(in)             :: in_buffer(:,:)
        integer, intent(in), optional   :: quality
        integer(c_int) , pointer               :: img_buffer(:,:)
        type(c_ptr)                     :: img
        integer                         :: w,h,c, quality_here
        integer                         :: status
        status       = 0
        c            = -1
        quality_here = 100
        w            = size(in_buffer,1)
        h            = size(in_buffer,2)
        if(w == 0 .or. h == 0) return
        if(present(quality)) quality_here = quality
        !        img_buffer = transfer(in_buffer,1_c_int8_t)
        img_buffer = in_buffer
        img=c_loc(img_buffer)
        ! if(c_associated(img, c_loc(img_buffer))) &
        !     stop 'save_jpeg_i4 img buffer and C pointer do not point to same target'
        fname = trim(fname)//c_null_char
        call write_jpeg(img, fname, w, h, quality_here, c, status)
        if(status /= 0 ) call simple_stop('save_jpeg_i4 call to write_jpeg returned status==1' )
        self%width = w
        self%height = h
    end function save_jpeg_i4

    function load_jpeg_r4(self, fname, out_buffer) result(status)
        class(jpg_img), intent(inout)    :: self
        character(len=*), intent(inout)  :: fname
        real, allocatable, intent(inout) :: out_buffer(:,:)
        integer(c_int) , pointer                :: img_buffer(:,:)
        type(c_ptr)                      :: img
        integer                          :: w,h,c, i,j
        integer                          :: status, shape(2)
        character(len=STDLEN)            :: fstr
        status = 1
        fstr   = trim(fname)//c_null_char
        call read_jpeg(fstr, img, w, h,c, status)
        if(status /=0) call simple_stop ("simple_jpg::load_jpeg_r4 read_jpeg failed ")
        w = size(img_buffer,1)
        h=  size(img_buffer,2)
        allocate(out_buffer(w,h))
        shape = (/ w, h /)
        call c_f_pointer(img,img_buffer, shape)
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
        class(jpg_img), intent(inout)       :: self
        character(len=*), intent(inout)     :: fname
        integer, allocatable, intent(inout) :: out_buffer(:,:)
        integer(c_int), pointer             :: img_buffer(:,:)
        type(c_ptr)                         :: img
        integer                             :: shape(2)
        integer                             :: i,j,w,h,c
        integer                             :: status
        character(len=STDLEN)               :: fstr
        status = 1
        fstr   = trim(fname)//c_null_char
        DebugPrint 'load_jpeg_i4 ', fstr
        call read_jpeg(fstr, img, w, h, c, status)
        if(status/=0) call simple_stop ("simple_jpg::load_jpeg_i4 read_jpeg failed ")
        !   w = size(img_buffer,1)
        !   h = size(img_buffer,2)
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
        DebugPrint 'load_jpeg_i4 done'
    end function load_jpeg_i4

 ! write a grayscale image as a color JPEG image

  subroutine write_jpeg_grayscale(name,image)
    character(*) name
    real :: image(:,:)
    real, allocatable :: temp(:,:,:)
    logical unnormalized
    allocate(temp(3,size(image,1),size(image,2)))
    temp(1,:,:) = image
    temp(2,:,:) = image
    temp(3,:,:) = image
    call write_jpeg(name,temp)
  end subroutine write_jpeg_grayscale

  subroutine rgb2channels(image)
    real, allocatable :: image(:,:,:)
    real, allocatable :: temp(:,:,:)
    integer i
    allocate(temp(size(image,2),size(image,3),size(image,1)))
    forall (i=1:size(image,1)) temp(:,:,i) = image(i,:,:)
    deallocate(image)
    allocate(image(size(temp,1),size(temp,2),size(temp,3)))
    image = temp
  end subroutine rgb2channels

  subroutine channels2rgb(image)
    real, allocatable :: image(:,:,:)
    real, allocatable :: temp(:,:,:)
    integer i
    allocate(temp(size(image,3),size(image,1),size(image,2)))
    forall (i=1:size(image,3)) temp(i,:,:) = image(:,:,i)
    deallocate(image)
    allocate(image(size(temp,1),size(temp,2),size(temp,3)))
    image = temp
  end subroutine channels2rgb

    subroutine test_jpg_export()
        character                          :: buf(1024)
        integer, allocatable               :: int32_buffer(:,:)
        real, allocatable                  :: real32_buffer(:,:)
        character(len=STDLEN), allocatable :: simple_path_str
        character(len=STDLEN), allocatable :: testimg
        integer status, bmp_size , width , height , pixel_size
        character(len=STDLEN)              :: cmd, fstr
        type(jpg_img)                      :: jpg
        debug           = .true.
        pixel_size      = 3
        width           = 227
        height          = 149
        status          = 1
        bmp_size        = width * height * pixel_size
        !        allocate(int32_buffer(bmp_size))
        simple_path_str = simple_getenv('SIMPLE_PATH')
        print *,"test_jpg_export: Starting"
        allocate(testimg, source=trim(simple_path_str)//"/bin/gui/ext/src/jpeg/testimg.jpg")
        cmd             = "convert " //trim(adjustl(testimg))//" -colorspace Gray  test.jpg"
        if(.not. file_exists("test.jpg")) then
            if(.not. file_exists(testimg)) call simple_stop("test_jpg_export needs the jpeglib testimg to run test")
            print *," Executing ", cmd
            call exec_cmdline(cmd)
        end if
        cmd = "identify  test.jpg"
        call exec_cmdline(cmd)
        print *,"test_jpg_export: Testing load_jpeg_i4 "
        fstr = "testimg.jpg"
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
