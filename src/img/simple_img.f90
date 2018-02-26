!*******************************************************************************
! Fortran module interfacing the GD Graphics Library http://www.boutell.com/gd/
! Version 1.0
! For documentation and all other relevant information:
! Mart.Rentmeester@nn-online.org    http://nn-online.org/code/
!*******************************************************************************
!
!! Modified by Michael Eager Feb 2018
module simple_img
    include 'simple_lib.f08'
    implicit none

    !*** cptr should be have the same size as a c pointer
    !*** It doesn't matter whether it is an integer or a real
    integer, parameter :: cptr = kind(5)
    integer, public, parameter :: max_colors = 256
    !buf_range(1)buf_range(1)im2 = gdImageScale(im, 1, 65535);

    type base_img
        private
        integer        :: width = 0
        integer        :: height = 0
        integer(cptr)  :: ptr = 0
        logical        :: fmode = .false.
    end type base_img

    type img_font
        private
        integer(cptr)  :: ptr = 0
    end type img_font


! interface
! subroutine cgd_image_create(x,y, ptr) bind(C, name="gdImageCreate")
! use, intrinsic :: iso_c_binding
! integer(c_int), value :: x,y
! integer(kind=kind(5)) :: ptr
! end subroutine CGD_IMAGE_CREATE
! end interface

public

contains

    subroutine read_png(file_name,buffer,status)
        character(*), intent(in)         :: file_name
        real, allocatable, intent(inout) :: buffer(:,:)
        integer, intent(out), optional  :: status
        integer             :: width,height, i,j, colorval,black,white
        type(base_img)      :: image

        if (present(status)) then
            call create_img_from_png(file_name,image,status)
        else
            call create_img_from_png(file_name,image)
        end if
        !  call allocate_color(image,0,0,0,black)
        !  call allocate_color(image,255,255,255,white)
        !  colorval=greyscale(image)
        !  print*,' Greyscale color ', colorval
        width  = get_width(image)
        height = get_height(image)
        allocate( buffer(height, width) )
        do i=1,width
            do j=1,height
                colorval = get_pixel(image,i,j)
                buffer(i,j) =REAL( red(image,colorval)*max_colors**2) + &
                    REAL(green(image,colorval)*max_colors) + &
                    REAL(blue(image,colorval))      ! REAL(MODULO(get_pixel(image,i,j),max_colors))
            end do
        end do

        print *," Input png file ", file_name
        print *,"  Number of colors ", number_of_colors(image)
        print *,"  Interlaced? ", get_interlaced(image)
        print *,"  Anti-aliased? ", anti_aliased()
        print *,"  Fmode ", get_fmode(image)
        call destroy_img(image)
    end subroutine read_png

    subroutine write_png(buffer,file_name,range_in,status)
        real, allocatable, intent(in) :: buffer(:,:)
        character(len=*), intent(in) :: file_name
        real, intent(in), optional   :: range_in(2)
        integer, intent(out), optional  :: status
        integer          :: width,height,i,j,k,colorval,ex, r,g,b, black,white ,colors(256,256,256)
        type(base_img)   :: image
        character(len=STDLEN) :: filename
        real                  :: buf_range(2), offset
        logical               :: no_range
        no_range = .true.
        if(.not. allocated(buffer)) then
            status=-1
            print*,"simple_img::write_png buffer not allocated"
            return
        end if
        write(filename,'(A)') file_name
        buf_range(1) = MINVAL(buffer)
        buf_range(2) = MAXVAL(buffer)
        if(present(range_in)) then
            buf_range=range_in
            no_range = .false.
        endif

        offset=0.0
        colors=0
        !colorval=greyscale(image)
        !print*,' Greyscale color ', colorval
        width  = size(buffer,1)
        height = size(buffer,2)

        print *,"Max value of buffer ", buf_range(2)
        if (maxval(buffer) >buf_range(2)) then
            buf_range(2) = maxval(buffer)
            print *,"Max value of buffer reset to ", buf_range(2)
        end if
        if (minval(buffer) < buf_range(1)) then
            buf_range(1) = minval(buffer)
            print *,"Min value of buffer reset to ", buf_range(1)
        end if
        if(buf_range(1) < 0.0 ) then
            offset = buf_range(1)
            print *,"Offset value of buffer reset to ", offset
        end if
        call create_img(width,height,image)
        !call save_alpha(image,.true.)

        ! call allocate_color(image,0,0,0,black)
        ! call allocate_color(image,255,255,255,white)
        ! call  allocate_color(image, 255, 0, 0,r)
        ! call allocate_color(image, 0, 255, 0,g)
        ! call  allocate_color(image, 0, 0, 255,b)
        ! do i=0,255
        !     do j=0,255
        !         do k=0,255
        !             call allocate_color(image,i,j,k,colors(i,j,k))
        !         enddo
        !     enddo
        ! enddo

        call set_interlace(image,.false.)
        do i=0,width-1
            do j=0,height-1
                if (no_range) then
                    ex= INT(buffer(i+1,j+1))
                else
                    colorval= INT( ( buffer(i+1,j+1) - offset / abs(buf_range(2)-buf_range(1)) )  * max_colors**3 )
                    if (colorval >= max_colors**3) colorval = max_colors**3 - 1
                    if (colorval < 0 ) colorval = 0
                    r = modulo(colorval ,  max_colors* max_colors) +1
                    g = modulo(colorval-r -1,  max_colors) +1
                    b = colorval - r - g + 1
                    ex = closest_color(image,r,g,b)
                    !  if(colors(r+1,g+1,b+1)==0)  call allocate_color(image,i,j,k,colors(r,g,b))
                end if
                call set_pixel(image, i, j, ex )
            end do
        end do

        !        colorval=greyscale(image)
        !        print*,' Greyscale color ', colorval
        call cgd_image_png(image%ptr, trim(adjustl(filename))//char(0) )
        call destroy_img(image)
    end subroutine write_png



    subroutine read_jpeg(file_name,buffer,status)
        character(*), intent(in)         :: file_name
        real, allocatable, intent(inout) :: buffer(:,:)
        integer, intent(out), optional  :: status
        integer             :: width,height, i,j, colorval,black,white
        type(base_img)      :: image

        if (present(status)) then
            call create_img_from_jpeg(file_name,image,status)
        else
            call create_img_from_jpeg(file_name,image)
        end if
        !  call allocate_color(image,0,0,0,black)
        !  call allocate_color(image,255,255,255,white)
        !  colorval=greyscale(image)
        !  print*,' Greyscale color ', colorval
        width  = get_width(image)
        height = get_height(image)
        allocate( buffer(height, width) )
        do i=1,width
            do j=1,height
                colorval = get_pixel(image,i,j)
                buffer(i,j) =REAL( red(image,colorval)*max_colors**2) + &
                    REAL(green(image,colorval)*max_colors) + &
                    REAL(blue(image,colorval))      ! REAL(MODULO(get_pixel(image,i,j),max_colors))
            end do
        end do

        print *," Input jpeg file ", file_name
        print *,"  Number of colors ", number_of_colors(image)
        print *,"  Interlaced? ", get_interlaced(image)
        print *,"  Anti-aliased? ", anti_aliased()
        print *,"  Fmode ", get_fmode(image)
        call destroy_img(image)
    end subroutine read_jpeg

    subroutine write_jpeg(buffer,file_name,range_in,status)
        real, allocatable, intent(in) :: buffer(:,:)
        character(len=*), intent(in) :: file_name
        real, intent(in), optional   :: range_in(2)
        integer, intent(out), optional  :: status
        integer          :: width,height,i,j,k,colorval,ex, r,g,b, black,white ,colors(256,256,256)
        type(base_img)   :: image
        character(len=STDLEN) :: filename
        real                  :: buf_range(2), offset
        logical               :: no_range
        integer, allocatable :: int_buffer(:,:)
        no_range = .true.
        if(.not. allocated(buffer)) then
            status=-1
            print*,"simple_img::write_jpeg buffer not allocated"
            return
        end if
        write(filename,'(A)') file_name
        buf_range(1) = MINVAL(buffer)
        buf_range(2) = MAXVAL(buffer)
        if(present(range_in)) then
            buf_range=range_in
            no_range = .false.
        endif

        offset=0.0
        colors=0

        width  = size(buffer,1)
        height = size(buffer,2)

        print *,"Max value of buffer ", buf_range(2)
        if (maxval(buffer) >buf_range(2)) then
            buf_range(2) = maxval(buffer)
            print *,"Max value of buffer reset to ", buf_range(2)
        end if
        if (minval(buffer) < buf_range(1)) then
            buf_range(1) = minval(buffer)
            print *,"Min value of buffer reset to ", buf_range(1)
        end if
        if(buf_range(1) < 0.0 ) then
            offset = buf_range(1)
            print *,"Offset value of buffer reset to ", offset
        end if
        allocate(int_buffer(width,height))
        int_buffer = INT(buffer)
        call cgd_image_create_truecolor(width,height,image)
        ! do i=0,width-1
        !     do j=0,height-1
        !         if (no_range) then
        !             ex= INT(buffer(i+1,j+1))
        !         else
        !             colorval= INT( ( buffer(i+1,j+1) - offset / abs(buf_range(2)-buf_range(1)) )  * max_colors**3 )
        !             if (colorval >= max_colors**3) colorval = max_colors**3 - 1
        !             if (colorval < 0 ) colorval = 0
        !             r = modulo(colorval ,  max_colors* max_colors) +1
        !             g = modulo(colorval-r -1,  max_colors) +1
        !             b = colorval - r - g + 1
        !             ex = closest_color(image,r,g,b)
        !             !  if(colors(r+1,g+1,b+1)==0)  call allocate_color(image,i,j,k,colors(r,g,b))
        !         end if
        !         call set_pixel(image, i, j, ex )
        !     end do
        ! end do
        call cgd_image_jpeg_buffer_put(image%ptr, width*height,int_buffer)
        !        colorval=greyscale(image)
        !        print*,' Greyscale color ', colorval
        call cgd_image_jpeg(image%ptr, trim(adjustl(filename))//char(0) )
        call destroy_img(image)
    end subroutine write_jpeg


    logical function compare_imgs(f1,f2)
        character(len=*), intent(in) :: f1,f2
        integer :: chk
        type(base_img)   :: image1, image2
        compare_imgs = .false.
        chk=0
        call create_img_from_png(f1,image1,chk)
        if(chk==0) then
            call create_img_from_png(f2,image2,chk)
            if(chk==0) then
                call cgd_image_compare(image1,image2,chk)
                compare_imgs = chk == 0
                call destroy_img(image1)
                call destroy_img(image2)
            else
                call destroy_img(image2)
            end if
        else
            call destroy_img(image1)
        end if

    end function compare_imgs

    subroutine create_jpeg(width,height,image,status)
        integer, intent(in)             :: width,height
        type(base_img), intent(out)   :: image
        integer, intent(out), optional  :: status

        if (present(status)) then
            call create_truecolor_img(width,height,image,status)
        else
            call create_truecolor_img(width,height,image)
        endif

    end subroutine create_jpeg

    subroutine create_img(width,height,image,status)
        integer, intent(in)             :: width,height
        type(base_img), intent(out)   :: image
        integer, intent(out), optional  :: status

        if (present(status)) then
            call create_palette_img(width,height,image,status)
        else
            call create_palette_img(width,height,image)
        endif

    end subroutine create_img

    subroutine create_palette_img(width,height,image,status)
        integer, intent(in)             :: width,height
        type(base_img), intent(out)     :: image
        integer, intent(out), optional  :: status

        if (image%ptr /= 0) call destroy_img(image)
        if (width <= 0 .or. height <= 0) then
            if (present(status)) status = 2
            return
        endif
        call cgd_image_create(width,height,image%ptr)
        if (present(status)) then
            if (image%ptr /= 0) then
                image%width = get_width(image)
                image%height = get_height(image)
                status = 0
            else
                status = 1
            endif
        endif

    end subroutine create_palette_img

    subroutine create_truecolor_img(width,height,image,status)
        integer, intent(in)             :: width,height
        type(base_img), intent(out)     :: image
        integer, intent(out), optional  :: status

        if (image%ptr /= 0) call destroy_img(image)

        if (width <= 0 .or. height <= 0) then
            if (present(status)) status = 2
            return
        endif
        call cgd_image_create_truecolor(width,height,image%ptr)
        if (present(status)) then
            if (image%ptr /= 0) then
                status = 0
                image%width = get_width(image)
                image%height = get_height(image)
            else
                status = 1
            endif
        endif
    end subroutine create_truecolor_img

    subroutine destroy_img(image)
        type(base_img), intent(inout)  :: image

        call cgd_image_destroy(image%ptr)
        image%ptr = 0
        image%width = 0
        image%height = 0
        image%fmode = .false.

    end subroutine destroy_img

    subroutine allocate_color(image,r,g,b,color)
        type(base_img), intent(in)  :: image
        integer, intent(in)         :: r,g,b
        integer, intent(out)        :: color

        call cgd_image_color_allocate(image%ptr,r,g,b,color)

    end subroutine allocate_color

    subroutine allocate_alpha_color(image,r,g,b,a,color)
        type(base_img), intent(in)  :: image
        integer, intent(in)         :: r,g,b,a
        integer, intent(out)        :: color

        call cgd_image_color_allocate_alpha(image%ptr,r,g,b,a,color)

    end subroutine allocate_alpha_color

    subroutine deallocate_color(image,color)
        type(base_img), intent(in)  :: image
        integer, intent(in)         :: color
        call cgd_image_color_deallocate(image%ptr,color)
    end subroutine deallocate_color

    subroutine set_pixel(image,x,y,c)
        type(base_img), intent(in)  :: image
        integer, intent(inout)         :: x,y,c
        call cgd_image_set_pixel(image%ptr,xm(image,x),ym(image,y),c)
    end subroutine set_pixel

    function get_pixel(image,x,y) result(color)
        type(base_img), intent(in)  :: image
        integer, intent(inout)         :: x,y
        integer  :: color
        call cgd_image_get_pixel(image%ptr,xm(image,x),ym(image,y),color)
    end function get_pixel

    subroutine draw_line(image,x1,y1,x2,y2,color)
        type(base_img), intent(in)  :: image
        integer, intent(in)         :: x1,y1,x2,y2,color

        call cgd_image_line(image%ptr,xm(image,x1),ym(image,y1),  &
            xm(image,x2),ym(image,y2),color)

    end subroutine draw_line

    subroutine draw_rectangle(image,x1,y1,x2,y2,color)
        type(base_img), intent(in)  :: image
        integer, intent(in)         :: x1,y1,x2,y2,color

        call cgd_image_rectangle(image%ptr,xm(image,x1),ym(image,y1),  &
            xm(image,x2),ym(image,y2),color)

    end subroutine draw_rectangle

    subroutine draw_polygon(image,n,x,y,color)
        type(base_img), intent(in)  :: image
        integer, intent(in)         :: n,x(n),y(n),color

        call cgd_image_polygon(image%ptr,n,xm(image,x),ym(image,y),color)

    end subroutine draw_polygon

    subroutine draw_arc(image,x,y,w,h,s,e,color)
        type(base_img), intent(in)  :: image
        integer, intent(in)         :: x,y,w,h,s,e,color

        call cgd_image_arc(image%ptr,xm(image,x),ym(image,y),w,h,s,e,color)

    end subroutine draw_arc

    subroutine draw_filled_polygon(image,n,x,y,color)
        type(base_img), intent(in)  :: image
        integer, intent(in)         :: n,x(n),y(n),color

        call cgd_image_filled_polygon(image%ptr,n,xm(image,x),ym(image,y),color)

    end subroutine draw_filled_polygon

    subroutine draw_filled_rectangle(image,x1,y1,x2,y2,color)
        type(base_img), intent(in)  :: image
        integer, intent(in)         :: x1,y1,x2,y2,color

        call cgd_image_filled_rectangle(image%ptr,xm(image,x1),ym(image,y1),  &
            xm(image,x2),ym(image,y2),color)

    end subroutine draw_filled_rectangle

    subroutine draw_filled_arc(image,x,y,w,h,s,e,color,style)
        type(base_img), intent(in)  :: image
        integer, intent(in)         :: x,y,w,h,s,e,color,style

        call cgd_image_filled_arc(image%ptr,xm(image,x),ym(image,y),w,h,s,e,color,style)

    end subroutine draw_filled_arc

    subroutine draw_filled_ellipse(image,x,y,w,h,color,style)
        type(base_img), intent(in)  :: image
        integer, intent(in)         :: x,y,w,h,color,style

        call cgd_image_filled_ellipse(image%ptr,xm(image,x),ym(image,y),w,h,color,style)

    end subroutine draw_filled_ellipse

    subroutine write_img_as_png(img,file)
        type(base_img), intent(in)  :: img
        character(*), intent(in)    :: file

        call cgd_image_png(img%ptr,trim(file)//char(0))

    end subroutine write_img_as_png

    subroutine create_img_from_png(file,image,status)
        character(*), intent(in)        :: file
        type(base_img), intent(out)     :: image
        integer, optional, intent(out)  :: status

        if (image%ptr /= 0) call destroy_img(image)
        call cgd_image_create_from_png(trim(file)//char(0),image%ptr)
        if (present(status)) then
            if (image%ptr /= 0) then
                status = 0
            else
                status = 1
            endif
        endif
        image%width = get_width(image)
        image%height = get_height(image)

    end subroutine create_img_from_png

#ifdef WITH_JPEG
    subroutine create_img_from_jpeg(file,image,status)
        character(*), intent(in)        :: file
        type(base_img), intent(out)     :: image
        integer, intent(out), optional  :: status

        if (image%ptr /= 0) call destroy_img(image)
        call cgd_image_create_from_jpeg(trim(file)//char(0),image%ptr)
        if (present(status)) then
            if (image%ptr /= 0) then
                status = 0
            else
                status = 1
            endif
        endif
        image%width = get_width(image)
        image%height = get_height(image)

    end subroutine create_img_from_jpeg

    subroutine write_img_as_jpeg(img,file,quality)
        type(base_img), intent(in)  :: img
        character(*), intent(in)    :: file
        integer, intent(in)         :: quality

        call cgd_image_jpeg(img%ptr,trim(file)//achar(0),quality)

    end subroutine write_img_as_jpeg
#endif

    subroutine fill_img(image,x,y,color)
        type(base_img), intent(in)  :: image
        integer, intent(in)         :: x,y,color

        call cgd_image_fill(image%ptr,xm(image,x),ym(image,y),color)
    end subroutine fill_img


    subroutine fill_to_border(image,x,y,bordercolor,color)
        type(base_img), intent(in)  :: image
        integer, intent(in)         :: x,y,bordercolor,color

        call cgd_image_fill_to_border(image%ptr,xm(image,x),ym(image,y),  &
            bordercolor,color)
    end subroutine fill_to_border

    function get_width(image)
        type(base_img), intent(in)  :: image
        integer                     :: get_width

        call cgd_width(image%ptr,get_width)
    end function get_width

    function get_height(image)
        type(base_img), intent(in)  :: image
        integer                     :: get_height
        call cgd_height(image%ptr,get_height)
    end function get_height

    function width(image)
        type(base_img), intent(in)  :: image
        integer                     :: width
        width = image%width
    end function width

    function height(image)
        type(base_img), intent(in)  :: image
        integer                     :: height
        height = image%height
    end function height

#ifdef WITH_XPM
    subroutine create_img_from_xpm(file,image,status)
        character(*), intent(in)        :: file
        type(base_img), intent(out)     :: image
        integer, optional, intent(out)  :: status

        if (image%ptr /= 0) call destroy_img(image)
        call cgd_image_create_from_xpm(trim(file)//char(0),image%ptr)
        if (present(status)) then
            if (image%ptr /= 0) then
                status = 0
            else
                status = 1
            endif
        endif
        image%width = get_width(image)
        image%height = get_height(image)
    end subroutine create_img_from_xpm
    subroutine write_img_as_xpm(img,file)
        type(base_img), intent(in)  :: img
        character(*), intent(in)    :: file

        call cgd_image_xpm(img%ptr,trim(file)//char(0))

    end subroutine write_img_as_xpm
    subroutine create_img_from_xbm(file,image,status)
        character(*), intent(in)        :: file
        type(base_img), intent(out)     :: image
        integer, optional, intent(out)  :: status

        if (image%ptr /= 0) call destroy_img(image)
        call cgd_image_create_from_xbm(trim(file)//char(0),image%ptr)
        if (present(status)) then
            if (image%ptr /= 0) then
                status = 0
            else
                status = 1
            endif
        endif
        image%width = get_width(image)
        image%height = get_height(image)
    end subroutine create_img_from_xbm
#endif

    subroutine set_brush(image,brush)
        type(base_img), intent(in)  :: image,brush

        call cgd_image_set_brush(image%ptr,brush%ptr)
    end subroutine set_brush

    subroutine set_tile(image,tile)
        type(base_img), intent(in)  :: image,tile

        call cgd_image_set_tile(image%ptr,tile%ptr)
    end subroutine set_tile

    subroutine set_style(image,style,length)
        type(base_img), intent(in)  :: image
        integer, intent(in)         :: length,style(length)

        call cgd_image_set_style(image%ptr,style,length)

    end subroutine set_style

    subroutine set_alpha_blending(image,b)
        type(base_img), intent(in)  :: image
        logical, intent(in)         :: b
        integer                     :: i

        i = 1
        if (b) i = 0
        call cgd_image_set_alpha_blending(image%ptr,i)

    end subroutine set_alpha_blending

    subroutine save_alpha(image,b)
        type(base_img), intent(in)  :: image
        logical, intent(in)         :: b
        integer                     :: i

        i = 0
        if (b) i = 1
        call cgd_image_save_alpha(image%ptr,i)

    end subroutine save_alpha

    function blue(image,color)
        type(base_img), intent(in)  :: image
        integer, intent(in)         :: color
        integer                     :: blue

        call cgd_blue(image%ptr,color,blue)

    end function blue

    function alpha(image,color)
        type(base_img), intent(in)  :: image
        integer, intent(in)         :: color
        integer                     :: alpha

        call cgd_alpha(image%ptr,color,alpha)

    end function alpha

    function greyscale(image) result(color)
        type(base_img), intent(in)  :: image
        integer         :: color
        call CGD_GREYSCALE(image%ptr, color)
    end function greyscale

    function pixel_color(image,x,y)
        type(base_img), intent(in)  :: image
        integer, intent(in)         :: x,y
        integer                     :: pixel_color

        call cgd_pixel_color(image,xm(image,x),ym(image,y),pixel_color)

    end function pixel_color

    function green(image,color)
        type(base_img), intent(in)  :: image
        integer, intent(in)         :: color
        integer                     :: green

        call cgd_green(image%ptr,color,green)

    end function green

    function red(image,color)
        type(base_img), intent(in)  :: image
        integer, intent(in)         :: color
        integer                     :: red

        call cgd_red(image%ptr,color,red)

    end function red

    function within_bounds(image,x,y)
        type(base_img), intent(in)  :: image
        integer, intent(in)         :: x,y
        logical                     :: within_bounds
        integer                     :: safe

        call cgd_image_bounds_safe(image,xm(image,x),ym(image,y),safe)
        within_bounds = (safe /= 0)

    end function within_bounds

    function closest_color(image,r,g,b)
        type(base_img), intent(in)  :: image
        integer, intent(in)         :: r,g,b
        integer                     :: closest_color

        call cgd_image_color_closest(image%ptr,r,g,b,closest_color)

    end function closest_color

    function hwb_closest_color(image,r,g,b)
        type(base_img), intent(in)  :: image
        integer, intent(in)         :: r,g,b
        integer                     :: hwb_closest_color

        call cgd_image_color_closest_hwb(image%ptr,r,g,b,hwb_closest_color)

    end function hwb_closest_color

    function closest_alpha_color(image,r,g,b,a)
        type(base_img), intent(in)  :: image
        integer, intent(in)         :: r,g,b,a
        integer                     :: closest_alpha_color

        call cgd_image_color_closest_alpha(image%ptr,r,g,b,a,closest_alpha_color)

    end function closest_alpha_color

    function exact_color(image,r,g,b)
        type(base_img), intent(in)  :: image
        integer, intent(in)         :: r,g,b
        integer                     :: exact_color

        call cgd_image_color_exact(image%ptr,r,g,b,exact_color)

    end function exact_color

    function resolve_alpha_color(image,r,g,b,a)
        type(base_img), intent(in)  :: image
        integer, intent(in)         :: r,g,b,a
        integer                     :: resolve_alpha_color

        call cgd_image_color_resolve_alpha(image%ptr,r,g,b,a,resolve_alpha_color)

    end function resolve_alpha_color

    function resolve_color(image,r,g,b)
        type(base_img), intent(in)  :: image
        integer, intent(in)         :: r,g,b
        integer                     :: resolve_color

        call cgd_image_color_resolve(image%ptr,r,g,b,resolve_color)

    end function resolve_color

    function get_interlaced(image)
        type(base_img), intent(in)  :: image
        logical                     :: get_interlaced
        integer                     :: i
        call cgd_image_get_interlaced(image%ptr,i)
        get_interlaced = (i /= 0)
    end function get_interlaced


    function transparent_color(image)
        type(base_img), intent(in)  :: image
        integer                     :: transparent_color

        call cgd_image_get_transparent(image%ptr,transparent_color)

    end function transparent_color

    subroutine set_transparent_color(image,color)
        type(base_img), intent(in)  :: image
        integer, intent(in)         :: color

        call cgd_image_transparent(image%ptr,color)

    end subroutine set_transparent_color

    subroutine set_interlace(image,interlace)
        type(base_img), intent(in)  :: image
        logical, intent(in)         :: interlace
        integer                     :: i = 0

        if (interlace) i = 1
        call cgd_image_interlace(image%ptr,i)

    end subroutine set_interlace

    integer function number_of_colors(image)
        type(base_img), intent(in)  :: image

        call cgd_image_colors_total(image%ptr,number_of_colors)
    end function number_of_colors

    integer function brushed()
        call cgd_brushed(brushed)
    end function brushed

    integer function styled()
        call cgd_styled(styled)
    end function styled

    integer function styled_brushed()
        call cgd_styled_brushed(styled_brushed)
    end function styled_brushed

    integer function tiled()
        call cgd_tiled(tiled)
    end function tiled

    integer function pie()
        call cgd_pie(pie)
    end function pie

    integer function chord()
        call cgd_chord(chord)
    end function chord

    integer function no_fill()
        call cgd_nofill(no_fill)
    end function no_fill

    integer function edged()
        call cgd_edged(edged)
    end function edged

    integer function transparent()
        call cgd_transparent(transparent)
    end function transparent

    integer function anti_aliased()
        call cgd_anti_aliased(anti_aliased)
    end function anti_aliased

    subroutine set_clip_rectangle(image,x1,y1,x2,y2)
        type(base_img), intent(in)  :: image
        integer, intent(in)         :: x1,y1,x2,y2
        integer                     :: xx1,yy1,xx2,yy2

        if (image%fmode) then
            xx1 = xm(image,x1)
            yy1 = ym(image,y2)
            xx2 = xm(image,x2)
            yy2 = ym(image,y1)
        else
            xx1 = x1
            yy1 = y1
            xx2 = x2
            yy2 = y2
        endif
        call cgd_image_set_clip(image%ptr,xx1,yy1,xx2,yy2)

    end subroutine set_clip_rectangle

    subroutine get_clip_rectangle(image,x1,y1,x2,y2)
        type(base_img), intent(in)  :: image
        integer, intent(out)        :: x1,y1,x2,y2
        integer                     :: s

        call cgd_image_get_clip(image%ptr,x1,y1,x2,y2)
        if (image%fmode) then
            x1 = xm(image,x1)
            s  = xm(image,y2)
            x2 = xm(image,x2)
            y2 = xm(image,y1)
            y1 = s
        endif

    end subroutine get_clip_rectangle

    subroutine set_anti_aliased(image,color)
        type(base_img), intent(in)  :: image
        integer, intent(in)         :: color

        call cgd_image_set_aa(image%ptr,color)

    end subroutine set_anti_aliased

    subroutine set_anti_aliased_dont_blend(image,color,bcolor)
        type(base_img), intent(in)  :: image
        integer, intent(in)         :: color,bcolor

        call cgd_image_set_aa_nb(image%ptr,color,bcolor)

    end subroutine set_anti_aliased_dont_blend

    subroutine set_fmode(image,mode)
        type(base_img), intent(inout)  :: image
        logical, intent(in)            :: mode
        image%fmode = mode
    end subroutine set_fmode

    function get_fmode(image)
        type(base_img), intent(in)  :: image
        logical                     :: get_fmode
        get_fmode = image%fmode
    end function get_fmode

    elemental function xm(image,x)
        type(base_img), intent(in)  :: image
        integer, intent(in)         :: x
        integer                     :: xm
        xm = x
        if (image%fmode) xm = x + 1
    end function xm

    elemental function ym(image,y)
        type(base_img), intent(in)  :: image
        integer, intent(in)         :: y
        integer                     :: ym
        ym = y
        if (image%fmode) ym = image%height - y
    end function ym

end module simple_img
