!*******************************************************************************
! Small test program for the fortran module interfacing the GD library.
! Make sure the variable character 'font' refers to an existing font file!
!*******************************************************************************

!! Modified by Michael Eager, Feb 2018
program simple_test_img_png
include 'simple_lib.f08'
use simple_img
implicit none
type(base_img)            :: image
integer                   :: c1,c2,c3,c4,c5
integer                   :: a1(4),a2(4)
real(rp)                  :: ptsize = 30.0_rp
real(rp)                  :: angle = 0.0_rp
integer                   :: w,h
integer                   :: r,x,y,status
real, allocatable         :: buffer(:,:)
character(len=256)        :: fout = 'gray-buffer.png'
a1(:) = (/ 5,47,12,55 /)
a2(:) = (/ 36,60,2,25 /)


 call create_img(64,64,image)
 call allocate_color(image,0,0,0,c1)
 call allocate_color(image,255,255,255,c2)
 call allocate_color(image,255,0,0,c3)
 call allocate_color(image,255,0,255,c4)
 call allocate_color(image,0,255,0,c5)

 call draw_line(image,1,1,64,64,c2)
 call draw_rectangle(image,10,10,53,53,c2)
 call draw_filled_rectangle(image,30,30,43,53,c3)
 call draw_filled_polygon(image,4,a1,a2,c4)
 call fill_img(image,1,5,c5)
 call write_img_as_png(image,"xxx.png")
 call destroy_img(image)
 write(*,*) 'test png'
 status=0
 call read_png("gray.png", buffer, status)
 if (status /= 0) write(*,*) "read_png failed"
 if (.not. allocated(buffer)) then
     write(*,*) "read_png failed, buffer not allocated"
 else
     print*, "read_png width ", size(buffer,1), " height ",  size(buffer,2)
     print*, "read_png max val ", MAXVAL(buffer)
     print *, buffer
 end if

 call write_png(buffer, fout, status)
 if (status /= 0) write(*,*) "write_png failed"
 deallocate(buffer)
 write(*,*) 'test read/write png from buffers'

 call create_img_from_png("bbb",image)
 call write_img_as_png(image,"bbb-from-png.png")
 call destroy_img(image)
 write(6,*) 'test png to png'

#ifndef WITHOUT_JPEG
 call create_img_from_jpeg("bbb.jpg",image)
 call write_img_as_png(image,"bbb-from-jpeg.png")
 call write_img_as_png(image,"bbb-from-jpeg.jpg")
 call destroy_img(image)
 write(*,*) 'test jpg'
#endif



#ifndef WITHOUT_XPM
 call create_img_from_xpm("bbb.xpm",image)
 call write_img_as_png(image,"bbb-from-xpm.png")
 call destroy_img(image)
 write(*,*) 'test xpm'
#endif

end program
