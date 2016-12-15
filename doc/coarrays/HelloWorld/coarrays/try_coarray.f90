program try_coarray
implicit none
real :: a[*]              ! Declare as coarray obj
real, codimension[*] :: b ! Or this way

! a and b below are local to the image
a = this_image()
b = this_image()*2

! Access a and b on other images
if( this_image() == 1 )then
    do img=1,num_images()
        print *, 'Image', this_image(), a[i], b[i]
    end do 
endif


end program try_coarray