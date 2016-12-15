program hello_image
implicit none

if( this_image() == 1 )then
    read(*,*) 
endif

sync all

write(*,*) "Hello from image ", this_image(), "out of ", num_images(), " total images"

if( this_image() == 1 )then
    read(*,*) 
endif

sync all

end program hello_image