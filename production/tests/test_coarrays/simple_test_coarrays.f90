program simple_test_coarrays
use simple_stream_module_api   
use simple_image,   only: image
use simple_imgfile, only: imgfile
implicit none
integer     :: ldim(3), i, j
real        :: smpd
type(image) :: img
integer     :: my_rank, num_procs
integer     :: a[*]
! Define a coarray variable
integer, codimension[*] :: counter
! Get the rank (ID) and number of processors
my_rank   = this_image()
num_procs = num_images()
if( this_image() == 1 )then
    ldim = [120,120,1]
    call img%new(ldim, smpd)
    call img%square(20)
    call img%write(string('squares_mrc.mrc'),1)
endif
! Write out a message from each rank
print *, "Hello from processor    ", my_rank,"out of",num_procs
! Synchronize to ensure all 'write' executions are done
sync all
! Increment the counter on the first rank
if( my_rank == 1 ) counter[1] = counter[1] + 1
! Again, ensure all processes are synchronized
sync all
! Print the incremented value from rank 1
if( my_rank == 1 )then
    print *, "The counter value is now", counter[1]
endif
a = 0
if( this_image() == 1 ) then
    a = 1
    print *, 'Image ', this_image(), 'has a value', a
    print *, 'Image ', this_image(), 'sending new value to image 2.'
    a[2] = 2 * a
endif
sync all
if( this_image() == 2 ) then
    a = 1
    print *, 'Image ', this_image(), 'now has a value', a
    print *, 'Image ', this_image(), 'sending new value to image 1.'
    a[1] = 2 * a
endif
sync all
print *, 'Image ', this_image(), &
  'sees that image 1 now has a value ', a[1]
sync all
call sleep(10)
sync all
end program simple_test_coarrays 
