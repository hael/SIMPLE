program simple_test_coarrays
implicit none
integer :: my_image, nimages, syncstat
my_image = this_image()
nimages  = num_images()
sync all(stat=syncstat)
if( syncstat /= 0 ) error stop 'simple_test_coarrays initial sync failed'
if( my_image == 1 ) write(*,'(A,I0,A)') 'simple_test_coarrays passed with ', nimages, ' images'
sync all(stat=syncstat)
if( syncstat /= 0 ) error stop 'simple_test_coarrays final sync failed'
end program simple_test_coarrays
