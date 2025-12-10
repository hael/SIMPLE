program simple_test_nice
include 'simple_lib.f08'
use simple_nice
implicit none
type(simple_nice_communicator) :: nice_communicator
call nice_communicator%init(1, "testserver")
call nice_communicator%cycle()
call sleep(5)
nice_communicator%view_micrographs%active         = .true.
nice_communicator%view_micrographs%thumbnail%path = "/tmp/cls2D_thumbnail.jpeg"
nice_communicator%view_micrographs%thumbnail%id   = 10 ! should be random
call nice_communicator%cycle()
call sleep(5)
nice_communicator%view_micrographs%active         = .true.
nice_communicator%view_micrographs%thumbnail%path = "/tmp/cls2D_thumnail.jpeg"
nice_communicator%view_micrographs%thumbnail%id   = 11! should be random
call nice_communicator%cycle()
call sleep(10)
call nice_communicator%terminate
end program simple_test_nice
