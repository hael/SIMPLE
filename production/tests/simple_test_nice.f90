program simple_test_nice
use simple_core_module_api
use simple_nice
implicit none
type(simple_nice_comm) :: nice_comm
call nice_comm%init(1, "testserver")
call nice_comm%cycle()
call sleep(5)
nice_comm%view_micrographs%active         = .true.
nice_comm%view_micrographs%thumbnail%path = "/tmp/cls2D_thumbnail.jpeg"
nice_comm%view_micrographs%thumbnail%id   = 10 ! should be random
call nice_comm%cycle()
call sleep(5)
nice_comm%view_micrographs%active         = .true.
nice_comm%view_micrographs%thumbnail%path = "/tmp/cls2D_thumnail.jpeg"
nice_comm%view_micrographs%thumbnail%id   = 11! should be random
call nice_comm%cycle()
call sleep(10)
call nice_comm%terminate
end program simple_test_nice
