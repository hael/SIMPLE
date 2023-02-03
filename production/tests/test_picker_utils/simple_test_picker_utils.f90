program simple_test_picker_utils
include 'simple_lib.f08'
use simple_picker_utils, only: picker_utils
use simple_image,        only: image
implicit none

character(len=*), parameter :: micname = '/home/elmlundho/cache/relion_tut/2_motion_correct/20170629_00021_frameImage_intg.mrc'
real,             parameter :: SMPD    = 0.885, MOLDIAM = 180.
type(image)        :: micimg
type(picker_utils) :: putils
integer            :: ldim(3), ifoo

call find_ldim_nptcls(micname, ldim, ifoo)
call micimg%new(ldim, SMPD)
call micimg%read(micname)
call putils%set_mics(micimg, SMPD)
call putils%bin_mic_shrink1([0.,20.])

end program simple_test_picker_utils
