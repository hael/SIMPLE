program simple_test_serialize
include 'simple_lib.f08'
use simple_image
use simple_stackops
use simple_ppca
implicit none
type(image)          :: img, img_msk, img_rev
real,    allocatable :: pcavec(:)
logical, allocatable :: l_mask(:,:,:)
integer, parameter   :: box = 256
call img%new([box,box,1], 1.0)
call img_rev%new([box,box,1], 1.0)
call img%square(60)
call img_msk%disc([box,box,1], 1.0, 80., l_mask)
call img_msk%vis
call img%vis
pcavec = img%serialize(l_mask)
call img_rev%unserialize(pcavec, l_mask)
call img_rev%vis
end program simple_test_serialize
