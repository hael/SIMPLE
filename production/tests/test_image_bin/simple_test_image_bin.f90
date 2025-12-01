program simple_test_image_bin
include 'simple_lib.f08'
use simple_image_bin
type(image_bin)      :: test_image_bin, test_image_bin_3D, ccimage
integer, allocatable :: imat_ccs(:,:,:)
integer              :: imat(4,4,1), imat_3D(3,3,2), nccs, i
real                 :: dist, dist_3D
! 2D test
print *, '2D TEST'
call test_image_bin%new_bimg(ldim=[4,4,1], smpd=1.)
imat(1,:,1) = [1,0,0,0]
imat(2,:,1) = [1,0,1,0]
imat(3,:,1) = [0,0,0,0]
imat(4,:,1) = [0,0,0,1]
call test_image_bin%set_imat(imat)
call test_image_bin%max_dist(dist)
print *, 'MAXIMUM DISTANCE FROM CENTER IS ', dist
call test_image_bin%find_ccs(ccimage)
call ccimage%get_nccs(nccs)
print *, nccs
call ccimage%get_imat(imat_ccs)
do i = 1,4
    print *, imat_ccs(i,:,1)
enddo
! extreme unit test case of 1
imat = 1
call test_image_bin%set_imat(imat)
call test_image_bin%max_dist(dist)
call test_image_bin%find_ccs(ccimage)
call ccimage%get_nccs(nccs)
print *, nccs
call ccimage%get_imat(imat_ccs)
do i = 1,4
    print *, imat_ccs(i,:,1)
enddo
! extreme unit test case of 0
imat = 0
call test_image_bin%set_imat(imat)
call test_image_bin%max_dist(dist)
call test_image_bin%find_ccs(ccimage)
call ccimage%get_nccs(nccs)
print *, nccs
call ccimage%get_imat(imat_ccs)
do i = 1,4
    print *, imat_ccs(i,:,1)
enddo
print *, ' '
! 3D test
print *, '3D TEST'
call test_image_bin_3D%new_bimg(ldim=[3,3,2],smpd=1.)
imat_3D = reshape( (/ 0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1 /), shape(imat_3D))
call test_image_bin_3D%set_imat(imat_3D)
call test_image_bin_3D%max_dist(dist_3D)
print *, 'MAXIMUM DISTANCE FROM CENTER IS ', dist_3D
end program simple_test_image_bin
