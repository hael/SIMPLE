program simple_test_avg_rotmat
include 'simple_lib.f08'
use simple_sym
use simple_ori
use simple_math
implicit none
integer,   parameter      :: N = 10
type(ori), allocatable    :: o(:)
real      :: rot_mats(3,3,N)
type(sym) :: pgrpsyms
integer   :: i, j
real      :: mat_avg(3,3), weights(N), total
call pgrpsyms%new('c1')
allocate(o(N))
do i = 1, N
    call o(i)%new(.true.)
enddo
! test with same rotation matrices
call pgrpsyms%rnd_euler(o(1))
do i = 2, N
    o(i) = o(1)
enddo
do i = 1, N
    rot_mats(:,:,i) = o(i)%get_mat()
    weights(i)      = 1.
enddo
call avg_rotmat(rot_mats, mat_avg, weights)
if( any(abs(mat_avg - rot_mats(1:3, 1:3, 1)) > epsilon(weights(1))) )then
    print *, 'FAIL! same rotation matrix case FAILED!'
    print *, rot_mats(:,:,1)
    print *, mat_avg
else
    print *, 'PASS! same rotation matrix case PASSED!'
endif
! test with random rotation matrices
call pgrpsyms%rnd_euler(o(1))
do i = 2, N
    call pgrpsyms%rnd_euler(o(i))
enddo
do i = 1, N
    rot_mats(:,:,i) = o(i)%get_mat()
    weights(i)      = 1.
enddo
call avg_rotmat(rot_mats, mat_avg, weights)
total = 0.
do i = 1, N
    total = total + rot_angle(mat_avg, rot_mats(:,:,i))
    print *, 'diff angle between "average" mat and mat ', int2str(i), ' is:', rot_angle(mat_avg, rot_mats(:,:,i))
enddo
print *, 'total difference between "average" mat and others is: ', total
print *, '-------'
do i = 1, N
    total = 0.
    do j = 1, N
        total = total + rot_angle(rot_mats(:,:,i), rot_mats(:,:,j))
        !print *, 'diff angle between mat ', int2str(i),' and mat ', int2str(j), ' is:', rot_angle(rot_mats(:,:,i), rot_mats(:,:,j))
    enddo
    print *, 'total difference between mat ', int2str(i),' and others is: ', total
enddo
end program simple_test_avg_rotmat
