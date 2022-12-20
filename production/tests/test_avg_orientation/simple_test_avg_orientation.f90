program simple_test_avg_orientation
include 'simple_lib.f08'
use simple_sym
use simple_ori
use simple_oris
implicit none
integer,   parameter      :: N = 10
type(ori), allocatable    :: o(:)
real      :: rot_mats(3,3,N)
type(sym) :: pgrpsyms
type(ori) :: o_avg
integer   :: i, j, geodesic_type
real      :: mat_avg(3,3), weights(N), total
call pgrpsyms%new('c1')
allocate(o(N))
do i = 1, N
    call o(i)%new(.true.)
enddo
call o_avg%new(.true.)
! test with same rotation matrices
call pgrpsyms%rnd_euler(o(1))
do i = 2, N
    o(i) = o(1)
enddo
do i = 1, N
    rot_mats(:,:,i) = o(i)%get_mat()
    weights(i)      = 1.
enddo
geodesic_type = 0
call geodesic_opt_ori(o, o_avg, weights, geodesic_type)
mat_avg = o_avg%get_mat()
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
call geodesic_opt_ori(o, o_avg, weights, geodesic_type)
total = 0.
do i = 1, N
    total = total + geodesic_dist_trace(o_avg, o(i))
    print *, 'diff angle between "average" mat and mat ', int2str(i), ' is:', geodesic_dist_trace(o_avg, o(i))
enddo
print *, 'total difference between "average" mat and others is: ', total/N
print *, '-------'
do i = 1, N
    total = 0.
    do j = 1, N
        if( j /= i ) total = total + geodesic_dist_trace(o(i), o(j))
    enddo
    print *, 'total difference between mat ', int2str(i),' and others is: ', total/(N-1.) ! distance of one against itself does not count
enddo
end program simple_test_avg_orientation