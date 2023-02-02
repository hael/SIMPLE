program simple_test_ori
include 'simple_lib.f08'
use simple_sym
use simple_ori
implicit none
type(ori) :: o_truth
type(sym) :: pgrpsyms
real :: vec(3), angle
call pgrpsyms%new('c1')
call o_truth%new(.true.)
call pgrpsyms%rnd_euler(o_truth)
print *, o_truth%get_mat()
call o_truth%get_axis_angle( vec, angle )
print *, vec
print *, angle/pi*180.
end program simple_test_ori
