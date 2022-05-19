program simple_test_bounds_from_mask3D
include 'simple_lib.f08'
use simple_image
implicit none
integer, parameter :: BOX=256, RAD=50
real,    parameter :: SMPD=1.
type(image) :: cube
logical, allocatable :: mask(:,:,:)
integer :: lb(3), ub(3)
call cube%new([BOX,BOX,BOX], SMPD)
call cube%square(RAD)
call cube%write('cube.mrc')
mask = cube%bin2logical()
call bounds_from_mask3D(mask, lb, ub)
print *, lb(1), ub(1)
print *, lb(2), ub(2)
print *, lb(3), ub(3)
end program simple_test_bounds_from_mask3D
