program simple_test_bounds_from_mask3D
use simple_core_module_api
use simple_math,  only: bounds_from_mask3D
use simple_image, only: image
implicit none
type(image)          :: cube
integer, parameter   :: BOX=256, RAD=50
real,    parameter   :: SMPD=1.
logical, allocatable :: mask(:,:,:)
integer              :: lb(3), ub(3)
call cube%new([BOX,BOX,BOX], SMPD)
call cube%square(RAD)
call cube%write(string('cube.mrc'))
mask = cube%bin2logical()
call bounds_from_mask3D(mask, lb, ub)
print *, lb(1), ub(1)
print *, lb(2), ub(2)
print *, lb(3), ub(3)
end program simple_test_bounds_from_mask3D
