program simple_test_rad_med
!$ use omp_lib
!$ use omp_lib_kinds
include 'simple_lib.f08'
use simple_image,          only: image
use simple_radial_medians, only: radial_medians
implicit none

integer, parameter   :: BOX  = 256
real,    parameter   :: SMPD = 1.0, SIG = 40.
type(image)          :: img
type(stats_struct)   :: stats
type(radial_medians) :: rad_med
real, allocatable    :: medians(:)
integer :: irad

call img%new([BOX,BOX,1], SMPD)
call img%gauimg2D(SIG, SIG)
call img%vis

call rad_med%new([BOX,BOX,1])
allocate(medians(rad_med%get_rad_max()), source=0.)
call rad_med%calc_radial_medians(img, stats, medians)
do irad = 1,rad_med%get_rad_max()
    print *, irad, medians(irad)
end do



end program simple_test_rad_med
