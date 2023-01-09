program simple_test_otsu
include 'simple_lib.f08'
implicit none
integer, parameter :: NSAMP_D1 = 200, NSAMP_D2 = 500, NSAMP_TOT = 700
real,    parameter :: MEAN_D1 = 2., MEAN_D2 = 2., SDEV_D1 = 2., SDEV_D2 = 2.
real,  allocatable :: dat1(:), dat2(:)
real    :: vec(NSAMP_TOT), ave, sdev, var, d, prob
integer :: cnt, i, ncorrect
logical :: err, mask_facit(NSAMP_TOT), mask(NSAMP_TOT)
cnt = 0
do i = 1, NSAMP_D1
    cnt = cnt + 1
    vec(cnt) = gasdev(MEAN_D1, SDEV_D1)
    mask_facit(cnt) = .true.
end do
do i = 1, NSAMP_D2
    cnt = cnt + 1
    vec(cnt) = gasdev(MEAN_D1, SDEV_D1)
    mask_facit(cnt) = .false.
end do
call otsu(NSAMP_TOT, vec, mask)
call moment(vec, ave, sdev, var, err, mask)
print *, 'distribution 1 ave/sdev: ', ave, sdev
call moment(vec, ave, sdev, var, err, .not. mask)
print *, 'distribution 2 ave/sdev: ', ave, sdev
ncorrect = max(count(mask.eqv.mask_facit),count(mask.eqv..not.mask_facit))
print *, '% correct assigned: ', 100. * (real(ncorrect) / real(NSAMP_TOT))
end program simple_test_otsu
