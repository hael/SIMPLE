program simple_test_otsu
include 'simple_lib.f08'
implicit none
integer, parameter :: NSAMP = 300, NPEAKS = 50
real,    parameter :: MEAN = 1, SDEV_P = 10
real,  allocatable :: dat1(:), dat2(:)
real    :: vec(NSAMP), ave, sdev, var, avg1, avg2, d
integer :: i, n1, n2
logical :: mask(NSAMP), err
do i = 1, NSAMP
    vec(i) = gasdev(MEAN, SDEV_P)
end do
call otsu(NSAMP, vec, mask)
dat1 = pack(vec, mask=     mask)
call moment(dat1, ave, sdev, var, err )
print *, 'dat1 stats avg/sdev: ', ave, sdev
dat2 = pack(vec, mask=.not.mask)
call moment(dat2, ave, sdev, var, err )
print *, 'dat2 stats avg/sdev: ', ave, sdev
n1   = size(dat1)
n2   = size(dat2)
avg1 = sum(dat1) / real(n1)
avg2 = sum(dat2) / real(n2)
d = abs(avg1 - avg2)
print *, 'd = ', d ! increases with variance


end program simple_test_otsu
