program simple_test_gaussian1D
include 'simple_lib.f08'
implicit none
integer, parameter :: N = 10000
real    :: x, dx, gauss(N), sig2, total_energy, x_max, g_max, g_min
integer :: cnt, i
sig2  = 1.
x_max = sqrt(sig2)
x     = 0.
dx    = 0.01
cnt   = 0
gauss = 0.
do while( x <= x_max )
    cnt        = cnt + 1
    gauss(cnt) = gaussian1D(x, avg=0., sigma_sq=sig2)
    x          = x + dx
enddo
gauss        = gauss / sum(gauss(1:cnt))
total_energy = sum(gauss(1:cnt))
print *, total_energy , gauss(1)
! testing multinomal
do i = 1, 10
    print *, multinomal( gauss(1:cnt) )
enddo
end program simple_test_gaussian1D
