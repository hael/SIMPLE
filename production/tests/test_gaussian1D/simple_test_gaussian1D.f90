program simple_test_gaussian1D
include 'simple_lib.f08'
use simple_eul_prob_tab, only: shift_sampling
implicit none
integer, parameter :: N = 10000
real    :: x, dx, gauss(N), sig2, total_energy, x_max, g_max, g_min, smpl_sh(2), cur_sh(2), sh_thres
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
! testing shift_sampling
cur_sh   = [0.6, -1.7]
sh_thres = 1.5
do i = 1, 10
    smpl_sh = shift_sampling( cur_sh, sh_thres )
    print *, smpl_sh
    if( any(smpl_sh < cur_sh-sh_thres) .or. any(smpl_sh > cur_sh+sh_thres) )then
        print *, 'UNIT TEST FAILED'
        stop
    endif
enddo
print *, 'UNIT TEST PASSED'
end program simple_test_gaussian1D
