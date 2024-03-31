program simple_test_gaussian1D
include 'simple_lib.f08'
implicit none
integer, parameter :: N = 10000
real    :: x, dx, gauss(N), sig2, total_energy, x_max, g_max, g_min, smpl_sh(2), cur_sh(2), sh_thres
integer :: cnt, i, which
sig2  = 1.
x_max = sqrt(sig2)
x     = 0.
dx    = 0.1
cnt   = 0
gauss = 0.
do while( x <= x_max )
    cnt        = cnt + 1
    gauss(cnt) = gaussian1D(x, avg=0., sigma_sq=sig2)
    x          = x + dx
enddo
gauss        = gauss / sum(gauss(1:cnt))
gauss(1:cnt) = (1. - gauss(1:cnt))
gauss        = gauss / sum(gauss(1:cnt))
total_energy = sum(gauss(1:cnt))
print *, total_energy
print *, gauss(1:cnt)
! testing multinomal
do i = 1, 10
    which = greedy_sampling( gauss(1:cnt), cnt )
    print *, which
enddo
! testing shift_sampling
cur_sh   = [0.6, -1.7]
sh_thres = 1.5
do i = 1, 10
    smpl_sh = shift_sampling( cur_sh, sh_thres )
    print *, 'cur_sh = ', cur_sh, '; sh_smpl = ', smpl_sh
    if( any(smpl_sh < cur_sh-sh_thres) .or. any(smpl_sh > cur_sh+sh_thres) )then
        print *, 'UNIT TEST FAILED'
        stop
    endif
enddo
print *, 'UNIT TEST PASSED'
contains
    ! shift multinomal sampling within a threshold, units are in Ang
    function shift_sampling( cur_sh, thres_bound ) result(sh)
        real,    intent(in)  :: cur_sh(2)
        real,    intent(in)  :: thres_bound
        integer, parameter   :: N_SMPL  = 100
        real,    parameter   :: SH_SAFE = .5
        real,    allocatable :: vals(:)
        integer :: i, sh_signs(2), which, dim, n_safe
        real    :: sh(2)
        real    :: d_sh, gauss_sh(N_SMPL), sh_vals(N_SMPL), sig2, d_thres, thres
        thres = thres_bound
        sh    = cur_sh
        if( thres < TINY ) return
        ! randomly pick the plus/minus for each x,y dimensions
        sh_signs = 1
        if( ran3() < 0.5 ) sh_signs(1) = -1
        if( ran3() < 0.5 ) sh_signs(2) = -1
        ! sampling for x
        d_sh   = thres / real(N_SMPL-1)
        n_safe = 1
        if( thres > SH_SAFE ) n_safe = 1 + floor(SH_SAFE / d_sh)
        sig2   = thres**2.
        do dim = 1, 2
            do i = 1, N_SMPL
                d_thres     = d_sh * real(i - 1)
                gauss_sh(i) = gaussian1D(d_thres, avg=0., sigma_sq=sig2)
                sh_vals(i)  = d_thres
            enddo
            gauss_sh = (1. - gauss_sh)
            gauss_sh = gauss_sh / sum(gauss_sh)
            which    = greedy_sampling( gauss_sh, N_SMPL, n_safe )
            sh(dim)  = cur_sh(dim) + real(sh_signs(dim)) * sh_vals(which)
        enddo
    end function shift_sampling
end program simple_test_gaussian1D
