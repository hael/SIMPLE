! % time 2D:    72.097917319251053     
!  Pearson correlation vs. fast 2D:    1.00000000    
!  Absolute difference vs. fast 2D:    2.45080600E-09
!  % time 3D:    59.664005833145232     
!  Pearson correlation vs. fast 3D:    1.00000000    
!  Absolute difference vs. fast 3D:    1.23470834E-09
program simple_test_kbinterpol_fast
use simple_core_module_api
use simple_test_utils, only: assert_true, assert_int
use simple_kbinterpol, only: kbinterpol
use simple_timer
implicit none

type(kbinterpol)      :: kbwin
real,       parameter :: ALPHA = 2.0_sp    !< oversampling factor
integer,    parameter :: BOX   = 256
real,    allocatable  :: kbw(:,:), kbw_fast(:,:), kbw3D(:,:,:), kbw3D_fast(:,:,:)
integer :: h, k, iwinsz, wdim, i, win(2,2), cnt, win3D(2,3)
real    :: loc(2), d(2), kbw_diff, kbw_corr, rmat(2,2), rmat3D(3,3), dists(3), loc3D(3), w(2)
real(timer_int_kind)    :: rt_kb, rt_kb_fast
integer(timer_int_kind) :: t_kb, t_kb_fast

call seed_rnd

! intepolation windows
kbwin = kbinterpol(RECWINSZ,ALPHA)
wdim  = kbwin%get_wdim()

allocate(kbw(wdim,wdim), kbw_fast(wdim,wdim), kbw3D(wdim,wdim,wdim), kbw3D_fast(wdim,wdim,wdim), source=1.)

! rotation & scale
call rotmat2D(ran3()*180., rmat)
rmat   = ALPHA * rmat
iwinsz = ceiling(RECWINSZ - 0.5)

! time original approach
t_kb = tic()
do k = -BOX,0
    do h = -BOX,BOX
        loc = matmul(real([h,k]), rmat)
        ! window
        win(1,:) = nint(loc)
        win(2,:) = win(1,:) + iwinsz
        win(1,:) = win(1,:) - iwinsz                   
        ! kernel
        kbw  = 1.
        do i = 1,wdim
            d = real(win(1,:) + i - 1, sp) - loc
            kbw(i,:) = kbw(i,:) * kbwin%apod(d(1))
            kbw(:,i) = kbw(:,i) * kbwin%apod(d(2))
        enddo
        kbw   = kbw / sum(kbw)
    end do
end do
rt_kb = toc(t_kb)

! time the new appraoch
t_kb_fast = tic()
do k = -BOX,0
    do h = -BOX,BOX
        loc = matmul(real([h,k]), rmat)          
        ! fast kernel
        call kbwin%apod_mat_2d(loc, iwinsz, wdim, kbw_fast)
    enddo
enddo
rt_kb_fast = toc(t_kb_fast)

print *, '% time 2D: ', 100. * rt_kb / rt_kb_fast 

! validate correctness
kbw_diff = 0.
kbw_corr = 0.
cnt      = 0
do k = -BOX,0
    do h = -BOX,BOX
        loc = matmul(real([h,k]), rmat)
        ! window
        win(1,:) = nint(loc)
        win(2,:) = win(1,:) + iwinsz
        win(1,:) = win(1,:) - iwinsz                   
        ! kernel
        kbw  = 1.
        do i = 1,wdim
            d = real(win(1,:) + i - 1, sp) - loc
            kbw(i,:) = kbw(i,:) * kbwin%apod(d(1))
            kbw(:,i) = kbw(:,i) * kbwin%apod(d(2))
        enddo
        kbw   = kbw / sum(kbw)
        ! fast kernel
        call kbwin%apod_mat_2d(loc, iwinsz, wdim, kbw_fast)
        ! accumulate absolute differences
        kbw_diff = kbw_diff + sum(abs(kbw - kbw_fast)) / real(wdim * wdim, sp)
        kbw_corr = kbw_corr + pearsn(kbw, kbw_fast)
        cnt = cnt + 1
    enddo
enddo

print *, 'Pearson correlation vs. fast 2D: ', kbw_corr / real(cnt, sp)

print *, 'Absolute difference vs. fast 2D: ', kbw_diff / real(cnt, sp)

! rotation & scale
call rnd_romat(rmat3D)
rmat3D = ALPHA * rmat3D
iwinsz = ceiling(RECWINSZ - 0.5)

! time original approach
t_kb  = tic()
do k = -BOX,0
    do h = -BOX,BOX
        loc3D = matmul(real([h,k,0]), rmat3D)
        ! window
        win3D(1,:) = nint(loc3D)
        win3D(2,:) = win3D(1,:) + iwinsz
        win3D(1,:) = win3D(1,:) - iwinsz
        ! kernel
        kbw3D = 1.
        do i=1,wdim
            dists    = real(win3D(1,:) + i - 1) - loc3D
            kbw3D(i,:,:) = kbw3D(i,:,:) * kbwin%apod(dists(1))
            kbw3D(:,i,:) = kbw3D(:,i,:) * kbwin%apod(dists(2))
            kbw3D(:,:,i) = kbw3D(:,:,i) * kbwin%apod(dists(3))
        enddo
        kbw3D = kbw3D / sum(kbw3D)
    end do
end do
rt_kb = toc(t_kb)

! time the fast approach
t_kb_fast  = tic()
do k = -BOX,0
    do h = -BOX,BOX
        loc3D = matmul(real([h,k,0]), rmat3D)
        ! fast kernel
        call kbwin%apod_mat_3d(loc3D, iwinsz, wdim, kbw3D_fast)
    end do
end do
rt_kb_fast = toc(t_kb_fast)

print *, '% time 3D: ', 100. * rt_kb / rt_kb_fast 

! validate correctness
kbw_diff   = 0.
kbw_corr   = 0.
cnt        = 0
do k = -BOX,0
    do h = -BOX,BOX
        loc3D = matmul(real([h,k,0]), rmat3D)
        ! window
        win3D(1,:) = nint(loc3D)
        win3D(2,:) = win3D(1,:) + iwinsz
        win3D(1,:) = win3D(1,:) - iwinsz
        ! kernel
        kbw3D = 1.
        do i=1,wdim
            dists    = real(win3D(1,:) + i - 1) - loc3D
            kbw3D(i,:,:) = kbw3D(i,:,:) * kbwin%apod(dists(1))
            kbw3D(:,i,:) = kbw3D(:,i,:) * kbwin%apod(dists(2))
            kbw3D(:,:,i) = kbw3D(:,:,i) * kbwin%apod(dists(3))
        enddo
        kbw3D = kbw3D / sum(kbw3D)
        ! fast kernel
        call kbwin%apod_mat_3d(loc3D, iwinsz, wdim, kbw3D_fast)
        ! accumulate absolute differences
        kbw_diff = kbw_diff + sum(abs(kbw3D - kbw3D_fast)) / real(wdim * wdim * wdim, sp)
        ! accumulate corrs
        kbw_corr = kbw_corr + pearsn(kbw3D, kbw3D_fast)
        cnt = cnt + 1
    end do
end do

print *, 'Pearson correlation vs. fast 3D: ', kbw_corr / real(cnt, sp)

print *, 'Absolute difference vs. fast 3D: ', kbw_diff / real(cnt, sp)


contains

    !>  \brief  for generating a random rotation matrix
    subroutine rnd_romat( rmat )
        ! Fast Random Rotation Matrices, Arvo, Graphics Gems III, 1992
        real, intent(out) :: rmat(3,3)
        real :: theta, phi, z, vx, vy, vz
        real :: r, st, ct, sx, sy
        ! init
        theta = ran3()*TWOPI
        phi   = ran3()*TWOPI
        z     = ran3()*2.
        ! V
        r  = sqrt( z )
        vx = r*sin( phi )
        vy = r*cos( phi )
        vz = sqrt( 2.-z )
        ! S=Vt*R; sz=vz
        st = sin( theta )
        ct = cos( theta )
        sx = vx*ct - vy*st
        sy = vx*st + vy*ct
        ! M
        rmat(1,1) = vx * sx - ct
        rmat(1,2) = vx * sy - st
        rmat(1,3) = vx * vz
        rmat(2,1) = vy * sx + st
        rmat(2,2) = vy * sy - ct
        rmat(2,3) = vy * vz
        rmat(3,1) = vz * sx
        rmat(3,2) = vz * sy
        rmat(3,3) = 1. - z
    end subroutine rnd_romat

end program simple_test_kbinterpol_fast