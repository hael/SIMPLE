!@descr: for all numerics tests
module simple_commanders_test_numerics
use simple_commanders_api
implicit none
#include "simple_local_flags.inc"

type, extends(commander_base) :: commander_test_eigh_test
  contains
    procedure :: execute      => exec_test_eigh_test
end type commander_test_eigh_test

type, extends(commander_base) :: commander_test_kbinterpol_fast
  contains
    procedure :: execute      => exec_test_kbinterpol_fast
end type commander_test_kbinterpol_fast

type, extends(commander_base) :: commander_test_maxnloc_test
  contains
    procedure :: execute      => exec_test_maxnloc_test
end type commander_test_maxnloc_test

type, extends(commander_base) :: commander_test_neigh
  contains
    procedure :: execute      => exec_test_neigh
end type commander_test_neigh

contains

subroutine exec_test_eigh_test( self, cline )
    use simple_linalg
    class(commander_test_eigh_test),    intent(inout) :: self
    class(cmdline),                     intent(inout) :: cline
    integer, parameter :: N = 5, N_EIGS = 3, N_L = 15000, N_EIGS_L = 12
    integer            :: i
    real               :: mat(N, N), eigvals(N_EIGS), eigvecs(N,N_EIGS), tmp(N,N), eig_vals(N), mat_ori(N, N)
    mat(:,1) = [ 0.67,-0.20, 0.19,-1.06, 0.46]
    mat(:,2) = [-0.20, 3.82,-0.13, 1.06,-0.48]
    mat(:,3) = [ 0.19,-0.13, 3.27, 0.11, 1.10]
    mat(:,4) = [-1.06, 1.06, 0.11, 5.86,-0.98]
    mat(:,5) = [ 0.46,-0.48, 1.10,-0.98, 3.54]
    mat_ori  = mat
    print *, 'SVDCMP result:'
    call svdcmp(mat, eig_vals, tmp)
    print *, eig_vals
    print *, 'EIGH result:'
    call eigh( N, mat_ori, N_EIGS, eigvals, eigvecs )
    print *, 'Selected eigenvalues', eigvals
    print *, 'Selected eigenvectors (stored columnwise)'
    do i = 1, N
        print *, eigvecs(i,:)
    enddo
    call test_eigh(N_L, N_EIGS_L)
    call simple_end('**** SIMPLE_TEST_EIGH_TEST_WORKFLOW NORMAL STOP ****')
end subroutine exec_test_eigh_test

subroutine exec_test_kbinterpol_fast( self, cline )
    use simple_test_utils, only: assert_true, assert_int
    use simple_kbinterpol, only: kbinterpol
    use simple_timer
    class(commander_test_kbinterpol_fast),    intent(inout) :: self
    class(cmdline),                     intent(inout) :: cline
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
    kbwin = kbinterpol(KBWINSZ,ALPHA)
    wdim  = kbwin%get_wdim()
    allocate(kbw(wdim,wdim), kbw_fast(wdim,wdim), kbw3D(wdim,wdim,wdim), kbw3D_fast(wdim,wdim,wdim), source=1.)
    ! rotation & scale
    call rotmat2D(ran3()*180., rmat)
    rmat   = ALPHA * rmat
    iwinsz = ceiling(KBWINSZ - 0.5)
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
    iwinsz = ceiling(KBWINSZ - 0.5)
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
    call simple_end('**** SIMPLE_TEST_KBINTERPOL_FAST_WORKFLOW NORMAL STOP ****')

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

end subroutine exec_test_kbinterpol_fast

subroutine exec_test_maxnloc_test( self, cline )
    class(commander_test_maxnloc_test), intent(inout) :: self
    class(cmdline),                     intent(inout) :: cline
    integer, parameter :: NNRS = 1000, NSEL = 10, NTST=100000
    real               :: arr(NNRS), arr_copy(NNRS)
    integer            :: indxarr(NNRS), i, loc(NSEL)
    type(ran_tabu)     :: rt
    rt = ran_tabu(NNRS)
    do i = 1, NNRS
        arr(i) = real(i)
    end do
    print *, 'testing maxnnloc'
    call rt%shuffle(arr)
    loc = maxnloc(arr, NSEL)
    arr_copy = arr
    indxarr = (/(i,i=1,NNRS)/)
    call hpsort(arr, indxarr)
    call reverse(indxarr)
    do i=1,NSEL
        print *, i, arr_copy(indxarr(i)), arr_copy(loc(i))
    end do
    print *, ''
    print *, 'testing minnnloc'
    call rt%shuffle(arr)
    loc = minnloc(arr, NSEL)
    arr_copy = arr
    indxarr = (/(i,i=1,NNRS)/)
    call hpsort(arr, indxarr)
    do i=1,NSEL
        print *, i, arr_copy(indxarr(i)), arr_copy(loc(i))
    end do
    call simple_end('**** SIMPLE_TEST_MAXNLOC_TEST_WORKFLOW NORMAL STOP ****')
end subroutine exec_test_maxnloc_test

subroutine exec_test_neigh( self, cline )
    class(commander_test_neigh),    intent(inout) :: self
    class(cmdline),                 intent(inout) :: cline
    character(len=*), parameter :: PGRP       = 'c2'
    integer,          parameter :: NSPACE     = 20000
    integer,          parameter :: NSPACE_SUB = 500
    logical            :: lnns(NSPACE), lclosest(NSPACE)
    real               :: nnn(NSPACE_SUB), angres, tot_nnn, euldists(NSPACE_SUB), inplrotdist
    integer            :: i, closest, subspace_inds(NSPACE_SUB)
    type(sym)          :: pgrpsym 
    type(ori)          :: o, o_close, osym
    type(stats_struct) :: nnn_stats, euldist_stats
    type(oris)         :: eulspace, eulspace_sub
    call pgrpsym%new(PGRP)
    call eulspace    %new(NSPACE,     is_ptcl=.false.)
    call eulspace_sub%new(NSPACE_SUB, is_ptcl=.false.)
    call pgrpsym%build_refspiral(eulspace)
    call pgrpsym%build_refspiral(eulspace_sub)
    lclosest = .false.
    write(*,'(a)') '>>> OBTAINING STATISTICS ABOUT THE SUBSPACE'
    do i = 1,NSPACE_SUB
        call progress_gfortran(i, NSPACE_SUB)
        call eulspace_sub%get_ori(i, o)
        closest           = pgrpsym%find_closest_proj(eulspace, o)
        call eulspace%get_ori(closest, o_close)
        call pgrpsym%sym_dists(o, o_close, osym, euldists(i), inplrotdist)
        lclosest(closest) = .true.
        subspace_inds(i)  = closest
    end do
    call calc_stats(euldists, euldist_stats)
    write(*,'(a,1x,f15.6)') '# UNIQUE PROJECTIONS IDENTIFIED :', real(count(lclosest))
    write(*,'(a,1x,f15.6)') '% UNIQUE PROJECTIONS IDENTIFIED :', (real(count(lclosest)) / real(NSPACE_SUB)) * 100.
    write(*,'(a)') '>>> DISTANCE STATISTICS FOR CORRESPONDING PROJECTION DIRECTIONS'
    write(*,'(a,1x,f15.6)') 'AVERAGE EULDIST (IN DEGREES) :', euldist_stats%avg
    write(*,'(a,1x,f15.6)') 'SDEV    EULDIST (IN DEGREES) :', euldist_stats%sdev
    write(*,'(a,1x,f15.6)') 'MINIMUM EULDIST (IN DEGREES) :', euldist_stats%minv
    write(*,'(a,1x,f15.6)') 'MAXIMUM EULDIST (IN DEGREES) :', euldist_stats%maxv
    write(*,'(a)') '>>> OBTAINING STATISTICS ABOUT THE NEAREST NEIGHBORS'
    angres = eulspace_sub%find_angres()
    write(*,'(a,1x,f15.6)') 'ANGULAR RESOLUTION OF SUBSPACE (IN DEGREES) :', angres
    tot_nnn = 0.
    do i = 1,NSPACE_SUB
        call progress_gfortran(i, NSPACE_SUB)
        call eulspace_sub%get_ori(i, o)
        call pgrpsym%nearest_proj_neighbors(eulspace, o, angres, lnns)
        nnn(i) = real(count(lnns))
        tot_nnn = tot_nnn + nnn(i)
    end do
    call calc_stats(nnn, nnn_stats)
    write(*,'(a,1x,f15.6)') 'TOTAL   # NN :', tot_nnn
    write(*,'(a,1x,f15.6)') 'AVERAGE # NN :', nnn_stats%avg
    write(*,'(a,1x,f15.6)') 'SDEV    # NN :', nnn_stats%sdev
    write(*,'(a,1x,f15.6)') 'MINIMUM # NN :', nnn_stats%minv
    write(*,'(a,1x,f15.6)') 'MAXIMUM # NN :', nnn_stats%maxv
    call simple_end('**** SIMPLE_TEST_NEIGH_WORKFLOW NORMAL STOP ****')
end subroutine exec_test_neigh

end module simple_commanders_test_numerics
