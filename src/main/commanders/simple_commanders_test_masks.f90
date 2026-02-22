!@descr: for all masks tests
module simple_commanders_test_masks
use simple_commanders_api
implicit none
#include "simple_local_flags.inc"

type, extends(commander_base) :: commander_test_bounds_from_mask3D_test
  contains
    procedure :: execute      => exec_test_bounds_from_mask3D_test
end type commander_test_bounds_from_mask3D_test

type, extends(commander_base) :: commander_test_graphene_mask
  contains
    procedure :: execute      => exec_test_graphene_mask
end type commander_test_graphene_mask

type, extends(commander_base) :: commander_test_image_bin
  contains
    procedure :: execute      => exec_test_image_bin
end type commander_test_image_bin

type, extends(commander_base) :: commander_test_mask
  contains
    procedure :: execute      => exec_test_mask
end type commander_test_mask

type, extends(commander_base) :: commander_test_msk_routines
  contains
    procedure :: execute      => exec_test_msk_routines
end type commander_test_msk_routines

type, extends(commander_base) :: commander_test_nano_mask
  contains
    procedure :: execute      => exec_test_nano_mask
end type commander_test_nano_mask

type, extends(commander_base) :: commander_test_otsu_test
  contains
    procedure :: execute      => exec_test_otsu_test
end type commander_test_otsu_test

type, extends(commander_base) :: commander_test_ptcl_center
  contains
    procedure :: execute      => exec_test_ptcl_center
end type commander_test_ptcl_center

contains

subroutine exec_test_bounds_from_mask3D_test( self, cline )
    use simple_math,  only: bounds_from_mask3D
    use simple_image, only: image
    class(commander_test_bounds_from_mask3D_test), intent(inout) :: self
    class(cmdline),                                intent(inout) :: cline
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
    call simple_end('**** SIMPLE_TEST_BOUNDS_FROM_MASK3D_TEST_WORKFLOW NORMAL STOP ****')
end subroutine exec_test_bounds_from_mask3D_test

subroutine exec_test_graphene_mask( self, cline )
    class(commander_test_graphene_mask), intent(inout) :: self
    class(cmdline),                      intent(inout) :: cline
    integer, parameter   :: box  = 160
    real,    parameter   :: smpd = 0.358
    real,    allocatable :: res(:)
    logical, allocatable :: graphene_mask(:)
    integer              :: i
    res = get_resarr( BOX, SMPD )
    graphene_mask = calc_graphene_mask( BOX, SMPD )
    write(*,*) 'RES OF GRAPHENE BAND 1: ', GRAPHENE_BAND1
    write(*,*) 'RES OF GRAPHENE BAND 2: ', GRAPHENE_BAND2
    do i=1,size(res)
        write(*,*) '>>> RESOLUTION:', res(i), graphene_mask(i)
    end do
    call simple_end('**** SIMPLE_TEST_GRAPHENE_MASK_WORKFLOW NORMAL STOP ****')
end subroutine exec_test_graphene_mask

subroutine exec_test_image_bin( self, cline )
    !use simple_image_bin
    class(commander_test_image_bin),    intent(inout) :: self
    class(cmdline),                     intent(inout) :: cline
    type(image_bin)      :: test_image_bin, test_image_bin_3D, ccimage
    integer, allocatable :: imat_ccs(:,:,:)
    integer              :: imat(4,4,1), imat_3D(3,3,2), nccs, i
    real                 :: dist, dist_3D
    ! 2D test
    print *, '2D TEST'
    call test_image_bin%new_bimg(ldim=[4,4,1], smpd=1.)
    imat(1,:,1) = [1,0,0,0]
    imat(2,:,1) = [1,0,1,0]
    imat(3,:,1) = [0,0,0,0]
    imat(4,:,1) = [0,0,0,1]
    call test_image_bin%set_imat(imat)
    call test_image_bin%max_dist(dist)
    print *, 'MAXIMUM DISTANCE FROM CENTER IS ', dist
    call test_image_bin%find_ccs(ccimage)
    call ccimage%get_nccs(nccs)
    print *, nccs
    call ccimage%get_imat(imat_ccs)
    do i = 1,4
        print *, imat_ccs(i,:,1)
    enddo
    ! extreme unit test case of 1
    imat = 1
    call test_image_bin%set_imat(imat)
    call test_image_bin%max_dist(dist)
    call test_image_bin%find_ccs(ccimage)
    call ccimage%get_nccs(nccs)
    print *, nccs
    call ccimage%get_imat(imat_ccs)
    do i = 1,4
        print *, imat_ccs(i,:,1)
    enddo
    ! extreme unit test case of 0
    imat = 0
    call test_image_bin%set_imat(imat)
    call test_image_bin%max_dist(dist)
    call test_image_bin%find_ccs(ccimage)
    call ccimage%get_nccs(nccs)
    print *, nccs
    call ccimage%get_imat(imat_ccs)
    do i = 1,4
        print *, imat_ccs(i,:,1)
    enddo
    print *, ' '
    ! 3D test
    print *, '3D TEST'
    call test_image_bin_3D%new_bimg(ldim=[3,3,2],smpd=1.)
    imat_3D = reshape( (/ 0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1 /), shape(imat_3D))
    call test_image_bin_3D%set_imat(imat_3D)
    call test_image_bin_3D%max_dist(dist_3D)
    print *, 'MAXIMUM DISTANCE FROM CENTER IS ', dist_3D
    call simple_end('**** SIMPLE_TEST_IMAGE_BIN_WORKFLOW NORMAL STOP ****')
end subroutine exec_test_image_bin

subroutine exec_test_mask( self, cline )
    class(commander_test_mask),         intent(inout) :: self
    class(cmdline),                     intent(inout) :: cline
    real,    parameter :: SMPD=0.356, MSK=56.
    integer, parameter :: EDGE=6, BOX=160
    type(image)        :: spher_msk
    type(image_bin)    :: spher_msk_bin
    integer            :: npix
    ! make a spherical mask
    call spher_msk%disc([BOX,BOX,BOX], SMPD, MSK, npix)
    ! transfer to binimg instance
    call spher_msk_bin%transfer2bimg(spher_msk)
    ! apply soft cosine edge
    call spher_msk_bin%cos_edge(EDGE)
    ! write
    call spher_msk_bin%write_bimg(string('spherical_mask.mrc'))
    call simple_end('**** SIMPLE_TEST_MASK_WORKFLOW NORMAL STOP ****')
end subroutine exec_test_mask

subroutine exec_test_msk_routines( self, cline )
    use simple_image
    use simple_imgarr_utils, only: write_imgarr
    class(commander_test_msk_routines), intent(inout) :: self
    class(cmdline),                     intent(inout) :: cline
    type(image)              :: img
    type(image), allocatable :: stk(:)
    integer :: i, nimgs
    integer :: ldim2(3), ldim3(3)
    real    :: smpd
    smpd  = 3.0
    nimgs = 8
    ! -----------------------------
    ! POSITIVE TESTS (2D serial)
    ! -----------------------------
    ldim2 = [256,256,1]
    call img%new(ldim2, smpd)
    call img%ran()
    call img%memoize_mask_coords
    call img%mask2D_soft(ldim2(1)/3.0)
    call img%mask2D_softavg(ldim2(1)/3.0)
    call img%mask2D_hard(ldim2(1)/3.0)
    write(*,*) "OK: 2D serial"
    ! -----------------------------
    ! POSITIVE TESTS (2D parallel)
    ! -----------------------------
    call unmemoize_mask_coords()
    allocate(stk(nimgs))
    do i = 1, nimgs
        call stk(i)%new(ldim2, smpd)
        call stk(i)%ran()
    end do
    ! Ensure memoization is established outside parallel
    call stk(1)%memoize_mask_coords
    !$omp parallel do default(shared) private(i) proc_bind(close) schedule(static)
    do i = 1, nimgs
        call stk(i)%mask2D_soft(ldim2(1)/3.0)
    end do
    !$omp end parallel do
    write(*,*) "OK: 2D parallel soft"
    call write_imgarr(stk, string('mask2D_soft_openmp.mrcs'))

    do i = 1, nimgs
        call stk(i)%new(ldim2, smpd)
        call stk(i)%ran()
    end do
    call stk(1)%memoize_mask_coords
    !$omp parallel do default(shared) private(i) proc_bind(close) schedule(static)
    do i = 1, nimgs
        call stk(i)%mask2D_softavg(ldim2(1)/3.0)
    end do
    !$omp end parallel do
    write(*,*) "OK: 2D parallel softavg"
    call write_imgarr(stk, string('mask2D_softavg_openmp.mrc'))

    do i = 1, nimgs
        call stk(i)%new(ldim2, smpd)
        call stk(i)%ran()
    end do
    call stk(1)%memoize_mask_coords
    !$omp parallel do default(shared) private(i) proc_bind(close) schedule(static)
    do i = 1, nimgs
        call stk(i)%mask2D_hard(ldim2(1)/3.0)
    end do
    !$omp end parallel do
    write(*,*) "OK: 2D parallel hard"
    call write_imgarr(stk, string('mask2D_hard_openmp.mrc'))
    ! -----------------------------
    ! POSITIVE TESTS (3D serial)
    ! -----------------------------
    ldim3 = [64,64,64]
    call img%new(ldim3, smpd)
    call img%ran()
    call img%memoize_mask_coords
    call img%mask3D_soft(ldim3(1)/3.0)
    call img%mask3D_softavg(ldim3(1)/3.0)
    call img%mask3D_hard(ldim3(1)/3.0)
    write(*,*) "OK: 3D serial"

    ! -----------------------------
    ! POSITIVE TESTS (3D parallel)
    ! -----------------------------
    call unmemoize_mask_coords()
    do i = 1, nimgs
        call stk(i)%new(ldim3, smpd)
        call stk(i)%ran()
    end do
    ! Pre-memoize outside parallel so 3D calls don't attempt memoize inside parallel
    call stk(1)%memoize_mask_coords
    !$omp parallel do default(shared) private(i) proc_bind(close) schedule(static)
    do i = 1, nimgs
        call stk(i)%mask3D_soft(ldim3(1)/3.0)
    end do
    !$omp end parallel do
    write(*,*) "OK: 3D parallel soft"
    do i = 1, nimgs
        call stk(i)%new(ldim3, smpd)
        call stk(i)%ran()
    end do
    call stk(1)%memoize_mask_coords
    !$omp parallel do default(shared) private(i) proc_bind(close) schedule(static)
    do i = 1, nimgs
        call stk(i)%mask3D_softavg(ldim3(1)/3.0)
    end do
    !$omp end parallel do
    write(*,*) "OK: 3D parallel softavg"
    do i = 1, nimgs
        call stk(i)%new(ldim3, smpd)
        call stk(i)%ran()
    end do
    call stk(1)%memoize_mask_coords
    !$omp parallel do default(shared) private(i) proc_bind(close) schedule(static)
    do i = 1, nimgs
        call stk(i)%mask3D_hard(ldim3(1)/3.0)
    end do
    !$omp end parallel do
    write(*,*) "OK: 3D parallel hard"
    write(*,*) "ALL MASK TESTS PASSED"
    call simple_end('**** SIMPLE_TEST_MSK_ROUTINES_WORKFLOW NORMAL STOP ****')
end subroutine exec_test_msk_routines

subroutine exec_test_nano_mask( self, cline )
    use simple_image,     only: image
    use simple_image_msk, only: automask2D
    use simple_parameters, only: parameters
    class(commander_test_nano_mask),    intent(inout) :: self
    class(cmdline),                     intent(inout) :: cline
    ! constants
    character(len=*), parameter :: STK='selected.spi'
    real,             parameter :: SMPD=0.358
    integer,          parameter :: NGROW=3, WINSZ=1, EDGE=12
    ! variables
    type(parameters)            :: params
    type(image),    allocatable :: imgs(:)
    real,           allocatable :: diams(:), shifts(:,:)
    integer                     ::  n, i, ldim(3)
    ! setup parameters
    call cline%set('smpd',    SMPD)
    call cline%set('amsklp',  20.)
    call cline%set('msk',     100.)
    call cline%set('automsk', 'no')
    call cline%set('part',    1.)
    call params%new(cline)
    ! read images
    call find_ldim_nptcls(string(STK), ldim, n)
    allocate(imgs(n))
    do i = 1, n
        call imgs(i)%new(ldim, SMPD)
        call imgs(i)%read(string(STK), i)
    end do
    ! mask
    call automask2D(params, imgs, NGROW, WINSZ, EDGE, diams, shifts)
    call simple_end('**** SIMPLE_TEST_NANO_MASK_WORKFLOW NORMAL STOP ****')
end subroutine exec_test_nano_mask

subroutine exec_test_otsu_test( self, cline )
    class(commander_test_otsu_test),    intent(inout) :: self
    class(cmdline),                     intent(inout) :: cline
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
    d    = abs(avg1 - avg2)
    print *, 'd = ', d ! increases with variance
    call simple_end('**** SIMPLE_TEST_OTSU_TEST_WORKFLOW NORMAL STOP ****')
end subroutine exec_test_otsu_test

subroutine exec_test_ptcl_center( self, cline )
    use simple_parameters,       only: parameters
    use simple_image,            only: image
    use simple_projector,        only: projector
    use simple_cmdline,          only: cmdline
    use simple_commanders_atoms, only: commander_pdb2mrc
    class(commander_test_ptcl_center),  intent(inout) :: self
    class(cmdline),                     intent(inout) :: cline
    integer,          parameter   :: ORI_IND = 15, NPLANES = 100, MAX_R = 45, CENTER_RAD = 40
    character(len=:), allocatable :: cmd
    real,             allocatable :: pspec(:)
    type(commander_pdb2mrc) :: xpdb2mrc
    type(parameters)        :: p
    type(cmdline)           :: cline_pdb2mrc
    type(image)             :: vol, noise, ptcl, ptcl_pad, roavg
    type(projector)         :: vol_pad
    type(oris)              :: spiral
    type(ori)               :: o1
    integer                 :: rc, ifoo, iind
    real                    :: ave, sdev, maxv, minv, masscen(3), sh(2)
    logical                 :: mrc_exists
    if( command_argument_count() < 4 )then
        write(logfhandle,'(a)') 'ERROR! Usage: simple_test_ptcl_center smpd=xx nthr=yy vol1=volume.mrc mskdiam=zz'
        write(logfhandle,'(a)') 'Example: https://www.rcsb.org/structure/1jyx with smpd=1. mskdiam=180'
        write(logfhandle,'(a)') 'DEFAULT TEST (example above) is running now...'
        inquire(file="1JYX.mrc", exist=mrc_exists)
        if( .not. mrc_exists )then
            write(*, *) 'Downloading the example dataset...'
            cmd = 'curl -s -o 1JYX.pdb https://files.rcsb.org/download/1JYX.pdb'
            call execute_command_line(cmd, exitstat=rc)
            write(*, *) 'Converting .pdb to .mrc...'
            call cline_pdb2mrc%set('smpd',                1.)
            call cline_pdb2mrc%set('pdbfile',     '1JYX.pdb')
            call cline_pdb2mrc%checkvar('smpd',            1)
            call cline_pdb2mrc%checkvar('pdbfile',         2)
            call cline_pdb2mrc%check()
            call xpdb2mrc%execute(cline_pdb2mrc)
            call cline_pdb2mrc%kill()
            cmd = 'rm 1JYX.pdb'
            call execute_command_line(cmd, exitstat=rc)
        endif
        call cline%set('smpd',           1.)
        call cline%set('nthr',          16.)
        call cline%set('vol1',   '1JYX.mrc')
        call cline%set('mskdiam',      180.)
        call cline%set('lp',             3.)
    else
        call cline%parse_oldschool
    endif
    call cline%checkvar('smpd',    1)
    call cline%checkvar('nthr',    2)
    call cline%checkvar('vol1',    3)
    call cline%checkvar('mskdiam', 4)
    call cline%checkvar('lp',      5)
    call cline%check
    call p%new(cline)
    print *, 'box = ', p%box
    print *, 'msk = ', p%msk
    call find_ldim_nptcls(p%vols(1), p%ldim, ifoo)
    call vol%new(p%ldim, p%smpd)
    call vol%read(p%vols(1))
    call vol%stats('foreground', ave, sdev, maxv, minv)
    call spiral%new(NPLANES, is_ptcl=.false.)
    call spiral%spiral
    call ptcl%new(    [p%box,   p%box,   1],       p%smpd)
    call noise%new(   [p%box,   p%box ,  1],       p%smpd)
    call vol_pad%new( [p%boxpd, p%boxpd, p%boxpd], p%smpd)
    call ptcl_pad%new([p%boxpd, p%boxpd, 1],       p%smpd)
    call vol%pad(vol_pad)
    call vol_pad%fft
    call vol_pad%expand_cmat(p%box)
    call spiral%get_ori(ORI_IND, o1)
    call vol_pad%fproject(o1,ptcl_pad)
    call ptcl_pad%ifft
    ! add noise in a small center region of the vol
    call noise%gauran(0., 0.01 * sdev)
    call ptcl%zero_and_unflag_ft
    call ptcl_pad%clip(ptcl)
    call ptcl%add(noise)
    iind = 1
    call ptcl%write(string('noisy_ptcl.mrc'), iind)
    allocate(pspec(ptcl%get_nyq()),source=0.)
    call ptcl%masscen(masscen)
    print *, 'masscen (original) = ', masscen
    call ptcl%remove_neg
    call test_center(ptcl)
    call ptcl%fft
    call ptcl%power_spectrum(pspec)
    print *, 'power spectrum = ', pspec(1:5)
    call ptcl%shift2Dserial(masscen(1:2))
    call ptcl%ifft
    iind = iind + 1
    call ptcl%write(string('noisy_ptcl.mrc'), iind)
    call ptcl%masscen(masscen)
    print *, 'masscen (after ori shifted) = ', masscen
    call ptcl%roavg(10, roavg)
    iind = iind + 1
    call roavg%write(string('noisy_ptcl.mrc'), iind)
    call roavg%fft
    call roavg%power_spectrum(pspec)
    print *, '(roavg) power spectrum = ', sum(pspec)
    sh = [10., -15.]
    call ptcl%fft
    call ptcl%shift2Dserial(sh)
    call ptcl%ifft
    iind = iind + 1
    call ptcl%write(string('noisy_ptcl.mrc'), iind)
    call ptcl%masscen(masscen)
    call ptcl%fft
    call ptcl%power_spectrum(pspec)
    print *, 'masscen (after fixed shifted) = ', masscen
    print *, 'power spectrum = ', pspec(1:5)
    call ptcl%ifft
    call ptcl%remove_neg
    call test_center(ptcl)
    call ptcl%roavg(10, roavg)
    iind = iind + 1
    call roavg%write(string('noisy_ptcl.mrc'), iind)
    call roavg%fft
    call roavg%power_spectrum(pspec)
    print *, '(roavg) power spectrum = ', sum(pspec)
    call simple_end('**** SIMPLE_TEST_PTCL_CENTER_WORKFLOW NORMAL STOP ****')

    contains

        subroutine test_center( img )
            type(image), intent(in) :: img
            type(image)   :: win_img, win_avg
            real, pointer :: rmat_ptr(:,:,:)
            integer       :: origin, center(2), i, ic, jc, cnt, cens(2,(CENTER_RAD+1)**2), img_ind
            real(dp)      :: up(MAX_R, (CENTER_RAD+1)**2), max_E, lip((CENTER_RAD+1)**2), zp(MAX_R, (CENTER_RAD+1)**2)
            logical       :: outside
            origin  = p%box/2
            cnt     = 0
            img_ind = 0
            call win_img%new([MAX_R*2, MAX_R*2, 1], 1.)
            print *, origin-CENTER_RAD/2, origin+CENTER_RAD/2
            do ic = origin-CENTER_RAD/2, origin+CENTER_RAD/2
                do jc = origin-CENTER_RAD/2, origin+CENTER_RAD/2
                    cnt         = cnt + 1
                    center      = [ic, jc]
                    cens(:,cnt) = center
                    call img%window_center(center, MAX_R, win_img, outside)
                    if( outside ) print *, 'Window outside the image boundary'
                    call win_img%roavg(10, win_avg)
                    call win_avg%get_rmat_ptr(rmat_ptr)
                    do i = 1, MAX_R
                        zp(i,cnt) = sum(rmat_ptr(MAX_R+1:MAX_R+i, MAX_R+1, 1))
                    enddo
                    do i = 1, MAX_R
                        up(i,cnt) = sum(zp(1:i,cnt))
                    enddo
                    call win_avg%kill
                    if( cnt == 1 .or. cnt == 21 )then
                        img_ind = img_ind + 1
                        call win_img%write(string('window_imgs.mrc'), img_ind)
                    endif
                enddo
            enddo
            max_E = maxval(up(MAX_R,:))
            do i = 1, cnt
                lip(i) = sum(dabs(up(:,i)/max_E - 1._dp))
            enddo
            print *, 'center = ', cens(:,minloc(lip))
        end subroutine test_center

end subroutine exec_test_ptcl_center

end module simple_commanders_test_masks
