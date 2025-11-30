program simple_test_ptcl_center
include 'simple_lib.f08'
use simple_cmdline,    only: cmdline
use simple_parameters, only: parameters
use simple_image,      only: image
use simple_projector,  only: projector
use simple_oris
use simple_ori
implicit none
integer,          parameter   :: ORI_IND = 15, NPLANES = 100, MAX_R = 45, CENTER_RAD = 40
character(len=:), allocatable :: cmd
real,             allocatable :: pspec(:)
type(parameters) :: p
type(cmdline)    :: cline
type(image)      :: vol, noise, ptcl, ptcl_pad, roavg
type(projector)  :: vol_pad
type(oris)       :: spiral
type(ori)        :: o1
integer          :: rc, ifoo, iind
real             :: ave, sdev, maxv, minv, masscen(3), sh(2)
logical          :: mrc_exists
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
        cmd = 'e2pdb2mrc.py 1JYX.pdb 1JYX.mrc'
        call execute_command_line(cmd, exitstat=rc)
        cmd = 'rm 1JYX.pdb'
        call execute_command_line(cmd, exitstat=rc)
    endif
    call cline%set('smpd'   , 1.)
    call cline%set('nthr'   , 16.)
    call cline%set('vol1'   , '1JYX.mrc')
    call cline%set('mskdiam', 180.)
    call cline%set('lp'   ,   3.)
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
call vol_pad%expand_cmat(p%alpha)
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

end program simple_test_ptcl_center