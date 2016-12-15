!==Program simple_delete
!
! <delete/begin> is a program for deleting images based on correlation matching. <delete/end> 
!
! The code is distributed with the hope that it will be useful, but _WITHOUT_ _ANY_ _WARRANTY_.
! Redistribution and modification is regulated by the GNU General Public License.
! *Author:* Hans Elmlund X-mas 2015
!
program simple_delete
use simple_cmdline  ! singleton
use simple_corrmat  ! singleton
use simple_jiffys,  only: simple_end, alloc_err
use simple_params,  only: params
use simple_build,   only: build
use simple_imgfile, only: imgfile
use simple_image,   only: image
implicit none
type(params)             :: p
type(build)              :: b
type(imgfile)            :: imghandle
type(image), allocatable :: imgs_del(:), imgs_all(:)
integer                  :: iimg, idel, nall, ndel, loc(1), cnt
integer, allocatable     :: selected4del(:)
real, allocatable        :: correlations(:,:)
logical                  :: debug=.false.
if( command_argument_count() < 2 )then
    write(*,'(a)', advance='no') 'SIMPLE_DELETE stk=<all_imgs.ext> stk2=<imgs2del.ext> outstk=<keep.ext>'
    write(*,'(a)') ' [nthr=<nr of OpenMP threads{1}>]'
    stop
endif
call parse_cmdline
call cmdcheckvar('stk',    1)
call cmdcheckvar('stk2',   2)
call cmdcheckvar('outstk', 3)
call cmdcheck
p = params()                 ! parameters generated
call b%build_general_tbox(p) ! general objects built
! find number of images to delete
call imghandle%open(p%stk2)
ndel = imghandle%getStackSz()
call imghandle%close
! find number of original images
call imghandle%open(p%stk)
nall = imghandle%getStackSz()
call imghandle%close
! read images
allocate(imgs_del(ndel), imgs_all(nall))
do idel=1,ndel
    call imgs_del(idel)%new([p%box,p%box,1], p%smpd)
    call imgs_del(idel)%read(p%stk2, idel)
end do
do iimg=1,nall
    call imgs_all(iimg)%new([p%box,p%box,1], p%smpd)
    call imgs_all(iimg)%read(p%stk, iimg)
end do
write(*,'(a)') '>>> CALCULATING CORRELATIONS'
call calc_cartesian_corrmat(imgs_del, imgs_all, correlations)
! find selected4del selected4del
allocate(selected4del(ndel))
do idel=1,ndel
    loc = maxloc(correlations(idel,:))
    selected4del(idel) = loc(1)
    if( debug ) print *, 'selected4del: ', loc(1), ' with corr: ', correlations(idel,loc(1))
end do
! write outoput
cnt = 0
do iimg=1,nall
    if( any(selected4del == iimg) )then
        cycle
    else
        cnt = cnt+1
        call b%img%read(p%stk, iimg)
        call b%img%write(p%outstk, cnt)
    endif
end do
call simple_end('**** SIMPLE_DELETE NORMAL STOP ****')
end program simple_delete