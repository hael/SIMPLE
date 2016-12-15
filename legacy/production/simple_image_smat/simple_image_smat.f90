!==Program simple_image_smat
!
! <image_smat/begin> is a program for creating a similarity matrix based on common line correlation. The idea
! being that it should be possible to cluster images based on their 3D similarity witout having a 3D model
! by only operating on class averages and find averages that fit well together in 3D. <image_smat/end>
!
! The code is distributed with the hope that it will be useful, but _WITHOUT_ _ANY_ _WARRANTY_.
! Redistribution and modification is regulated by the GNU General Public License.
! *Author:* Hans Elmlund 2016
!
program simple_image_smat
use simple_defs     ! singleton
use simple_corrmat  ! singleton
use simple_jiffys,  ! singleton
use simple_params,  only: params
use simple_build,   only: build
use simple_ori,     only: ori
use simple_imgfile, only: imgfile
use simple_image,   only: image
use simple_cmdline, only: cmdline
implicit none
type(params)         :: p
type(build)          :: b
type(cmdline)        :: cline
integer              :: iptcl, alloc_stat, funit, io_stat
real, allocatable    :: corrmat(:,:)
logical              :: debug=.false.
if( command_argument_count() < 2 )then
    write(*,'(a)', advance='no') 'SIMPLE_IMAGE_SMAT stk=<stack.ext> smpd=<sampling distance(in A)>'
    write(*,'(a)', advance='no') ' [lp=<low-pass limit(in A)>] [msk=<mask radius(in pixels)>]'
    write(*,'(a)') ' [hp=<high-pass limit(in A)>] [nthr=<nr of OpenMP threads{1}>]'
    stop
endif
call cline%parse
call cline%checkvar('stk',  1)
call cline%checkvar('smpd', 2)
call cline%check                              ! checks args and prints cmdline to cmdline.dat
call cline%set('prg', 'image_smat')
p = params(cline, .false.)                    ! constants & derived constants produced
call b%build_general_tbox(p, cline, .false., .true.) ! general objects built (no oritab reading)
allocate(b%imgs_sym(p%nptcls), stat=alloc_stat)
call alloc_err('In: simple_image_smat, 1', alloc_stat)
do iptcl=1,p%nptcls
    call b%imgs_sym(iptcl)%new([p%box,p%box,1], p%smpd, p%imgkind)
    call b%imgs_sym(iptcl)%read(p%stk, iptcl)
end do
write(*,'(a)') '>>> CALCULATING CORRELATIONS'
if( cline%defined('lp') )then
    if( .not. cline%defined('msk') ) stop 'need mask radius (msk) 4 Fourier corr calc!'
    call calc_cartesian_corrmat(b%imgs_sym, corrmat, p%msk, p%lp)
else
    if( cline%defined('msk') )then
        call calc_cartesian_corrmat(b%imgs_sym, corrmat, p%msk)
    else
        call calc_cartesian_corrmat(b%imgs_sym, corrmat)
    endif
endif
funit = get_fileunit()
open(unit=funit, status='REPLACE', action='WRITE', file='img_smat.bin', access='STREAM')
write(unit=funit,pos=1,iostat=io_stat) corrmat
if( io_stat .ne. 0 )then
    write(*,'(a,i0,a)') 'I/O error ', io_stat, ' when writing to image_smat.bin'
    stop 'I/O error; simple_image_smat'
endif
close(funit)
call simple_end('**** SIMPLE_IMAGE_SMAT NORMAL STOP ****')
end program simple_image_smat
