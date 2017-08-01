!==Program simple_iminfo
!
! <iminfo/begin> is a program for printing header information in MRC and SPIDER stacks and volumes <iminfo/end>
!
! The code is distributed with the hope that it will be useful, but _WITHOUT_ _ANY_ _WARRANTY_.
! Redistribution and modification is regulated by the GNU General Public License.
! *Author:* Hans Elmlund 2015
!
program simple_iminfo
use simple_defs     ! singleton
use simple_jiffys   ! singleton
use simple_cmdline, only: cmdline
use simple_params,  only: params
use simple_build,   only: build
use simple_jiffys,  only: simple_end, fname2format
use simple_image,   only: image
use simple_imgfile
implicit none
type(params)      :: p
type(image)       :: img
type(cmdline)     :: cline
character(len=20) :: conv
character(len=1)  :: form
integer           :: ldim(3), maxim, i, iform, n_nans, mode
real              :: smpd, sdev, ave, minv, maxv
if( command_argument_count() < 1 )then
    write(*,'(a)', advance='no') 'SIMPLE_IMINFO [fname=<filename.ext>] [box=<box size(pixels)>]'
    write(*,'(a)') ' [smpd=<sampling distance(in A)>] [stats=<yes|no|print{no}>]'
    stop
endif
call cline%parse
call cline%set('prg', 'iminfo')
p = params(cline) ! constants & derived constants produced
if( cline%defined('fname') )then
    call find_ldim_nptcls(p%fname, ldim, maxim, doprint=.true.)
endif
p%box  = ldim(1)
p%smpd = smpd
call img%new([ldim(1),ldim(2),1],p%smpd)
if( p%vis .eq. 'yes' .or. p%stats .ne. 'no' )then
    do i=1,maxim
        call img%read(p%fname, i)
        if( p%stats .ne. 'no' )then
            call img%cure(maxv, minv, ave, sdev, n_nans)
            if( p%stats .eq. 'print' .or. n_nans > 0 )then
                write(*,*) '*********IMAGE********', i, '*******'
                write(*,*) 'maxv = ',   maxv
                write(*,*) 'minv = ',   minv
                write(*,*) 'ave = ',    ave
                write(*,*) 'sdev = ',   sdev
                write(*,*) 'n_nans = ', n_nans
            endif
        endif
        if( p%vis .eq. 'yes' ) call img%vis
    end do
endif
call simple_end('**** SIMPLE_IMINFO NORMAL STOP ****')
end program simple_iminfo