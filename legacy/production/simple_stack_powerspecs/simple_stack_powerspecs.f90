!==Program simple_stack_powerspecs
!
!<stack_powerspecs/begin> is a program for stacking powerspectra that unblur generated <stack_powerspecs/end>
!
! The code is distributed with the hope that it will be useful, but _WITHOUT_ _ANY_ _WARRANTY_.
! Redistribution and modification is regulated by the GNU General Public License.
! *Author:* Hans Elmlund 2011-08-16.
!
program simple_stack_powerspecs
use simple_build,     only: build
use simple_params,    only: params
use simple_imgfile,   only: imgfile
use simple_image,     only: image
use simple_cmdline,   only: cmdline
use simple_defs       ! singleton
use simple_jiffys     ! singleton
implicit none
type(params)                       :: p
type(build)                        :: b
type(cmdline)                      :: cline
integer                            :: nspecs, ldim(3), ispec, ifoo
character(len=STDLEN), allocatable :: specnames(:)
type(image)                        :: mask, tmp
real                               :: mm(2)
logical, parameter                 :: debug = .false.
if( command_argument_count() < 3 )then
    write(*,'(a)',advance='no') 'SIMPLE_STACK_POWERSPECS filetab=<specs.txt> smpd=<sampling distance(in A)>'
    write(*,'(a)',advance='no') ' outstk=<powerspecs.ext> [lp=low-pass limit(in A){6}]'
    write(*,'(a)') ' [clip=<clip2box{256}>]'
    stop
endif
call cline%parse
call cline%checkvar('filetab', 1)
call cline%checkvar('smpd',    2)
call cline%checkvar('outstk',  3)
call cline%check
if( .not. cline%defined('lp') )   call cline%set('lp',    6. )
if( .not. cline%defined('clip') ) call cline%set('clip', 256.)
p = params(cline,checkpara=.false.)             ! constants & derived constants produced
call b%build_general_tbox(p,cline,do3d=.false.) ! general stuff built
call read_filetable(p%filetab, specnames)
nspecs = size(specnames)
if( debug ) write(*,*) 'read the spec filenames'
call find_ldim_nptcls(specnames(1),ldim,ifoo)
ldim(3) = 1 ! to correct for the stupide 3:d dim of mrc stacks
if( debug ) write(*,*) 'logical dimension: ', ldim
p%box    = ldim(1)
! create mask
call tmp%new([p%clip,p%clip,1], p%smpd)
tmp = cmplx(1.,0.)
call tmp%bp(0.,p%lp,0.)
call tmp%ft2img('real', mask)
call mask%write('resolution_mask.mrc', 1)
! prepare b%img and tmp for reading
call b%img%new([p%box,p%box,1], p%smpd)
tmp = 0.0
! LOOP OVER SPECTRA
do ispec=1,nspecs
    if( .not. file_exists(specnames(ispec)) )then
        write(*,*) 'inputted spec file does not exist: ', trim(adjustl(specnames(ispec)))
    endif
    call b%img%read(specnames(ispec))
    call b%img%clip(tmp)  
    mm = tmp%minmax()
    if( debug ) print *, 'min/max: ', mm(1), mm(2)
    call tmp%write(p%outstk, ispec)
    call progress(ispec, nspecs)
end do
call simple_end('**** SIMPLE_STACK_POWERSPECS NORMAL STOP ****')
end program simple_stack_powerspecs