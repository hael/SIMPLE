!==Program simple_npeaks
!
! <npeaks/begin> is a program for checking the number of nonzero orientation weights (number of correlation peaks included in the weighted reconstruction) <npeaks/end>
!
! The code is distributed with the hope that it will be useful, but _WITHOUT_ _ANY_ _WARRANTY_.
! Redistribution and modification is regulated by the GNU General Public License.
! *Author:* Hans Elmlund 2015
!
program simple_npeaks
use simple_build,   only: build
use simple_oris,    only: oris
use simple_params,  only: params
use simple_jiffys,  only: simple_end
use simple_cmdline, only: cmdline
implicit none
type(build)   :: b
type(params)  :: p
type(cmdline) :: cline
if( command_argument_count() < 2 )then
    write(*,'(a)',advance='no') 'SIMPLE_NPEAKS smpd=<sampling distance(in A)> box=<image size(in pixels)>'
    write(*,'(a)',advance='no') ' lp=<low-pass limit(A){20}> [nspace=<nr of projection directions{1000}>] '
    write(*,'(a)') ' [moldiam=<molecular diameter(A)>] [pgrp=<cn|dn|t|o|i{c1}>] '
    stop
endif
call cline%parse
call cline%checkvar('smpd', 1)
call cline%checkvar('box',  2)
call cline%checkvar('lp',   3)
if( .not. cline%defined('nspace') ) call cline%set('nspace', 1000.)
if( .not. cline%defined('lp') )     call cline%set('lp',       20.)
call cline%check
call cline%set('prg', 'npeaks')
p = params(cline) ! parameters generated
call b%build_general_tbox(p, cline)
p%npeaks = min(10,b%e%find_npeaks(p%lp, p%moldiam))
write(*,'(A,1X,I4)') '>>> NPEAKS:', p%npeaks
! END GRACEFULLY
call simple_end('**** SIMPLE_NPEAKS NORMAL STOP ****')
end program