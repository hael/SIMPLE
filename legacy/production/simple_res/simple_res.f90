!==Program simple_res
!
! <res/begin> is a program for checking the low-pass resolution limit for a given Fourier index. <res/end>
!
! The code is distributed with the hope that it will be useful, but _WITHOUT_ _ANY_ _WARRANTY_.
! Redistribution and modification is regulated by the GNU General Public License.
! *Author:* Hans Elmlund 2015
!
program simple_res
use simple_params,  only: params
use simple_cmdline, only: cmdline
use simple_timing
implicit none
type(params)  :: p
type(cmdline) :: cline
real          :: lp
if( command_argument_count() < 2 )then
    write(*,'(a)', advance='no') 'SIMPLE_RES smpd=<sampling distance(in A)>'
    write(*,'(a)') ' find=<Fourier index> box=<box size (in pixels)>'
    stop
endif
call cline%parse
call cline%checkvar('smpd', 1)
call cline%checkvar('find', 2)
call cline%checkvar('box',  3)
call cline%check
call cline%set('prg', 'res')
p = params(cline)
lp = (real(p%box-1)*p%smpd)/real(p%find)
write(*,'(A,1X,f7.2)') '>>> LOW-PASS LIMIT:', lp
end program
