!==Program simple_find
!
! <find/begin> is a program for calculating the Fourier index wen given low-pass limit in \AA{} (\texttt{lp}),
! box size in pixels (\texttt{box}) and sampling distance in \AA{} (\texttt{box}). <find/end>
!
! The code is distributed with the hope that it will be useful, but _WITHOUT_ _ANY_ _WARRANTY_.
! Redistribution and modification is regulated by the GNU General Public License.
! *Author:* Hans Elmlund 2015
!
program simple_find
use simple_params, only: params
use simple_cmdline
implicit none
type(params) :: p
integer      :: find
if( command_argument_count() < 3 )then
    write(*,'(a)', advance='no') 'SIMPLE_FIND lp=<low-pass limit(in A)> box=<box size (in pixels)>'
    write(*,'(a)') ' smpd=<sampling distance(in A)>'
    stop
endif
call parse_cmdline
call cmdcheckvar('lp',   1)
call cmdcheckvar('box',  2)
call cmdcheckvar('smpd', 3)
call cmdcheck
p = params()
find = int((real(p%box-1)*p%smpd)/p%lp)
write(*,'(A,1X,I4)') '>>> FOURIER INDEX:', find
end program