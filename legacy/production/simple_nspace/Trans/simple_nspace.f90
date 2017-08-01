!==Program simple_nspace
!
! <nspace/begin> is a program for checking the theoretical resolution limit for different numbers 
! of discrete projection directions <nspace/end>
!
! The code is distributed with the hope that it will be useful, but _WITHOUT_ _ANY_ _WARRANTY_.
! Redistribution and modification is regulated by the GNU General Public License.
! *Author:* Hans Elmlund 2015
!
program simple_nspace
use simple_cmdline, only: cmdline
use simple_jiffys,  only: simple_end
use simple_math,    only: resang
use simple_oris,    only: oris
use simple_params,  only: params
use simple_timing
implicit none
type(oris)    :: o
type(params)  :: p
type(cmdline) :: cline
real          :: ares
integer       :: i
if( command_argument_count() < 1 )then
    write(*,'(a)') 'SIMPLE_NSPACE moldiam=<molecular diameter (in A)>'
    stop
endif
call cline%parse
call cline%checkvar('moldiam', 1)
call cline%check
call cline%set('prg', 'nspace')
p = params(cline) ! parameters generated
do i=500,5000,500
    o = oris(i)
    call o%spiral
    ares = o%find_angres()
    write(*,'(A,1X,I7,1X,A,1X,F5.2)') 'NR OF PROJDIRS:', i, 'RESOLUTION:', resang(ares, p%moldiam)
end do
call simple_end('**** SIMPLE_NSPACE NORMAL STOP ****')
end program
