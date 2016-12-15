!==Program simple_check_box
!
! <check_box/begin> is a program for checking the image dimensions of MRC and SPIDER stacks and volumes <check_box/end>
!
! The code is distributed with the hope that it will be useful, but _WITHOUT_ _ANY_ _WARRANTY_.
! Redistribution and modification is regulated by the GNU General Public License.
! *Author:* Hans Elmlund 2015
!
program simple_check_box
use simple_cmdline, only: cmdline
use simple_jiffys,  only: simple_end
use simple_params,  only: params
use simple_timing
implicit none
type(params)  :: p
type(cmdline) :: cline
if( command_argument_count() < 1 )then
    write(*,'(a)') 'SIMPLE_CHECK_BOX [stk=<ptcls.ext>] [vol1=<vol.ext>]'
    stop
endif
call cline%parse
p = params(cline) ! parameters generated
write(*,'(A,1X,I7)') '>>> BOX:', p%box
! END GRACEFULLY
call simple_end('**** SIMPLE_CHECK_BOX NORMAL STOP ****')
end program
