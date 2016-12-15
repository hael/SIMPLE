!==Program simple_check_nptcls
!
! <check_nptcls/begin> is a program for checking the number of images in MRC and SPIDER stacks <check_nptcls/end>
!
! The code is distributed with the hope that it will be useful, but _WITHOUT_ _ANY_ _WARRANTY_.
! Redistribution and modification is regulated by the GNU General Public License.
! *Author:* Hans Elmlund 2015
!
program simple_check_nptcls
use simple_cmdline, only: cmdline
use simple_jiffys,  only: simple_end
use simple_params,  only: params
use simple_timing
implicit none
type(params)  :: p
type(cmdline) :: cline
if( command_argument_count() < 1 )then
    write(*,'(a)') 'SIMPLE_CHECK_NPTCLS stk=<ptcls.ext>'
    stop
endif
call cline%parse
call cline%checkvar('stk', 1)
call cline%check
call cline%set('prg', 'check_nptcls')
p = params(cline) ! parameters generated
write(*,'(A,1X,I7)') '>>> NPTCLS:', p%nptcls
! END GRACEFULLY
call simple_end('**** SIMPLE_CHECK_NPTCLS NORMAL STOP ****')
end program
