!==Program simple_check2d_conv
!
! <check2d_conv/begin> is a program for checking if a PRIME2D run has converged. The statistics outputted
! include (1) the overlap between the distribution of parameters for succesive runs. (2) The percentage of search space 
! scanned, i.e. how many reference images are evaluated on average. (3) The average correlation between the images 
! and their corresponding best matching reference section. If convergence to a local optimum is achieved, the fraction
! increases. Convergence is achieved if the parameter distribution overlap is larger than $0.95$ and more than 99\%
! of the reference sections need to be searched to find an improving solution. <check2d_conv/end>
!
! The code is distributed with the hope that it will be useful, but _WITHOUT_ _ANY_ _WARRANTY_.
! Redistribution and modification is regulated by the GNU General Public License.
! *Author:* Hans Elmlund X-mas 2015
!
program simple_check2D_conv
use simple_cmdline, only: cmdline
use simple_build,   only: build
use simple_params,  only: params
use simple_jiffys,  only: simple_end
implicit none
type(params)  :: p
type(build)   :: b
type(cmdline) :: cline
logical       :: converged
if( command_argument_count() < 4 )then
    write(*,'(a)',advance='no') 'SIMPLE_CHECK2D_CONV smpd=<sampling distance(in A)> box=<image size(in pixels)>'
    write(*,'(a)',advance='no') ' oritab=<clustering doc> nptcls=<nr particle images>'
    write(*,'(a)') ' [lp=<low-pass limit{20}>]'
    stop
endif
call cline%parse
if( .not. cline%defined('lp') ) call cline%set('lp', 20.)
call cline%checkvar('smpd',   1)
call cline%checkvar('oritab', 2)
call cline%checkvar('nptcls', 3)
call cline%check
call cline%set('prg', 'check2D_conv')
p = params(cline)                   ! parameters generated
call b%build_general_tbox(p, cline) ! general objects built
p%ncls = b%a%get_ncls()
converged = b%conv%check_conv2D()   ! convergence check
call simple_end('**** SIMPLE_CHECK2D_CONV STOP ****')    
end program 