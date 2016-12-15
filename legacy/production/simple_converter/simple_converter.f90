!==Program simple_converter
!
! <converter/begin> is a program for converting between SPIDER and MRC formats. 
! <converter/end> 
!
! The code is distributed with the hope that it will be useful, but _WITHOUT_ _ANY_ _WARRANTY_.
! Redistribution and modification is regulated by the GNU General Public License.
! *Author:* Hans Elmlund X-mas 2015
!
program simple_converter
use simple_cmdline, only: cmdline
use simple_jiffys,  only: simple_end, progress
use simple_params,  only: params
use simple_build,   only: build
implicit none
type(params)  :: p
type(build)   :: b
type(cmdline) :: cline
integer       :: iptcl
if( command_argument_count() < 2 )then
    write(*,'(a)',advance='no') 'SIMPLE_CONVERTER [stk=<input particle stack>] [vol1=<invol.ext>]'
    write(*,'(a)') ' [outstk=<output particle stack>] [outvol=<outvol.ext>]'
    stop
endif
call cline%parse
p = params(cline, allow_mix=.true.) ! parameters generated
call b%build_general_tbox(p, cline) ! general objects built
if( cline%defined('stk') )then
    do iptcl=1,p%nptcls
        call progress(iptcl, p%nptcls)
        call b%img%read(p%stk, iptcl)
        call b%img%write(p%outstk, iptcl)
    end do 
else if( cline%defined('vol1') )then
    call b%vol%read(p%vols(1))
    call b%img%write(p%outvol)
else
    stop 'not enough arguments to execute simple_converter'
endif
call simple_end('**** SIMPLE_CONVERTER NORMAL STOP ****')
end program simple_converter