!==Program cavgassemble
!
! <cavgassemble/begin> is a program that assembles class averages when the clustering program (\prgname{simple\_prime2D}) 
! has been executed in distributed mode. <cavgassemble/end>
!
! The code is distributed with the hope that it will be useful, but _WITHOUT_ _ANY_ _WARRANTY_.
! Redistribution and modification is regulated by the GNU General Public License.
! *Author:* Hans Elmlund 2016
!
program simple_cavgassemble
use simple_defs                ! singleton
use simple_cmdline,            only: cmdline
use simple_build,              only: build
use simple_params,             only: params
use simple_hadamard2D_matcher, only: prime2D_assemble_sums_from_parts, prime2D_write_sums
use simple_jiffys,             only: simple_end
type(params)  :: p
type(build)   :: b
type(cmdline) :: cline
if( command_argument_count() < 7 )then
    write(*,'(a)', advance='no') 'SIMPLE_CAVGASSEMBLE stk=<stack.ext> smpd=<sampling distance(in A)>'
    write(*,'(a)', advance='no') ' msk=<mask radius(in pixels)> ncls=<nr clusters> oritab=<previous clustering doc>'
    write(*,'(a)', advance='no') ' which_iter=<iteration nr> npart=<nr partitions> [ctf=<yes|no|flip|mul{no}>]' 
    write(*,'(a)', advance='no') ' [nthr=<nr openMP threads{1}>] [inner=<inner mask radius(in pixels)>]'
    write(*,'(a)') ' [width=<pixels falloff inner mask{10}>]'
    stop
endif
call cline%parse
call cline%checkvar('stk',        1)
call cline%checkvar('smpd',       2)
call cline%checkvar('msk',        3)
call cline%checkvar('ncls',       4)
call cline%checkvar('oritab',     5)
call cline%checkvar('which_iter', 6)
call cline%checkvar('npart',      7)
call cline%check
call cline%set('prg', 'cavgassemble')
p = params(cline) ! constants & derived constants produced
call b%build_general_tbox(p,cline) ! general objects built
call b%build_hadamard_prime2D_tbox(p)
call prime2D_assemble_sums_from_parts(b, p)
call prime2D_write_sums( b, p, p%which_iter)
call simple_end('**** SIMPLE_CAVGASSEMBLE NORMAL STOP ****')
end program simple_cavgassemble
