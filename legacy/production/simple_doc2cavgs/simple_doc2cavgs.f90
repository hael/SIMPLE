!==Program simple_doc2cavgs
!
! <doc2cavgs/begin> is a program for generating class averages. We use it to
! to re-generate class averages when \prgname{simple\_prime2d} has been run on downscaled images. If the images processed
! with PRIME2D were downscaled from 200x200 to 100x00, set \texttt{mul=2}.
! <comment/begin>The program assumes that the images have been phase-flipped, as no CTF correction by Wiener restoration
! is implemented yet. <comment/end> <doc2cavgs/end> 
!
! The code is distributed with the hope that it will be useful, but _WITHOUT_ _ANY_ _WARRANTY_.
! Redistribution and modification is regulated by the GNU General Public License.
! *Author:* Hans Elmlund X-mas 2015
!
program simple_doc2cavgs
use simple_cmdline,            only: cmdline
use simple_jiffys,             only: simple_end
use simple_hadamard2D_matcher, only: prime2D_assemble_sums, prime2D_write_sums
use simple_params,             only: params
use simple_build,              only: build
implicit none
type(params)  :: p
type(build)   :: b
type(cmdline) :: cline
if( command_argument_count() < 3 )then
    write(*,'(a)',advance='no') 'SIMPLE_DOC2CAVGS stk=<stack.ext> smpd=<sampling distance(in A)>'
    write(*,'(a)',advance='no') ' msk=<mask radius(in pixels)> ncls=<nr of clusters> oritab=<previous'
    write(*,'(a)',advance='no') ' clustering doc> which_iter=<iteration nr> [mul=<shift multiplication'
    write(*,'(a)',advance='no') ' factor{1}>] [nthr=<nr of OpenMP threads{1}>] [ctf=<yes|no|flip|mul{no}>]'
    write(*,'(a)') ' [deftab=<defocus info file>] [remap=<yes|no{no}>]'
    stop
endif
call cline%parse
call cline%checkvar('stk',        1)
call cline%checkvar('smpd',       2)
call cline%checkvar('msk',        3)
call cline%checkvar('ncls',       4)
call cline%checkvar('oritab',     5)
call cline%checkvar('which_iter', 6)
call cline%check
call cline%set('prg', 'doc2cavgs')
p = params(cline)                   ! parameters generated
call b%build_general_tbox(p, cline) ! general objects built
call b%build_hadamard_prime2D_tbox(p)
if( p%srch_inpl .eq. 'yes' )then
    write(*,'(a)') '>>> GENERATING CLUSTERS FOR SEARCH WITH IN-PLANE ALIGNMENT'
else
    write(*,'(a)') '>>> GENERATING CLUSTERS FOR SEARCH WITHOUT IN-PLANE ALIGNMENT'
endif
if( p%remap .eq. 'yes' )then
    call b%a%remap_classes
    p%ncls = b%a%get_ncls()
endif
if( cline%defined('which_iter') )then
    call prime2D_assemble_sums(b, p, mul=p%mul)
    call prime2D_write_sums(b, p, p%which_iter)
else
    call prime2D_assemble_sums(b, p, mul=p%mul)
    call prime2D_write_sums(b, p)
endif
call simple_end('**** SIMPLE_DOC2CAVGS NORMAL STOP ****')
end program simple_doc2cavgs
