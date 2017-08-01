!==Program simple_prime2D_init
!
! <prime2D_init/begin> is a program for initialising prime2D. We use it to produce the initial random
! references when executing \prgname{simple\_prime2D}. 
! <comment/begin>The program assumes that the images have been phase-flipped, as no CTF correction by 
! Wiener restoration is implemented yet. The random clustering and in-plane alignment will be printed in 
! the file \texttt{prime2D\_startdoc.txt} produced by the program. This file is used together with the 
! initial references (\texttt{startcavgsmsk.ext}) to execute \prgname{simple\_prime2D}.
! <comment/end><prime2D_init/end> 
!
! The code is distributed with the hope that it will be useful, but _WITHOUT_ _ANY_ _WARRANTY_.
! Redistribution and modification is regulated by the GNU General Public License.
! *Author:* Hans Elmlund X-mas 2015
!
program simple_prime2D_init
use simple_cmdline,            only: cmdline
use simple_jiffys,             only: simple_end
use simple_hadamard2D_matcher, only: prime2D_assemble_sums, prime2D_norm_sums, prime2D_write_sums
use simple_params,             only: params
use simple_build,              only: build
implicit none
type(params)  :: p
type(build)   :: b
type(cmdline) :: cline
integer       :: ncls_in_oritab, icls
if( command_argument_count() < 3 )then
    write(*,'(a)',advance='no') 'SIMPLE_PRIME2D_INIT stk=<stack.ext> smpd=<sampling distance(in A)>'
    write(*,'(a)',advance='no') ' ncls=<nr of clusters> [nthr=<nr of OpenMP threads{1}>]'
    write(*,'(a)',advance='no') ' [oritab=<input doc>] [filwidth=<filament width(in A)>]'
    write(*,'(a)',advance='no') ' [ctf=<yes|no|flip|mul{no}>] [kv=<acceleration voltage(in kV){300.}>]'
    write(*,'(a)',advance='no') ' [cs=<spherical aberration constant(in mm){2.7}>]'
    write(*,'(a)') ' [fraca=<frac amp contrast{0.07}>] [deftab=<defocus info file>]'
    stop
endif
call cline%parse
call cline%checkvar('stk',  1)
call cline%checkvar('smpd', 2)
call cline%checkvar('ncls', 3)
if( .not. cline%defined('eo') )then
    call cline%set('eo', 'no')
endif
call cline%check
p = params(cline)                   ! parameters generated
call b%build_general_tbox(p, cline) ! general objects built
call b%build_hadamard_prime2D_tbox(p)
write(*,'(a)') '>>> GENERATING INITIAL CLUSTER CENTERS'
if( cline%defined('oritab') )then
    call b%a%remap_classes
    ncls_in_oritab = b%a%get_ncls()
    if( p%ncls < ncls_in_oritab ) stop 'Inputted ncls < ncls_in_oritab; not allowed!'
    if( p%ncls > ncls_in_oritab )then
        call b%a%expand_classes(p%ncls)
    endif
else
    if( p%srch_inpl .eq. 'yes' )then
        call b%a%rnd_cls(p%ncls)
    else
        call b%a%rnd_cls(p%ncls, srch_inpl=.false.)
    endif
endif
p%oritab = 'prime2D_startdoc.txt'
call b%a%write(p%oritab)
if( cline%defined('filwidth') )then
    do icls=1,p%ncls
        call b%cavgs(icls)%bin_filament(p%filwidth)
    end do
else
    call prime2D_assemble_sums(b, p, mul=p%mul)
endif
call prime2D_write_sums(b, p)
call simple_end('**** SIMPLE_PRIME2D_INIT NORMAL STOP ****')
end program simple_prime2D_init
