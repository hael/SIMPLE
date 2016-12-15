!==Program simple_recvol
!
! <recvol/begin> is a program for reconstructing volumes from MRC and SPIDER stacks, given input orientations and state 
! assignments (obtained by program \prgname{simple\_prime3D}). The algorithm is based on direct 
! Fourier inversion with a Kaiser-Bessel (KB) interpolation kernel. This window function reduces the real-space ripple 
! artifacts associated with direct moving windowed-sinc interpolation. The feature sought when implementing this algorithm 
! was to enable quick, reliable reconstruction from aligned individual particle images. <comment/begin> \texttt{mul} is 
! used to scale the origin shifts if down-sampled were used for alignment and the original images are used for reconstruction. 
! This program can be run in distributed mode using \prgname{distr\_simple.pl}.  
! \texttt{ctf}, \texttt{kv}, \texttt{fraca}, \texttt{cs} and \texttt{deftab} are used to communicate CTF information to the program. 
! \texttt{ctf=yes}, \texttt{ctf=flip} or \texttt{ctf=mul} turns on the Wiener restoration. If the images were pre-multiplied with 
! CTF set \texttt{ctf=mul} or if the images were phase-flipped set \texttt{ctf=flip}. \texttt{amsklp}, \texttt{mw}, and \texttt{edge} 
! are parameters that control the solvent mask: the low-pass limit used to generate the envelope; the molecular weight of the molecule 
! (protein assumed but it works reasonably well also for RNA; slight modification of \texttt{mw} might be needed). The \texttt{inner} 
! parameter controls the radius of the soft-edged mask used to remove the unordered DNA/RNA core of spherical icosahedral viruses. 
! The \texttt{even} and \texttt{odd} parameters allow you to reconstruct either the even or the odd pair. <comment/end> <recvol/end>
!
! The code is distributed with the hope that it will be useful, but _WITHOUT_ _ANY_ _WARRANTY_.
! Redistribution and modification is regulated by the GNU General Public License.
! *Author:* Hans Elmlund 2015
!

program simple_recvol
use simple_rec_master, only: exec_rec
use simple_build,      only: build
use simple_params,     only: params
use simple_jiffys,     only: simple_end
use simple_cmdline,    only: cmdline
use simple_defs        ! singleton
implicit none
type(params)  :: p
type(build)   :: b
type(cmdline) :: cline
if( command_argument_count() < 3 )then
    write(*,'(a)', advance='no') 'SIMPLE_RECVOL stk=<ptcls.ext> smpd=<sampling distance(in A)>'
    write(*,'(a)', advance='no') ' oritab=<algndoc.txt>'
    write(*,'(a)', advance='no') ' [frac=<fraction of ptcls to include{1.}>] [nthr=<nr of openMP threads{1}>]  '
    write(*,'(a)') '  [pgrp=<cn|dn|t|o|i{c1}>]'
    write(*,'(a)') '**less commonly used**'
    write(*,'(a)', advance='no') ' [mul=<shift multiplication factor{1}>] [ctf=<yes|no|flip|mul{no}>]'
    write(*,'(a)', advance='no') ' [kv=<acceleration voltage(in kV){300.}>] [fraca=<frac amp contrast{0.07}>]'
    write(*,'(a)', advance='no') ' [cs=<spherical aberration constant(in mm){2.7}>] [deftab=<text file defocus values>]'
    write(*,'(a)', advance='no') ' [state=<state to reconstruct{all}>] [inner=<inner mask radius(in pixels)>]'
    write(*,'(a)', advance='no') ' [width=<pixels falloff inner mask{10}>] [even=<yes|no{no}>] [odd=<yes|no{no}>]'
    write(*,'(a)') ' '
    stop
endif
call cline%parse
call cline%checkvar('stk',    1)
call cline%checkvar('smpd',   2)
call cline%checkvar('oritab', 3)
call cline%set('trs', 5.) ! to assure that shifts are being used
call cline%check
call cline%set('prg', 'recvol')
p = params(cline)         ! constants & derived constants produced
call b%build_general_tbox(p, cline)
call b%build_rec_tbox(p)
call exec_rec(b, p, cline)
call simple_end('**** SIMPLE_RECVOL NORMAL STOP ****')    
end program simple_recvol
