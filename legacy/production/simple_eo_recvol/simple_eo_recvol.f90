!==Program simple_eo_recvol
!
! <eo_recvol/begin> is a program for reconstructing volumes from MRC or SPIDER stacks, given input orientations and state 
! assignments (obtained by program \prgname{simple\_prime3D}). The algorithm is based on direct 
! Fourier inversion with a Kaiser-Bessel (KB) interpolation kernel. This window function reduces the real-space ripple 
! artefacts associated with direct moving windowed-sinc interpolation. The feature sought when implementing this algorithm
! was to enable quick, reliable reconstruction from aligned individual particle images. The even and odd pairs are 
! automatically reconstructed, the FSC calculated, and the Wiener filter formalism used for image restoration (CTF 
! correction). Use \prgname{distr\_simple.pl} for distributed execution. <comment/begin> 
! \texttt{mul} is used to scale the origin shifts if down-sampled images were used for alignment 
! and the original images are used for reconstruction. \texttt{ctf}, \texttt{kv}, \texttt{fraca},
! \texttt{cs} and \texttt{deftab} are used to communicate CTF information to the program. \texttt{ctf=yes}, \texttt{ctf=flip} or
! \texttt{ctf=mul} turns on the Wiener restoration. If you input CTF info to the program, please ensure that the correct kV,
! Cs and fraca (fraction of amplitude contrast) parameters are inputted as well. If the images were pre-multiplied with the CTF, 
! set \texttt{ctf=mul} or if the images were phase-flipped set \texttt{ctf=flip}. \texttt{amsklp} and \texttt{mw} 
! parameters control the solvent mask: the low-pass limit used to generate the envelope; the molecular weight of the 
! molecule (protein assumed but it works reasonably well also for RNA; slight modification of \texttt{mw} might be needed). 
! The \texttt{inner} parameter controls the radius of the soft-edged mask used to remove the unordered DNA/RNA core of spherical 
! icosahedral viruses. <comment/end> <eo_recvol/end>
!
! The code is distributed with the hope that it will be useful, but _WITHOUT_ _ANY_ _WARRANTY_.
! Redistribution and modification is regulated by the GNU General Public License.
! *Author:* Hans Elmlund 2011
!
program simple_eo_recvol
use simple_build,      only: build
use simple_params,     only: params
use simple_rec_master, only: exec_eorec
use simple_jiffys,     only: simple_end
use simple_cmdline,    only: cmdline
use simple_defs        ! singleton
implicit none
type(params)  :: p
type(build)   :: b
type(cmdline) :: cline
if( command_argument_count() < 4 )then
    write(*,'(a)', advance='no') 'SIMPLE_EO_RECVOL stk=<ptcls.ext>'
    write(*,'(a)', advance='no') ' smpd=<sampling distance(in A)> oritab=<algndoc.txt> [frac=<fraction ptcls to'
    write(*,'(a)') ' include{1.}>] [mw=<molecular weight(in kD)>] [nthr=<nr openMP threads{1}>] [pgrp=<cn|dn|t|o|i{c1}>]'
    write(*,'(a)') '**less commonly used**'
    write(*,'(a)', advance='no') '[mul=<shift multiplication factor{1}>]'
    write(*,'(a)', advance='no') ' [ctf=<yes|no|flip|mul{no}>] [kv=<acceleration voltage(in kV){300.}>]'
    write(*,'(a)', advance='no') ' [fraca=<frac amp contrast{0.07}>] [cs=<spherical aberration constant(in mm){2.7}>]'
    write(*,'(a)') ' [deftab=<text file defocus values>] [state=<state to reconstruct{all}>]'
    stop
endif
call cline%parse
call cline%checkvar('stk',    1)
call cline%checkvar('smpd',   2)
call cline%checkvar('oritab', 3)
call cline%set('trs',    5.) ! to assure that shifts are being used
call cline%check
call cline%set('prg', 'eo_recvol')
p = params(cline)            ! constants & derived constants produced
call b%build_general_tbox(p, cline)
call b%build_eo_rec_tbox(p)
call exec_eorec(b, p, cline)
call simple_end('**** SIMPLE_EO_RECVOL NORMAL STOP ****')    
end program simple_eo_recvol