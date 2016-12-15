!==Program simple_symsrch
!
! <symsrch/begin> is a program for searching for the principal symmetry axis of a volume 
! reconstructed without assuming any point-group symmetry. Our philosophy is to start off without assuming any 
! symmetry, and then analyzse the reconstructed volume to identify the correct point-group symmetry. 
! \prgname{simple\_symsrch} can be used for this purpose. The program takes as input the asymmetrical reconstruction 
! (obtained with \prgname{simple\_prime3D}), the alignment document for all the particle images that 
! have gone into the reconstruction, and the desired point-group symmetry. It then projects the reconstruction in 20 
! (default option) even directions, uses common lines-based optimisation to identify the principal symmetry 
! axis, applies the rotational transformation to the inputted orientations, and produces a new alignment document. 
! Input this document to \texttt{simple\_recvol} or \texttt{simple\_eo\_recvol} together with the images and the 
! point-group symmetry to generate a symmetrised map. If you are unsure about the point-group, you should of course 
! test many different point-groups and compare the asymmetric map with the symmetrised maps. SIMPLE now implements 
! most point-groups: c- and d-groups, as well as tetrahedral, octahedral, and icosahedral groups. <comment/begin>
! The \texttt{state} parameter allows you to apply symmetry for the given state. <comment/end> <symsrch/end>
!
! The code is distributed with the hope that it will be useful, but _WITHOUT_ _ANY_ _WARRANTY_.
! Redistribution and modification is regulated by the GNU General Public License.
! *Author:* Hans Elmlund 2015
!
program simple_symsrch
use simple_jiffys   ! singleton
use simple_cmdline,   only: cmdline
use simple_oris,      only: oris
use simple_params,    only: params
use simple_build,     only: build
use simple_symsrcher, only: symsrch_master
implicit none
type(params), target :: p
type(build), target  :: b
type(cmdline)        :: cline
type(oris)           :: o
if( command_argument_count() < 4 )then
    write(*,'(a)', advance='no') 'SIMPLE_SYMSRCH [vol1=<vol.ext>] [stk=<ptcls.ext>] smpd=<sampling distance(in A)>'
    write(*,'(a)', advance='no') ' [oritab=<input alignment doc>] pgrp=<cn|dn|t|o|i{c1}> outfile=<output alignment doc>'
    write(*,'(a)', advance='no') ' lp=<low-pass limit(in A)> [amsklp=<low-pass limit for centering mask(in A){50}>]'
    write(*,'(a)', advance='no') ' [hp=<high-pass limit(in A)>] [nthr=<nr openMP threads{1}>] [nspace=<nr of projs{20}>]'
    write(*,'(a)') ' [compare=<yes|no{no}>]'
    stop
endif
call cline%parse
if( cline%defined('vol1') .or. cline%defined('stk') )then
    ! all ok
else
    stop 'Either vol1 or stk needs to be defined on the command line!'
endif
if( cline%defined('vol1') )then
    if( .not. cline%defined('oritab') ) stop 'Need oritab to go with the volume!'
endif
call cline%checkvar('smpd',    1)
call cline%checkvar('pgrp',    2)
call cline%checkvar('outfile', 3)
call cline%checkvar('lp',      4)
if( .not. cline%defined('nspace') )then
    if( cline%defined('vol1') )then
        call cline%set('nptcls', 20.) ! 20 projections 4 symsrch
        call cline%set('nspace', 20.) ! 20 projections 4 symsrch
    else
        call cline%set('nptcls', 1.) ! single-particle molecule search
        call cline%set('nspace', 1.) ! single-particle molecule search
    endif
else
    call cline%set('nptcls', cline%get_rarg('nspace'))
endif
if( .not. cline%defined('amsklp')  ) call cline%set('amsklp',   50.)
if( .not. cline%defined('compare') ) call cline%set('compare', 'no')
call cline%check ! checks args and prints cmdline to cmdline.dat
p = params(cline)  ! constants & derived constants produced
if( cline%defined('stk') )then
    p%nptcls = 1
    p%nspace = 1
endif
call b%build_general_tbox(p, cline, .true., .true.) ! general objects built (no oritab reading)
call b%build_comlin_tbox(p)                         ! objects for common lines based alignment built
call symsrch_master( cline, p, b, o )
call o%write(p%outfile)
call simple_end('**** SIMPLE_SYMSRCH NORMAL STOP ****')
end program simple_symsrch