!==Program simple_orisops
!
! <orisops/begin> is a program for analyzing SIMPLE orientation/parameter files (text files containing input parameters and/or 
! parameters estimated by \prgname{simple\_prime2D} or \prgname{simple\_prime3D}. <comment/begin> If two orientation tables 
! (\texttt{oritab} and \texttt{oritab2}) are inputted, the program provides statistics of the distances between the orientations 
! in the two documents. These statistics include the sum of angular distances between the orientations, the average angular 
! distance between the orientations, the standard deviation of angular distances, the minimum angular distance, and the maximum 
! angular distance. If only \texttt{oritab} is inputted, there are a few options available. If \texttt{errify=yes}, then the 
! program introduces uniform random angular errors $\in{[\texttt{-angerr,angerr}]}$, and uniform origin shift errors 
! $\in{[\texttt{-sherr,sherr}]}$, and uniform random defocus errors $\in{[\texttt{-deferr,deferr}]}$. If $\texttt{nstates}>1$ 
! then random states are assigned $\in{[1,\texttt{nstates}]}$. If \texttt{mirr=2d}, then the Euler angles in \texttt{oritab} 
! are mirrored according to the relation \texttt{e1=e1, e2=180.+e2, e3=-e3}. If \texttt{mirr=3d}, then the Euler angles in 
! \texttt{oritab} are mirrored according to the relation $R=M(M\cdot{}R)$, where $R$ is the rotation matrix calculated from 
! the Euler angle triplet and $M$ is a 3D reflection matrix (like a unit matrix but with the 3,3-component sign swapped). If 
! \texttt{e1}, \texttt{e2}, or \texttt{e3} is inputted, the orientations in \texttt{oritab} are rotated correspondingly. If 
! you input \texttt{state} as well, you rotate \textit{only} the orientations assigned to state \texttt{state}. If \texttt{mul} 
! is defined, you multiply the origin shifts with \texttt{mul}. If \texttt{zero=yes}, then the shifts are zeroed. If none of 
! the above described parameter are defined, and \texttt{oritab} is still defined, the program projects the 3D orientation into 
! the xy-plane and plots the resulting vector (this is useful for checking orientation coverage). If \texttt{oritab} is not 
! defined, the program generates random Euler angles $e1\in{[0,360]}$, $e2\in{[0,180]}$, and $e3\in{[0,360]}$ and random origin 
! shifts $x\in{[\texttt{-trs,yrs}]}$ and $y\in{[\texttt{-trs,yrs}]}$. If \texttt{ndiscrete} is set to an integer number > 0, the 
! orientations produced are randomly sampled from the set of \texttt{ndiscrete} quasi-even projection directions, and the in-plane 
! parameters are assigned randomly, as described above. If \texttt{even=yes}, then all \texttt{nptcls} orientations are assigned 
! quasi-even projection directions, and random in-plane parameters. If \texttt{nstates} is set to some integer number > 0, then 
! states are assigned randomly $\in{[1,\texttt{nstates}]}$. If \texttt{zero=yes} in this mode of execution, the projection 
! directions are zeroed and only the in-plane parameters are kept intact. If \texttt{errify=yes} and \texttt{astigerr} defined, 
! then uniform random astigmatism errors are introduced $\in{[\texttt{-astigerr,astigerr}]}$. <comment/end> <orisops/end>
!
! The code is distributed with the hope that it will be useful, but _WITHOUT_ _ANY_ _WARRANTY_.
! Redistribution and modification is regulated by the GNU General Public License.
! *Author:* Hans Elmlund 2015
!
program simple_orisops
use simple_cmdline, only: cmdline
use simple_build,   only: build
use simple_params,  only: params
use simple_ori,     only: ori
use simple_oris,    only: oris
use simple_jiffys,  only: simple_end, nlines
use simple_math,    only: normvec
implicit none
type(build)   :: b
type(ori)     :: orientation
type(oris)    :: o
type(params)  :: p
type(cmdline) :: cline
real          :: normal(3), mind, maxd, avgd, sdevd, sumd, vard
real          :: mind2, maxd2, avgd2, sdevd2, vard2, homo_cnt, homo_avg
integer       :: s, i, j, cnt
logical       :: err 
if( command_argument_count() < 1 )then
    write(*,'(a)', advance='no') 'SIMPLE_ORISOPS [nptcls=<number of oris>] [oritab=<input alignment doc>]'
    write(*,'(a)', advance='no') ' [oritab2=<2nd input alignment doc>] [plot=<yes|no{no}>]'
    write(*,'(a)', advance='no') ' [outfile=<output alignment doc>] [e1=<1st Euler{0}>] [e2=<2nd Euler{0}>]'
    write(*,'(a)', advance='no') ' [e3=<3d Euler{0}>] [trs=<origin shift(in pixels){0}>] [nstates=<nr of states{1}>]'
    write(*,'(a)', advance='no') ' [pgrp=<cn|dn|t|o|i{c1}>] [ctf=<yes|no{no}>] [kv=<acceleration voltage(in kV){300.}>]'
    write(*,'(a)', advance='no') ' [fraca=<frac amp contrast{0.07}>] [cs=<spherical aberration constant(in mm){2.7}>]'
    write(*,'(a)', advance='no') ' [defocus=<defocus(in microns){3.0}>] [deftab=<text file defocus values>]'
    write(*,'(a)', advance='no') ' [angerr=<angular error(in degrees){0}>] [sherr=<shift error(in pixels){0}>]'
    write(*,'(a)', advance='no') ' [deferr=<defocuserror(in microns){1.0}>] [astigerr=<astigmatism error(in microns){0.}>]'
    write(*,'(a)', advance='no') ' [even=<yes|no{no}>] [zero=<yes|no{no}>] [discrete=<yes|no{no}>] [ndiscrete=<nr of'
    write(*,'(a)', advance='no') ' discrete orientations>] [state=<state 2 process>] [errify=<yes|no{no}>] [deftab=<text'
    write(*,'(a)', advance='no') ' file defocus values>] [mul=<shift multiplication factor{1}>] [mirr=<no|2d|3d{no}>]'
    write(*,'(a)', advance='no') ' [xsh=<x shift(pixels){0}>] [ysh=<y shift(pixels){0}>]'
    write(*,'(a)', advance='no') ' [zsh=<z shift(pixels){0}>] [ctfstats=<yes|no{no}>] [trsstats=<yes|no{no}>]'
    write(*,'(a)', advance='no') ' [hist=<var2plot>] [ncls=<number of clusters>] [minp=<cluster population{10}>]'
    write(*,'(a)') ' [clustvalid=<yes|homo|no|{no}>] [thres=<homogeneity treshold{0.7}>]'
    stop
endif
call cline%parse
call cline%set('prg', 'orisops')
if( .not. cline%defined('thres') ) call cline%set('thres', 0.7)
p = params(cline)
call b%build_general_tbox(p, cline)
if( cline%defined('oritab2') )then
    if( .not. cline%defined('oritab') ) stop 'need oritab for comparison'
    if( nlines(p%oritab) .ne. nlines(p%oritab2) )then
        stop 'inconsistent number of lines in the two oritabs!'
    endif
    o = oris(p%nptcls)
    call o%read(p%oritab2)
    call b%a%diststat(o, sumd, avgd, sdevd, mind, maxd)
    write(*,'(a,1x,f15.6)') 'SUM OF ANGULAR DISTANCE BETWEEN ORIENTATIONS  :', sumd
    write(*,'(a,1x,f15.6)') 'AVERAGE ANGULAR DISTANCE BETWEEN ORIENTATIONS :', avgd
    write(*,'(a,1x,f15.6)') 'STANDARD DEVIATION OF ANGULAR DISTANCES       :', sdevd
    write(*,'(a,1x,f15.6)') 'MINIMUM ANGULAR DISTANCE                      :', mind
    write(*,'(a,1x,f15.6)') 'MAXIMUM ANGULAR DISTANCE                      :', maxd
else if( cline%defined('oritab') )then
    if( cline%defined('hist') )then
        call b%a%histogram(p%hist)
        goto 999
    endif
    if( p%ctfstats .eq. 'yes' )then
        call b%a%stats('ctfres', avgd, sdevd, vard, err )
        call b%a%minmax('ctfres', mind, maxd)
        write(*,'(a,1x,f8.2)') 'AVERAGE CTF RESOLUTION               :', avgd
        write(*,'(a,1x,f8.2)') 'STANDARD DEVIATION OF CTF RESOLUTION :', sdevd
        write(*,'(a,1x,f8.2)') 'MINIMUM CTF RESOLUTION (BEST)        :', mind
        write(*,'(a,1x,f8.2)') 'MAXIMUM CTF RESOLUTION (WORST)       :', maxd
        call b%a%stats('dfx', avgd, sdevd, vard, err )
        call b%a%minmax('dfx', mind, maxd)
        call b%a%stats('dfy', avgd2, sdevd2, vard2, err )
        call b%a%minmax('dfy', mind2, maxd2)
        write(*,'(a,1x,f8.2)') 'AVERAGE DF                           :', (avgd+avgd2)/2.
        write(*,'(a,1x,f8.2)') 'STANDARD DEVIATION OF DF             :', (sdevd+sdevd2)/2.
        write(*,'(a,1x,f8.2)') 'MINIMUM DF                           :', (mind+mind2)/2.
        write(*,'(a,1x,f8.2)') 'MAXIMUM DF                           :', (maxd+maxd2)/2.
        goto 999
    endif
    if( p%trsstats .eq. 'yes' )then
        call b%a%stats('x', avgd, sdevd, vard, err )
        call b%a%minmax('x', mind, maxd)
        call b%a%stats('y', avgd2, sdevd2, vard2, err )
        call b%a%minmax('y', mind2, maxd2)
        write(*,'(a,1x,f8.2)') 'AVERAGE TRS               :', (avgd+avgd2)/2.
        write(*,'(a,1x,f8.2)') 'STANDARD DEVIATION OF TRS :', (sdevd+sdevd2)/2.
        write(*,'(a,1x,f8.2)') 'MINIMUM TRS               :', (mind+mind2)/2.
        write(*,'(a,1x,f8.2)') 'MAXIMUM TRS               :', (maxd+maxd2)/2.
        goto 999
    endif
    if( p%errify .eq. 'yes' )then   ! introduce error in input orientations
        if( cline%defined('angerr') .or. cline%defined('sherr') ) call b%a%introd_alig_err(p%angerr, p%sherr)
        if( p%ctf .eq. 'yes' ) call b%a%introd_ctf_err(p%deferr)
    endif
    if( p%mirr .eq. '2d' )then      ! mirror input Eulers
        call b%a%mirror2d
    endif
    if( p%mirr .eq. '3d' )then      ! mirror input Eulers
        call b%a%mirror3d
    endif
    if( cline%defined('e1') )then ! rotate input Eulers
        call orientation%new
        call orientation%e1set(p%e1)
        call orientation%e2set(p%e2)
        call orientation%e3set(p%e3) 
        if( cline%defined('state') )then
            do i=1,b%a%get_noris()
                s = nint(b%a%get(i, 'state'))
                if( s == p%state )then
                    call b%a%rot(i,orientation)
                endif
            end do
        else
            call b%a%rot(orientation)
        endif
    endif
    if( cline%defined('mul') )then
        call b%a%mul_shifts(p%mul)
    endif
    if( p%zero  .eq. 'yes' ) call b%a%zero_shifts
    if( p%plot  .eq. 'yes' )then ! plot polar vectors                          
        do i=1,b%a%get_noris()
            normal = b%a%get_normal(i)
            write(*,'(1x,f7.2,3x,f7.2)') normal(1), normal(2)
        end do
    endif
    if( p%discrete .eq. 'yes' )then
        if( cline%defined('ndiscrete') )then
            call b%a%discretize(p%ndiscrete)
        else
            stop 'need ndiscrete to be defined!'
        endif
    endif
    if( cline%defined('xsh') )then
        call b%a%map3dshift22d([p%xsh,p%ysh,p%zsh])
    endif
    if( p%clustvalid .eq. 'yes' )then
        if( cline%defined('ncls') )then
            write(*,'(a,3x,f5.1)') '>>> COHESION: ',   b%a%cohesion_norm('class',p%ncls)*100.
            write(*,'(a,1x,f5.1)') '>>> SEPARATION: ', b%a%separation_norm('class',p%ncls)*100.
        else if( cline%defined('nstates') )then
            write(*,'(a,3x,f5.1)') '>>> COHESION: ',   b%a%cohesion_norm('state',p%nstates)*100.
            write(*,'(a,1x,f5.1)') '>>> SEPARATION: ', b%a%separation_norm('state',p%nstates)*100.
        else
            stop 'need ncls/nstates as input for clustvalid'
        endif
    else if( p%clustvalid .eq. 'homo' )then
        if( cline%defined('ncls') )then
            call b%a%homogeneity('class', p%minp, p%thres, homo_cnt, homo_avg)
            write(*,'(a,1x,f5.1)') '>>> THIS % OF CLUSTERS CONSIDERED HOMOGENEOUS: ', homo_cnt*100.
            write(*,'(a,1x,f5.1)') '>>> AVERAGE HOMOGENEITY:                       ', homo_avg*100.
        else if( cline%defined('nstates') )then
            call b%a%homogeneity('state', p%minp, p%thres, homo_cnt, homo_avg)
            write(*,'(a,1x,f5.1)') '>>> THIS % OF CLUSTERS CONSIDERED HOMOGENEOUS: ', homo_cnt*100.
            write(*,'(a,13x,f5.1)') '>>> AVERAGE HOMOGENEITY:                      ', homo_avg*100.
        else
            stop 'need ncls/nstates as input for clustvalid'
        endif
    endif
    if( cline%defined('nstates') )then
        call b%a%rnd_states(p%nstates)
    endif
else ! make orientations
    if( cline%defined('ncls') )then
        o = oris(p%ncls)
        call o%spiral(p%nsym, p%eullims)
        call b%a%new(p%ncls*p%minp)
        cnt = 0
        do i=1,p%ncls
            orientation = o%get_ori(i)
            do j=1,p%minp
                cnt = cnt+1
                call b%a%set_ori(cnt, orientation)
            end do
        end do
        if( p%zero .ne. 'yes' ) call b%a%rnd_inpls(p%trs)   
    else if( cline%defined('ndiscrete') )then
        if( p%ndiscrete > 0 )then
            call b%a%rnd_oris_discrete(p%ndiscrete, p%nsym, p%eullims)
        endif
        call b%a%rnd_inpls(p%trs)
    else if( p%even .eq. 'yes' )then
        call b%a%spiral(p%nsym, p%eullims)
        call b%a%rnd_inpls(p%trs)
    else
        call b%a%rnd_oris(p%trs) 
    endif
    if( p%nstates > 1 ) call b%a%rnd_states(p%nstates)
    if( cline%defined('astigerr') )then
        if( p%ctf .eq. 'yes' ) call b%a%rnd_ctf(p%defocus, p%deferr, p%astigerr)
    else
        if( p%ctf .eq. 'yes' ) call b%a%rnd_ctf(p%defocus, p%deferr)
    endif
endif
call b%a%write(p%outfile)
999 call simple_end('**** SIMPLE_ORISOPS NORMAL STOP ****')
end program simple_orisops