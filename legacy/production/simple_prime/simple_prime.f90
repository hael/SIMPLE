!==Program simple_prime
!
! <prime/begin> is an \textit{ab inito} reconstruction/low-resolution refinement program based on probabilistic projection matching. 
! PRIME is shorthand for PRobabilistic Initial 3D Model Generation for Single-Particle Cryo-Electron Microscopy. If you process 
! images of a molecule with a diameter of 200 \AA{} it should be possible to obtain an 8 \AA map using 1000 reference sections 
! (the default setting in PRIME) provided that the images are of sufficient quality, they are many enough, and they are sampled 
! finely enough. We do \textit{not} recommend using more than 1000 reference sections in attempt to push for high-resolution. 
! Then it is more efficient to use our new continuous probabilistic refinement code, implemented by the \prgname{simple\_oasis} 
! program, described below. \prgname{simple\_oasis} is at least ten times faster than PRIME for high-resolution refinement. We do 
! not recommend using PRIME for heterogeneity analysis (starting off with many random blobs), because there are more effective ways 
! to deal with the heterogeneity problem. We will address the heterogeneity problem in the next SIMPLE 2.2 release. However, if you 
! suspect that you have heterogeneity of the kind where totally different species co-exist it could be worth trying initialisation 
! with a few blobs. You should use phase-flipped images for initial model production with PRIME (phase flipping can be done with 
! \prgname{simple\_stackops}). Do \textit{not} search the origin shifts initially, when the model is of very low quality. If your 
! images are far off centre, use \prgname{simple\_stackops} with option \texttt{shalgn=yes} instead to shiftalign the images 
! beforehand (the algorithm implemented is the same as EMAN's \texttt{cenalignint} program). We recommend running the first round 
! of PRIME with dynamic resolution stepping \texttt{dynlp=yes}. The \texttt{dynlp} option implements a heuristic resolution 
! weighting/update scheme. The initial low-pass limit is set so that each image receives ten nonzero orientation weights. When 
! quasi-convergence has been reached, the limit is updated one Fourier index at the time until PRIME reaches the condition where 
! six nonzero orientation weights is assigned to each image. FSC-based filtering is unfortunately not possible to do in the 
! \textit{ab initio} reconstruction step, because when the orientations are mostly random, the FSC overestimates the resolution.
! Once the initial model has converged, we recommend start searching the shifts (by setting \texttt{trs} to some nonzero value), 
! apply the FSC for resolution-weighting (by setting \texttt{eo=yes}), and turn on the Wiener restoration by setting 
! \texttt{ctf=yes|flip|mul} where \texttt{yes} instructs PRIME to take care of all CTF correction, \texttt{flip} indicates that 
! the images have been phase-flipped beforehand and \texttt{mul} indicates that the images have been multiplied with the CTF 
! beforehand. To use Wiener restoration you also need to input CTF parameters, for example via \texttt{deftab=defocus\_values.txt}. 
! Remember that the defocus values should be given in microns and the astigmatism angle in degrees (one row of the file 
! \texttt{defocus\_values.txt} may look like: \texttt{dfx=3.5  dfy=3.3  angast=20.0}). There are many ways of using (and probably 
! also abusing) \prgname{simple\_prime}. We will walk you through an examples (see section \ref{recsym}, below). Use 
! \prgname{distr\_simple.pl}, described below, for distributed PRIME execution. <comment/begin> Since the search is probabilistic, 
! we figured that an elegant convergence criterion could be formulated based on the variance of the distribution of orientations 
! assigned to each image. This works well for asymmetrical reconstructions, but for symmetrical reconstructions the variance 
! increases towards the end of the run, when the shape most consistent with the point group is being established. Note that we do 
! not assume any point-group symmetry in the initial runs. However, the \prgname{simple\_symsrch} program (described below) can 
! be used to align the reconstruction to its symmetry axis so that future searches can be restricted to the asymmetric unit in 
! refinement (see section \ref{recsym}, below). Less commonly used, and less obvious input parameters are \texttt{nspace}, which 
! controls the number of reference projections, \texttt{amsklp}, which controls the low-pass limit used in the automask routine, 
! \texttt{maxits}, which controls the maximum number of iterations executed, \texttt{pgrp}, which controls the point-group symmetry, 
! assuming that the starting volume is aligned to its principal symmetry axis, \texttt{edge}, which controls the size of the 
! softening edge in the automask routine. Setting \texttt{diversify=no} removes the stochastic component of the search and 
! evaluates the entire neighbourhood in every round. This would correspond to standard projection matching with weighted 
! orientation assignment. It is not a recommended option, but implemented for testing purposes. Setting \texttt{noise=yes} 
! produces a random noise starting volume instead of a random blob, if no input volume is given.<comment/end> <prime/end>
!
! The code is distributed with the hope that it will be useful, but _WITHOUT_ _ANY_ _WARRANTY_.
! Redistribution and modification is regulated by the GNU General Public License.
! *Author:* Hans Elmlund 2015
!
program simple_prime
use simple_defs     ! singleton
use simple_cmdline  ! singleton
use simple_jiffys,  only: simple_end
use simple_matcher, only: prime_exec, prime_find_resrange
use simple_params,  only: params
use simple_build,   only: build
implicit none
type(params) :: p
type(build)  :: b
integer      :: i, startit
logical      :: update_res=.false., converged=.false.
real         :: lpstart, lpstop
if( command_argument_count() < 3 )then
    write(*,'(a)',advance='no') 'SIMPLE_PRIME stk=<stack.ext> [vol1=<invol.ext>] [vol2=<refvol_2.ext> etc.]'
    write(*,'(a)',advance='no') ' smpd=<sampling distance(in A)> msk=<mask radius(in pixels)>'
    write(*,'(a)',advance='no') ' [oritab=<previous alignment doc>] [trs=<origin shift(in pixels){0}>]'
    write(*,'(a)',advance='no') ' [lp=<low-pass limit{20}>] [dynlp=<yes|no{yes}>]'
    write(*,'(a)',advance='no') ' [nstates=<nstates to reconstruct>] [frac=<fraction of ptcls to include{1}>]'
    write(*,'(a)',advance='no') ' [mw=<molecular weight (in kD)>] [nthr=<nr of OpenMP threads{1}>]'
    write(*,'(a)',advance='no') ' [startit=<start iteration>]'
    write(*,'(a)',advance='no') ' [lpstop=<stay at this low-pass limit (in A)>] [deftab=<text file defocus values>]'
    write(*,'(a)',advance='no') ' [nspace=<nr reference sections{1000}>] [eo=<yes|no{no}>]' 
    write(*,'(a)',advance='no') ' [amsklp=<automask low-pass limit(in A)>] [pgrp=<cn|dn|t|o|i{c1}>] [ctf=<yes|no|'
    write(*,'(a)',advance='no') 'flip|mul{no}>] [kv=<acceleration voltage(in kV){300.}>] [cs=<spherical'
    write(*,'(a)') ' aberration constant(in mm){2.7}>] [fraca=<frac amp contrast{0.07}>] [hp=<high-pass limit(in A)>]'
    write(*,'(a)') '**less commonly used**'
    write(*,'(a)',advance='no') ' [maxits=<max iterations{100}>] [fracvalid=<fraction of particles 4 validation{0.}>] '
    write(*,'(a)',advance='no') ' [lpvalid=<low-pass limit validptcls{20}>] [edge=<edge size for softening molecular '
    write(*,'(a)',advance='no') 'envelope(in pixels){3}>] [find=<Fourier index>] [fstep=<Fourier step size>]'
    write(*,'(a)',advance='no') ' [diversify=<yes|no{yes}>] [noise=<yes|no{no}>]'  
    write(*,'(a)',advance='no') ' [time_per_image=<{100}>] [dens=<density(e.g. 9.368 Da/A3 4 gold clusters){0.}>]'
    write(*,'(a)',advance='no') ' [trsstep=<origin shift stepsize{1}>] [outfile=<output alignment doc 4 parallell jobs>]'
    write(*,'(a)',advance='no') ' [nvox=<nr of voxels in mask{0}>] [inner=<inner mask'
    write(*,'(a)') ' radius(in pixels)>] [width=<pixels falloff inner mask{10}>] [norec=<yes|no{no}>]'
    stop
endif
call parse_cmdline
call cmdcheckvar('stk',  1)
call cmdcheckvar('smpd', 2)
call cmdcheckvar('msk',  3)
if( .not. defined_cmd_arg('nspace') )then
    call set_cmdline('nspace', 1000.)
endif
if( defined_cmd_arg('lp') .or. defined_cmd_arg('find')) call set_cmdline('dynlp', 'no')
call set_cmdline('refine', 'no')
call cmdcheck
p = params()                 ! parameters generated
call b%build_general_tbox(p) ! general objects built
if( .not. defined_cmd_arg('eo') ) p%eo = 'no' ! default
if( p%eo .eq. 'yes' ) p%dynlp = 'no'
if( defined_cmd_arg('lp') .or. (defined_cmd_arg('find') .or. p%eo .eq. 'yes') )then
    ! alles ok!
else
    stop 'need a starting low-pass limit (set lp or find)!'
endif
if( p%oritab .eq. '' .and. p%vols(1) .ne. '' )then ! set default deterministic
    p%intensify = 'yes' ! if no oris but vol do exhaustive search
    p%diversify = 'no'
    if( defined_cmd_arg('lp') .or. (defined_cmd_arg('find') .or. p%eo .eq. 'yes') )then
        ! ok
    else
        write(*,'(a)') 'WARNING! You are using intensifying search with automatic'
        write(*,'(a)') 'resolution stepping, probably not what you want!'
    endif
endif
! build prime objects
call b%build_prime_tbox(p) 
if( defined_cmd_arg('part') )then
    if( .not. defined_cmd_arg('outfile') ) stop 'need unique output file for parallel jobs'
    if( defined_cmd_arg('find') ) p%lp = b%img%get_lp(p%find)
    call prime_exec(b, p, 0, update_res, converged) ! partition or not, depending on 'part'       
else
    if( p%dynlp .eq. 'yes' )then
        call prime_find_resrange( b, p, lpstart, lpstop ) ! determine resolution range
        if( .not. defined_cmd_arg('lpstop') )  p%lpstop = lpstop
        if( defined_cmd_arg('lpstart') )then
            p%lp = p%lpstart
        else
            p%lp = lpstart
        endif
        if( p%lpstop >= p%lp ) stop 'lpstop of lower resolution than the limit (lp)!'
    endif
    p%find = int((real(p%box-1)*p%smpd)/p%lp)
    startit = 1
    if( defined_cmd_arg('startit') ) startit = p%startit
    do i=startit,p%maxits
        call prime_exec(b, p, i, update_res, converged)
        if( update_res )then
            p%find = p%find+p%fstep
            p%lp = max(p%lpstop,b%img%get_lp(p%find))
        endif
        if(converged) exit
    end do
endif
call simple_end('**** SIMPLE_PRIME NORMAL STOP ****')
end program simple_prime