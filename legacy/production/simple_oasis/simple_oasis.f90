!==Program simple_oasis
!
! <oasis/begin> is our new refinement code, implementing continuous probabilistic orientation search. The OASIS (for Optimization 
! by Adaptive Sampling around Important Solutions) implements a probabilistic model of all parameters subject to estimation as the 
! basis for the orientation search. The stochastic optimiser, used to estimate the values of the latent variables of the statistical 
! model, displays much better general function optimisation performance than the BFGS algorithm, Powell's direction set method, or 
! the downhill simplex method due to Nelder and Mead. The OASIS search algorithm requires orders of magnitude less computations 
! than traditional global search techniques, such as simulated annealing or genetic algorithms. The OASIS code was developed for 
! high-resolution refinement and heterogeneity analysis. Although heterogeneity analysis is implemented, it is considered experimental 
! at this stage, since it has not been fully tested. <oasis/end>
!
! The code is distributed with the hope that it will be useful, but _WITHOUT_ _ANY_ _WARRANTY_.
! Redistribution and modification is regulated by the GNU General Public License.
! *Author:* Hans Elmlund 2015
!
program simple_oasis
use simple_defs     ! singleton
use simple_cmdline  ! singleton
use simple_jiffys,  only: simple_end
use simple_hadamard_matcher, only: prime_exec
use simple_params,  only: params
use simple_build,   only: build
implicit none
type(params) :: p
type(build)  :: b
integer      :: i, startit
logical      :: update_res=.false., converged=.false.
if( command_argument_count() < 3 )then
    write(*,'(a)',advance='no') 'SIMPLE_OASIS stk=<stack.ext> vol1=<invol.ext> [vol2=<refvol_2.ext> etc.]'
    write(*,'(a)',advance='no') ' smpd=<sampling distance(in A)> msk=<mask radius(in pixels)>'
    write(*,'(a)',advance='no') ' oritab=<previous alignment doc> [refine=<soft|het|stoch{soft}>]'
    write(*,'(a)',advance='no') ' [lp=<low-pass limit{20}>] [frac=<fraction of ptcls to include{1}>]'
    write(*,'(a)',advance='no') ' [mw=<molecular weight (in kD)>] [nthr=<nr of OpenMP threads{1}>]'
    write(*,'(a)',advance='no') ' [startit=<start iteration>] [lpstart=<start low-pass limit(in A)>]'
    write(*,'(a)',advance='no') ' [lpstop=<stay at this low-pass limit (in A)>] [deftab=<text file defocus values>]'
    write(*,'(a)',advance='no') ' [eo=<yes|no{yes}>] [amsklp=<automask low-pass limit(in A){20}>]' 
    write(*,'(a)',advance='no') ' [pgrp=<cn|dn|t|o|i{c1}>] [ctf=<yes|no|flip|mul{no}>]'
    write(*,'(a)',advance='no') ' [kv=<acceleration voltage(in kV){300.}>] [cs=<spherical'
    write(*,'(a)') ' aberration constant(in mm){2.7}>] [fraca=<frac amp contrast{0.07}>] [hp=<high-pass limit(in A)>]'
    write(*,'(a)') '**less commonly used**'
    write(*,'(a)',advance='no') '[nsample=<nr of stochastic search samples{100}>]'
    write(*,'(a)',advance='no') ' [maxits=<max iterations{100}>] [edge=<edge size for softening molecular'
    write(*,'(a)',advance='no') ' envelope(in pixels){3}>] [dens=<density(e.g. 9.368 Da/A3 4 gold clusters){0.}>]'
    write(*,'(a)',advance='no') ' [nvox=<nr of voxelsin mask{0}>] [inner=<inner mask radius(in pixels)>]'
    write(*,'(a)') ' [width=<pixels falloff inner mask{10}>] [norec=<yes|no{no}>]'
    stop
endif
call parse_cmdline
call cmdcheckvar('stk',    1)
call cmdcheckvar('vol1',   2)
call cmdcheckvar('smpd',   3)
call cmdcheckvar('msk',    4)
call cmdcheckvar('oritab', 5)
if( .not. defined_cmd_arg('refine')) call set_cmdline('refine', 'soft')
call set_cmdline('dynlp', 'no')
call cmdcheck
p = params()                 ! parameters generated
call b%build_general_tbox(p) ! general objects built
if( .not. defined_cmd_arg('eo') ) p%eo = 'yes' ! default
if( defined_cmd_arg('lp') .or. defined_cmd_arg('find') .or. p%eo .eq. 'yes' )then
    ! ok
else
    stop 'Need starting low-pass limit (lp) 4 refinement!'
endif
call b%build_hadamard_prime_tbox(p) ! prime objects built
startit = 1
if( defined_cmd_arg('startit') ) startit = p%startit
if( defined_cmd_arg('part') )then
    if( .not. defined_cmd_arg('outfile') ) stop 'need unique output file for parallel jobs'
    call prime_exec(b, p, 0, update_res, converged) ! partition or not, depending on 'part'       
else
    do i=startit,p%maxits
        call prime_exec(b, p, i, update_res, converged)
        if(converged) exit
    end do
endif
call simple_end('**** SIMPLE_OASIS NORMAL STOP ****')
end program simple_oasis