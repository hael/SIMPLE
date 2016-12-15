!==Program simple_resrange
!
! <resrange/begin> is a program for estimating the resolution range used in the heuristic resolution-stepping
! scheme in the PRIME3D initial model production procedure. The initial low-pass limit is set so that each image 
! receives ten nonzero orientation weights. When quasi-convergence has been reached, the limit is updated one 
! Fourier index at the time, until PRIME reaches the condition where six nonzero orientation weights are assigned 
! to each image. FSC-based filtering is unfortunately not possible to do in the \textit{ab initio} reconstruction 
! step, because when the orientations are mostly random, the FSC overestimates the resolution. This program is 
! used internally when executing PRIME in distributed mode. We advise you to check 
! the starting and stopping low-pass limits before executing PRIME3D using this program. <comment/begin>
! The resolution range estimate depends on the molecular diameter, which is estimated based on the box size. 
! If you want to override this estimate, set \texttt{moldiam} to the desired value (in \AA{}). This may be necessary 
! if your images have a lot of background "padding". However, for 
! starting model generation it is probably better to clip the images snugly around the particle, because smaller 
! images equal less computation. <comment/end> <resrange/end>
!
! The code is distributed with the hope that it will be useful, but _WITHOUT_ _ANY_ _WARRANTY_.
! Redistribution and modification is regulated by the GNU General Public License.
! *Author:* Hans Elmlund 2015
!
program simple_resrange
use simple_cmdline,            only: cmdline
use simple_hadamard3D_matcher, only: prime3D_find_resrange
use simple_params,             only: params
use simple_build,              only: build
use simple_jiffys,             only: simple_end
use simple_timing
implicit none
type(params)  :: p
type(build)   :: b
type(cmdline) :: cline
!start of the execution commands
call timestamp()
call start_Alltimers_cpu()
if( command_argument_count() < 1 )then
    write(*,'(a)',advance='no') 'SIMPLE_RESRANGE smpd=<sampling distance(in A)>'
    write(*,'(a)',advance='no') ' [nspace=<nr of reference sections{1000}>] [pgrp=<cn|dn|t|o|i{c1}>]'
    write(*,'(a)') ' [box=<image size(in pixels)>] [moldiam=<molecular diameter(in A))>]'
    stop
endif
call cline%parse
call cline%checkvar('smpd', 1)
call cline%check
call cline%set('prg', 'resrange')
if( .not. cline%defined('nspace') ) call cline%set('nspace', 1000.)
if( cline%defined('box') .or. cline%defined('moldiam') )then
    p = params(cline)                     ! parameters generated
    call b%build_general_tbox(p, cline)   ! general objects built
    call b%build_hadamard_prime3D_tbox(p) ! prime objects built
    call prime3D_find_resrange( b, p, p%lp, p%lpstop )
    write(*,'(A,1X,F5.1)') '>>> LP START:', p%lp
    write(*,'(A,2X,F5.1)') '>>> LP STOP:', p%lpstop
    write(*,'(A,2X,F5.1)') '>>> HP:', p%hp
else
    stop 'need either box size or moldiam to estimate resrange; simple_resrange'
endif
! END GRACEFULLY
call simple_end('**** SIMPLE_RESRANGE NORMAL STOP ****')
!*******************************************************************************
!    Environment shutdown
!
!*******************************************************************************
!shutting down the timers
call stop_Alltimers_cpu()

end program
