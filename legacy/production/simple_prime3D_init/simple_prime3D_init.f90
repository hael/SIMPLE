!==Program simple_prime3D_init
!
! <prime3D_init/begin> is a program for generating a random initial model for initialisation of PRIME3D when 
! executed in distributed mode. It is assumed that the images have been phase flipped. If the 
! data set is large (>5000 images), generating a random model can be quite slow. To speedup, set \texttt{nran} 
! to some smaller number, resulting in \texttt{nran} images selected randomly for reconstruction. <prime3D_init/end> 
!
! The code is distributed with the hope that it will be useful, but _WITHOUT_ _ANY_ _WARRANTY_.
! Redistribution and modification is regulated by the GNU General Public License.
! *Author:* Hans Elmlund 2015
!
program simple_prime3D_init
use simple_cmdline,            only: cmdline
use simple_jiffys,             only: simple_end
use simple_hadamard3D_matcher, only: gen_random_model, prime3D_find_resrange
use simple_params,             only: params
use simple_build,              only: build
use simple_timing
implicit none
type(params)       :: p
type(build)        :: b
type(cmdline)      :: cline
integer, parameter :: MAXIMGS=1000
call timestamp()
call start_Alltimers_cpu()
if( command_argument_count() < 3 )then
    write(*,'(a)',advance='no') 'SIMPLE_PRIME3D_INIT stk=<stack.ext> smpd=<sampling distance(in A)>'
    write(*,'(a)',advance='no') ' msk=<mask radius(in pixels)> [nspace=<nr reference sections{1000}>]'
    write(*,'(a)',advance='no') ' [nran=<size of random sample>] [lp=<low-pass limit(in A)>]'
    write(*,'(a)') ' [ctf=<yes|no|flip|mul{no}>] [deftab=<text file defocus values>] [nthr=<nr OpenMP threads{1}>]'
    write(*,'(a)') '**less commonly used**'
    write(*,'(a)',advance='no') '[pgrp=<cn|dn|t|o|i{c1}>] [npeaks=<nr nonzero orientation weights{1}>]'
    write(*,'(a)',advance='no') ' [ctf=<yes|no|flip|mul{no}>] [kv=<acceleration voltage(in kV){300.}>]'
    write(*,'(a)',advance='no') ' [fraca=<frac amp contrast{0.07}>] [cs=<spherical aberration constant(in mm){2.7}>]'
    write(*,'(a)',advance='no') ' [deftab=<text file with defocus values>] [inner=<inner mask radius(in pixels)>]'
    write(*,'(a)') ' [width=<pixels falloff inner mask{10}>] [xfel=<yes|no{no}>]'
    stop
endif
call cline%parse
call cline%checkvar('stk',  1)
call cline%checkvar('smpd', 2)
call cline%checkvar('msk',  3)
call cline%check
if( .not. cline%defined('eo') )     call cline%set('eo', 'no')
if( .not. cline%defined('nspace') ) call cline%set('nspace', 1000.)
p = params(cline)                    ! parameters generated
if( p%l_xfel )then
    if( cline%defined('msk') .or. cline%defined('mw') .or.&
    cline%defined('nvox') .or. cline%defined('mskfile') )then
        stop 'no mask allowed when processing XFEL patterns; simple_prime2'
    endif
endif
if( p%ctf .ne. 'no')then
    if( .not. cline%defined('deftab') )&
    &stop 'need texfile with defocus/astigmatism values for ctf .ne. no mode exec'
endif
call b%build_general_tbox(p, cline)   ! general objects built
call b%build_hadamard_prime3D_tbox(p) ! prime3D objects built
! determine resolution range
if( cline%defined('lp') ) call prime3D_find_resrange( b, p, p%lp, p%lpstop )
! determine the number of peaks
if( .not. cline%defined('npeaks') ) p%npeaks = min(10,b%e%find_npeaks(p%lp, p%moldiam))
! generate the random model
if( cline%defined('nran') )then
    call gen_random_model( b, p, p%nran )
else
    if( p%nptcls > MAXIMGS )then
         call gen_random_model( b, p, MAXIMGS )
    else
        call gen_random_model( b, p )
    endif
endif
call simple_end('**** SIMPLE_PRIME3D_INIT NORMAL STOP ****')
!*******************************************************************************
!    Environment shutdown
!
!*******************************************************************************
! shutting down the timers
call stop_Alltimers_cpu()
end program simple_prime3D_init
