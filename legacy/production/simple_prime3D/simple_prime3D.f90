!==Program simple_prime3D
!
! <prime3D/begin> is an \textit{ab inito} reconstruction/refinement program based on probabilistic projection matching. 
! PRIME is shorthand for PRobabilistic Initial 3D Model Generation for Single-Particle Cryo-Electron Microscopy. You should 
! use phase-flipped images for initial model production with PRIME3D (phase flipping can be done with \prgname{simple\_stackops}). 
! Do not search the origin shifts initially, when the model is of very low quality. If your 
! images are far off centre, use \prgname{simple\_stackops} with option \texttt{shalgn=yes} instead to shiftalign the images 
! beforehand (the algorithm implemented is the same as EMAN's \texttt{cenalignint} program). We recommend running the first round 
! of PRIME with the default dynamic resolution stepping \texttt{dynlp=yes}. The \texttt{dynlp} option implements a heuristic resolution 
! weighting/update scheme. The initial low-pass limit is set so that each image receives ten nonzero orientation weights. When 
! quasi-convergence has been reached, the limit is updated one Fourier index at the time until PRIME reaches the condition where 
! six nonzero orientation weights are assigned to each image. FSC-based filtering is unfortunately not possible to do in the 
! ab initio reconstruction step, because when the orientations are mostly random, the FSC overestimates the resolution.
! Once the initial model has converged, we recommend start searching the shifts (by setting \texttt{trs} to some nonzero value), 
! applying the FSC for resolution-weighting (by setting \texttt{eo=yes}). You should NOT give \texttt{ctf=flip} on the command line
! unless the model has converged. Giving \texttt{ctf=flip} on the command lines signal to PRIME that tou have obtain a reconstruction
! of decent resolution and you want to take it further by applying Wiener restoration by resolution-weighting the reocnstructed volume
! more accurately. In order to be able to use Wiener restoration you also need to input CTF parameters, for example via 
! \texttt{deftab=defocus\_values.txt}. Remember that the defocus values should be given in microns and the astigmatism angle in 
! degrees (one row of the file \texttt{defocus\_values.txt} may look like: \texttt{dfx=3.5} \texttt{dfy=3.3} \texttt{angast=20.0}). 
! <comment/begin>Note that we do not assume any point-group symmetry in the initial runs. However, the \prgname{simple\_symsrch} 
! program can be used to align the reconstruction to its symmetry axis so that future searches can be restricted to the asymmetric unit in 
! refinement. Less commonly used and less obvious input parameters are \texttt{nspace}, which  controls the number of reference projections, 
! \texttt{amsklp}, which controls the low-pass limit used in the automask routine, \texttt{maxits}, which controls the maximum number of 
! iterations executed, \texttt{pgrp}, which controls the point-group symmetry, assuming that the starting volume is aligned to its principal 
! symmetry axis, \texttt{edge}, which controls the size of the softening edge in the automask routine.<comment/end> <prime3D/end>
!
! The code is distributed with the hope that it will be useful, but _WITHOUT_ _ANY_ _WARRANTY_.
! Redistribution and modification is regulated by the GNU General Public License.
! *Authors:* Frederic Bonnet, Cyril Reboul & Hans Elmlund 2015
!
program simple_prime3D
use simple_defs                ! singleton
use simple_jiffys,             ! singleton
use simple_cmdline,            only: cmdline
use simple_hadamard3D_matcher, only: prime3D_exec, prime3D_find_resrange
use simple_params,             only: params
use simple_build,              only: build
use simple_timing
use simple_cuda_defs
use simple_cuda
use simple_file_utils
use simple_file_defs
!temporary
!use simple_err_defs
use simple_file_highlev
!use simple_eglossary
!use simple_error_handling
!use simple_dynamic_memory
!use simple_systemQuery_cpu
!use simple_deviceQuery_gpu
implicit none
type(params)  :: p
type(build)   :: b
type(cmdline) :: cline
integer       :: i, startit
logical       :: update_res=.false., converged=.false.
real          :: lpstart, lpstop
! benchmark control logicals data structure
type(t_bench):: s_bench
! filename string
character(len=10) :: char_out
character(len=80) :: tmr_name
! file handlers
integer      :: unt
integer      :: length
! CUDA err variable for the return function calls
integer      :: err
!function calls
integer      :: get_length_of_string_c
integer      :: convert_int2char_pos_c
integer      :: convert_int2char_indexed_c
integer      :: strlen
call timestamp()
call start_Alltimers_cpu()
call simple_file_lib_initialise()
! starting the cuda environment
call simple_cuda_init(err)
if (err .ne. RC_SUCCESS ) write(*,*) 'cublas init failed'
if( command_argument_count() < 3 )then
    write(*,'(a)',advance='no') 'SIMPLE_PRIME3D stk=<stack.ext> vol1=<invol.ext> [vol2=<refvol_2.ext> etc.]'
    write(*,'(a)',advance='no') ' smpd=<sampling distance(in A)> msk=<mask radius(in pixels)>'
    write(*,'(a)',advance='no') ' [oritab=<previous alignment doc>] [trs=<origin shift(in pixels){0}>]'
    write(*,'(a)',advance='no') ' [lp=<low-pass limit{20}>] [dynlp=<yes|no{yes}>]'
    write(*,'(a)',advance='no') ' [nstates=<nstates to reconstruct>] [frac=<fraction of ptcls to include{1}>]'
    write(*,'(a)',advance='no') ' [automsk=<yes|no{no}>] [mw=<molecular weight(in kD)>] [amsklp=<automask low-pass'
    write(*,'(a)',advance='no') ' limit(in A){20}>] [nthr=<nr of OpenMP threads{1}>]'
    write(*,'(a)',advance='no') ' [use_gpu=<yes|no{no}>] [fix_gpu=<yes|no{no}>]'
    write(*,'(a)',advance='no') ' [set_gpu=<(1-MAX_N_GPU){0}>]'
    write(*,'(a)',advance='no') ' [startit=<start iteration>] [refine=<no|shc|neigh|shcneigh|qcont|qcontneigh{no}>]'
    write(*,'(a)',advance='no') ' [lpstop=<stay at this low-pass limit (in A)>] [deftab=<text file defocus values>]'
    write(*,'(a)',advance='no') ' [nspace=<nr reference sections{1000}>] [eo=<yes|no{no}>]' 
    write(*,'(a)',advance='no') '  [pgrp=<cn|dn|t|o|i{c1}>]'
    write(*,'(a)',advance='no') ' [ctf=<yes|no|flip|mul{no}>] [kv=<acceleration voltage(in kV){300.}>] [cs=<spherical'
    write(*,'(a)',advance='no') ' aberration constant(in mm){2.7}>] [fraca=<frac amp contrast{0.07}>] [hp=<high-pass'
    write(*,'(a)') ' limit(in A)>] [xfel=<yes|no{no}>]'
    write(*,'(a)') '**less commonly used**'
    write(*,'(a)',advance='no') ' [maxits=<max iterations{500}>] [shbarrier=<yes|no{yes}>]'
    write(*,'(a)',advance='no') ' [noise=<yes|no{no}>] [npeaks=<number of nonzero orientation weights>]'
    write(*,'(a)',advance='no') ' [edge=<edge size for softening molecular envelope(in pixels){15}>]'
    write(*,'(a)',advance='no') ' [binwidth=<binary layers grown for molecular envelope(in pixels){1}>]'
    write(*,'(a)') ' [inner=<inner mask radius(in pixels)>] [width=<pixels falloff inner mask{10}>]'
    stop
endif
call cline%parse
call cline%checkvar('stk',  1)
call cline%checkvar('vol1', 2)
call cline%checkvar('smpd', 3)
call cline%checkvar('msk',  4)
if( .not. cline%defined('nspace') )then
    call cline%set('nspace', 1000.)
endif
if( cline%defined('lp') .or. cline%defined('find')) call cline%set('dynlp', 'no')
if( .not. cline%defined('refine') )then
    call cline%set('refine', 'no')
endif
call cline%check
p = params(cline)                 ! parameters generated
if( p%l_xfel )then
    if( cline%defined('msk') .or. cline%defined('mw') .or.&
    cline%defined('nvox') .or. cline%defined('automsk') )then
        stop 'no mask allowed when processing XFEL patterns; simple_prime3D'
    endif
endif
if( str_has_substr(p%refine,'neigh') .or. str_has_substr(p%refine,'qcont'))then
    if( .not. cline%defined('oritab') )then
        stop 'need oritab input for execution of prime3D with refine mode'
    endif
endif
call b%build_general_tbox(p, cline)         ! general objects built
if( .not. cline%defined('eo') ) p%eo = 'no' ! default
if( p%eo .eq. 'yes' ) p%dynlp = 'no'    
if( cline%defined('lp') .or. cline%defined('find')&
.or. p%eo .eq. 'yes' .or. p%dynlp .eq. 'yes' )then
    ! alles ok!
else
   stop 'need a starting low-pass limit (set lp or find)!'
endif
! build prime objects
call b%build_hadamard_prime3D_tbox(p)
if( cline%defined('part') )then
   if( .not. cline%defined('outfile') ) stop 'need unique output file for parallel jobs'
   if( cline%defined('find') ) p%lp = b%img%get_lp(p%find)
   call prime3D_exec(b, p, cline, 0, update_res, converged) ! partition or not, depending on 'part'
else
   if( p%dynlp .eq. 'yes' )then
      call prime3D_find_resrange( b, p, lpstart, lpstop ) ! determine resolution range
      if( cline%defined('lpstart') )then
         p%lp = p%lpstart
      else
         p%lp = lpstart
      endif
      if( cline%defined('lpstop') )then
         ! defined aleady
      else
         p%lpstop = lpstop
      endif
   endif
   p%find = int((real(p%box-1)*p%smpd)/p%lp)
   startit = 1
   if( cline%defined('startit') ) startit = p%startit
   call start_timer_cpu('z_prime3D_exec')
   length = get_length_of_string_c(p%maxits)
   unt = 1
   do i=startit,p%maxits
      unt = 100 + i
      !err = convert_int2char_pos_c(char_out,i)
      err = convert_int2char_indexed_c(char_out,i,startit,p%maxits)
      tmr_name = 'z_prime3D_iter'
      tmr_name = tmr_name(1:strlen(tmr_name))//char_out(1:strlen(char_out))
      tmr_name = tmr_name(1:strlen(tmr_name))//".asc"
      !write(*,*) "tmr_name: ",tmr_name
      if (b%s_bench%bench_write_i==1 .or. ibench_write .eqv. .true. ) call file_open(tmr_name,unt,'unknown','asis','readwrite')
      if (b%s_bench%bench_i==1 .or. ibench .eqv. .true. ) call start_timer_cpu(tmr_name)
      call prime3D_exec(b, p, cline, i, update_res, converged)
      if (b%s_bench%bench_i==1 .or. ibench .eqv. .true. ) call stop_timer_cpu(tmr_name)
      if( update_res )then
         p%find = p%find+p%fstep
         p%lp = max(p%lpstop,b%img%get_lp(p%find))
      endif
      if(converged) exit
   end do
   call stop_timer_cpu('z_prime3D_exec')
endif
call simple_end('**** SIMPLE_PRIME3D NORMAL STOP ****')
!******************************************************
!    Environment shutdown
!
!******************************************************
!shutting down the environment
call simple_cuda_shutdown()
!shutting down the timers
call stop_Alltimers_cpu()
end program simple_prime3D
