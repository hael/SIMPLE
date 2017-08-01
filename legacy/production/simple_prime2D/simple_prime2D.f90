!==Program simple_prime2D
!
! <prime2D/begin> is a reference-free 2D alignment/clustering algorithm adopted from the prime3D probabilistic \textit{ab initio}
! 3D reconstruction algorithm. It is assumed that the images are phase-flipped (phase flipping can 
! be done with \prgname{simple\_stackops}). Do \textit{not} search the origin shifts initially, when the cluster centers are of low 
! quality. If your images are far off centre, use \prgname{simple\_stackops} with option \texttt{shalgn=yes} instead to shiftalign the images 
! beforehand (the algorithm implemented is the same as EMAN's \texttt{cenalignint} program).
! Use \prgname{distr\_simple.pl} for distributed PRIME2D execution. <prime2D/end>
!
! The code is distributed with the hope that it will be useful, but _WITHOUT_ _ANY_ _WARRANTY_.
! Redistribution and modification is regulated by the GNU General Public License.
! *Authors:* Hans Elmlund X-mas 2015
!
program simple_prime2D
use simple_defs                ! singleton
use simple_jiffys              ! singleton
use simple_cmdline,            only: cmdline
use simple_hadamard2D_matcher, only: prime2D_exec
use simple_params,             only: params
use simple_build,              only: build
use simple_timing
use simple_cuda_defs
use simple_cuda
implicit none
type(params)  :: p
type(build)   :: b
type(cmdline) :: cline
integer       :: i, startit, ncls_from_refs, lfoo(3)
logical       :: converged=.false.
integer       :: err ! CUDA err variable for the return function calls
call timestamp()
! call start_Alltimers_cpu()
! starting the cuda environment
call simple_cuda_init(err)
if( err .ne. RC_SUCCESS ) write(*,*) 'cublas init failed'
if( command_argument_count() < 3 )then
    write(*,'(a)',advance='no') 'SIMPLE_PRIME2D stk=<stack.ext> smpd=<sampling distance(in A)>'
    write(*,'(a)',advance='no') ' msk=<mask radius(in pixels)> ncls=<nr of clusters>'
    write(*,'(a)',advance='no') ' refs=<initial_references.ext> [oritab=<previous clustering doc>]'
    write(*,'(a)',advance='no') ' [lp=<low-pass limit(in A){20}>] [trs=<origin shift(in pixels){0}>]'
    write(*,'(a)',advance='no') ' [use_gpu=<yes|no{no}>] [ctf=<yes|no|flip|mul{no}>] [kv=<acceleration '
    write(*,'(a)',advance='no') ' voltage(in kV){300.}>] [cs=<spherical aberration constant(in mm){2.7}>] '
    write(*,'(a)',advance='no') ' [fraca=<frac amp contrast{0.07}>] [deftab=<defocus info file>]'
    write(*,'(a)',advance='no') ' [nthr=<nr of OpenMP threads{1}>] [startit=<start iteration>]'
    write(*,'(a)',advance='no') ' [hp=<high-pass limit(in A)>] [srch_inpl=<yes|no{yes}>] [maxits=<max'
    write(*,'(a)',advance='no') ' iterations{500}>] [automsk=<yes|no{no}>] [amsklp=<automask low-pass limit(in A){25}>]'
    write(*,'(a)') ' [inner=<inner mask radius(in pixels)>] [width=<pixels falloff inner mask(in pixels){10}>]'    
    stop
endif
call cline%parse
call cline%checkvar('stk',    1)
call cline%checkvar('smpd',   2)
call cline%checkvar('msk',    3)
call cline%checkvar('ncls',   4)
call cline%checkvar('refs',   5)
if( .not. cline%defined('lp') )then
    call cline%set('lp', 20.)
endif
if( .not. cline%defined('amsklp') )then
    call cline%set('amsklp', 25.)
endif
if( .not. cline%defined('edge') )then
    call cline%set('edge', 20.)
endif
if( .not. cline%defined('eo') )then
    call cline%set('eo', 'no')
endif
call cline%check
p = params(cline)                     ! parameters generated
p%boxmatch = p%box                    !!!!!!!!!!!!!!!!!! 4 NOW
call b%build_general_tbox(p, cline)   ! general objects built
call b%build_hadamard_prime2D_tbox(p) ! 2D Hadamard matcher built
if( p%srch_inpl .eq. 'no' )then
    if( .not. cline%defined('oritab') )then
        stop 'need oritab for this mode (srch_inpl=no) of execution!'
    endif
endif
if( cline%defined('refs') )then
    call find_ldim_nptcls(p%refs, lfoo, ncls_from_refs)
    ! consistency check
    if( p%ncls /=  ncls_from_refs ) stop 'nrefs /= inputted ncls'
endif
! execute
if( cline%defined('part') )then
    if( .not. cline%defined('outfile') ) stop 'need unique output file for parallel jobs'
    call prime2D_exec(b, p, cline, 0, converged) ! partition or not, depending on 'part'       
else
    startit = 1
    if( cline%defined('startit') ) startit = p%startit
    do i=startit,p%maxits
        call prime2D_exec(b, p, cline, i, converged)
        if(converged) exit
    end do
endif
call simple_end('**** SIMPLE_PRIME2D NORMAL STOP ****')
!******************************************************
!    Environment shutdown
!******************************************************
! shutting down CUDA
call simple_cuda_shutdown()
! shutting down timers
! call stop_Alltimers_cpu()
end program simple_prime2D
