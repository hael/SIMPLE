program simple_cont3D
use simple_defs            ! singleton
use simple_cmdline,        only: cmdline
use simple_jiffys,         only: simple_end
use simple_cont3D_matcher, only: cont3D_exec
use simple_params,         only: params
use simple_build,          only: build
implicit none
type(params)  :: p
type(build)   :: b
type(cmdline) :: cline
integer       :: i, startit
logical       :: converged=.false.
real          :: lpstart, lpstop
if( command_argument_count() < 3 )then
    write(*,'(a)',advance='no') 'SIMPLE_CONT3D stk=<stack.ext> vol1=<invol.ext> [vol2=<refvol_2.ext> etc.]'
    write(*,'(a)',advance='no') ' smpd=<sampling distance(in A)> msk=<mask radius(in pixels)>'
    write(*,'(a)',advance='no') ' oritab=<previous alignment doc> trs=<origin shift(in pixels){0}>'
    write(*,'(a)',advance='no') ' [frac=<fraction of ptcls to include{1}>] [automsk=<yes|no{no}>]'
    write(*,'(a)',advance='no') ' [mw=<molecular weight(in kD)>] [nthr=<nr of OpenMP threads{1}>]'
    write(*,'(a)',advance='no') ' [startit=<start iteration>] [lpstop=<stay at this low-pass limit (in A)>]'
    write(*,'(a)',advance='no') ' [deftab=<text file defocus values>]'
    write(*,'(a)',advance='no') ' [amsklp=<automask low-pass limit(in A){20}>] [pgrp=<cn|dn|t|o|i{c1}>]'
    write(*,'(a)',advance='no') ' [ctf=<yes|no|flip|mul{no}>] [kv=<acceleration voltage(in kV){300.}>]'
    write(*,'(a)',advance='no') ' [cs=<spherical aberration constant(in mm){2.7}>] [fraca=<frac amp'
    write(*,'(a)') ' contrast{0.07}>] [hp=<high-passlimit(in A)>] [xfel=<yes|no{no}>]'
    write(*,'(a)') '**less commonly used**'
    write(*,'(a)',advance='no') ' [maxits=<max iterations{500}>]'
    write(*,'(a)',advance='no') ' [edge=<edge size for softening molecular envelope(in pixels){15}>]'
    write(*,'(a)',advance='no') ' [inner=<inner mask radius(in pixels)>] [width=<pixels falloff inner mask{10}>]'
    stop
endif
call cline%parse
call cline%checkvar('stk',        1)
call cline%checkvar('vol1',       2)
call cline%checkvar('smpd',       3)
call cline%checkvar('msk',        4)
call cline%checkvar('oritab',     5)
call cline%checkvar('trs',        6)
call cline%set('eo',     'yes')
call cline%set('refine', 'yes')
call cline%check
call cline%set('prg', 'cont3D')
p = params(cline)
if( p%xfel .eq. 'yes' )then
    if( cline%defined('msk') .or. cline%defined('mw') .or.&
    cline%defined('nvox') .or. cline%defined('automsk') )then
        stop 'no mask allowed when processing XFEL patterns; simple_prime3D'
    endif
endif
call b%build_general_tbox(p, cline)
call b%build_cont3D_tbox(p)
if( cline%defined('part') )then
    if( .not. cline%defined('outfile') ) stop 'need unique output file for parallel jobs'
    call cont3D_exec(b, p, cline, 0, converged) ! partition or not, depending on 'part'
else
    startit = 1
    if( cline%defined('startit') ) startit = p%startit
    do i=startit,p%maxits
        call cont3D_exec(b, p, cline, i, converged)
        if(converged) exit
    end do
endif
end program simple_cont3D