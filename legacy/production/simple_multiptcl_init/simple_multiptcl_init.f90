program simple_multiptcl_init
use simple_rec_master, only: exec_rec, exec_eorec
use simple_build,      only: build
use simple_params,     only: params
use simple_jiffys,     only: simple_end
use simple_cmdline,    only: cmdline
use simple_defs        ! singleton
implicit none
type(params)      :: p
type(build)       :: b
type(cmdline)     :: cline
real, allocatable :: optlp
if( command_argument_count() < 3 )then
    write(*,'(a)', advance='no') 'SIMPLE_MULTIPTCL_INIT stk=<ptcls.ext> smpd=<sampling distance(in A)>'
    write(*,'(a)', advance='no') ' oritab=<PRIME3D doc> [msk=<mask radius(in pixels)>] [nstates=<nr states>]'
    write(*,'(a)', advance='no') ' [lp=<low-pass limit{20}>] [eo=<yes|no{yes}>] [frac=<fraction of ptcls to include{1.}>]'
    write(*,'(a)', advance='no') ' [nthr=<nr of openMP threads{1}>] [pgrp=<cn|dn|t|o|i{c1}>] [norec=<yes|no{no}>]'
    write(*,'(a)') ' [state2split=<state group 2 split>]'
    write(*,'(a)') '**less commonly used**'
    write(*,'(a)', advance='no') '[mul=<shift multiplication factor{1}>] [ctf=<yes|no|flip|mul{no}>] [kv=<acceleration'
    write(*,'(a)', advance='no') ' voltage(in kV){300.}>] [fraca=<frac amp contrast{0.07}>] [cs=<spherical'
    write(*,'(a)', advance='no') ' aberration constant(in mm){2.7}>] [deftab=<text file defocus values>]'
    write(*,'(a)', advance='no') ' [errify=<yes|no{no}>] [inner=<inner mask radius(in pixels)>]'
    write(*,'(a)') ' [width=<pixels falloff inner mask{10}>] [zero=<yes|no>{no}]'
    stop
endif
call cline%parse
call cline%checkvar('stk',    1)
call cline%checkvar('smpd',   2)
call cline%checkvar('oritab', 3)
call cline%set('trs', 3.) ! to assure that shifts are being used
call cline%check
call cline%set('prg', 'multiptcl_init')
p = params(cline)         ! constants & derived constants produced
call b%build_general_tbox(p, cline)
if( cline%defined('state2split') )then
    call b%a%split_state(p%state2split)
    p%nstates = p%nstates+1
else
    call b%a%rnd_states(p%nstates)
endif
if( p%errify.eq.'yes' )call errify_oris
if( p%norec .eq. 'no' )then
    if( cline%defined('lp') )then
        call b%build_rec_tbox(p)
        call exec_rec(b, p, cline, 'startvol')
    else
        call b%build_eo_rec_tbox(p)
        call exec_eorec(b, p, cline, 'startvol')
    endif
endif
if( p%zero.eq.'yes' )call b%a%set_all('corr',0.)
call b%a%write('multiptcl_startdoc.txt')
call simple_end('**** SIMPLE_MULTIPTCL_INIT NORMAL STOP ****')    

contains
    
    subroutine errify_oris()
        real :: sherr,angerr
        if( cline%defined('trs') )then
            sherr = p%trs
        else
            sherr = 3.
        endif
        angerr = 15.
        write(*,'(A,F6.2,A)')'>>> IN-PLANE SHIFT   ERROR INTRODUCED: ',sherr,' PIXELS'
        write(*,'(A,F6.2,A)')'>>> IN-PLANE ANGULAR ERROR INTRODUCED: ',angerr,' DEGREES'
        call b%a%introd_alig_err( angerr,sherr )
    end subroutine errify_oris
    
end program simple_multiptcl_init
