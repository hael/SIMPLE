!==Program volassemble
!
! <volassemble/begin> is a program that assembles volume(s) when the reconstruction program (\prgname{simple\_recvol}) 
! has been executed in distributed mode.  <comment/begin>  \texttt{lp} is used to low-pass filter the assembled volume 
! to the given resolution. \texttt{find} allows for using the integer Fourier index instead of \texttt{lp}. The \texttt{find} 
! option is used by \prgname{distr\_simple.pl}, when PRIME is executed in initial model production mode. \texttt{inner} 
! applies a soft-edged inner mask. An inner mask is used for icosahedral virus reconstruction, because the DNA or RNA core 
! is often unordered, and if not removed it may negatively impact the alignment. The \texttt{width} parameter controls the 
! fall-off of the edge of the \texttt{inner} mask. \texttt{even} is used to assemble the even reconstruction. \texttt{odd} 
! is used to assemble the odd reconstruction, and \texttt{eo} is used to assemble both the even and the odd reconstruction. 
! Normally, you don't fiddle with these parameters, but they are used internally by \prgname{distr\_simple.pl}. <comment/end>
! <volassemble/end>
!
! The code is distributed with the hope that it will be useful, but _WITHOUT_ _ANY_ _WARRANTY_.
! Redistribution and modification is regulated by the GNU General Public License.
! *Author:* Hans Elmlund 2015
!
program simple_volassemble
use simple_defs           ! singleton
use simple_jiffys,        ! singleton
use simple_cmdline,       only: cmdline
use simple_build,         only: build
use simple_params,        only: params
use simple_reconstructor, only: reconstructor
implicit none
type(params)                  :: p
type(build)                   :: b
type(cmdline)                 :: cline
character(len=:), allocatable :: fbody
integer                       :: part, s, ss, endit, i, cnt, state4name
type(reconstructor)           :: recvol_read
logical                       :: here(2), debug=.false.
if( command_argument_count() < 4 )then
    write(*,'(a)', advance='no') 'SIMPLE_VOLASSEMBLE stk=<ptcls.ext> npart=<nr partitions>'
    write(*,'(a)', advance='no') ' smpd=<sampling distance(in A)> oritab=<algndoc.txt>'
    write(*,'(a)', advance='no') ' [nthr=<nr openMP threads{1}>]'
    write(*,'(a)') ' [even=<yes|no{no}>] [odd=<yes|no{no}>] [eo=<yes|no{no}>] [xfel=<yes|no{no}>]'
    stop
endif
call cline%parse
call cline%checkvar('stk',   1)
call cline%checkvar('npart', 2)
call cline%checkvar('smpd',  3)
if( cline%defined('state') ) call cline%set('nstates', 1.)
if( .not. cline%defined('eo') )then
    call cline%set('eo', 'no')
endif
call cline%check
p = params(cline) ! constants & derived constants produced
call b%build_general_tbox(p,cline) ! general objects built
call b%build_rec_tbox(p)     ! reconstruction toolbox built
! rebuild b%vol according to box size (beacuse it is otherwise boxmatch)
call b%vol%new([p%box,p%box,p%box], p%smpd)
if( cline%defined('find') )then
    p%lp = b%img%get_lp(p%find)
endif
call recvol_read%new([p%boxpd,p%boxpd,p%boxpd], p%smpd)
call recvol_read%alloc_rho(p)
endit = 1
if( p%eo .eq. 'yes' ) endit = 2
cnt = 0
do ss=1,p%nstates
    if( cline%defined('state') )then
        s = 1
    else
        s = ss
    endif
    if( debug ) write(*,*) 'processing state: ', s
    call b%recvol%reset
    do part=1,p%npart
        cnt = cnt+1
        call progress(cnt,p%nstates*p%npart)
        if( cline%defined('state') )then
            state4name = p%state
        else
            state4name = s
        endif
        allocate(fbody, source='recvol'//'_state'//int2str_pad(state4name,2)//'_part'//int2str_pad(part,p%numlen))
        if( debug ) write(*,*) 'processing fbody: ', fbody
        do i=1,endit
            if( cline%defined('even') .or. cline%defined('odd') )then
                if( p%even .eq. 'yes' .and. p%odd .eq. 'no' )then
                    p%vols(s) = fbody//'_even'//p%ext
                    p%masks(s) = 'rho_'//fbody//'_even'//p%ext
                else if( p%odd .eq. 'yes' .and. p%even .eq. 'no' )then
                    p%vols(s) = fbody//'_odd'//p%ext
                    p%masks(s) = 'rho_'//fbody//'_odd'//p%ext
                else if( p%odd .eq. 'yes' .and. p%even .eq. 'yes' )then
                    stop 'ERROR! Cannot have even=yes and odd=yes simultaneously'
                endif
            else
                if( p%eo .eq. 'yes' )then
                    if( i == 1 )then
                        p%vols(s) = fbody//'_odd'//p%ext
                        p%masks(s) = 'rho_'//fbody//'_odd'//p%ext
                    else
                        p%vols(s) = fbody//'_even'//p%ext
                        p%masks(s) = 'rho_'//fbody//'_even'//p%ext
                    endif   
                else
                    p%vols(s)  = fbody//p%ext
                    p%masks(s) = 'rho_'//fbody//p%ext
                endif
            endif
            call assemble(p%vols(s), p%masks(s))
        end do
        deallocate(fbody)
    end do
    if( p%even .eq. 'yes' .and. p%odd .eq. 'no' )then
        call normalize('recvol_state'//int2str_pad(s,2)//'_even'//p%ext)
    else if( p%odd .eq. 'yes' .and. p%even .eq. 'no' )then
        call normalize('recvol_state'//int2str_pad(s,2)//'_odd'//p%ext)
    else if( p%odd .eq. 'yes' .and. p%even .eq. 'yes' )then
        stop 'ERROR! Cannot have even=yes and odd=yes simultaneously'
    else
        call normalize('recvol_state'//int2str_pad(s,2)//p%ext)
    endif
end do

call simple_end('**** SIMPLE_VOLASSEMBLE NORMAL STOP ****')

contains

    subroutine assemble( recnam, kernam )
        character(len=*), intent(in) :: recnam
        character(len=*), intent(in) :: kernam
        inquire(FILE=recnam, EXIST=here(1))
        inquire(FILE=kernam, EXIST=here(2))
        if( all(here) )then     
            call recvol_read%read(recnam)
            if( debug )then
                if( recvol_read%contains_nans() )then
                    write(*,*) 'WARNING! recvol: ', recnam, 'contains NaN:s'
                endif
            endif
            call recvol_read%read_rho(kernam)
            if( debug )then
                if( recvol_read%rho_contains_nans() )then
                    write(*,*) 'WARNING! kernel: ', kernam, 'contains NaN:s'
                endif
            endif
            call b%recvol%sum(recvol_read)
            if( debug )then
                if( b%recvol%contains_nans() )then
                    write(*,*) 'WARRNING! summed image part contains NaN:s'
                endif 
                if( b%recvol%rho_contains_nans() )then
                    write(*,*) 'WARNING! summed kernel part contains NaN:s'
                endif
            endif
        else
            if( .not. here(1) ) write(*,'(A,A,A)') 'WARNING! ', adjustl(trim(recnam)), ' missing'
            if( .not. here(2) ) write(*,'(A,A,A)') 'WARNING! ', adjustl(trim(kernam)), ' missing'
            return
        endif
    end subroutine
    
    subroutine normalize( recnam )
        character(len=*), intent(in) :: recnam
        call b%recvol%sampl_dens_correct
        call b%recvol%bwd_ft
        call b%recvol%clip(b%vol)
        call b%vol%write(recnam, del_if_exists=.true.)
    end subroutine 
    
end program simple_volassemble