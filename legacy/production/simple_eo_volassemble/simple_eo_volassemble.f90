!==Program simple_eo_volassemble
!
! <eo_volassemble/begin> is a program that assembles volume(s) when the reconstruction program (\prgname{simple\_eo\_recvol}) 
! has been executed in distributed mode using \prgname{distr\_simple.pl}. <comment/begin> \texttt{inner} applies a soft-edged 
! inner mask. An inner mask is used for icosahedral virus reconstruction, because the DNA or RNA core is often unordered and 
! if not removed it may negatively impact the alignment. The \texttt{width} parameter controls the fall-off of the edge of the 
! \texttt{inner} mask. <comment/end> <eo_volassemble/end>
!
! The code is distributed with the hope that it will be useful, but _WITHOUT_ _ANY_ _WARRANTY_.
! Redistribution and modification is regulated by the GNU General Public License.
! *Author:* Hans Elmlund 2011-08-16.
!
program simple_eo_volassemble
use simple_defs              ! singleton
use simple_jiffys,           ! singleton
use simple_cmdline,          only: cmdline
use simple_build,            only: build
use simple_params,           only: params
use simple_eo_reconstructor, only: eo_reconstructor
implicit none
type(params)                  :: p
type(build)                   :: b
type(cmdline)                 :: cline
integer                       :: part, s, alloc_stat, cnt, n, ss, state4name
type(eo_reconstructor)        :: eorecvol_read
character(len=:), allocatable :: fname
real, allocatable             :: res05s(:), res0143s(:)
real                          :: res
logical                       :: debug=.false.
if( command_argument_count() < 4 )then
    write(*,'(a)', advance='no') 'SIMPLE_EO_VOLASSEMBLE stk=<ptcls.ext> npart=<nr partitions>'
    write(*,'(a)', advance='no') ' msk=<mask radius(in pixels)> smpd=<sampling distance(in A)> oritab=<algndoc.txt> '
    write(*,'(a)', advance='no') ' [nstates=<number of states>] [nthr=<nr openMP threads{1}>]'
    write(*,'(a)', advance='no') ' [lpstop=<stay at this low-pass limit(in A)>]'
    write(*,'(a)') ' [state=<state to reconstruct{all}>] [ctf=<yes|no|flip|mul{no}>]'
    stop
endif
call cline%parse
call cline%checkvar('stk',   1)
call cline%checkvar('npart', 2)
call cline%checkvar('msk',   3)
call cline%checkvar('smpd',  4)
if( cline%defined('state') ) call cline%set('nstates', 1.)
call cline%check
call cline%set('prg', 'eo_volassemble')
p = params(cline) ! constants & derived constants produced
call b%build_general_tbox(p, cline)  ! general objects built
if( cline%defined('nstates') )then
    if( p%nstates /= b%a%get_nstates() ) stop 'Inconsistent number of states between command-line and oritab'
endif
call b%build_eo_rec_tbox(p)   ! reconstruction toolbox built
allocate(res05s(p%nstates), res0143s(p%nstates), stat=alloc_stat)
call alloc_err("In: simple_eo_volassemble", alloc_stat)
! rebuild b%vol according to box size (beacuse it is otherwise boxmatch)
call b%vol%new([p%box,p%box,p%box], p%smpd)
call eorecvol_read%new(p)
n = p%nstates*p%npart
cnt = 0
do ss=1,p%nstates
    if( cline%defined('state') )then
        s = 1
    else
        s = ss
    endif
    if( debug ) write(*,*) 'processing state: ', s
    call b%eorecvol%reset_all
    do part=1,p%npart
        cnt = cnt+1
        call progress(cnt,n)
        if( cline%defined('state') )then
            state4name = p%state
        else
            state4name = s
        endif
        allocate(fname, source='recvol'//'_state'//int2str_pad(state4name,2)//'_part'//int2str_pad(part,p%numlen))
        if( debug ) write(*,*) 'processing file: ', fname
        call assemble(fname)
        deallocate(fname)
    end do
    call normalize('recvol_state'//int2str_pad(state4name,2))
end do
! set the resolution limit according to the worst resolved model
res  = maxval(res0143s)
p%lp = max( p%lpstop,res )
write(*,'(a,1x,F6.2)') '>>> LOW-PASS LIMIT:', p%lp
write(0,'(a)') "GENERATED VOLUMES: recvol*.ext"
call simple_end('**** SIMPLE_EO_VOLASSEMBLE NORMAL STOP ****')

contains

    subroutine assemble( fbody )
        character(len=*), intent(in) :: fbody
        call eorecvol_read%read_eos(fbody)
        ! sum the Fourier coefficients
        call b%eorecvol%sum(eorecvol_read)
    end subroutine
    
    subroutine normalize( recnam )
        use simple_image, only: image
        character(len=*), intent(in)  :: recnam
        integer :: ldim(3)
        real    :: smpd
        call b%eorecvol%sum_eos
        call b%eorecvol%sampl_dens_correct_eos(s)
        call b%eorecvol%get_res(res05s(s), res0143s(s))
        call b%eorecvol%sampl_dens_correct_sum(b%vol)
        call b%eorecvol%write_eos(recnam)
        call b%vol%write(recnam//p%ext, del_if_exists=.true.)
    end subroutine
    
end program simple_eo_volassemble