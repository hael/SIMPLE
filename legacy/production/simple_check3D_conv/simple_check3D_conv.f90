!==Program simple_check3D_conv
!
! <check3D_conv/begin> is a program for checking if a PRIME3D run has converged. The statistics outputted
! include (1) angle of feasible region, which is proportional to the angular resolution of the set of discrete
! projection directions being searched. (2) The average angular distance between orientations in the present and
! previous iteration. In the early iterations, the distance is large because a diverse set of orientations is explored.
! If convergence to a local optimum is achieved, the distance decreases. (3) The percentage of search space 
! scanned, i.e. how many reference images are evaluated on average. (4) The average correlation between the images 
! and their corresponding best matching reference sections. (5) The average standard deviation of the Euler angles. 
! Convergence is achieved if the angular distance between the orientations in successive iterations falls significantly 
! below the angular resolution of the search space and more than 99\% of the reference sections need to be matched on
! average. <check3D_conv/end>
!
! The code is distributed with the hope that it will be useful, but _WITHOUT_ _ANY_ _WARRANTY_.
! Redistribution and modification is regulated by the GNU General Public License.
! *Author:* Hans Elmlund 2011
!
program simple_check3D_conv
use simple_cmdline, only: cmdline
use simple_build,   only: build
use simple_params,  only: params
use simple_math,    only: rad2deg, get_lplim
use simple_jiffys,  only: file2rarr, simple_end,int2str_pad
implicit none
type(params)      :: p
type(build)       :: b
type(cmdline)     :: cline
real, allocatable :: maplp(:)
integer           :: istate, loc(1)
logical           :: here, limset, converged
if( command_argument_count() < 3 )then
    write(*,'(a)',advance='no') 'SIMPLE_CHECK3D_CONV smpd=<sampling distance(in A)> box=<image size(in pixels)>'
    write(*,'(a)',advance='no') ' oritab=<alignment doc> nptcls=<nr particle images> [nstates=<number of states>]'
    write(*,'(a)',advance='no') ' [lp=<low-pass limit{20}>] [eo=<yes|no{no}>] [pgrp=<cn|dn|t|o|i{c1}>]'
    write(*,'(a)',advance='no') ' [nspace=<nr reference sections{1000}>] [find=<Fourier index>]'
    write(*,'(a)') ' [refine=<no|neigh|shc|shcneigh|qcont|qcontneigh{no}>]'
    stop
endif
call cline%parse
if( .not. cline%defined('nspace') ) call cline%set('nspace', 1000.)
if( .not. cline%defined('lp') )     call cline%set('lp',       20.)
call cline%checkvar('smpd',   1)
call cline%checkvar('oritab', 2)
call cline%checkvar('nptcls', 3)
call cline%check
call cline%set('prg', 'check3D_conv')
p = params(cline)                   ! parameters generated
call b%build_general_tbox(p, cline) ! general objects built

! nstates consistency check
if( cline%defined('nstates') )then
    if( p%nstates /= b%a%get_nstates() ) stop 'Inconsistent number of states between command-line and oritab'
endif

limset = .false.
if( p%eo .eq. 'yes' )then
    allocate( maplp(p%nstates) )
    do istate=1,p%nstates
        p%fsc = 'fsc_state'//int2str_pad(istate,2)//'.bin'
        inquire(file=p%fsc, exist=here)
        if( here )then
            b%fsc(istate,:) = file2rarr(p%fsc)
            maplp(istate)   = max(b%img%get_lp(get_lplim(b%fsc(istate,:))),2.*p%smpd)
        else
            write(*,*) 'Tried to open the fsc file: ', trim(p%fsc)
            stop 'but it does not exist!'
        endif
    enddo
    loc     = maxloc( maplp )
    p%state = loc(1)            ! state with worst low-pass
    p%lp    = maplp( p%state )  ! worst lp
    p%fsc   =  'fsc_state'//int2str_pad(p%state,2)//'.bin'
    deallocate(maplp)
    limset = .true.
endif

! Let find override the command line input lp (if given)
if( .not. limset .and. cline%defined('find') )then
    p%lp = b%img%get_lp(p%find)
    limset = .true.
endif

! Method for setting lp with lowest priority is lp on the command line
if( cline%defined('lp') ) limset = .true.

! If we arrived here and the limit wasn't set: fall over
if( limset )then
    ! we are happy
else
    ! we fall over
    stop 'No method available to set low-pass limit! ABORTING...'
endif

! calculate angular threshold
p%thres = rad2deg(atan(p%lp/(p%moldiam/2.)))

! check convergence
converged = b%conv%check_conv3D()
call simple_end('**** SIMPLE_CHECK3D_CONV STOP ****')    
end program 