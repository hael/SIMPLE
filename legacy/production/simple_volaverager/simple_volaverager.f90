!==Program simple_volaverager
!
! <simple\_volaverager/begin> is a program for averaging volumes accordin to state label in oritab.
! <simple\_volaverager/end> 
!
! The code is distributed with the hope that it will be useful, but _WITHOUT_ _ANY_ _WARRANTY_.
! Redistribution and modification is regulated by the GNU General Public License.
! *Author:* Hans Elmlund Feb 2016
!
program simple_volaverager
use simple_jiffys,  ! singleton
use simple_cmdline, only: cmdline
use simple_params,  only: params
use simple_build,   only: build
use simple_image,   only: image
use simple_defs
implicit none
type(params)                       :: p
type(build)                        :: b
type(cmdline)                      :: cline
integer, allocatable               :: ptcls(:)
character(len=STDLEN), allocatable :: volnames(:)
type(image)                        :: vol_avg
integer                            :: istate, ivol, nvols, funit_vols, numlen, ifoo
character(len=:), allocatable      :: fname
character(len=1)                   :: fformat
logical                            :: debug=.true.
if( command_argument_count() < 2 )then
    write(*,'(a)') 'SIMPLE_VOLAVERAGER vollist=<list of volumes> oritab=oritab=<state labels>'
    stop
endif
call cline%parse
call cline%checkvar('vollist', 1)
call cline%checkvar('oritab',  2)
call cline%check
p = params(cline) ! parameters generated
! READ THE VOLNAMES
nvols = nlines(p%vollist)
if( debug ) print *, 'number of volumes: ', nvols
allocate(volnames(nvols))
funit_vols = get_fileunit()
open(unit=funit_vols, status='old', file=p%vollist)
do ivol=1,nvols
    read(funit_vols,'(a256)') volnames(ivol)
    if( debug ) print *, 'read volname: ', volnames(ivol)
end do
close(funit_vols)
! FIND LOGICAL DIMENSION
call find_ldim_nptcls(volnames(1), p%ldim, ifoo)
p%box  = p%ldim(1)
! BUILD GENERAL TOOLBOX
call b%build_general_tbox(p, cline) ! general objects built
! FIGURE OUT THE FILE EXTENSION
fformat = fname2format(volnames(1))
select case(fformat)
    case('M')
        p%ext = '.mrc'
    case('S')
        p%ext = '.spi'
    case('D')
        p%ext = '.mrc'
    case('B')
        p%ext = '.mrc'
    case DEFAULT
        stop 'This file format is not supported by SIMPLE; simple_volaverager'
end select
if( debug ) print *, 'file extension: ', p%ext
! AVERAGE THE STATES
call vol_avg%copy(b%vol)
p%nstates = b%a%get_nstates()
if( debug ) print *, 'number of states: ', p%nstates
numlen = len(int2str(p%nstates))
do istate=1,p%nstates
    if( debug ) print *, 'processing state: ', istate
    ptcls = nint(b%a%get_ptcls_in_state(istate))
    vol_avg = 0.
    do ivol=1,size(ptcls)
        call b%vol%read(volnames(ptcls(ivol)))
        if( debug ) print *, 'read volume: ', volnames(ptcls(ivol))
        call vol_avg%add(b%vol)
    end do
    call vol_avg%div(real(size(ptcls)))
    allocate(fname, source='sumvol_state'//int2str_pad(istate, numlen)//p%ext)
    if( debug ) print *, 'trying to write volume to file: ', fname
    call vol_avg%write(fname)
    deallocate(ptcls,fname)
end do
call simple_end('**** SIMPLE_VOLAVERAGER NORMAL STOP ****')
end program simple_volaverager
