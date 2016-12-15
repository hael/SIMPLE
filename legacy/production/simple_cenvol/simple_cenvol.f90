!==Program simple_cenvol
!
! <cenvol/begin> is a program for centering a volume and mapping the shift parameters back to the particle images
! <comment/begin> Often, when class averages are used for 3D processing and the paramaters are mapped back to the
! particles, the reconstructed volume is off-centre. This program is useful for making sure that the mask doesn't
! cut off-centre volumes <comment/end> <cenvol/end>
!
! The code is distributed with the hope that it will be useful, but _WITHOUT_ _ANY_ _WARRANTY_.
! Redistribution and modification is regulated by the GNU General Public License.
! *Author:* Cyril Reboul & Hans Elmlund 2016
!
program simple_cenvol
use simple_defs     ! singleton
use simple_cmdline, only: cmdline
use simple_jiffys,  only: simple_end, int2str
use simple_params,  only: params
use simple_build,   only: build
implicit none
type(params), target :: p
type(build), target  :: b
type(cmdline)        :: cline
integer              :: istate
real, allocatable    :: shvec(:,:)
logical, parameter   :: debug=.false.
if( command_argument_count() < 4 )then
    write(*,'(a)', advance='no') 'SIMPLE_CENVOL vol1=<vol1.ext> [vol2=<vol2.ext> etc.]'
    write(*,'(a)', advance='no') ' smpd=<sampling distance(in A)> oritab=<input alignment doc>'
    write(*,'(a)') ' outfile=<output alignment doc> [amsklp=<low-pass limit for centering mask(in A){50}> '
    stop
endif
call cline%parse
call cline%checkvar('vol1',    1)
call cline%checkvar('smpd',    2)
call cline%checkvar('oritab',  3)
call cline%checkvar('outfile', 4)
if( .not. cline%defined('amsklp') ) call cline%set('amsklp', 50.)
call cline%check ! checks args and prints cmdline to cmdline.dat
call cline%set('prg', 'cenvol')
p = params(cline)  ! constants & derived constants produced
call b%build_general_tbox(p, cline, .true.)
! center volume(s)
allocate(shvec(p%nstates,3))
do istate=1,p%nstates
    call b%vol%read(p%vols(istate))
    shvec(istate,:) = b%vol%center(p%amsklp,p%msk)
    if( debug )then
        call b%vol%shift(-shvec(istate,1), -shvec(istate,2), -shvec(istate,3))
        call b%vol%write('shifted_vol_state'//int2str(istate)//p%ext)
    endif
    ! transfer the 3D shifts to 2D
    call b%a%map3dshift22d(-shvec(istate,:), state=istate)
end do
call b%a%write(p%outfile)
call simple_end('**** SIMPLE_CENVOL NORMAL STOP ****')
end program
