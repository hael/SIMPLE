program simple_test_install
include 'simple_lib.f08'

use simple_testfuns          ! use all in there
use simple_image,            only: image

implicit none
type( image )         :: cube, img
real                  :: smpd
integer               :: box, nspace, msk
character(len=8)      :: datestr
character(len=STDLEN) :: folder, cmd
character(len=300)    :: command
call seed_rnd
call date_and_time(date=datestr)
folder = trim('./SIMPLE_TEST_INSTALL_'//datestr)
call simple_mkdir(trim( folder ))
call simple_chdir( folder )
! dummy data
box    = 96
smpd   = 2.
nspace = 64
msk    = nint(real(box)/3.)
! volume
call img%new( [box,box,box], smpd )
call img%square( nint(real(box)/12.) )
call cube%new( [box,box,box], smpd )
call cube%square( nint(real(box)/16.) )
call cube%shift([16.,16.,16.])
call img%add( cube )
call cube%new( [box,box,box], smpd )
call cube%square( nint(real(box)/10.) )
call cube%shift([4.,-16.,0.])
call img%add( cube )
call cube%kill
call img%write( 'cubes.mrc' )
call img%kill
write(logfhandle,*)'>>> WROTE TEST VOLUME cubes.mrc'
! test units
command = 'simple_test_units'
call exec_cmdline( trim(command) )
! test search
command = 'simple_test_srch vol1=cubes.mrc msk='//int2str(msk)//&
    & ' smpd='//real2str(smpd)//' verbose=no'
call exec_cmdline( trim(command) )
! end
call simple_chdir(PATH_PARENT)
call simple_end('**** SIMPLE_TEST_INSTALL NORMAL STOP ****')
end program simple_test_install
