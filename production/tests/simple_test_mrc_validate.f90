program simple_test_mrc_validate
use simple_core_module_api
use simple_cmdline, only: cmdline
use simple_image,   only: image
implicit none
#include "simple_local_flags.inc"
type(cmdline) :: cline
type(image)   :: vol
type(string)  :: vol_file
real          :: smpd
integer       :: ldim(3), ifoo
call cline%parse_oldschool
if( .not. cline%defined('vol') ) THROW_HARD('The vol keyword is required; expected vol=volume.mrc')
if( .not. cline%defined('smpd') ) THROW_HARD('The smpd keyword is required; expected smpd=<Angstrom per voxel>')
vol_file = cline%get_carg('vol')
smpd     = cline%get_rarg('smpd')
call find_ldim_nptcls(vol_file, ldim, ifoo)
write(logfhandle,'(a,3(1x,i0))') 'Input volume dimensions:', ldim
call vol%new(ldim, smpd)
call vol%read(vol_file)
call vol%write(string('vol_simple.mrc'))
call vol%kill
call cline%kill
end program simple_test_mrc_validate
