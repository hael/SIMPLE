program simple_test_nano_detect_atoms

include 'simple_lib.f08'
use simple_nano_detect_atoms
use simple_nano_picker_utils
use simple_nanoparticle
use simple_nanoparticle_utils
use simple_image
use simple_parameters
use simple_strings, only: int2str
implicit none
#include "simple_local_flags.inc"
real                     :: smpd, mskdiam
character(len=2)         :: element
character(len=100)       :: filename_exp
logical                  :: use_valid_corr
type(nano_picker)        :: test_exp2
type(nanoparticle)       :: nano
type(image)              :: simatms, raw_img
real                     :: a(3), corr
type(parameters), target :: params
integer                  :: ldim(3)
! Inputs
filename_exp     = 'rec_merged.mrc'
element          = 'PT'
smpd             = 0.358
mskdiam          = 28.4
!Henry's method (for comparison)
print *, 'OLD METHOD: '
params_glob            => params
params_glob%element    = element
params_glob%smpd       = smpd
call nano%new(filename_exp,msk=(mskdiam / smpd)/2)
call nano%identify_atomic_pos(a, l_fit_lattice=.true.)
call nano%simulate_atoms(simatms)
call simatms%write('simatms_henry.mrc')
call nano%get_img_raw(raw_img)
corr = simatms%real_corr(raw_img)
print *, 'correlation to rec_merged.mrc is ', corr
call nano%write_centers('old_method_positions_in_angstroms',which='valid_corr')
call nano%kill
call simatms%kill
call raw_img%kill
print *, ' '
end program simple_test_nano_detect_atoms