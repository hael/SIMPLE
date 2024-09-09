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
real                     :: smpd, mskdiam, corr_thres, cs_thres
character(len=2)         :: element
character(len=100)       :: filename_exp
logical                  :: use_valid_corr, use_cs_thres
type(nano_picker)        :: test_exp2
type(nanoparticle)       :: nano
type(image)              :: simatms, raw_img
real                     :: a(3), corr
logical                  :: use_auto_corr_thres
type(parameters), target :: params
integer                  :: cs_thres_int, ldim(3)

! Inputs
filename_exp     = 'rec_merged.mrc'
element          = 'PT'
smpd             = 0.358
mskdiam          = 28.4
corr_thres       = 0.3
cs_thres         = 4.0 ! setting higher to demonstrate new method


!Henry's method (for comparison)
print *, 'OLD METHOD: '
params_glob            => params
params_glob%element    = element
params_glob%smpd       = smpd
params_glob%corr_thres = corr_thres
cs_thres_int           = anint(cs_thres)
use_auto_corr_thres    = .false.
call nano%new(filename_exp,msk=(mskdiam / smpd)/2)
call nano%identify_atomic_pos(a, l_fit_lattice=.true., use_cs_thres=use_cs_thres,&
               &use_auto_corr_thres=use_auto_corr_thres, cs_thres=cs_thres_int)
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