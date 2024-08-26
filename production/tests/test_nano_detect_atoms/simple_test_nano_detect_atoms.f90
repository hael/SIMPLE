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
type(nano_picker)        :: test_exp
real                     :: smpd, dist_thres, mskdiam, corr_thres, cs_thres
character(len=2)         :: element
character(len=100)       :: filename_exp, pdbfile_ref
integer                  :: offset, peak_thres_level, intensity_level
logical                  :: circle, denoise
logical                  :: use_valid_corr, use_cs_thres
! for Henry's test
type(nano_picker)        :: test_exp2
type(nanoparticle)       :: nano
type(image)              :: simatms, raw_img
real                     :: a(3), corr
logical                  :: use_auto_corr_thres
type(parameters), target :: params
integer                  :: cs_thres_int, ldim(3)
! for boximg examination
type(image)              :: boximg
real                     :: center_pos(3)

! Inputs
filename_exp     = 'rec_merged.mrc'
pdbfile_ref      = 'reference.pdb'
element          = 'PT'
smpd             = 0.358
offset           = 2
peak_thres_level = 2
intensity_level  = 2
circle           = .true.  ! whether to use cube or sphere for correlation calculation (.true. = sphere, default value is .false.)
use_valid_corr   = .true.  ! whether to discard boxes based on a low valid correlation
use_cs_thres     = .true.  ! whether to discard boxes based on a low contact score
mskdiam          = 28.4
corr_thres       = 0.3
cs_thres         = 2.0

print *, 'NEW METHOD: '
call test_exp%new(smpd, element, filename_exp, peak_thres_level, offset, mskdiam=mskdiam, intensity_level=intensity_level, circle=circle)
call test_exp%exec_nano_picker(corr_thres=corr_thres,cs_thres=cs_thres)
call test_exp%calc_per_atom_corr('final_positions_and_corrs')

!OUTPUT FILES
call test_exp%write_pdb('experimental_centers_TEST')
call test_exp%write_positions_and_scores('pos_and_scores_centers_refined.csv', 'centers')
corr = test_exp%whole_map_correlation()
print *, 'Correlation to rec_merged.mrc is: ', corr
call test_exp%kill
print *, ' '

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
! call nano%identify_atomic_pos(a, l_fit_lattice=.true., use_cs_thres=use_cs_thres,&
!                 &use_auto_corr_thres=use_auto_corr_thres)
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

! sanity checks
print *, 'SANITY CHECKS'
call test_exp2%new(smpd, element, filename_exp, peak_thres_level, offset, mskdiam=mskdiam, intensity_level=intensity_level, circle=circle)
call test_exp2%set_pdb('rec_merged_ATMS.pdb')
call test_exp2%calc_per_atom_corr('old_method_sanity_check')
call test_exp2%kill
call compare_pick('old_method_positions_in_angstroms.pdb', 'old_method_sanity_check.pdb')

end program simple_test_nano_detect_atoms