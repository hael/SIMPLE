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
type(nano_picker)        :: test_exp4
real                     :: smpd, dist_thres, mskdiam
character(len=2)         :: element
character(len=100)       :: filename_exp, filename_sim, pdbfile_ref
integer                  :: offset, peak_thres_level
logical                  :: circle, denoise, debug, use_euclids, use_zscores
! for Henry's test
type(nanoparticle)       :: nano
real                     :: a(3)
logical                  :: use_cs_thres, use_auto_corr_thres
type(parameters), target :: params

! Inputs
filename_exp     = 'rec_merged.mrc'
filename_sim     = 'simulated_NP.mrc'
pdbfile_ref      = 'reference.pdb'
element          = 'PT'
smpd             = 0.358
offset           = 2
peak_thres_level = 2
dist_thres       = 3.
circle           = .true.
denoise          = .false.
debug            = .true.
use_euclids      = .false.
use_zscores      = .false.
mskdiam          = 27.

print *, 'NEW METHOD: '
call test_exp4%new(smpd, element, filename_exp, peak_thres_level, offset, denoise, use_euclids, mskdiam)
call test_exp4%simulate_atom()
call test_exp4%setup_iterators()
call test_exp4%match_boxes(circle=circle)
if (debug) call test_exp4%write_dist(  'corr_dist_before_high_filter.csv','corr'   )
if (debug) call test_exp4%write_dist(   'int_dist_before_high_filter.csv','avg_int')
if (debug) call test_exp4%write_dist('euclid_dist_before_high_filter.csv','euclid' )
if (debug) then
    call test_exp4%find_centers()
    call test_exp4%write_positions_and_scores('pos_and_scores_centers_before_high_filter.csv','centers')
    call test_exp4%write_positions_and_scores('pos_and_intensities_before_high_filter.csv','intensities')
    call test_exp4%write_positions_and_scores('pos_and_euclids_before_high_filter.csv','euclid')
end if
!call test_exp4%identify_high_scores(use_zscores=use_zscores)
call test_exp4%identify_threshold()
call test_exp4%apply_threshold
if (debug) call test_exp4%write_dist('corr_dist_after_high_filter.csv','corr'   )
if (debug) call test_exp4%write_dist( 'int_dist_after_high_filter.csv','avg_int')
call test_exp4%distance_filter(dist_thres)
if (debug) call test_exp4%write_dist(  'corr_dist_after_dist_filter_high.csv','corr'   )
if (debug) call test_exp4%write_dist(   'int_dist_after_dist_filter_high.csv','avg_int')
if (debug) call test_exp4%write_dist('euclid_dist_after_dist_filter_high.csv','euclid' )
call test_exp4%find_centers()
call test_exp4%write_positions_and_scores('pos_and_scores_centers_pre_discard.csv','centers')
call test_exp4%discard_atoms()
call test_exp4%calc_per_atom_corr
! OUTPUT FILES
call test_exp4%write_positions_and_scores('pos_and_scores.csv','pixels')
call test_exp4%write_positions_and_scores('pos_and_scores_centers.csv','centers')
call test_exp4%write_positions_and_scores('pos_and_intensities.csv','intensities')
call test_exp4%write_positions_and_scores('pos_and_euclids.csv','euclid')
call test_exp4%write_pdb('experimental_centers_discard')
call test_exp4%compare_pick('experimental_centers_discard.pdb',trim(pdbfile_ref))
call test_exp4%write_NP_image('result.mrc')
call test_exp4%kill
print *, ' '

! Henry's method (for comparison)
print *, 'OLD METHOD: '
params_glob         => params
params_glob%element = element
params_glob%smpd    = smpd
use_cs_thres        = .false.
use_auto_corr_thres = .false.
call nano%new(filename_exp,msk=mskdiam)
call nano%identify_atomic_pos(a, l_fit_lattice=.true., use_cs_thres=use_cs_thres,&
                &use_auto_corr_thres=use_auto_corr_thres)
call nano%kill

end program simple_test_nano_detect_atoms