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
character(len=100)       :: filename_exp, pdbfile_ref
integer                  :: offset, peak_thres_level, intensity_level
logical                  :: circle, denoise, debug, use_euclids, use_zscores
logical                  :: use_valid_corr, use_cs_thres
! for Henry's test
type(nanoparticle)       :: nano
real                     :: a(3), corr
logical                  :: use_auto_corr_thres
type(parameters), target :: params

! Inputs
filename_exp     = 'rec_merged.mrc'
pdbfile_ref      = 'reference.pdb'
element          = 'PT'
smpd             = 0.358
offset           = 2
peak_thres_level = 2
dist_thres       = 3.
intensity_level  = 2
circle           = .true.  ! whether to use cube or sphere for correlation calculation
denoise          = .false. ! whether to denoise the experimental volume before picking
debug            = .false. ! whether to write additional output files for the purpose of debugging
use_euclids      = .false. ! whether to use euclidean distance between boximg and simulated atom instead of correlation
use_zscores      = .false. ! whether to use z-scores when determining outlier correlation scores
use_valid_corr   = .true.  ! whether to discard boxes based on a low valid correlation
use_cs_thres     = .true.  ! whether to discard boxes based on a low contact score
mskdiam          = 39.3

print *, 'NEW METHOD: '
call test_exp4%new(smpd, element, filename_exp, peak_thres_level, offset, denoise, use_euclids, mskdiam, intensity_level)
call test_exp4%simulate_atom()
call test_exp4%setup_iterators()
call test_exp4%match_boxes(circle=circle)
if (debug) call test_exp4%write_dist(  'corr_dist_before_high_filter.csv','corr'   )
if (debug) call test_exp4%write_dist(   'int_dist_before_high_filter.csv','avg_int')
if (debug) call test_exp4%write_dist('euclid_dist_before_high_filter.csv','euclid' )
if (debug) then
    call test_exp4%find_centers()
    call test_exp4%write_positions_and_scores('pos_and_scores_centers_before_high_filter.csv','centers'    )
    call test_exp4%write_positions_and_scores('pos_and_intensities_before_high_filter.csv'   ,'intensities')
    call test_exp4%write_positions_and_scores('pos_and_euclids_before_high_filter.csv'       ,'euclid'     )
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
call test_exp4%refine_threshold(20,max_thres=0.7)
!call test_exp4%write_positions_and_scores('pos_and_scores_centers.csv','centers')
call test_exp4%discard_atoms(use_valid_corr=use_valid_corr, use_cs_thres=use_cs_thres)
!call test_exp4%write_positions_and_scores('pos_and_scores_centers_discarded.csv','centers')
call test_exp4%calc_per_atom_corr
! OUTPUT FILES
if (debug) call test_exp4%write_positions_and_scores('pos_and_scores.csv','pixels')
if (debug) call test_exp4%write_positions_and_scores('pos_and_scores_centers.csv','centers')
if (debug) call test_exp4%write_positions_and_scores('pos_and_intensities.csv','intensities')
if (debug) call test_exp4%write_positions_and_scores('pos_and_euclids.csv','euclid')
call test_exp4%write_pdb('experimental_centers_TEST')
call test_exp4%compare_pick('experimental_centers_TEST.pdb',trim(pdbfile_ref))
if (debug) call test_exp4%write_NP_image('result.mrc')
corr = test_exp4%whole_map_correlation('reference.pdb')
print *, 'Correlation to reference.pdb is: ', corr
corr = test_exp4%whole_map_correlation('startvol.mrc.pdb')
print *, 'Correlation to startvol.mrc.pdb is: ', corr
corr = test_exp4%whole_map_correlation()
print *, 'Correlation to rec_merged.mrc is: ', corr
call test_exp4%kill
print *, ' '

! Henry's method (for comparison)
! print *, 'OLD METHOD: '
! params_glob         => params
! params_glob%element = element
! params_glob%smpd    = smpd
! use_auto_corr_thres = .true.
! call nano%new(filename_exp,msk=mskdiam)
! call nano%identify_atomic_pos(a, l_fit_lattice=.true., use_cs_thres=use_cs_thres,&
!                 &use_auto_corr_thres=use_auto_corr_thres)
! call nano%kill

end program simple_test_nano_detect_atoms