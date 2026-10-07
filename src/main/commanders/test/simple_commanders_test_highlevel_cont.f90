!@descr: continuous 3D pose refinement gate of simple_commanders_test_highlevel
submodule (simple_commanders_test_highlevel) simple_commanders_test_highlevel_cont
implicit none
#include "simple_local_flags.inc"

contains

module subroutine exec_test_cont_refine3D_1jxy( self, cline )
    class(commander_test_cont_refine3D_1jxy), intent(inout) :: self
    class(cmdline),                           intent(inout) :: cline
    type(parameters) :: params
    type(string)     :: cwd_saved, fixture_root
    integer          :: status
    logical          :: all_ok
    call simple_getcwd(cwd_saved)
    fixture_root = filepath(cwd_saved, 'test_cont_refine3D_1jxy_'//int2str(get_process_id()))
    if( dir_exists(fixture_root) ) call simple_rmdir(fixture_root)
    call simple_mkdir(fixture_root)
    call simple_chdir(fixture_root, status)
    if( status /= 0 ) THROW_HARD('TEST_CONT_REFINE3D_1JXY FAILED: could not enter fixture directory')
    call params%new(cline)
    all_ok = .true.
    call run_cont_refine3D_gate(params%nthr, all_ok)
    call simple_chdir(cwd_saved, status)
    if( status /= 0 ) THROW_HARD('TEST_CONT_REFINE3D_1JXY FAILED: could not restore original directory')
    if( all_ok )then
        call simple_rmdir(fixture_root)
        write(logfhandle,'(a)') 'PASS: refine=cont validated against the simulation truth'
        call simple_end('**** SIMPLE_TEST_CONT_REFINE3D_1JXY NORMAL STOP ****')
    else
        THROW_HARD('TEST_CONT_REFINE3D_1JXY FAILED')
    endif
end subroutine exec_test_cont_refine3D_1jxy

!> The end-to-end gate of continuous Cartesian pose refinement (N20 of the pose_cont
!  refactoring plan, stages a-e and h; f and g follow with the polish in Phase 9). Fixture: the
!  recipe of E25 (the embedded 1JYX model, box 144, 1.3 A, 1 000 particles, the CTF spread of
!  E25, SNR 1, no simulated shift), simulated by simulate_particles under the fixed test seed,
!  in a project whose ptcl3D poses are the truth perturbed by exactly 15 degrees about a random
!  axis and 2 pixels in a random direction (the seed of a discrete workflow: state, half-set,
!  projection direction 1). Correctness is measured against the truth orientations (ruling R5):
!  per-particle rotation and shift error, fraction inside the basin width, and the frame- and
!  hand-independent pair metric where a run builds its own reference.
!  a  refine3D refine=cont objfun=euclid, shared memory, one pass against the truth reference
!     at lp 8 A with the E25 bounds (trs 5, athres_cont 15); no polar in-plane route runs
!  b  the same with objfun=cc
!  c  stage a distributed over two parts: poses agree with stage a
!  d  a-c (and e): no polar reprojection model or probability table in the run directory, and
!     the projection direction no discrete pass would leave alone is the seed's in every particle
!  e  refine3D_auto pose_cont=only from the perturbed project: runs to maxits, samples nsample
!     particles per iteration, runs no pose initialization or registration pass, and reduces
!     the pose error
!  f  refine3D refine=neigh against the truth reference with pose_cont=yes and with no, two
!     trailing iterations over 30% samples: with the polish the error over the updated particles
!     is no worse; the polish refines exactly the sample (corr_cart on the updated particles only,
!     the others keep their seeds) after the discrete pass's polar in-plane route; the discrete
!     iteration after a polish scores finite; the last assembly represents the final project's N
!  g  refine3D_auto pose_cont=yes from the perturbed project (ref_pose_init=cc as in h):
!     completes, the polish follows the last main-loop iteration over its whole sample, and the
!     error is no worse than h's polar run without the polish
!  h  the normal entry: a polar refine3D_auto from the perturbed project, its poses first
!     initialized against the independent truth map (ref_pose_init=cc, the polar workflow's path
!     for poor input poses), then refine3D_auto pose_cont=only on its output: the polar run
!     improves on the seeds (the precondition of the stage), and the continuation accepts its
!     poses and leaves the pose error no worse. 1JYX is D2-symmetric and the gate runs in c1, so
!     the global searches of the polar workflow assign any of the four symmetry-related poses
!     (Phase 8 finding); stage h measures the error up to them (sym_pose_error)
!  Every metric goes to metrics.tsv.
subroutine run_cont_refine3D_gate( nthr, all_ok )
    use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
    use simple_atoms,              only: atoms
    use simple_molecule_data,      only: molecule_data, betagal_1jyx
    use simple_ui,                 only: make_ui
    use simple_image,              only: image
    use simple_commanders_refine3D, only: commander_refine3D, commander_refine3D_auto
    use simple_cartft_pose_opt,    only: right_increment_rotation
    use simple_refine3D_fnames,    only: refine3D_reproj_model_fname, refine3D_trail_manifest_fname
    use simple_trail_chain_manifest, only: trail_chain_manifest, TRAIL_MANIFEST_OK
    use simple_test_gate,          only: test_gate, NO_FLOOR
    use simple_test_truth_metrics, only: pair_pose_error
    integer, intent(in)    :: nthr
    logical, intent(inout) :: all_ok
    character(len=*), parameter :: GATE_DIR   = 'cont_refine3D_gate'
    character(len=*), parameter :: ATOM_VOL   = '1JYX_atoms.mrc'
    character(len=*), parameter :: TRUTH_VOL  = '1JYX.mrc'
    character(len=*), parameter :: PTCL_STK   = 'simulated_particles.mrc'
    character(len=*), parameter :: TRUTH_ORIS = 'simulated_oris.txt'
    character(len=*), parameter :: PROJNAME   = 'cont_gate'
    integer, parameter :: BOX        = 144
    integer, parameter :: NPTCLS     = 1000
    real,    parameter :: SMPD       = 1.3
    real,    parameter :: MSKDIAM    = 120.
    real,    parameter :: MSKRAD     = 46.          !< pixels, the soft mask of the truth reference (E25)
    real,    parameter :: SNR        = 1.
    real,    parameter :: LP         = 8.
    real,    parameter :: TRS        = 5., ATHRES_CONT = 15.   !< the bounds of E25 and of Phase 0 (O6)
    real,    parameter :: ROT_ERR    = 15.          !< injected rotation (deg)
    real,    parameter :: SH_ERR     = 2.           !< injected shift (pixels)
    integer, parameter :: PERTURBATION_SEED = 20260915
    integer, parameter :: NPAIRS     = 20000
    integer, parameter :: GATE_SEED  = 20260930
    integer, parameter :: AUTO_MAXITS = 4, AUTO_NSAMPLE = 500, CONT_MAXITS = 2
    integer, parameter :: POLISH_MAXITS = 2
    real,    parameter :: POLISH_UPDATE_FRAC = 0.3
    ! Floors, declared before the first run (section 10.1; ruling R5: no map or FSC floor).
    ! Stages a-c run E25 through the production path: E25's floors at 8 A apply, F1 a tenth of
    ! the injected error (median rotation <= 1.5 deg, median shift <= 0.2 px) and F2 the
    ! fraction inside the basin width lp/(mskdiam/2) = 7.64 deg >= 0.95.
    real,    parameter :: MAX_MED_ROT = 0.1*ROT_ERR, MAX_MED_SH = 0.1*SH_ERR, MIN_BASIN = 0.95
    ! Stage c: the distributed run prepares the same references and particles; it differs from
    ! stage a only at the rounding level of threaded evaluation (E25 differs between threaded
    ! runs in the third decimal of a degree): median difference <= 0.01 deg, at most 1% of the
    ! particles differing by more than 0.1 deg.
    real,    parameter :: MAX_MED_DIFF = 0.01, DIFF_LIM = 0.1, MAX_FRAC_DIFF = 0.01
    ! Stage e refines from a reference reconstructed from the perturbed seeds, each particle in
    ! 2 of the 4 iterations: the floor is the reduction rule of 10.1 at half the injected error
    ! (median rotation <= 7.5 deg, median shift <= 1 px), and below the seed's error.
    real,    parameter :: MAX_MED_ROT_E = 0.5*ROT_ERR, MAX_MED_SH_E = 0.5*SH_ERR
    ! Stage h: the continuation leaves the pair pose error no worse than the polar run's, within
    ! 0.1 deg (the rounding-level variation of a 20 000-pair median)
    real,    parameter :: MAX_PAIR_GAIN_LOSS = 0.1
    ! Stages f and g (Phase 9; declared before their first run): the polish leaves the error
    ! (symmetry-aware, over the particles the passes updated) no worse than the same run without
    ! it, within the 0.1 deg of stage h
    integer, parameter :: NSYM_1JYX = 4             !< order of the D2 point group of 1JYX
    real,    parameter :: SYM_CLUSTER_RAD = 3.      !< deg, the radius of a symmetry-mate cluster
    type(commander_simulate_particles) :: xsim
    type(commander_new_project)        :: xnew_project
    type(commander_refine3D)           :: xrefine3D
    type(commander_refine3D_auto)      :: xauto
    type(cmdline)       :: cl
    type(atoms)         :: molecule
    type(molecule_data) :: mol
    type(image)         :: vol
    type(sp_project)    :: seeded, out_a, out_c, run_proj
    type(oris)          :: truth, stats
    type(ctfparams)     :: ctfvars
    type(test_gate)     :: gate
    type(string)        :: root, gate_root, stk_abs, truth_abs, proj_dir
    integer, allocatable :: inds(:)
    real    :: med_rot, med_sh, basin, seed_rot, seed_sh, pair_seed, pair_before, pair_after, frac5
    real    :: sym_seed, sym_before, sym_after, centers(3,3,NSYM_1JYX), sym_nopolish, sym_polish
    real    :: corr_mean
    integer :: sym_convention, conv_tmp
    logical :: updated(NPTCLS), l_sample_ok, l_finite
    type(sp_project) :: out_f
    real    :: diffs(NPTCLS)
    integer :: i, status
    call make_ui
    write(logfhandle,'(a)') '>>> TEST_CONT_REFINE3D_1JXY: refine=cont gate'
    call simple_getcwd(root)
    if( file_exists(GATE_DIR) ) call simple_rmdir(GATE_DIR, status)
    call simple_mkdir(GATE_DIR)
    call simple_chdir(GATE_DIR, status)
    if( status /= 0 ) THROW_HARD('Could not enter '//GATE_DIR)
    call simple_getcwd(gate_root)
    call gate%new(string('metrics.tsv'))
    ! ---- the truth: 1JYX centred; the reference is its soft-masked copy (E25) ----
    mol = betagal_1jyx()
    call molecule%pdb2mrc(volfile=string(ATOM_VOL), smpd=SMPD, mol=mol, center_pdb=.true., vol_dim=[BOX,BOX,BOX])
    call molecule%kill
    call vol%new([BOX,BOX,BOX], SMPD)
    call vol%read(string(ATOM_VOL))
    call vol%mask3D_soft(MSKRAD, backgr=0.)
    call vol%write(string(TRUTH_VOL))
    ! the half-maps of an external reference of independent provenance (ref_pose_init=cc, stage h)
    call vol%write(string('1JYX_even.mrc'))
    call vol%write(string('1JYX_odd.mrc'))
    call vol%kill
    truth_abs = simple_abspath(string(TRUTH_VOL))
    ! ---- particles: the CTF spread and noise of E25, no simulated shift ----
    call cl%set('prg',      'simulate_particles')
    call cl%set('vol1',     ATOM_VOL)
    call cl%set('smpd',     SMPD)
    call cl%set('mskdiam',  MSKDIAM)
    call cl%set('nthr',     nthr)
    call cl%set('nptcls',   NPTCLS)
    call cl%set('pgrp',     'c1')
    call cl%set('snr',      SNR)
    call cl%set('ctf',      'yes')
    call cl%set('kv',       300.)
    call cl%set('cs',       2.7)
    call cl%set('fraca',    0.1)
    call cl%set('defocus',  1.5)
    call cl%set('dferr',    0.5)
    call cl%set('astigerr', 0.1)
    call cl%set('bfac',     0.)
    call cl%set('sherr',    0.)
    call xsim%execute(cl)
    call cl%kill
    stk_abs = simple_abspath(string(PTCL_STK))
    call truth%new(NPTCLS, is_ptcl=.true.)
    call truth%read(string(TRUTH_ORIS), [1,NPTCLS])
    ! ---- the seeded project: truth CTF, perturbed poses, the fields of a discrete workflow ----
    call cl%set('projname',  PROJNAME)
    call cl%set('qsys_name', 'local')
    call xnew_project%execute(cl)             ! creates and enters PROJNAME/
    call cl%kill
    call simple_getcwd(proj_dir)
    call seeded%read(simple_abspath(string(PROJNAME//'.simple')))
    ctfvars%smpd    = SMPD
    ctfvars%kv      = 300.
    ctfvars%cs      = 2.7
    ctfvars%fraca   = 0.1
    ctfvars%ctfflag = CTFFLAG_YES
    call seeded%add_stk(stk_abs, ctfvars)
    call seeded%os_cls2D%new(1, is_ptcl=.false.)
    call seeded%os_cls2D%set_all2single('state', 1.)
    call set_fixed_seed(PERTURBATION_SEED)
    do i = 1, NPTCLS
        call seeded%os_ptcl3D%set(i, 'dfx',    truth%get(i, 'dfx'))
        call seeded%os_ptcl3D%set(i, 'dfy',    truth%get(i, 'dfy'))
        call seeded%os_ptcl3D%set(i, 'angast', truth%get(i, 'angast'))
        call seeded%os_ptcl3D%set_state(i, 1)
        call seeded%os_ptcl3D%set(i, 'eo',   real(mod(i,2)))
        call seeded%os_ptcl3D%set(i, 'proj', 1.)
        call perturb(i)
        call seeded%os_ptcl2D%set(i, 'class', 1.)
        call seeded%os_ptcl2D%set(i, 'corr',  0.5)
        call seeded%os_ptcl2D%set(i, 'dfx',    truth%get(i, 'dfx'))
        call seeded%os_ptcl2D%set(i, 'dfy',    truth%get(i, 'dfy'))
        call seeded%os_ptcl2D%set(i, 'angast', truth%get(i, 'angast'))
        call seeded%os_ptcl2D%set_state(i, 1)
    enddo
    allocate(inds(NPTCLS))
    inds = [(i, i=1,NPTCLS)]
    call pose_errors(seeded%os_ptcl3D, seed_rot, seed_sh, basin)
    pair_seed = pair_pose_error(seeded%os_ptcl3D, truth, inds, inds, NPAIRS, GATE_SEED, frac5)
    call gate%report('seed_pair_pose_error_deg',       pair_seed)
    call sym_pose_error(seeded%os_ptcl3D, .true., sym_convention, centers, sym_seed)
    call gate%report('seed_sym_pose_error_deg',        sym_seed)
    call gate%report('seed_median_rotation_error_deg', seed_rot)
    call gate%report('seed_median_shift_error_px',     seed_sh)
    call simple_chdir(gate_root, status)
    ! ---- a, b, c: refine3D refine=cont against the truth reference ----
    call run_refine3D_cont('a_euclid', 'euclid', 0, out_a)
    call gate_e25('a', out_a)
    ! a pure Cartesian pass runs no separate polar in-plane refinement (inpl_cont=no)
    call gate%check('a_no_polar_inplane_route', .not. any(out_a%os_ptcl3D%get_all('cont_inpl_attempted') > 0.5))
    call run_refine3D_cont('b_cc', 'cc', 0, run_proj)
    call gate_e25('b', run_proj)
    call run_refine3D_cont('c_distr', 'euclid', 2, out_c)
    call gate_e25('c', out_c)
    do i = 1, NPTCLS
        diffs(i) = rad2deg(geodesic_angle(out_a%os_ptcl3D, out_c%os_ptcl3D, i))
    enddo
    call gate%metric('c_median_rotation_difference_to_a_deg', median_nocopy(diffs), MAX_MED_DIFF, &
        &median_nocopy(diffs) <= MAX_MED_DIFF)
    call gate%metric('c_fraction_differing_by_more_than_0.1deg', real(count(diffs > DIFF_LIM))/real(NPTCLS), &
        &MAX_FRAC_DIFF, real(count(diffs > DIFF_LIM))/real(NPTCLS) <= MAX_FRAC_DIFF)
    ! ---- e: refine3D_auto pose_cont=only from the perturbed project ----
    call enter_stage('e_auto')
    call cl%set('prg',       'refine3D_auto')
    call cl%set('projfile',  simple_abspath(string(PROJNAME//'.simple')))
    call cl%set('mkdir',     'no')
    call cl%set('pose_cont', 'only')
    call cl%set('pgrp',      'c1')
    call cl%set('mskdiam',   MSKDIAM)
    call cl%set('nsample',   AUTO_NSAMPLE)
    call cl%set('maxits',    AUTO_MAXITS)
    call cl%set('nthr',      nthr)
    call xauto%execute(cl)
    call cl%kill
    call run_proj%read(simple_abspath(string(PROJNAME//'.simple')))
    call stats%new(1, is_ptcl=.false.)
    call stats%read(string(STATS_FILE))
    call gate%check('e_ran_to_maxits', nint(stats%get(1, 'ITERATION')) == AUTO_MAXITS)
    call gate%check('e_nsample_particles_per_iteration', &
        &abs(stats%get(1, 'PERCEN_PARTICLES_SAMPLED') - 100.*real(AUTO_NSAMPLE)/real(NPTCLS)) < 0.01)
    ! the refine=cont convergence rule's stable fraction of the last iteration (maxits sets minits here)
    call gate%report('e_cont_stable_pct', stats%get(1, 'CONT_STABLE_PCT'))
    call stats%kill
    call check_no_polar('e', run_proj)
    call pose_errors(run_proj%os_ptcl3D, med_rot, med_sh, basin)
    call gate%metric('e_median_rotation_error_deg', med_rot, MAX_MED_ROT_E, med_rot <= MAX_MED_ROT_E .and. med_rot < seed_rot)
    call gate%metric('e_median_shift_error_px',     med_sh,  MAX_MED_SH_E,  med_sh  <= MAX_MED_SH_E  .and. med_sh  < seed_sh)
    call gate%report('e_fraction_in_basin_8A', basin)
    call simple_chdir(gate_root, status)
    ! ---- f: a discrete refine3D with and without the polish ----
    write(logfhandle,'(a)') '>>> CONT_REFINE3D_1JXY STAGE f: pose_cont=no'
    call run_refine3D_discrete('f_nopolish', 'no', out_f)
    updated = out_f%os_ptcl3D%get_all('updatecnt') > 0.5
    call sym_pose_error(out_f%os_ptcl3D, .true., conv_tmp, centers, sym_nopolish, updated)
    call gate%report('f_nopolish_sym_pose_error_deg', sym_nopolish)
    write(logfhandle,'(a)') '>>> CONT_REFINE3D_1JXY STAGE f: pose_cont=yes'
    call run_refine3D_discrete('f_polish', 'yes', out_f)
    updated = out_f%os_ptcl3D%get_all('updatecnt') > 0.5
    call gate%report('f_updated_particles', real(count(updated)))
    call sym_pose_error(out_f%os_ptcl3D, .true., conv_tmp, centers, sym_polish, updated)
    call gate%metric('f_polish_sym_pose_error_deg', sym_polish, sym_nopolish + MAX_PAIR_GAIN_LOSS, &
        &sym_polish <= sym_nopolish + MAX_PAIR_GAIN_LOSS)
    ! the polish refines the sample of its discrete pass and nothing else
    l_sample_ok = count(updated) > 0 .and. count(updated) < NPTCLS
    do i = 1, NPTCLS
        if( updated(i) )then
            if( .not. out_f%os_ptcl3D%get(i, 'corr_cart') > 0. ) l_sample_ok = .false.
        else
            if( out_f%os_ptcl3D%get(i, 'corr_cart') /= 0. ) l_sample_ok = .false.
            if( any(out_f%os_ptcl3D%get_euler(i) /= seeded%os_ptcl3D%get_euler(i)) ) l_sample_ok = .false.
        endif
    enddo
    call gate%check('f_polish_refines_exactly_the_sample', l_sample_ok)
    ! the discrete pass before the polish ran its polar in-plane route (inpl_cont=yes); the
    ! polish runs with inpl_cont=no and leaves cont_inpl_attempted as the discrete pass set it
    call gate%check('f_discrete_pass_ran_the_polar_inplane_route', &
        &any(out_f%os_ptcl3D%get_all('cont_inpl_attempted') > 0.5 .and. updated))
    ! the discrete iteration after a polish ran and scored its particles
    l_finite  = .true.
    corr_mean = 0.
    do i = 1, NPTCLS
        if( .not. updated(i) ) cycle
        if( .not. ieee_is_finite(out_f%os_ptcl3D%get(i, 'corr')) ) l_finite = .false.
        corr_mean = corr_mean + out_f%os_ptcl3D%get(i, 'corr')
    enddo
    corr_mean = corr_mean / real(max(1, count(updated)))
    call gate%check('f_discrete_scores_finite_after_polish', l_finite .and. corr_mean > 0.)
    call gate%report('f_polish_mean_corr', corr_mean)
    ! ---- h: a polar refine3D_auto, then the continuation on its output ----
    call enter_stage('h_continuation')
    call cl%set('prg',       'refine3D_auto')
    call cl%set('projfile',  simple_abspath(string(PROJNAME//'.simple')))
    call cl%set('mkdir',     'no')
    call cl%set('vol1',      truth_abs)
    call cl%set('ref_pose_init', 'cc')
    call cl%set('pgrp',      'c1')
    call cl%set('mskdiam',   MSKDIAM)
    call cl%set('nthr',      nthr)
    call xauto%execute(cl)
    call cl%kill
    call run_proj%read(simple_abspath(string(PROJNAME//'.simple')))
    pair_before = pair_pose_error(run_proj%os_ptcl3D, truth, inds, inds, NPAIRS, GATE_SEED, frac5)
    call gate%report('h_polar_pair_pose_error_deg', pair_before)
    call sym_pose_error(run_proj%os_ptcl3D, .true., sym_convention, centers, sym_before)
    call gate%metric('h_polar_sym_pose_error_deg', sym_before, sym_seed, sym_before < sym_seed)
    call pose_errors(run_proj%os_ptcl3D, med_rot, med_sh, basin)
    call gate%report('h_polar_median_rotation_error_deg', med_rot)
    call cl%set('prg',       'refine3D_auto')
    call cl%set('projfile',  simple_abspath(string(PROJNAME//'.simple')))
    call cl%set('mkdir',     'no')
    call cl%set('pose_cont', 'only')
    call cl%set('pgrp',      'c1')
    call cl%set('mskdiam',   MSKDIAM)
    call cl%set('maxits',    CONT_MAXITS)
    call cl%set('nthr',      nthr)
    call xauto%execute(cl)
    call cl%kill
    call run_proj%read(simple_abspath(string(PROJNAME//'.simple')))
    pair_after = pair_pose_error(run_proj%os_ptcl3D, truth, inds, inds, NPAIRS, GATE_SEED, frac5)
    call gate%report('h_continued_pair_pose_error_deg', pair_after)
    ! against the symmetry mates and frame of the polar result
    call sym_pose_error(run_proj%os_ptcl3D, .false., sym_convention, centers, sym_after)
    call gate%metric('h_continued_sym_pose_error_deg', sym_after, sym_before + MAX_PAIR_GAIN_LOSS, &
        &sym_after <= sym_before + MAX_PAIR_GAIN_LOSS)
    call pose_errors(run_proj%os_ptcl3D, med_rot, med_sh, basin)
    call gate%report('h_continued_median_rotation_error_deg', med_rot)
    call gate%report('h_continued_median_shift_error_px',     med_sh)
    call simple_chdir(gate_root, status)
    ! ---- g: refine3D_auto pose_cont=yes, as h's polar run with the polish ----
    write(logfhandle,'(a)') '>>> CONT_REFINE3D_1JXY STAGE g'
    call enter_stage('g_auto_polish')
    call cl%set('prg',           'refine3D_auto')
    call cl%set('projfile',      simple_abspath(string(PROJNAME//'.simple')))
    call cl%set('mkdir',         'no')
    call cl%set('vol1',          truth_abs)
    call cl%set('ref_pose_init', 'cc')
    call cl%set('pose_cont',     'yes')
    call cl%set('pgrp',          'c1')
    call cl%set('mskdiam',       MSKDIAM)
    call cl%set('nthr',          nthr)
    call xauto%execute(cl)
    call cl%kill
    call run_proj%read(simple_abspath(string(PROJNAME//'.simple')))
    call stats%new(1, is_ptcl=.false.)
    call stats%read(string(STATS_FILE))
    call gate%check('g_polish_followed_the_last_iteration_over_its_sample', &
        &nint(stats%get(1, 'POSE_CONT_ATTEMPTS')) == nint(stats%get(1, 'PERCEN_PARTICLES_SAMPLED')*real(NPTCLS)/100.))
    call gate%report('g_last_polish_improved_pct', stats%get(1, 'POSE_CONT_IMPROVED_PCT'))
    call stats%kill
    call gate%check('g_every_particle_polished', all(run_proj%os_ptcl3D%get_all('corr_cart') > 0.))
    call sym_pose_error(run_proj%os_ptcl3D, .true., conv_tmp, centers, sym_polish)
    call gate%metric('g_sym_pose_error_deg', sym_polish, sym_before + MAX_PAIR_GAIN_LOSS, &
        &sym_polish <= sym_before + MAX_PAIR_GAIN_LOSS)
    call simple_chdir(gate_root, status)
    all_ok = all_ok .and. gate%passed()
    call gate%kill
    call seeded%kill
    call out_f%kill
    call out_a%kill
    call out_c%kill
    call run_proj%kill
    call truth%kill
    call simple_chdir(root, status)

  contains

    !> the truth pose of particle i perturbed by exactly ROT_ERR about a random axis and
    !! SH_ERR in a random direction (as E25 perturbs its seeds)
    subroutine perturb( iptcl )
        integer, intent(in) :: iptcl
        real(dp) :: rotmat(3,3), axis(3), z, azimuth, radial
        real     :: uniform(3)
        call random_number(uniform)
        z       = 2._dp*real(uniform(1),dp) - 1._dp
        azimuth = 2._dp*DPI*real(uniform(2),dp)
        radial  = sqrt(max(0._dp, 1._dp - z*z))
        axis    = [radial*cos(azimuth), radial*sin(azimuth), z]
        rotmat  = right_increment_rotation(real(truth%get_mat(iptcl),dp), real(deg2rad(ROT_ERR),dp)*axis)
        call seeded%os_ptcl3D%set_euler(iptcl, real(dm2euler(rotmat)))
        azimuth = 2._dp*DPI*real(uniform(3),dp)
        call seeded%os_ptcl3D%set_shift(iptcl, truth%get_2Dshift(iptcl) + SH_ERR*real([cos(azimuth), sin(azimuth)]))
    end subroutine perturb

    !> the rotation angle (radians) between the poses of particle i in two fields
    real function geodesic_angle( os1, os2, iptcl )
        class(oris), intent(in) :: os1, os2
        integer,     intent(in) :: iptcl
        real :: m1(3,3), m2(3,3)
        m1 = os1%get_mat(iptcl)
        m2 = os2%get_mat(iptcl)
        geodesic_angle = acos(max(-1., min(1., (sum(m1*m2) - 1.)/2.)))
    end function geodesic_angle

    !> median rotation (deg) and shift (pixels) errors against the truth, and the fraction of
    !! particles inside the basin width at LP
    subroutine pose_errors( os, med_rot_out, med_sh_out, basin_out )
        class(oris), intent(in)  :: os
        real,        intent(out) :: med_rot_out, med_sh_out, basin_out
        real :: rot(NPTCLS), sh(NPTCLS)
        integer :: j
        do j = 1, NPTCLS
            rot(j) = geodesic_angle(os, truth, j)
            sh(j)  = norm2(os%get_2Dshift(j) - truth%get_2Dshift(j))
        enddo
        basin_out   = real(count(rot <= LP/(MSKDIAM/2.)))/real(NPTCLS)
        med_rot_out = rad2deg(median_nocopy(rot))
        med_sh_out  = median_nocopy(sh)
    end subroutine pose_errors

    !> The median rotation error (deg) against the truth up to the symmetry mates of 1JYX and a
    !! global frame: the frame difference of particle i, D = R_truth^T R or R R_truth^T (both
    !! conventions are tried and the tighter is kept), takes one of NSYM_1JYX values for a correct
    !! assignment; the values are found as the centres of the densest clusters of radius
    !! SYM_CLUSTER_RAD, and the error of particle i is its distance to the nearest centre. With
    !! l_new the centres and the convention are found and returned, otherwise those given are used.
    subroutine sym_pose_error( os, l_new, convention, cent, err, sel )
        class(oris),       intent(in)    :: os
        logical,           intent(in)    :: l_new
        integer,           intent(inout) :: convention
        real,              intent(inout) :: cent(3,3,NSYM_1JYX)
        real,              intent(out)   :: err
        logical, optional, intent(in)    :: sel(NPTCLS)   !< the particles measured (all by default)
        real    :: d(3,3,NPTCLS,2), c(3,3,NSYM_1JYX,2), med(2), e(NPTCLS), rt(3,3), re(3,3)
        logical :: claimed(NPTCLS), lsel(NPTCLS)
        integer :: j, k, jj, iconv, best, nbest, nn
        lsel = .true.
        if( present(sel) ) lsel = sel
        do j = 1, NPTCLS
            rt = truth%get_mat(j)
            re = os%get_mat(j)
            d(:,:,j,1) = matmul(transpose(rt), re)
            d(:,:,j,2) = matmul(re, transpose(rt))
        enddo
        do iconv = 1, 2
            if( l_new )then
                claimed = .not. lsel
                do k = 1, NSYM_1JYX
                    best  = 0
                    nbest = -1
                    do j = 1, NPTCLS
                        if( claimed(j) ) cycle
                        nn = 0
                        do jj = 1, NPTCLS
                            if( .not. claimed(jj) .and. mat_angle(d(:,:,j,iconv), d(:,:,jj,iconv)) <= SYM_CLUSTER_RAD ) nn = nn + 1
                        enddo
                        if( nn > nbest )then
                            nbest = nn
                            best  = j
                        endif
                    enddo
                    if( best == 0 ) best = 1
                    c(:,:,k,iconv) = d(:,:,best,iconv)
                    do jj = 1, NPTCLS
                        if( mat_angle(d(:,:,best,iconv), d(:,:,jj,iconv)) <= SYM_CLUSTER_RAD ) claimed(jj) = .true.
                    enddo
                enddo
            else
                c(:,:,:,iconv) = cent
                if( iconv /= convention ) cycle
            endif
            do j = 1, NPTCLS
                e(j) = huge(1.)
                do k = 1, NSYM_1JYX
                    e(j) = min(e(j), mat_angle(d(:,:,j,iconv), c(:,:,k,iconv)))
                enddo
            enddo
            med(iconv) = median(pack(e, lsel))
        enddo
        if( l_new )then
            convention = merge(1, 2, med(1) <= med(2))
            cent       = c(:,:,:,convention)
        endif
        err = med(convention)
    end subroutine sym_pose_error

    !> the rotation angle (deg) between two rotation matrices
    real function mat_angle( a, b )
        real, intent(in) :: a(3,3), b(3,3)
        mat_angle = rad2deg(acos(max(-1., min(1., (sum(a*b) - 1.)/2.))))
    end function mat_angle

    !> a run directory under the gate holding a copy of the seeded project; entered
    subroutine enter_stage( dir )
        character(len=*), intent(in) :: dir
        call simple_mkdir(dir)
        call simple_chdir(string(dir), status)
        if( status /= 0 ) THROW_HARD('Could not enter '//dir)
        call seeded%write(simple_abspath(string(PROJNAME//'.simple'), check_exists=.false.))
    end subroutine enter_stage

    !> one refine3D refine=cont pass against the truth reference on the shift-then-joint route of
    !! the E25 floors (production runs the joint stage alone); nparts > 0 distributes it
    subroutine run_refine3D_cont( dir, objfun, nparts, out )
        character(len=*), intent(in)    :: dir, objfun
        integer,          intent(in)    :: nparts
        type(sp_project), intent(inout) :: out
        call enter_stage(dir)
        call cl%set('prg',      'refine3D')
        call cl%set('projfile', simple_abspath(string(PROJNAME//'.simple')))
        call cl%set('mkdir',    'no')
        call cl%set('vol1',     truth_abs)
        call cl%set('refine',   'cont')
        call cl%set('objfun',   objfun)
        call cl%set('pgrp',     'c1')
        call cl%set('mskdiam',  MSKDIAM)
        call cl%set('lp',       LP)
        call cl%set('trs',      TRS)
        call cl%set('athres_cont', ATHRES_CONT)
        call cl%set('cont_route',  'shift_then_joint')
        call cl%set('maxits',   1)
        if( nparts > 0 )then
            call cl%set('nparts', nparts)
            call cl%set('nthr',   max(1, nthr/nparts))
        else
            call cl%set('nthr',   nthr)
        endif
        call xrefine3D%execute(cl)
        call cl%kill
        call out%read(simple_abspath(string(PROJNAME//'.simple')))
        call check_no_polar(dir(1:1), out)
        call simple_chdir(gate_root, status)
    end subroutine run_refine3D_cont

    !> two trailing refine3D refine=neigh iterations over 30% samples, with or without the polish
    !! (stage f); the last assembly must represent the N of the final project (C1 of the aftermath)
    subroutine run_refine3D_discrete( dir, pose_cont, out )
        character(len=*), intent(in)    :: dir, pose_cont
        type(sp_project), intent(inout) :: out
        type(trail_chain_manifest) :: manifest
        integer, allocatable :: nrep(:), nsmp(:)
        integer :: mstatus
        call enter_stage(dir)
        call cl%set('prg',         'refine3D')
        call cl%set('projfile',    simple_abspath(string(PROJNAME//'.simple')))
        call cl%set('mkdir',       'no')
        call cl%set('vol1',        truth_abs)
        call cl%set('refine',      'neigh')
        call cl%set('pose_cont',   pose_cont)
        call cl%set('objfun',      'euclid')
        call cl%set('pgrp',        'c1')
        call cl%set('mskdiam',     MSKDIAM)
        call cl%set('lp',          LP)
        call cl%set('trs',         TRS)
        call cl%set('athres_cont', ATHRES_CONT)
        call cl%set('update_frac', POLISH_UPDATE_FRAC)
        call cl%set('trail_rec',   'yes')
        call cl%set('maxits',      POLISH_MAXITS)
        call cl%set('nthr',        nthr)
        call xrefine3D%execute(cl)
        call cl%kill
        call out%read(simple_abspath(string(PROJNAME//'.simple')))
        call manifest%read(refine3D_trail_manifest_fname(1), mstatus)
        call out%os_ptcl3D%get_group_update_counts('state', 1, nrep, nsmp)
        call gate%check(dir//'_trailing_population_of_the_current_sample', &
            &mstatus == TRAIL_MANIFEST_OK .and. nint(manifest%get_mrep()) == nrep(1))
        call gate%report(dir//'_trailing_represented_population', manifest%get_mrep())
        call manifest%kill
        call simple_chdir(gate_root, status)
    end subroutine run_refine3D_discrete

    !> the floors of E25 at 8 A on the poses of a refine3D refine=cont pass
    subroutine gate_e25( tag, out )
        character(len=*), intent(in)    :: tag
        type(sp_project), intent(inout) :: out
        call pose_errors(out%os_ptcl3D, med_rot, med_sh, basin)
        call gate%metric(tag//'_median_rotation_error_deg', med_rot, MAX_MED_ROT, med_rot <= MAX_MED_ROT)
        call gate%metric(tag//'_median_shift_error_px',     med_sh,  MAX_MED_SH,  med_sh  <= MAX_MED_SH)
        call gate%metric(tag//'_fraction_in_basin_8A',      basin,   MIN_BASIN,   basin   >= MIN_BASIN)
        call gate%report(tag//'_mean_corr_cart', out%os_ptcl3D%get_avg('corr_cart'))
    end subroutine gate_e25

    !> stage d: no polar reprojection model, no probability table, and the seed's projection
    !! direction in every particle (a discrete pass would assign its own)
    subroutine check_no_polar( tag, out )
        character(len=*), intent(in) :: tag
        type(sp_project), intent(in) :: out
        type(string), allocatable :: files(:)
        integer :: j, npolar
        logical :: l_seed_proj
        npolar = 0
        if( file_exists(refine3D_reproj_model_fname('even')) ) npolar = npolar + 1
        if( file_exists(refine3D_reproj_model_fname('odd'))  ) npolar = npolar + 1
        call simple_list_files(DIST_FBODY//'*', files)
        npolar = npolar + size(files)
        call simple_list_files(ASSIGNMENT_FBODY//'*', files)
        npolar = npolar + size(files)
        call gate%check('d_'//tag//'_no_polar_model_or_probability_table', npolar == 0)
        l_seed_proj = .true.
        do j = 1, NPTCLS
            if( out%os_ptcl3D%get_proj(j) /= 1 ) l_seed_proj = .false.
        enddo
        call gate%check('d_'//tag//'_no_discrete_pass_assigned_a_projection', l_seed_proj)
        call simple_list_files('*', files)
        call gate%report('d_'//tag//'_files_in_run_directory', real(size(files)))
    end subroutine check_no_polar

end subroutine run_cont_refine3D_gate

end submodule simple_commanders_test_highlevel_cont
