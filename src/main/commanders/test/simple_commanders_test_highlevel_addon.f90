!@descr: solve3D add-on gate and snapshot generator of simple_commanders_test_highlevel
submodule (simple_commanders_test_highlevel) simple_commanders_test_highlevel_addon
implicit none
#include "simple_local_flags.inc"

contains

!> solve3D_addon end to end: solve3D on a seeded selection of a first
!! set of simulated particles, the add-on on a larger project that appends a
!! second set, checked against the simulation truth and the base run. The
!! fixture directory (a few hundred MB) is removed on success and kept for
!! inspection on failure.
module subroutine exec_test_solve3D_addon( self, cline )
    class(commander_test_solve3D_addon), intent(inout) :: self
    class(cmdline),                         intent(inout) :: cline
    type(parameters) :: params
    type(string)     :: cwd_saved, fixture_root
    integer          :: status
    logical          :: all_ok
    call simple_getcwd(cwd_saved)
    fixture_root = filepath(cwd_saved, 'test_solve3D_addon_'//int2str(get_process_id()))
    if( dir_exists(fixture_root) ) call simple_rmdir(fixture_root)
    call simple_mkdir(fixture_root)
    call simple_chdir(fixture_root, status)
    if( status /= 0 ) THROW_HARD('TEST_SOLVE3D_ADDON FAILED: could not enter fixture directory')
    call params%new(cline)
    all_ok = .true.
    call run_solve3D_addon_gate(params%nthr, all_ok)
    call simple_chdir(cwd_saved, status)
    if( status /= 0 ) THROW_HARD('TEST_SOLVE3D_ADDON FAILED: could not restore original directory')
    if( all_ok )then
        call simple_rmdir(fixture_root)
        write(logfhandle,'(a)') 'PASS: solve3D_addon validated against the simulation truth and the base run'
        call simple_end('**** SIMPLE_TEST_SOLVE3D_ADDON NORMAL STOP ****')
    else
        THROW_HARD('TEST_SOLVE3D_ADDON FAILED')
    endif
end subroutine exec_test_solve3D_addon

!> Generates append-only cumulative projects for stream/add-on integration tests.
module subroutine exec_generate_solve3D_addon_snapshots( self, cline )
    class(commander_generate_solve3D_addon_snapshots), intent(inout) :: self
    class(cmdline),                                        intent(inout) :: cline
    type(parameters) :: params
    type(sp_project) :: source, generated, snapshot
    type(image)      :: img
    type(string)     :: source_stk, stack_fname, project_fname, project_name
    integer, allocatable :: chunk_first(:), chunk_last(:), chunk_stkind(:), chunk_snapshot(:), snapshot_last(:)
    integer :: nptcls, nremaining, naddons, naddon_base, naddon_extra
    integer :: nchunks, nchunks_snapshot, numlen
    integer :: iptcl, ind_in_stk, stkind, previous_stkind, isnapshot, previous_snapshot
    integer :: ichunk, local_ind
    call cline%set('mkdir', 'no')
    call params%new(cline)
    call source%read(params%projfile)
    nptcls = source%os_ptcl3D%get_noris()
    if( source%os_ptcl2D%get_noris() /= nptcls ) THROW_HARD('ptcl2D and ptcl3D lengths differ')
    if( nptcls < 1 ) THROW_HARD('input project has no particles')
    if( params%nsnapshots < 2 ) THROW_HARD('nsnapshots must be at least 2')
    if( params%nptcls_base < 1 .or. params%nptcls_base >= nptcls ) THROW_HARD('nptcls_base must lie within the input particle range')
    nremaining = nptcls - params%nptcls_base
    naddons     = params%nsnapshots - 1
    if( nremaining < naddons ) THROW_HARD('each addon snapshot must contribute at least one particle')
    naddon_base  = nremaining / naddons
    naddon_extra = mod(nremaining, naddons)
    allocate(snapshot_last(params%nsnapshots))
    snapshot_last(1) = params%nptcls_base
    do isnapshot = 2, params%nsnapshots
        snapshot_last(isnapshot) = snapshot_last(isnapshot - 1) + naddon_base
        if( isnapshot - 1 <= naddon_extra ) snapshot_last(isnapshot) = snapshot_last(isnapshot) + 1
    enddo
    allocate(chunk_first(nptcls), chunk_last(nptcls), chunk_stkind(nptcls), chunk_snapshot(nptcls))
    nchunks          = 0
    previous_stkind = 0
    previous_snapshot = 0
    isnapshot = 1
    do iptcl = 1, nptcls
        if( iptcl > snapshot_last(isnapshot) ) isnapshot = isnapshot + 1
        call source%map_ptcl_ind2stk_ind('ptcl2D', iptcl, stkind, ind_in_stk)
        if( isnapshot /= previous_snapshot .or. stkind /= previous_stkind )then
            nchunks = nchunks + 1
            chunk_first(nchunks)    = iptcl
            chunk_stkind(nchunks)   = stkind
            chunk_snapshot(nchunks) = isnapshot
            if( nchunks > 1 ) chunk_last(nchunks - 1) = iptcl - 1
            previous_snapshot = isnapshot
            previous_stkind   = stkind
        endif
    enddo
    chunk_last(nchunks) = nptcls
    generated = source
    call generated%os_mic%kill
    call generated%os_stk%kill
    call generated%os_ptcl2D%kill
    call generated%os_ptcl3D%kill
    call generated%os_cls2D%kill
    call generated%os_cls3D%kill
    call generated%os_out%kill
    call generated%jobproc%kill
    if( generated%projinfo%isthere(1, 'sigma2_state') ) call generated%projinfo%delete_entry('sigma2_state')
    if( generated%projinfo%isthere(1, 'solve3D_manifest') ) call generated%projinfo%delete_entry('solve3D_manifest')
    if( generated%projinfo%isthere(1, 'solve3D_run_id') ) call generated%projinfo%delete_entry('solve3D_run_id')
    call generated%os_stk%new(nchunks, is_ptcl=.false.)
    call generated%os_ptcl2D%new(nptcls, is_ptcl=.true.)
    call generated%os_ptcl3D%new(nptcls, is_ptcl=.true.)
    call img%new([source%get_box(), source%get_box(), 1], source%get_smpd(), wthreads=.false.)
    numlen = len(int2str(nchunks))
    do ichunk = 1, nchunks
        stack_fname = simple_abspath(string('snapshot_stack'//int2str_pad(ichunk, numlen)//STK_EXT), check_exists=.false.)
        call generated%os_stk%transfer_ori(ichunk, source%os_stk, chunk_stkind(ichunk))
        call generated%os_stk%set(ichunk, 'stk',        stack_fname)
        call generated%os_stk%set(ichunk, 'fromp',      chunk_first(ichunk))
        call generated%os_stk%set(ichunk, 'top',        chunk_last(ichunk))
        call generated%os_stk%set(ichunk, 'nptcls',     chunk_last(ichunk) - chunk_first(ichunk) + 1)
        call generated%os_stk%set(ichunk, 'nptcls_stk', chunk_last(ichunk) - chunk_first(ichunk) + 1)
        local_ind = 0
        do iptcl = chunk_first(ichunk), chunk_last(ichunk)
            local_ind = local_ind + 1
            call source%get_stkname_and_ind('ptcl2D', iptcl, source_stk, ind_in_stk)
            call img%read(source_stk, ind_in_stk)
            call img%write(stack_fname, local_ind, del_if_exists=(local_ind == 1))
            call generated%os_ptcl2D%transfer_ori(iptcl, source%os_ptcl2D, iptcl)
            call generated%os_ptcl3D%transfer_ori(iptcl, source%os_ptcl3D, iptcl)
            call generated%os_ptcl2D%set(iptcl, 'stkind', ichunk)
            call generated%os_ptcl3D%set(iptcl, 'stkind', ichunk)
            call generated%os_ptcl2D%set(iptcl, 'indstk', local_ind)
            call generated%os_ptcl3D%set(iptcl, 'indstk', local_ind)
            call generated%os_ptcl2D%set(iptcl, 'pind', iptcl)
            call generated%os_ptcl3D%set(iptcl, 'pind', iptcl)
        enddo
    enddo
    call img%kill
    do isnapshot = 1, params%nsnapshots
        nchunks_snapshot = count(chunk_snapshot(1:nchunks) <= isnapshot)
        snapshot = generated
        snapshot%os_stk    = generated%os_stk%extract_subset(1, nchunks_snapshot)
        snapshot%os_ptcl2D = generated%os_ptcl2D%extract_subset(1, snapshot_last(isnapshot))
        snapshot%os_ptcl3D = generated%os_ptcl3D%extract_subset(1, snapshot_last(isnapshot))
        project_name  = 'snapshot'//int2str(isnapshot)
        project_fname = simple_abspath(project_name//METADATA_EXT, check_exists=.false.)
        call snapshot%projinfo%set(1, 'projname', project_name)
        call snapshot%projinfo%set(1, 'projfile', project_fname)
        call snapshot%write(project_fname)
        write(logfhandle,'(a,i0,a,i0,a,a)') '>>> SNAPSHOT ', isnapshot, ': ', snapshot_last(isnapshot), &
            &' particles; ', project_fname%to_char()
        call snapshot%kill
    enddo
    write(logfhandle,'(a,i0,a,i0,a,i0)') '>>> BASE PARTICLES: ', params%nptcls_base, &
        &'; ADDON PARTICLES: ', naddon_base, ' OR ', naddon_base + min(1, naddon_extra)
    call generated%kill
    call source%kill
end subroutine exec_generate_solve3D_addon_snapshots

!> solve3D_addon gate on symmetry-broken 6VXX particles: solve3D on a seeded selection of the
!  first NBASE rows, then the add-on on all NPTCLS rows (same project basename). Gates frozen-input
!  integrity, manifests, cohort coverage and poses, union map vs truth and base; metrics.tsv.
subroutine run_solve3D_addon_gate( nthr, all_ok )
    use simple_atoms,                   only: atoms
    use simple_molecule_data,           only: molecule_data, sars_cov2_spkgp_6vxx
    use simple_ui,                      only: make_ui
    use simple_commanders_solve3D,      only: commander_solve3D_addon
    use simple_solve3D_manifest,        only: solve3D_manifest
    use simple_sigma2_state_file,       only: sigma2_state_digest_file
    use simple_sigma2_files,            only: canonical_sigma2_consumable
    use simple_refine3D_fnames,         only: refine3D_state_vol_fname
    use simple_test_gate,               only: test_gate
    use simple_test_truth_metrics,      only: dock_both_hands, compare_to_truth, pair_pose_error, add_gaussian_blob
    use simple_solve3D_addon_report, only: solve3D_addon_report, ADDON_REPORT_FNAME
    use simple_image,                   only: image
    use, intrinsic :: iso_fortran_env, only: int64
    integer, intent(in)    :: nthr
    logical, intent(inout) :: all_ok
    character(len=*), parameter :: GATE_DIR    = 'solve3D_addon_gate'
    character(len=*), parameter :: TRUTH_VOL   = 'truth_6VXX_blob.mrc'
    character(len=*), parameter :: PTCL_STK    = 'simulated_particles.mrc'
    character(len=*), parameter :: STK_BASE    = 'particles_base.mrc'   !< the base set: particles 1-NBASE
    character(len=*), parameter :: STK_NEW     = 'particles_new.mrc'    !< the set appended after it
    character(len=*), parameter :: TRUTH_ORIS  = 'simulated_oris.txt'
    character(len=*), parameter :: PROJNAME    = 'addon_gate'
    character(len=*), parameter :: STRICT_DIR  = 'strict'
    real,    parameter :: SMPD        = 2.2
    integer, parameter :: BOX         = 112
    real,    parameter :: MSKDIAM     = 180.
    integer, parameter :: NPTCLS      = 3000   ! rows of the current project: both sets
    integer, parameter :: NBASE       = 2000   ! rows of the frozen project: the base set
    real,    parameter :: BASE_SELECT = 0.75   ! seeded selection of the base set (about 1500 frozen particles)
    real,    parameter :: SNR         = 0.2
    integer, parameter :: NSAMPLE     = 500    ! sampled in the base run and in the add-on (trailing from stage 5)
    real,    parameter :: GATE_LPSTART = 20.
    real,    parameter :: GATE_LPSTOP  = 8.
    integer, parameter :: NCLS2D      = 10     ! class labels for class-balanced sampling
    integer, parameter :: NPAIRS      = 20000
    integer, parameter :: GATE_SEED   = 20260926
    real,    parameter :: DOCK_HP     = 100.
    real,    parameter :: DOCK_LP     = 20.
    ! blob: amplitude relative to the map maximum, width and position in A
    real,    parameter :: BLOB_AMP    = 1.5, BLOB_SIGMA = 11., BLOB_POS(3) = [45., 25., 30.]
    ! floors: worst observed run plus a margin (~30% pose median, 5 deg excess, 0.05/0.03 corr, 1 A FSC)
    real,    parameter :: MAX_COHORT_POSE_ERR  = 15.  !< cohort-frozen pair median; random poses give ~40
    real,    parameter :: MAX_POSE_ERR_EXCESS  = 5.   !< cohort-frozen above frozen-frozen
    real,    parameter :: MIN_COVERAGE        = 0.9   !< cohort particles with updatecnt > 0
    real,    parameter :: MIN_UNION_CORR      = 0.9   !< docked union map vs truth
    real,    parameter :: MAX_CORR_LOSS       = 0.03  !< union map correlation below the base map's
    real,    parameter :: MAX_FSC_LOSS        = 1.0   !< union FSC=0.143 above the base map's (A)
    real,    parameter :: MIN_REPORT_CORR     = 0.9   !< report: union vs base map up to the base FSC=0.143 resolution
    type(commander_solve3D)           :: xsolve3D
    type(commander_solve3D_addon)     :: xaddon
    type(solve3D_addon_report)        :: addon_report
    type(commander_simulate_particles):: xsim
    type(commander_new_project)       :: xnew_project
    type(cmdline)       :: cl
    type(atoms)         :: molecule
    type(molecule_data) :: mol
    type(sp_project)    :: spproj, strict_proj, out_proj, frz_proj
    type(oris)          :: truth
    type(ctfparams)     :: ctfvars
    type(solve3D_manifest) :: man_base, man_out, man_pub
    type(sp_project)          :: pub_proj
    type(string)        :: root, stk_abs, stk_base_abs, stk_new_abs, full_proj, strict_proj_fname
    type(string)        :: frozen_run_proj, out_run_proj
    type(image)         :: img
    type(string)        :: frozen_sigma, frozen_vol, addon_vol, base_vol, cwd_here, truth_abs
    character(len=STDLEN) :: msg
    integer(int64)      :: dig_proj0
    integer, allocatable :: frozen_inds(:), cohort_inds(:)
    logical, allocatable :: l_frozen(:)
    type(test_gate)     :: gate
    integer :: i, status, nf, nc, nposed, box_vol
    real    :: r, smpd_vol, corr_addon, corr_base, fsc_addon, fsc_base, err_cf, err_ff, coverage
    real    :: frac_ff, frac_cf
    logical :: found, l_same, l_pub, l_frz_row
    call make_ui
    write(logfhandle,'(a)') '>>> TEST_SOLVE3D_ADDON: solve3D_addon gate'
    call simple_getcwd(root)
    if( file_exists(GATE_DIR) )then
        call simple_rmdir(GATE_DIR, status)
        if( status /= 0 ) THROW_HARD('Could not reset '//GATE_DIR)
    endif
    call simple_mkdir(GATE_DIR)
    call simple_chdir(GATE_DIR, status)
    if( status /= 0 ) THROW_HARD('Could not enter '//GATE_DIR)
    call gate%new(string('metrics.tsv'))
    ! ---- the truth: 6VXX with an off-axis blob ----
    mol = sars_cov2_spkgp_6vxx()
    call molecule%pdb2mrc(smpd=SMPD, volfile=string(TRUTH_VOL), mol=mol, center_pdb=.true., vol_dim=[BOX,BOX,BOX])
    call molecule%kill()
    call add_gaussian_blob(string(TRUTH_VOL), BLOB_POS, BLOB_SIGMA, BLOB_AMP)
    truth_abs = simple_abspath(string(TRUTH_VOL))
    ! ---- particles ----
    call cl%set('prg',     'simulate_particles')
    call cl%set('vol1',    TRUTH_VOL)
    call cl%set('smpd',    SMPD)
    call cl%set('mskdiam', MSKDIAM)
    call cl%set('nthr',    nthr)
    call cl%set('nptcls',  NPTCLS)
    call cl%set('pgrp',    'c1')
    call cl%set('snr',     SNR)
    call cl%set('ctf',     'yes')
    call cl%set('sherr',   0.0)
    call xsim%execute(cl)
    call cl%kill
    stk_abs = simple_abspath(string(PTCL_STK))
    ! the base set and the set appended after it, each its own stack
    call img%new([BOX,BOX,1], SMPD)
    do i = 1, NPTCLS
        call img%read(stk_abs, i)
        if( i <= NBASE )then
            call img%write(string(STK_BASE), i)
        else
            call img%write(string(STK_NEW), i - NBASE)
        endif
    enddo
    call img%kill
    stk_base_abs = simple_abspath(string(STK_BASE))
    stk_new_abs  = simple_abspath(string(STK_NEW))
    call truth%new(NPTCLS, is_ptcl=.true.)
    call truth%read(string(TRUTH_ORIS), [1,NPTCLS])
    ! ---- the current project: both sets ----
    call cl%set('projname',  PROJNAME)
    call cl%set('qsys_name', 'local')
    call xnew_project%execute(cl)             ! creates and enters PROJNAME/
    call cl%kill
    full_proj = simple_abspath(string(PROJNAME//'.simple'))
    call spproj%read(full_proj)
    ctfvars%smpd    = SMPD
    ctfvars%kv      = 300.
    ctfvars%cs      = 2.7
    ctfvars%fraca   = 0.1
    ctfvars%ctfflag = CTFFLAG_YES
    call spproj%add_stk(stk_base_abs, ctfvars)
    call spproj%add_stk(stk_new_abs,  ctfvars)
    call set_fixed_seed(GATE_SEED)
    call spproj%os_cls2D%new(NCLS2D, is_ptcl=.false.)
    call spproj%os_cls2D%set_all2single('state', 1.)
    do i = 1, NPTCLS
        call spproj%os_ptcl3D%set(i, 'dfx',    truth%get(i, 'dfx'))
        call spproj%os_ptcl3D%set(i, 'dfy',    truth%get(i, 'dfy'))
        call spproj%os_ptcl3D%set(i, 'angast', truth%get(i, 'angast'))
        call random_number(r)
        ! a 2D classification stand-in: the 3D workflows need a searched ptcl2D
        ! field and class labels for class-balanced sampling, nothing else
        call spproj%os_ptcl2D%set(i, 'class', 1 + int(r*real(NCLS2D)))
        call spproj%os_ptcl2D%set(i, 'corr',  0.5)
        call spproj%os_ptcl2D%set(i, 'dfx',    truth%get(i, 'dfx'))
        call spproj%os_ptcl2D%set(i, 'dfy',    truth%get(i, 'dfy'))
        call spproj%os_ptcl2D%set(i, 'angast', truth%get(i, 'angast'))
        call spproj%os_ptcl2D%set_state(i, 1)
        call spproj%os_ptcl3D%set_state(i, 1)
    enddo
    call spproj%write(full_proj)
    ! ---- the frozen project: the base set alone (its import gives the current
    ! project's first stack and rows), a seeded selection, same basename, own
    ! directory ----
    allocate(l_frozen(NPTCLS), source=.false.)
    do i = 1, NBASE
        call random_number(r)
        l_frozen(i) = r < BASE_SELECT
    enddo
    strict_proj = spproj
    strict_proj%os_stk    = spproj%os_stk%extract_subset(1, 1)
    strict_proj%os_ptcl2D = spproj%os_ptcl2D%extract_subset(1, NBASE)
    strict_proj%os_ptcl3D = spproj%os_ptcl3D%extract_subset(1, NBASE)
    do i = 1, NBASE
        if( .not. l_frozen(i) )then
            call strict_proj%os_ptcl2D%set_state(i, 0)
            call strict_proj%os_ptcl3D%set_state(i, 0)
        endif
    enddo
    call simple_mkdir(STRICT_DIR)
    strict_proj_fname = simple_abspath(string(STRICT_DIR//'/'//PROJNAME//'.simple'), check_exists=.false.)
    call strict_proj%write(strict_proj_fname)
    call strict_proj%kill
    call spproj%kill
    ! ---- the base run on the frozen project ----
    call simple_getcwd(cwd_here)
    call simple_chdir(string(STRICT_DIR), status)
    call cl%set('prg',            'solve3D')
    call cl%set('projfile',       strict_proj_fname)
    call cl%set('mkdir',          'yes')
    call cl%set('pgrp',           'c1')
    call cl%set('mskdiam',        MSKDIAM)
    call cl%set('nthr',           nthr)
    call cl%set('nsample',        NSAMPLE)
    call cl%set('force_lp_range', 'yes')
    call cl%set('lpstart',        GATE_LPSTART)
    call cl%set('lpstop',         GATE_LPSTOP)
    call xsolve3D%execute(cl)
    call cl%kill
    call simple_getcwd(frozen_run_proj)
    frozen_run_proj = frozen_run_proj//'/'//PROJNAME//'.simple'
    call simple_chdir(cwd_here, status)
    ! the frozen inputs, before the add-on
    call frz_proj%read(frozen_run_proj)
    call man_base%read_registered(frz_proj, frozen_run_proj, status, msg)
    call gate%check('base_manifest_registered', status == 0)
    call man_base%validate_frozen(frz_proj, status, msg)
    call gate%check('base_manifest_valid_frozen_input', status == 0)
    call man_base%get_artifact('sigma2_state', 0, frozen_sigma, found)
    call frz_proj%get_vol('vol', 1, frozen_vol, smpd_vol, box_vol)
    base_vol   = frozen_vol
    dig_proj0  = sigma2_state_digest_file(frozen_run_proj)
    ! ---- the add-on on the current project ----
    call cl%set('prg',             'solve3D_addon')
    call cl%set('projfile',        full_proj)
    call cl%set('projfile_frozen', frozen_run_proj)
    call cl%set('nthr',            nthr)
    call cl%set('addon_diag',      'yes')
    call xaddon%execute(cl)
    call cl%kill
    call simple_getcwd(out_run_proj)
    out_run_proj = out_run_proj//'/'//PROJNAME//'.simple'
    call simple_chdir(cwd_here, status)
    ! ---- provenance and isolation ----
    call gate%check('frozen_project_unchanged', sigma2_state_digest_file(frozen_run_proj) == dig_proj0)
    call gate%check('frozen_sigma2_unchanged',  man_base%matches_artifact('sigma2_state', 0, frozen_sigma))
    call gate%check('frozen_map_unchanged',     man_base%matches_artifact('vol', 1, frozen_vol))
    ! the frozen term was weighted by the base run's committed residual sigma2
    ! state (its copy in the add-on run is byte-equal after every accumulation)
    call gate%check('frozen_sigma2_consumed_as_committed', &
        &man_base%matches_artifact('sigma2_state', 0, string('1_solve3D_addon/frozen/frozen_sigma2_state.bin')))
    call man_base%kill
    call out_proj%read(out_run_proj)
    ! the final reconstruction bootstrapped the union's sigma2 state over
    ! every particle at native sampling, so the output is the frozen input of
    ! a next add-on
    call gate%check('output_registers_the_union_sigma2_state', &
        &canonical_sigma2_consumable(out_proj, out_proj%os_ptcl3D, BOX, SMPD, .true., msg))
    call man_out%read_registered(out_proj, out_run_proj, status, msg)
    call gate%check('addon_manifest_registered', status == 0)
    call man_out%validate_frozen(out_proj, status, msg)
    call gate%check('addon_manifest_valid_frozen_input', status == 0)
    ! all done: the finished project replaced the original current project
    ! file, registering the add-on's manifest by absolute path
    call pub_proj%read(full_proj)
    l_pub = pub_proj%os_ptcl3D%get_noris() == out_proj%os_ptcl3D%get_noris()
    if( l_pub )then
        do i = 1, out_proj%os_ptcl3D%get_noris()
            l_pub = l_pub .and. pub_proj%os_ptcl3D%get_state(i) == out_proj%os_ptcl3D%get_state(i) .and. &
                &all(abs(pub_proj%os_ptcl3D%get_euler(i) - out_proj%os_ptcl3D%get_euler(i)) < 1.e-4)
        enddo
    endif
    call gate%check('original_project_replaced_by_the_output', l_pub)
    call man_pub%read_registered(pub_proj, full_proj, status, msg)
    call gate%check('published_project_registers_the_addon_manifest', status == 0 .and. &
        &man_pub%get_run_id() == man_out%get_run_id())
    call man_pub%kill
    call pub_proj%kill
    call man_out%kill
    call gate%check('addon_diag_map_written', file_exists(string('1_solve3D_addon/addon_diag/')// &
        &refine3D_state_vol_fname(1)))
    ! ---- the validation report against the base solution ----
    if( file_exists(string('1_solve3D_addon/'//ADDON_REPORT_FNAME)) )then
        call addon_report%read(string('1_solve3D_addon/'//ADDON_REPORT_FNAME))
        call gate%check('addon_report_no_regression', .not. addon_report%any_regressed())
        call gate%metric('addon_report_union_base_corr', addon_report%get_corr(1), MIN_REPORT_CORR, &
            &addon_report%get_corr(1) >= MIN_REPORT_CORR)
        call gate%report('addon_report_dshell0143', real(addon_report%get_dshell(1)))
        call gate%report('addon_report_cohort_base_fsc0143_A', addon_report%get_cohort_res0143(1))
        call addon_report%kill
    else
        call gate%check('addon_report_written', .false.)
    endif
    ! frozen rows identical to the frozen project; cohort rows (the appended
    ! set included) posed
    nf = count(l_frozen)
    nc = NPTCLS - nf
    allocate(frozen_inds(0), cohort_inds(0))
    l_same = .true.
    nposed = 0
    do i = 1, NPTCLS
        l_frz_row = i <= NBASE
        if( l_frz_row ) l_frz_row = frz_proj%os_ptcl3D%get_state(i) > 0 .and. frz_proj%os_ptcl3D%get_updatecnt(i) > 0
        if( l_frz_row )then
            frozen_inds = [frozen_inds, i]
            l_same = l_same .and. all(abs(out_proj%os_ptcl3D%get_euler(i) - frz_proj%os_ptcl3D%get_euler(i)) < 1.e-3) &
                &.and. all(abs(out_proj%os_ptcl3D%get_2Dshift(i) - frz_proj%os_ptcl3D%get_2Dshift(i)) < 1.e-4) &
                &.and. out_proj%os_ptcl3D%get_state(i) == frz_proj%os_ptcl3D%get_state(i) &
                &.and. out_proj%os_ptcl3D%get_eo(i) == frz_proj%os_ptcl3D%get_eo(i) &
                &.and. out_proj%os_ptcl3D%get_updatecnt(i) == frz_proj%os_ptcl3D%get_updatecnt(i)
        else
            cohort_inds = [cohort_inds, i]
            if( out_proj%os_ptcl3D%get_state(i) > 0 .and. out_proj%os_ptcl3D%get_updatecnt(i) > 0 ) nposed = nposed + 1
        endif
    enddo
    call gate%check('frozen_rows_restored_exactly', l_same)
    call gate%check('every_particle_active', out_proj%count_state_gt_zero() == NPTCLS)
    coverage = real(nposed) / real(max(1,size(cohort_inds)))
    call gate%metric('cohort_coverage', coverage, MIN_COVERAGE, coverage >= MIN_COVERAGE)
    ! ---- poses against the truth ----
    err_ff = pair_pose_error(out_proj%os_ptcl3D, truth, frozen_inds, frozen_inds, NPAIRS, GATE_SEED + 1, frac_ff)
    err_cf = pair_pose_error(out_proj%os_ptcl3D, truth, cohort_inds, frozen_inds, NPAIRS, GATE_SEED + 1, frac_cf)
    call gate%report('frozen_frozen_pairs_within_5deg', frac_ff)
    call gate%report('cohort_frozen_pairs_within_5deg', frac_cf)
    call gate%report('frozen_frozen_pair_pose_error_deg', err_ff)
    call gate%metric('cohort_frozen_pair_pose_error_deg', err_cf, MAX_COHORT_POSE_ERR, err_cf <= MAX_COHORT_POSE_ERR)
    call gate%metric('cohort_pose_error_excess_deg', err_cf - err_ff, MAX_POSE_ERR_EXCESS, &
        &err_cf - err_ff <= MAX_POSE_ERR_EXCESS)
    ! ---- maps against the truth ----
    call out_proj%get_vol('vol', 1, addon_vol, smpd_vol, box_vol)
    call dock_and_compare(addon_vol, 'union', corr_addon, fsc_addon)
    call dock_and_compare(base_vol,  'base',  corr_base,  fsc_base)
    call gate%metric('union_map_truth_corr', corr_addon, MIN_UNION_CORR, corr_addon >= MIN_UNION_CORR)
    call gate%report('base_map_truth_corr', corr_base)
    call gate%metric('union_minus_base_corr', corr_addon - corr_base, -MAX_CORR_LOSS, &
        &corr_addon - corr_base >= -MAX_CORR_LOSS)
    call gate%metric('union_truth_fsc0143_A', fsc_addon, fsc_base + MAX_FSC_LOSS, &
        &fsc_addon > 0. .and. fsc_addon <= fsc_base + MAX_FSC_LOSS)
    call gate%report('base_truth_fsc0143_A', fsc_base)
    all_ok = all_ok .and. gate%passed()
    call gate%kill
    call out_proj%kill
    call frz_proj%kill
    call truth%kill
    call simple_chdir(root, status)

contains

    !> dock a map onto the truth in both hands and score it there; a missing
    !! map scores correlation 0 and resolution -1 (fails its floors)
    subroutine dock_and_compare( fname, tag, corr, fsc0143 )
        class(string),    intent(in)  :: fname
        character(len=*), intent(in)  :: tag
        real,             intent(out) :: corr, fsc0143
        real :: cc_direct, cc_mirror, fsc05
        corr    = 0.
        fsc0143 = -1.
        if( .not. file_exists(fname) ) return
        call dock_both_hands(truth_abs, fname, MSKDIAM, DOCK_HP, DOCK_LP, 'gate_'//tag, &
            &string('gate_'//tag//'_docked.mrc'), cc_direct, cc_mirror)
        call compare_to_truth(truth_abs, string('gate_'//tag//'_docked.mrc'), MSKDIAM, corr, fsc05, fsc0143)
    end subroutine dock_and_compare

end subroutine run_solve3D_addon_gate

end submodule simple_commanders_test_highlevel_addon
