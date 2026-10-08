!@descr: FLEX PCA workflow gates: shared memory against distributed, and the per-state refinement
submodule (simple_commanders_test_highlevel) simple_commanders_test_highlevel_flex
use, intrinsic :: iso_fortran_env,      only: real32, real64
use simple_commanders_flex_pca,         only: commander_flex_pca
use simple_commanders_refine3D,         only: commander_refine3D_auto
use simple_state_weight_set,            only: state_weight_set
use simple_sigma2_state_file,           only: sigma2_state_digest_file
use simple_flex_pca_application_tester, only: create_flex_pca_phantom_fixture, &
    &build_flex_pca_phantom_project, configure_flex_pca_phantom, read_state_weight_table, FLEX_PHANTOM_BOX, &
    &FLEX_PHANTOM_NPTCLS, FLEX_PHANTOM_NCOMP, FLEX_PHANTOM_MAX_STATES, &
    &FLEX_PHANTOM_SMPD, FLEX_PHANTOM_MSKDIAM, FLEX_PHANTOM_LABEL_ACCURACY_MIN
use simple_flex_pca_embedding_io, only: read_deconv_block
use simple_image,                 only: image
use simple_oris,                  only: oris
use simple_sp_project,            only: sp_project
use simple_test_gate,             only: test_gate
use simple_test_truth_metrics,    only: compare_to_truth
use simple_ui,                    only: make_ui
implicit none
#include "simple_local_flags.inc"

integer, parameter :: FLEX_GATE_SEED = 202610063
!> the per-state refinement: its seed, its iteration budget and how much the weighted map may trail the
!! hard one in truth correlation (plan section 9); it runs with twice the FLEX modes' threads, the
!! processors this serial entry reserves
integer, parameter :: STATE_REFINE_SEED   = 202610071
integer, parameter :: STATE_REFINE_MAXITS = 5
real,    parameter :: STATE_REFINE_MARGIN = 0.01
real,    parameter :: MAP_CORR_LP = 8.0

type :: flex_gate_result
    logical :: valid = .false.
    logical :: embedding_found = .false.
    integer :: nstates = 0
    integer :: matched_a = 0
    integer :: matched_b = 0
    real    :: label_accuracy = 0.0
    real    :: pc_corr = 0.0
    real(real64) :: noise_scale = 1.0_real64
    integer, allocatable :: predicted_truth(:)
    real(real64), allocatable :: z(:,:)
    type(string) :: map_a
    type(string) :: map_b
end type flex_gate_result

contains

module subroutine exec_test_flex_pca_blobs( self, cline )
    class(commander_test_flex_pca_blobs), intent(inout) :: self
    class(cmdline),                        intent(inout) :: cline
    class(parameters), allocatable :: params
    type(string) :: cwd_saved, fixture_root
    integer :: nthr, status
    logical :: all_ok
    call simple_getcwd(cwd_saved)
    fixture_root = filepath(cwd_saved, 'test_flex_pca_blobs_'//int2str(get_process_id()))
    if( dir_exists(fixture_root) ) call simple_rmdir(fixture_root)
    call simple_mkdir(fixture_root)
    call simple_chdir(fixture_root, status)
    if( status /= 0 ) THROW_HARD('FLEX_PCA_BLOBS FAILED: could not enter fixture directory')
    allocate(params)
    call params%new(cline)
    nthr = params%nthr
    deallocate(params)
    all_ok = .true.
    call run_flex_pca_blobs_gate(nthr, all_ok)
    call simple_chdir(cwd_saved, status)
    if( status /= 0 ) THROW_HARD('FLEX_PCA_BLOBS FAILED: could not restore original directory')
    CWD_GLOB = cwd_saved%to_char()
    if( all_ok )then
        call simple_rmdir(fixture_root)
        write(logfhandle,'(A)') 'PASS: FLEX PCA shared/distributed workflow agrees with simulation truth'
        call simple_end('**** SIMPLE_TEST_FLEX_PCA_BLOBS NORMAL STOP ****')
    else
        THROW_HARD('FLEX_PCA_BLOBS FAILED')
    endif
end subroutine exec_test_flex_pca_blobs

subroutine run_flex_pca_blobs_gate( nthr, all_ok )
    integer, intent(in)    :: nthr
    logical, intent(inout) :: all_ok
    type(test_gate) :: gate
    type(flex_gate_result) :: shared, distributed
    type(string) :: root, truth_a, truth_b, truth_mean, truth_diff, oritab, stack
    real :: latent_corr, truth_label_agreement, map_a_corr, map_b_corr
    call make_ui
    call simple_getcwd(root)
    CWD_GLOB = root%to_char()
    call gate%new(string('metrics.tsv'))
    call create_flex_pca_phantom_fixture(root, nthr, truth_a, truth_b, truth_mean, truth_diff, oritab, stack)
    call run_flex_mode('shared', 1, nthr, root, truth_a, truth_b, truth_mean, truth_diff, &
        &oritab, stack, gate, shared)
    call run_flex_mode('distributed', 2, nthr, root, truth_a, truth_b, truth_mean, truth_diff, &
        &oritab, stack, gate, distributed)
    call gate%check('shared and distributed runs publish valid result sets', shared%valid .and. distributed%valid)
    if( shared%valid .and. distributed%valid )then
        call gate%report('published_state_count_difference', real(abs(shared%nstates - distributed%nstates)))
        truth_label_agreement = real(count(shared%predicted_truth == distributed%predicted_truth)) / &
            &real(FLEX_PHANTOM_NPTCLS)
        call gate%report('shared_distributed_truth_label_agreement', truth_label_agreement)
        if( shared%embedding_found .and. distributed%embedding_found )then
            latent_corr = best_latent_correlation(shared%z, distributed%z)
            call gate%report('shared_distributed_latent_correlation', latent_corr)
        endif
        call map_pair_correlation(shared%map_a, distributed%map_a, map_a_corr)
        call map_pair_correlation(shared%map_b, distributed%map_b, map_b_corr)
        call gate%report('shared_distributed_mode_a_map_correlation', map_a_corr)
        call gate%report('shared_distributed_mode_b_map_correlation', map_b_corr)
    endif
    ! the per-state refinement of the shared run's states, after the FLEX comparisons it does not touch
    if( shared%valid ) call run_state_refinements(root, 2*nthr, truth_a, truth_b, shared, gate)
    all_ok = all_ok .and. gate%passed()
    call gate%kill
    call kill_result(shared)
    call kill_result(distributed)
end subroutine run_flex_pca_blobs_gate

subroutine run_flex_mode( tag, nparts, nthr, root, truth_a, truth_b, truth_mean, truth_diff, &
    &oritab, stack, gate, result )
    character(len=*), intent(in) :: tag
    integer,          intent(in) :: nparts, nthr
    class(string),    intent(in) :: root, truth_a, truth_b, truth_mean, truth_diff, oritab, stack
    type(test_gate),  intent(inout) :: gate
    type(flex_gate_result), intent(out) :: result
    class(commander_flex_pca), allocatable :: flex_pca
    class(cmdline), allocatable :: cline
    type(string) :: run_dir
    integer :: status
    run_dir = root//'/'//trim(tag)
    call simple_mkdir(run_dir)
    call build_flex_pca_phantom_project(run_dir//'/phantom.simple', stack, oritab)
    call simple_chdir(run_dir, status)
    if( status /= 0 ) THROW_HARD('FLEX_PCA_BLOBS FAILED: could not enter run directory')
    CWD_GLOB = run_dir%to_char()
    allocate(flex_pca, cline)
    call configure_flex_pca_phantom(cline, truth_mean, 'pcg', 'pcg', nparts, nthr)
    call set_fixed_seed(FLEX_GATE_SEED, propagate=.true.)
    call flex_pca%execute(cline)
    call cline%kill
    deallocate(cline, flex_pca)
    call inspect_flex_result(tag, truth_a, truth_b, truth_diff, gate, result)
    call simple_chdir(root, status)
    if( status /= 0 ) THROW_HARD('FLEX_PCA_BLOBS FAILED: could not restore fixture directory')
    CWD_GLOB = root%to_char()
end subroutine run_flex_mode

subroutine inspect_flex_result( tag, truth_a, truth_b, truth_diff, gate, result )
    character(len=*), intent(in) :: tag
    class(string),    intent(in) :: truth_a, truth_b, truth_diff
    type(test_gate),  intent(inout) :: gate
    type(flex_gate_result), intent(out) :: result
    class(sp_project), allocatable :: project
    class(oris), allocatable :: field
    real(real32), allocatable :: weights(:,:)
    real(real64), allocatable :: precision(:,:,:)
    integer, allocatable :: labels(:), deconv_labels(:), truth_counts(:,:), state_truth(:)
    real, allocatable :: corr(:,:)
    type(string) :: map
    character(len=STDLEN) :: message
    integer :: status, i, state, truth_label, predicted, nlabel_mismatch, nweight_mismatch, map_box
    real :: best_score, fsc05, fsc0143, map_smpd
    result = flex_gate_result()
    allocate(project, field)
    call project%read(string('phantom.simple'))
    field = project%os_ptcl3D
    call read_state_weight_table(project, field, result%nstates, weights, labels, status, message)
    call gate%check(trim(tag)//': delivered flex-weight set validates', status == 0)
    if( status /= 0 )then
        write(logfhandle,'(A,A)') '>>> FLEX_PCA_BLOBS weight-store failure: ', trim(message)
        call field%kill
        call project%kill
        deallocate(field, project)
        return
    endif
    call gate%check(trim(tag)//': publication retains both truth modes', result%nstates >= 2)
    call gate%check(trim(tag)//': publication respects the state ceiling', &
        &result%nstates <= FLEX_PHANTOM_MAX_STATES)
    allocate(truth_counts(result%nstates,2), source=0)
    allocate(state_truth(result%nstates), source=0)
    allocate(result%predicted_truth(FLEX_PHANTOM_NPTCLS), source=0)
    nlabel_mismatch = 0
    nweight_mismatch = 0
    do i = 1, FLEX_PHANTOM_NPTCLS
        truth_label = merge(1, 2, mod(i,2) == 1)
        if( labels(i) >= 1 .and. labels(i) <= result%nstates ) &
            &truth_counts(labels(i),truth_label) = truth_counts(labels(i),truth_label) + 1
        if( int(labels(i)) /= project%os_ptcl3D%get_state(i) ) nlabel_mismatch = nlabel_mismatch + 1
        if( labels(i) > 0 )then
            predicted = maxloc(weights(:,i), dim=1)
            if( int(labels(i)) /= predicted ) nweight_mismatch = nweight_mismatch + 1
        endif
    enddo
    do state = 1, result%nstates
        state_truth(state) = maxloc(truth_counts(state,:), dim=1)
    enddo
    do i = 1, FLEX_PHANTOM_NPTCLS
        if( labels(i) >= 1 .and. labels(i) <= result%nstates ) &
            &result%predicted_truth(i) = state_truth(labels(i))
    enddo
    result%label_accuracy = real(sum(maxval(truth_counts,dim=2))) / real(FLEX_PHANTOM_NPTCLS)
    call gate%metric(trim(tag)//': truth_label_accuracy', result%label_accuracy, &
        &FLEX_PHANTOM_LABEL_ACCURACY_MIN, result%label_accuracy >= FLEX_PHANTOM_LABEL_ACCURACY_MIN)
    call gate%check(trim(tag)//': weight flags and project labels agree', nlabel_mismatch == 0)
    call gate%check(trim(tag)//': hard labels select maximum kernel weight', nweight_mismatch == 0)
    call inspect_eigenvolume(tag, truth_diff, gate, result%pc_corr)
    allocate(corr(result%nstates,2), source=0.0)
    do state = 1, result%nstates
        call project%get_vol('vol', state, map, map_smpd, map_box)
        call gate%check(trim(tag)//': delivered state map '//int2str(state)//' exists', file_exists(map))
        if( file_exists(map) )then
            call compare_to_truth(truth_a, map, FLEX_PHANTOM_MSKDIAM, corr(state,1), fsc05, fsc0143, &
                &corr_lp=MAP_CORR_LP)
            call compare_to_truth(truth_b, map, FLEX_PHANTOM_MSKDIAM, corr(state,2), fsc05, fsc0143, &
                &corr_lp=MAP_CORR_LP)
        endif
        call map%kill
    enddo
    call match_truth_maps(corr, result%matched_a, result%matched_b, best_score)
    call gate%check(trim(tag)//': two distinct maps match the two truth modes', &
        &result%matched_a > 0 .and. result%matched_b > 0)
    if( result%matched_a > 0 .and. result%matched_b > 0 )then
        call gate%check(trim(tag)//': mode-A map prefers truth A', &
            &corr(result%matched_a,1) > corr(result%matched_a,2))
        call gate%check(trim(tag)//': mode-B map prefers truth B', &
            &corr(result%matched_b,2) > corr(result%matched_b,1))
        call project%get_vol('vol', result%matched_a, map, map_smpd, map_box)
        result%map_a = simple_abspath(map)
        call map%kill
        call project%get_vol('vol', result%matched_b, map, map_smpd, map_box)
        result%map_b = simple_abspath(map)
        call map%kill
    endif
    write(logfhandle,'(A,A,A)') '>>> FLEX_PCA_BLOBS ', trim(tag), &
        &' state truth counts/correlations: state nA nB corrA corrB'
    do state = 1, result%nstates
        write(logfhandle,'(A,I0,2(1X,I0),2(1X,F9.6))') '>>>   ', state, truth_counts(state,:), corr(state,:)
    enddo
    call read_deconv_block('flex_pca_embedding.bin', FLEX_PHANTOM_NPTCLS, FLEX_PHANTOM_NCOMP, &
        &result%z, precision, deconv_labels, result%noise_scale, result%embedding_found)
    call gate%check(trim(tag)//': deconvolved embedding is published', result%embedding_found)
    if( result%embedding_found ) call gate%report(trim(tag)//': deconvolution_noise_scale', real(result%noise_scale))
    result%valid = result%nstates >= 2 .and. result%nstates <= FLEX_PHANTOM_MAX_STATES .and. &
        &result%matched_a > 0 .and. result%matched_b > 0
    call field%kill
    call project%kill
    deallocate(field, project)
    deallocate(weights, labels, truth_counts, state_truth, corr)
    if( allocated(precision) ) deallocate(precision)
    if( allocated(deconv_labels) ) deallocate(deconv_labels)
end subroutine inspect_flex_result

subroutine inspect_eigenvolume( tag, truth_diff, gate, corr )
    character(len=*), intent(in) :: tag
    class(string),    intent(in) :: truth_diff
    type(test_gate),  intent(inout) :: gate
    real,             intent(out) :: corr
    class(image), allocatable :: truth, pc
    type(string) :: pcfile
    pcfile = 'flex_pca_polished_pc001.mrc'
    corr = 0.0
    call gate%check(trim(tag)//': polished leading eigenvolume is delivered', file_exists(pcfile))
    if( .not. file_exists(pcfile) ) return
    allocate(truth, pc)
    call truth%new([FLEX_PHANTOM_BOX,FLEX_PHANTOM_BOX,FLEX_PHANTOM_BOX], FLEX_PHANTOM_SMPD, wthreads=.false.)
    call pc%new([FLEX_PHANTOM_BOX,FLEX_PHANTOM_BOX,FLEX_PHANTOM_BOX], FLEX_PHANTOM_SMPD, wthreads=.false.)
    call truth%read(truth_diff)
    call pc%read(pcfile)
    corr = abs(truth%real_corr(pc))
    call gate%report(trim(tag)//': leading_eigenvolume_truth_correlation', corr)
    call truth%kill
    call pc%kill
    deallocate(truth, pc)
end subroutine inspect_eigenvolume

!> The per-state refinement (plan section 9), on the shared-memory run's published project: each
!! truth-matched state is refined by refine3D_auto state=X twice from that same parent, with the frozen
!! FLEX weights and with the hard labels. The weighted map must correlate with its truth (same frame,
!! low-pass 8 A) at least as well as the hard one, minus STATE_REFINE_MARGIN; the parent project must
!! stay byte-identical.
subroutine run_state_refinements( root, nthr, truth_a, truth_b, shared, gate )
    class(string),          intent(in)    :: root, truth_a, truth_b
    integer,                intent(in)    :: nthr
    type(flex_gate_result), intent(in)    :: shared
    type(test_gate),        intent(inout) :: gate
    type(string)   :: run_dir, parent
    integer(int64) :: checksum_before
    real    :: corr_flex, corr_hard
    integer :: status, k, state
    run_dir = root//'/shared'
    parent  = run_dir//'/phantom.simple'
    checksum_before = sigma2_state_digest_file(parent)
    do k = 1, 2
        state = merge(shared%matched_a, shared%matched_b, k == 1)
        if( k == 1 )then
            call refine_flex_state(run_dir, parent, state, 'flex', nthr, truth_a, gate, corr_flex)
            call refine_flex_state(run_dir, parent, state, 'no',   nthr, truth_a, gate, corr_hard)
        else
            call refine_flex_state(run_dir, parent, state, 'flex', nthr, truth_b, gate, corr_flex)
            call refine_flex_state(run_dir, parent, state, 'no',   nthr, truth_b, gate, corr_hard)
        endif
        call gate%report('state'//int2str(state)//'_hard_refined_truth_correlation', corr_hard)
        call gate%metric('state'//int2str(state)//'_weighted_refined_truth_correlation', corr_flex, &
            &corr_hard - STATE_REFINE_MARGIN, corr_flex >= corr_hard - STATE_REFINE_MARGIN)
    enddo
    call gate%check('the parent project is byte-identical after the per-state refinements', &
        &sigma2_state_digest_file(parent) == checksum_before)
    call simple_chdir(root, status)
    if( status /= 0 ) THROW_HARD('FLEX_PCA_BLOBS FAILED: could not restore fixture directory')
    CWD_GLOB = root%to_char()
end subroutine run_state_refinements

!> refine3D_auto state=X m_estimator=mode in a new run directory of run_dir, from projfile; corr is the
!! refined map's correlation with its truth
subroutine refine_flex_state( run_dir, projfile, state, mode, nthr, truth, gate, corr )
    class(string),    intent(in)    :: run_dir, projfile, truth
    integer,          intent(in)    :: state, nthr
    character(len=*), intent(in)    :: mode
    type(test_gate),  intent(inout) :: gate
    real,             intent(out)   :: corr
    class(commander_refine3D_auto), allocatable :: xrefine
    class(cmdline),                 allocatable :: cline
    class(sp_project),              allocatable :: work
    type(state_weight_set) :: wset
    type(string)           :: work_dir, work_projfile, map, imgkind
    character(len=STDLEN)  :: tag, message
    real    :: map_smpd, fsc05, fsc0143
    integer :: map_box, status, i, ind
    logical :: l_registered
    corr = 0.
    tag  = 'state'//int2str(state)//'_'//trim(mode)
    call simple_chdir(run_dir, status)
    if( status /= 0 ) THROW_HARD('FLEX_STATE_REFINE FAILED: could not enter run directory')
    CWD_GLOB = run_dir%to_char()
    allocate(xrefine, cline, work)
    call cline%set('prg',         'refine3D_auto')
    call cline%set('projfile',    projfile)
    call cline%set('mkdir',       'yes')
    call cline%set('state',       state)
    call cline%set('m_estimator', mode)
    call cline%set('pgrp',        'c1')
    call cline%set('mskdiam',     FLEX_PHANTOM_MSKDIAM)
    call cline%set('maxits',      STATE_REFINE_MAXITS)
    call cline%set('nthr',        nthr)
    call set_fixed_seed(STATE_REFINE_SEED, propagate=.true.)
    call xrefine%execute(cline)
    ! refine3D_auto ran in its own numbered directory
    call simple_getcwd(work_dir)
    work_projfile = work_dir//'/refine3D_auto_state'//int2str_pad(state,2)//'.simple'
    call gate%check(trim(tag)//': the work project is written', file_exists(work_projfile))
    if( file_exists(work_projfile) )then
        call work%read(work_projfile)
        call gate%check(trim(tag)//': the work project holds one state', work%os_ptcl3D%get_n('state') == 1)
        call work%get_vol('vol', 1, map, map_smpd, map_box)
        call gate%check(trim(tag)//': the refined map exists', file_exists(map))
        if( file_exists(map) ) call compare_to_truth(truth, map, FLEX_PHANTOM_MSKDIAM, corr, fsc05, fsc0143, &
            &corr_lp=MAP_CORR_LP)
        if( trim(mode) == 'flex' )then
            call wset%new(work, work%os_ptcl3D, status, message)
            call gate%check(trim(tag)//': the work project''s single-state weight set validates', status == 0)
            if( status == 0 ) call gate%check(trim(tag)//': the work set records the parent state', &
                &wset%get_parent_state() == state .and. wset%get_nstates() == 1)
            call wset%kill
            ind = 0
            do i = 1, work%os_out%get_noris()
                if( .not. work%os_out%isthere(i, 'imgkind') ) cycle
                imgkind = work%os_out%get_str(i, 'imgkind')
                if( imgkind%to_char() == 'vol' .and. work%os_out%get_state(i) == 1 ) ind = i
            enddo
            l_registered = .false.
            if( ind > 0 ) l_registered = work%os_out%isthere(ind, 'mass') .and. work%os_out%isthere(ind, 'ess')
            call gate%check(trim(tag)//': the final map is registered with applied mass and effective sample size', &
                &l_registered)
        endif
        call work%kill
    endif
    call gate%report(trim(tag)//'_truth_correlation', corr)
    call cline%kill
    deallocate(xrefine, cline, work)
end subroutine refine_flex_state

subroutine match_truth_maps( corr, matched_a, matched_b, best_score )
    real,    intent(in)  :: corr(:,:)
    integer, intent(out) :: matched_a, matched_b
    real,    intent(out) :: best_score
    integer :: state_a, state_b
    best_score = -huge(best_score)
    matched_a = 0
    matched_b = 0
    do state_a = 1, size(corr,1)
        do state_b = 1, size(corr,1)
            if( state_a == state_b ) cycle
            if( corr(state_a,1) + corr(state_b,2) > best_score )then
                best_score = corr(state_a,1) + corr(state_b,2)
                matched_a = state_a
                matched_b = state_b
            endif
        enddo
    enddo
end subroutine match_truth_maps

real function best_latent_correlation( a, b ) result( score )
    real(real64), intent(in) :: a(:,:), b(:,:)
    real :: corr(FLEX_PHANTOM_NCOMP,FLEX_PHANTOM_NCOMP)
    integer :: i, j
    do i = 1, FLEX_PHANTOM_NCOMP
        do j = 1, FLEX_PHANTOM_NCOMP
            corr(i,j) = abs(vector_correlation(a(:,i), b(:,j)))
        enddo
    enddo
    score = max(min(corr(1,1),corr(2,2)), min(corr(1,2),corr(2,1)))
end function best_latent_correlation

real function vector_correlation( a, b ) result( corr )
    real(real64), intent(in) :: a(:), b(:)
    real(real64) :: ac(size(a)), bc(size(b)), denom
    ac = a - sum(a) / real(size(a),real64)
    bc = b - sum(b) / real(size(b),real64)
    denom = sqrt(sum(ac*ac) * sum(bc*bc))
    if( denom > 0.0_real64 )then
        corr = real(sum(ac*bc) / denom)
    else
        corr = 0.0
    endif
end function vector_correlation

subroutine map_pair_correlation( a, b, corr )
    class(string), intent(in) :: a, b
    real,          intent(out) :: corr
    real :: fsc05, fsc0143
    call compare_to_truth(a, b, FLEX_PHANTOM_MSKDIAM, corr, fsc05, fsc0143, corr_lp=MAP_CORR_LP)
end subroutine map_pair_correlation

subroutine kill_result( result )
    type(flex_gate_result), intent(inout) :: result
    if( allocated(result%predicted_truth) ) deallocate(result%predicted_truth)
    if( allocated(result%z) ) deallocate(result%z)
    call result%map_a%kill
    call result%map_b%kill
    result = flex_gate_result()
end subroutine kill_result

module subroutine exec_test_write_state_weights_labels( self, cline )
    class(commander_test_write_state_weights_labels), intent(inout) :: self
    class(cmdline),                                   intent(inout) :: cline
    call write_state_weights_from_labels(cline, 1.0)
    call simple_end('**** SIMPLE_TEST_WRITE_STATE_WEIGHTS_LABELS NORMAL STOP ****')
end subroutine exec_test_write_state_weights_labels

module subroutine exec_test_write_state_weights_mixed( self, cline )
    class(commander_test_write_state_weights_mixed), intent(inout) :: self
    class(cmdline),                                  intent(inout) :: cline
    call write_state_weights_from_labels(cline, 0.7)
    call simple_end('**** SIMPLE_TEST_WRITE_STATE_WEIGHTS_MIXED NORMAL STOP ****')
end subroutine exec_test_write_state_weights_mixed

!> Publish a PARTITION state weight set over the project's labelled particles (state > 0): w_own on the
!! labelled state, the rest spread equally over the other states (w_own = 1 gives the hard labels).
!! Written in the working directory and registered in the project (nstates from the command line).
subroutine write_state_weights_from_labels( cline, w_own )
    use simple_state_weight_set, only: state_weight_set
    class(cmdline), intent(inout) :: cline
    real,           intent(in)    :: w_own
    class(sp_project), allocatable :: project
    type(state_weight_set) :: wset
    type(string)           :: projfile
    real,    allocatable   :: weights(:,:)
    integer, allocatable   :: pinds(:), labels(:)
    integer :: nstates, nsel, i, s
    if( .not. cline%defined('projfile') ) THROW_HARD('projfile= is required')
    if( .not. cline%defined('nstates')  ) THROW_HARD('nstates= is required')
    projfile = cline%get_carg('projfile')
    nstates  = cline%get_iarg('nstates')
    allocate(project)
    call project%read(projfile)
    nsel = 0
    do i = 1, project%os_ptcl3D%get_noris()
        s = project%os_ptcl3D%get_state(i)
        if( s >= 1 .and. s <= nstates ) nsel = nsel + 1
    enddo
    allocate(pinds(nsel), labels(nsel), weights(nsel,nstates))
    nsel = 0
    do i = 1, project%os_ptcl3D%get_noris()
        s = project%os_ptcl3D%get_state(i)
        if( s < 1 .or. s > nstates ) cycle
        nsel = nsel + 1
        pinds(nsel)  = i
        labels(nsel) = s
        if( nstates > 1 )then
            weights(nsel,:) = (1. - w_own) / real(nstates - 1)
        else
            weights(nsel,:) = 0.
        endif
        weights(nsel,s) = w_own
    enddo
    call wset%publish(project, project%os_ptcl3D, projfile, pinds, weights, labels, 'test_tool')
    write(logfhandle,'(A,I0,A,I0,A,F5.2)') '>>> STATE WEIGHT SET PUBLISHED FROM LABELS: ', nsel, ' particles x ', nstates, &
        &' states, own-state weight ', w_own
    call wset%kill
    call project%kill
    deallocate(project, pinds, labels, weights)
end subroutine write_state_weights_from_labels

end submodule simple_commanders_test_highlevel_flex
