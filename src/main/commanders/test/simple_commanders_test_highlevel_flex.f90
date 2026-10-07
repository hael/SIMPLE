!@descr: shared-memory and distributed FLEX PCA workflow gate
submodule (simple_commanders_test_highlevel) simple_commanders_test_highlevel_flex
use, intrinsic :: iso_fortran_env, only: int32, real32, real64
use simple_commanders_flex_pca, only: commander_flex_pca
use simple_flex_pca_application_tester, only: create_flex_pca_phantom_fixture, &
    &build_flex_pca_phantom_project, configure_flex_pca_phantom, FLEX_PHANTOM_BOX, &
    &FLEX_PHANTOM_NPTCLS, FLEX_PHANTOM_NCOMP, FLEX_PHANTOM_MAX_STATES, &
    &FLEX_PHANTOM_SMPD, FLEX_PHANTOM_MSKDIAM, FLEX_PHANTOM_LABEL_ACCURACY_MIN
use simple_flex_pca_embedding_io, only: read_deconv_block
use simple_flex_weights_state, only: flex_weights_store
use simple_image, only: image
use simple_oris, only: oris
use simple_sp_project, only: sp_project
use simple_test_gate, only: test_gate
use simple_test_truth_metrics, only: compare_to_truth
use simple_ui, only: make_ui
implicit none
#include "simple_local_flags.inc"

integer, parameter :: FLEX_GATE_SEED = 202610063
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
    type(flex_weights_store) :: store
    real(real32), allocatable :: weights(:,:)
    real(real64), allocatable :: scalars(:,:), precision(:,:,:)
    integer(int32), allocatable :: labels(:)
    integer, allocatable :: deconv_labels(:), truth_counts(:,:), state_truth(:)
    real, allocatable :: corr(:,:)
    type(string) :: map
    character(len=STDLEN) :: message
    integer :: status, i, state, truth_label, predicted, nlabel_mismatch, nweight_mismatch, map_box
    real :: best_score, fsc05, fsc0143, map_smpd
    result = flex_gate_result()
    allocate(project, field)
    call project%read(string('phantom.simple'))
    field = project%os_ptcl3D
    call store%new(project, field, FLEX_PHANTOM_BOX, FLEX_PHANTOM_SMPD, status, message)
    if( status == 0 ) call store%take(result%nstates, weights, labels, scalars)
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
        call project%get_vol('vol_flex', state, map, map_smpd, map_box)
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
        call project%get_vol('vol_flex', result%matched_a, map, map_smpd, map_box)
        result%map_a = simple_abspath(map)
        call map%kill
        call project%get_vol('vol_flex', result%matched_b, map, map_smpd, map_box)
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
    deallocate(weights, labels, scalars, truth_counts, state_truth, corr)
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

end submodule simple_commanders_test_highlevel_flex
