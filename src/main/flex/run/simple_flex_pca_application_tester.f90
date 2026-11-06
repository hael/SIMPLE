!@descr: library-tier end-to-end tests of the flex_pca application on a deterministic two-state phantom
!! The fixture exercises the production shared-memory commander path for every basis/state backend pair,
!! and once more on the gridding pair with the bandwidth cross-validation (nbins > 1) and a covariance box
!! below the native one (trial maps at box_crop, delivered maps at the native box).
!! It uses a project-backed simulation because publication is part of the application contract: state maps,
!! kernel weights and hard labels must agree. Truth-quality metrics are printed for phase-0 calibration;
!! state-count, chance-separation, same-frame map ordering and weight-set consistency are independent guards.
submodule (simple_flex_pca_application) simple_flex_pca_application_tester
use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
use, intrinsic :: iso_fortran_env, only: real32
use simple_defs,                      only: CWD_GLOB, STDLEN
use simple_syslib,                    only: file_exists, simple_chdir, simple_mkdir
!$ use omp_lib,                       only: omp_get_max_active_levels, omp_set_max_active_levels
use simple_commanders_flex_pca,       only: commander_flex_pca
use simple_oris,                      only: oris
use simple_sp_project,                only: sp_project
use simple_test_truth_metrics,        only: compare_to_truth
use simple_test_flex_phantom_fixture, only: create_flex_pca_phantom_fixture, build_flex_pca_phantom_project,&
    &read_state_weight_table, configure_flex_pca_phantom, FLEX_PHANTOM_BOX, FLEX_PHANTOM_NPTCLS, FLEX_PHANTOM_SMPD,&
    &FLEX_PHANTOM_MSKDIAM, FLEX_PHANTOM_LABEL_ACCURACY_MIN
use simple_test_utils,                only: assert_int, assert_true, enter_fixture, leave_fixture, set_fixed_seed, tests_failed
use simple_ui,                        only: make_ui
implicit none
#include "simple_local_flags.inc"


integer, parameter         :: BOX                             = FLEX_PHANTOM_BOX
integer, parameter         :: DEFAULT_NPTCLS                  = FLEX_PHANTOM_NPTCLS
real,    parameter         :: SMPD                            = FLEX_PHANTOM_SMPD
real,    parameter         :: MSKDIAM                         = FLEX_PHANTOM_MSKDIAM
real,    parameter         :: LABEL_ACCURACY_MIN              = FLEX_PHANTOM_LABEL_ACCURACY_MIN
integer, parameter         :: NARMS                           = 5
integer, parameter         :: NTHR                            = 4
integer, parameter         :: SEED_APP                        = 202610063
real,    parameter         :: MAP_CORR_LP                     = 8.0

character(len=8),  parameter :: BASIS_BACKEND(NARMS) = [character(len=8) :: &
    &'gridding', 'pcg',      'gridding', 'pcg',      'gridding']
character(len=8),  parameter :: STATE_BACKEND(NARMS) = [character(len=8) :: &
    &'gridding', 'gridding', 'pcg',      'pcg',      'gridding']
character(len=18), parameter :: ARM_TAG(NARMS)       = [character(len=18) :: &
    &'grid_basis_grid', 'pcg_basis_grid', 'grid_basis_pcg', 'pcg_basis_pcg', 'grid_cv_crop48']
!> bandwidth cross-validation bins of each arm; 1 skips the cross-validation
integer,           parameter :: ARM_NBINS(NARMS)     = [1, 1, 1, 1, 4]
!> covariance box of each arm (box_crop)
integer,           parameter :: ARM_BOX_CROP(NARMS)  = [BOX, BOX, BOX, BOX, 48]

type :: arm_result
    integer :: nstates = 0
    real    :: label_accuracy = 0.0
    real    :: pc_corr = 0.0
    real    :: map_corr = 0.0
    real    :: seconds = 0.0
end type arm_result

contains

    module subroutine run_all_flex_pca_application_tests()
        type(string) :: cwd_saved, root
        integer      :: nfail0
        write(*,'(A)') '**** running flex_pca application phantom ****'
        nfail0 = tests_failed
        call enter_fixture('flex_pca_two_state', cwd_saved, root)
        CWD_GLOB = root%to_char()
        call make_ui
        call test_two_state_backend_matrix(root)
        call leave_fixture(cwd_saved, root, nfail0)
        CWD_GLOB = cwd_saved%to_char()
    end subroutine run_all_flex_pca_application_tests

    !> The same simulated particles and truth are used in every arm. Backend changes therefore
    !! cannot hide behind a changed random draw. The application itself is re-seeded per arm.
    subroutine test_two_state_backend_matrix( root )
        class(string), intent(in) :: root
        type(arm_result) :: result(NARMS)
        type(string) :: truth_a, truth_b, truth_mean, truth_diff, oritab, stack, arm_dir
        integer :: arm
        call create_flex_pca_phantom_fixture(root, NTHR, truth_a, truth_b, truth_mean, truth_diff, oritab, stack)
        do arm = 1, NARMS
            arm_dir = root//'/'//trim(ARM_TAG(arm))
            call simple_mkdir(arm_dir)
            call build_flex_pca_phantom_project(arm_dir//'/phantom.simple', stack, oritab)
            call simple_chdir(arm_dir)
            CWD_GLOB = arm_dir%to_char()
            call run_application_arm(BASIS_BACKEND(arm), STATE_BACKEND(arm), ARM_NBINS(arm), ARM_BOX_CROP(arm), &
                &truth_mean, result(arm))
            call validate_arm(ARM_TAG(arm), truth_a, truth_b, truth_diff, result(arm))
            call simple_chdir(root)
            CWD_GLOB = root%to_char()
        enddo
        write(logfhandle,'(A)') '>>> FLEX_PCA TWO-STATE PHANTOM PHASE-0 OBSERVATIONS'
        do arm = 1, NARMS
            write(logfhandle,'(A,A,A,I0,A,F7.4,A,F7.4,A,F7.4,A,F9.1)') '>>>   ', trim(ARM_TAG(arm)), &
                &' states=', result(arm)%nstates, ' label_accuracy=', result(arm)%label_accuracy, &
                &' pc_corr=', result(arm)%pc_corr, ' min_matched_corr=', result(arm)%map_corr, &
                &' commander_seconds=', result(arm)%seconds
        enddo
    end subroutine test_two_state_backend_matrix

    subroutine run_application_arm( basis_backend, state_backend, nbins, box_crop, truth_mean, result )
        character(len=*), intent(in)  :: basis_backend, state_backend
        integer,          intent(in)  :: nbins, box_crop
        class(string),    intent(in)  :: truth_mean
        type(arm_result), intent(out) :: result
        class(cmdline),            allocatable :: cline
        class(commander_flex_pca), allocatable :: flex_pca
        integer(timer_int_kind) :: started
        integer :: levels_saved
        allocate(cline, flex_pca)
        result = arm_result()
        call configure_flex_pca_phantom(cline, truth_mean, basis_backend, state_backend, 1, NTHR, nbins=nbins, &
            &box_crop=box_crop)
        call set_fixed_seed(SEED_APP, propagate=.true.)
        ! the run leaves the process's OpenMP nesting level as it found it
        !$ levels_saved = omp_get_max_active_levels()
        !$ call omp_set_max_active_levels(2)
        started = tic()
        call flex_pca%execute(cline)
        result%seconds = real(toc(started))
        !$ call assert_int(2, omp_get_max_active_levels(), 'flex_pca leaves the OpenMP nesting level as found')
        !$ call omp_set_max_active_levels(levels_saved)
        call cline%kill
        deallocate(flex_pca, cline)
    end subroutine run_application_arm

    subroutine validate_arm( tag, truth_a, truth_b, truth_diff, result )
        character(len=*), intent(in)    :: tag
        class(string),    intent(in)    :: truth_a, truth_b, truth_diff
        type(arm_result), intent(inout) :: result
        class(sp_project), allocatable :: project
        class(oris),       allocatable :: field
        real(real32), allocatable :: weights(:,:)
        integer,      allocatable :: labels(:)
        integer, allocatable :: truth_counts(:,:)
        real,    allocatable :: corr(:,:)
        integer :: status, i, state, state_a, state_b, matched_a, matched_b, truth_label, predicted
        integer :: nlabel_mismatch, nweight_mismatch, nselection_mismatch
        character(len=STDLEN) :: message
        character(len=STDLEN) :: state_tag
        real :: best_score
        allocate(project, field)
        call project%read(string('phantom.simple'))
        field = project%os_ptcl3D
        call read_state_weight_table(project, field, result%nstates, weights, labels, status, message)
        call assert_int(0, status, trim(tag)//': delivered flex-weight set validates')
        if( status /= 0 )then
            write(logfhandle,'(A,A)') '>>> FLEX_PCA PHANTOM weight-store failure: ', trim(message)
            call field%kill
            call project%kill
            deallocate(field, project)
            if( allocated(weights) ) deallocate(weights)
            if( allocated(labels)  ) deallocate(labels)
            return
        endif
        call assert_true(result%nstates >= 2, trim(tag)//': publication retains representatives of both truth modes')
        call assert_true(result%nstates <= 3, trim(tag)//': automatic merge respects the requested state ceiling')
        allocate(truth_counts(result%nstates,2), source=0)
        nlabel_mismatch = 0
        nweight_mismatch = 0
        nselection_mismatch = 0
        do i = 1, DEFAULT_NPTCLS
            if( (sum(weights(:,i)) > 0.) .neqv. (project%os_ptcl3D%get_state(i) > 0) ) &
                &nselection_mismatch = nselection_mismatch + 1
            truth_label = merge(1, 2, mod(i,2) == 1)
            if( labels(i) >= 1 .and. labels(i) <= result%nstates ) &
                &truth_counts(labels(i),truth_label) = truth_counts(labels(i),truth_label) + 1
            if( int(labels(i)) /= project%os_ptcl3D%get_state(i) ) nlabel_mismatch = nlabel_mismatch + 1
            if( labels(i) > 0 )then
                predicted = maxloc(weights(:,i), dim=1)
                if( int(labels(i)) /= predicted ) nweight_mismatch = nweight_mismatch + 1
            endif
        enddo
        result%label_accuracy = real(sum(maxval(truth_counts,dim=2))) / real(DEFAULT_NPTCLS)
        call assert_true(result%label_accuracy >= LABEL_ACCURACY_MIN, &
            &trim(tag)//': published state labels separate the two truth modes from chance')
        call assert_int(0, nlabel_mismatch, trim(tag)//': weight flags and project labels agree')
        call assert_int(0, nweight_mismatch, trim(tag)//': hard labels select the maximum kernel weight')
        call assert_int(0, nselection_mismatch, trim(tag)//': the weight set selects exactly the particles with state > 0')
        result%pc_corr = leading_pc_correlation(tag, truth_diff)
        call assert_true(ieee_is_finite(result%pc_corr), trim(tag)//': leading eigenvolume correlation is finite')
        allocate(corr(result%nstates,2))
        do state = 1, result%nstates
            write(state_tag,'(A,"_s",I0)') trim(tag), state
            call state_truth_correlations(project, state, truth_a, truth_b, trim(state_tag), corr(state,:))
        enddo
        best_score = -huge(best_score)
        matched_a = 0
        matched_b = 0
        do state_a = 1, result%nstates
            do state_b = 1, result%nstates
                if( state_a == state_b ) cycle
                if( corr(state_a,1) + corr(state_b,2) > best_score )then
                    best_score = corr(state_a,1) + corr(state_b,2)
                    matched_a = state_a
                    matched_b = state_b
                endif
            enddo
        enddo
        if( matched_a > 0 .and. matched_b > 0 )then
            result%map_corr = min(corr(matched_a,1), corr(matched_b,2))
            call assert_true(corr(matched_a,1) > corr(matched_a,2), &
                &trim(tag)//': matched mode-A map prefers truth A over truth B')
            call assert_true(corr(matched_b,2) > corr(matched_b,1), &
                &trim(tag)//': matched mode-B map prefers truth B over truth A')
        else
            result%map_corr = 0.0
            call assert_true(.false., trim(tag)//': two distinct delivered maps can be matched to truth')
        endif
        write(logfhandle,'(A,A,A)') '>>> FLEX_PCA PHANTOM ', trim(tag), &
            &' state truth counts/correlations: state nA nB corrA corrB'
        do state = 1, result%nstates
            write(logfhandle,'(A,I0,2(1X,I0),2(1X,F9.6))') '>>>   ', state, truth_counts(state,:), corr(state,:)
        enddo
        write(logfhandle,'(A,I0,A,I0,A,F9.6)') '>>>   matched truth states: A=', matched_a, &
            &' B=', matched_b, ' min_own_corr=', result%map_corr
        call field%kill
        call project%kill
        deallocate(field, project)
        call ab_dump(tag, weights, labels)
        deallocate(weights, labels, truth_counts, corr)
    end subroutine validate_arm

    !> With SIMPLE_FLEX_AB_DIR set, the arm's embedding file and its published weights and hard labels
    !! are kept there as <tag>_embedding.bin and <tag>_weights.txt (row: label, then one weight per
    !! state), so two trees can be compared (the FLEX A/B of the flex_on_simple plan, section 8.2)
    subroutine ab_dump( tag, weights, labels )
        use simple_fileio, only: simple_copy_file
        character(len=*), intent(in) :: tag
        real(real32),     intent(in) :: weights(:,:)   !< (nstates, nptcls)
        integer,          intent(in) :: labels(:)
        character(len=4096) :: dir
        integer :: dirlen, status, u, i
        call get_environment_variable('SIMPLE_FLEX_AB_DIR', dir, dirlen, status)
        if( status /= 0 .or. dirlen == 0 ) return
        call simple_copy_file(string('flex_pca_embedding.bin'), string(trim(dir)//'/'//trim(tag)//'_embedding.bin'))
        open(newunit=u, file=trim(dir)//'/'//trim(tag)//'_weights.txt', status='replace', action='write')
        do i = 1, size(labels)
            write(u,'(I0,*(1X,ES16.9))') labels(i), weights(:,i)
        enddo
        close(u)
    end subroutine ab_dump

    !> on the eigenvolume's box (the covariance box, box_crop), the truth Fourier-cropped to it
    real function leading_pc_correlation( tag, truth_diff ) result( corr )
        use simple_imghead, only: find_ldim_nptcls
        character(len=*), intent(in) :: tag
        class(string),    intent(in) :: truth_diff
        type(image) :: truth, pc
        type(string) :: pcfile
        real    :: smpd_pc
        integer :: ldim_pc(3), ifoo
        corr = 0.0
        pcfile = 'flex_pca_polished_pc001.mrc'
        call assert_true(file_exists(pcfile), trim(tag)//': polished leading eigenvolume is delivered')
        if( .not. file_exists(pcfile) ) return
        call find_ldim_nptcls(pcfile, ldim_pc, ifoo)
        smpd_pc = SMPD * real(BOX) / real(ldim_pc(1))
        call truth%read_and_crop(truth_diff, SMPD, ldim_pc(1), smpd_pc)
        call pc%new(ldim_pc, smpd_pc, wthreads=.false.)
        call pc%read(pcfile)
        corr = abs(truth%real_corr(pc))
        call truth%kill
        call pc%kill
    end function leading_pc_correlation

    subroutine state_truth_correlations( project, state, truth_a, truth_b, tag, corr )
        class(sp_project), intent(in)  :: project
        integer,           intent(in)  :: state
        class(string),     intent(in)  :: truth_a, truth_b
        character(len=*),  intent(in)  :: tag
        real,              intent(out) :: corr(2)
        type(string) :: map
        real :: smpd, fsc05, fsc0143
        integer :: box
        call project%get_vol('vol', state, map, smpd, box)
        call assert_true(file_exists(map), trim(tag)//': delivered state map exists')
        if( .not. file_exists(map) )then
            corr = 0.0
            return
        endif
        call compare_to_truth(truth_a, map, MSKDIAM, corr(1), fsc05, fsc0143, corr_lp=MAP_CORR_LP)
        call compare_to_truth(truth_b, map, MSKDIAM, corr(2), fsc05, fsc0143, corr_lp=MAP_CORR_LP)
        call map%kill
    end subroutine state_truth_correlations

end submodule simple_flex_pca_application_tester
