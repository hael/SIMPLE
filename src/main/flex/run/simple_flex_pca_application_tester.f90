!@descr: library-tier end-to-end tests of the flex_pca application on a deterministic two-state phantom
!! The fixture exercises the production shared-memory commander path for every basis/state backend pair,
!! and once more on the gridding pair with the bandwidth cross-validation (nbins > 1) and a covariance box
!! below the native one (trial maps at box_crop, delivered maps at the native box).
!! It uses a project-backed simulation because publication is part of the application contract: state maps,
!! kernel weights and hard labels must agree. Truth-quality metrics are printed for phase-0 calibration;
!! state-count, chance-separation, same-frame map ordering and weight-set consistency are independent guards.
module simple_flex_pca_application_tester
use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
use, intrinsic :: iso_fortran_env, only: int32, real32, real64
use simple_core_module_api, only: ctfflag_yes, ctfparams, cwd_glob, file_exists, int32, logfhandle, oris, real32, &
    &real64, simple_abspath, simple_chdir, simple_exception, simple_mkdir, stdlen, string, tic, timer_int_kind, toc
use simple_cmdline,             only: cmdline
use simple_commanders_flex_pca, only: commander_flex_pca
use simple_commanders_sim,      only: commander_simulate_particles
use simple_state_weight_set,    only: state_weight_set
use simple_image,               only: image
use simple_oris,                only: oris
use simple_sp_project,          only: sp_project
use simple_test_truth_metrics,  only: add_gaussian_blob, compare_to_truth
use simple_test_utils
use simple_type_defs,           only: CTFFLAG_YES, ctfparams
use simple_ui,                  only: make_ui
implicit none
private

public :: run_all_flex_pca_application_tests
public :: create_flex_pca_phantom_fixture, build_flex_pca_phantom_project, read_state_weight_table
public :: configure_flex_pca_phantom

integer, parameter         :: BOX                             = 64
integer, parameter         :: NPTCLS                          = 2000
integer, parameter         :: NHALF                           = NPTCLS / 2
integer, parameter         :: NARMS                           = 5
integer, parameter         :: MIN_NEFF                        = 100
integer, parameter         :: NTHR                            = 4
integer, parameter         :: SEED_A                          = 202610061
integer, parameter         :: SEED_B                          = 202610062
integer, parameter         :: SEED_APP                        = 202610063
real,    parameter         :: SMPD                            = 2.0
real,    parameter         :: MSKDIAM                         = 100.0
real,    parameter         :: KV                              = 300.0
real,    parameter         :: CS                              = 2.7
real,    parameter         :: FRACA                           = 0.1
real,    parameter         :: DEFOCUS                         = 1.5
real,    parameter         :: SIM_SNR                         = 1.0
real,    parameter         :: MAP_CORR_LP                     = 8.0
! The delivered state ids are recoded independently to their majority truth class. With at most
! three state ids there are only 2^3 binary recodings; under random truth labels Hoeffding's union
! bound puts P(best accuracy >= 0.60) below 8*exp(-2*NPTCLS*0.10^2) < 4e-17. This is an independent
! detection floor, not the phase-0 production observation used later as a regression floor.
real,    parameter         :: LABEL_ACCURACY_MIN              = 0.60
integer, parameter, public :: FLEX_PHANTOM_BOX                = BOX
integer, parameter, public :: FLEX_PHANTOM_NPTCLS             = NPTCLS
integer, parameter, public :: FLEX_PHANTOM_NCOMP              = 2
integer, parameter, public :: FLEX_PHANTOM_MAX_STATES         = 3
real,    parameter, public :: FLEX_PHANTOM_SMPD               = SMPD
real,    parameter, public :: FLEX_PHANTOM_MSKDIAM            = MSKDIAM
real,    parameter, public :: FLEX_PHANTOM_LABEL_ACCURACY_MIN = LABEL_ACCURACY_MIN

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

    subroutine run_all_flex_pca_application_tests()
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

    !> Build the deterministic truth maps, shared orientations and interleaved two-state stack.
    !! Both the library fixture and the high-level shared/distributed gate use this exact recipe.
    subroutine create_flex_pca_phantom_fixture( root, nthr, truth_a, truth_b, truth_mean, truth_diff, oritab, stack )
        class(string), intent(in)  :: root
        integer,       intent(in)  :: nthr
        type(string),  intent(out) :: truth_a, truth_b, truth_mean, truth_diff, oritab, stack
        type(string) :: stk_a, stk_b
        call write_phantom_truths(root, truth_a, truth_b, truth_mean, truth_diff)
        oritab = root//'/truth_orientations.txt'
        call write_uniform_orientations(oritab)
        call simple_chdir(root)
        stk_a = simulate_stack(truth_a, oritab, string('particles_state_a.mrcs'), SEED_A, nthr)
        stk_b = simulate_stack(truth_b, oritab, string('particles_state_b.mrcs'), SEED_B, nthr)
        stack = interleave_stacks(stk_a, stk_b, string('particles_two_state.mrcs'))
        call stk_a%kill
        call stk_b%kill
    end subroutine create_flex_pca_phantom_fixture

    !> Common asymmetric body plus one narrow lobe moved by exactly eight pixels between A and B.
    subroutine write_phantom_truths( root, truth_a, truth_b, truth_mean, truth_diff )
        class(string), intent(in)  :: root
        type(string),  intent(out) :: truth_a, truth_b, truth_mean, truth_diff
        integer, parameter :: NBLOBS        = 6
        real,    parameter :: POS(3,NBLOBS) = reshape([ &
            &-18., -10.,   4.,   14.,   8.,  -8.,   -6.,  18.,  12., &
            &  4., -16., -14.,   20., -12.,  16.,  -20.,  13., -11.], [3,NBLOBS])
        real,    parameter :: SIGMA(NBLOBS) = [7.0, 6.0, 5.5, 6.5, 4.5, 5.0]
        real,    parameter :: AMP(NBLOBS)   = [1.0, 0.85, 0.75, 0.70, 0.55, 0.50]
        type(image) :: vol_a, vol_b, vol_mean, vol_diff
        real, allocatable :: base(:,:,:), a(:,:,:), b(:,:,:)
        real :: ctr(3), d2
        integer :: blob, ix, iy, iz
        truth_a    = root//'/truth_state_a.mrc'
        truth_b    = root//'/truth_state_b.mrc'
        truth_mean = root//'/truth_mean.mrc'
        truth_diff = root//'/truth_a_minus_b.mrc'
        allocate(base(BOX,BOX,BOX), source=0.0)
        ctr = real([BOX,BOX,BOX]) / 2.0 + 1.0
        do iz = 1, BOX
            do iy = 1, BOX
                do ix = 1, BOX
                    do blob = 1, NBLOBS
                        d2 = sum((([real(ix),real(iy),real(iz)] - ctr) * SMPD - POS(:,blob))**2)
                        base(ix,iy,iz) = base(ix,iy,iz) + AMP(blob) * exp(-0.5 * d2 / SIGMA(blob)**2)
                    enddo
                enddo
            enddo
        enddo
        call vol_a%new([BOX,BOX,BOX], SMPD, wthreads=.false.)
        call vol_a%set_rmat(base, .false.)
        call vol_a%write(truth_a)
        call vol_a%write(truth_b)
        call add_gaussian_blob(truth_a, [-8.0,-2.0,0.0], 4.0, 0.85)
        call add_gaussian_blob(truth_b, [ 8.0,-2.0,0.0], 4.0, 0.85)
        call vol_a%read(truth_a)
        call vol_b%new([BOX,BOX,BOX], SMPD, wthreads=.false.)
        call vol_b%read(truth_b)
        a = vol_a%get_rmat()
        b = vol_b%get_rmat()
        call vol_mean%new([BOX,BOX,BOX], SMPD, wthreads=.false.)
        call vol_mean%set_rmat(0.5 * (a + b), .false.)
        call vol_mean%write(truth_mean)
        call vol_diff%new([BOX,BOX,BOX], SMPD, wthreads=.false.)
        call vol_diff%set_rmat(a - b, .false.)
        call vol_diff%write(truth_diff)
        call vol_a%kill
        call vol_b%kill
        call vol_mean%kill
        call vol_diff%kill
        deallocate(base, a, b)
    end subroutine write_phantom_truths

    subroutine write_uniform_orientations( fname )
        class(string), intent(in) :: fname
        class(oris), allocatable :: poses
        allocate(poses)
        call poses%new(NHALF, is_ptcl=.true.)
        call poses%spiral()
        call poses%write(fname, [1,NHALF])
        call poses%kill
        deallocate(poses)
    end subroutine write_uniform_orientations

    function simulate_stack( volume, oritab, outstk, seed, nthr ) result( stack )
        class(string), intent(in) :: volume, oritab, outstk
        integer,       intent(in) :: seed, nthr
        type(string) :: stack
        class(commander_simulate_particles), allocatable :: simulate
        class(cmdline), allocatable :: cline
        allocate(simulate, cline)
        call set_fixed_seed(seed, propagate=.true.)
        call cline%set('prg',      'simulate_particles')
        call cline%set('mkdir',    'no')
        call cline%set('vol1',     volume)
        call cline%set('oritab',   oritab)
        call cline%set('outstk',   outstk)
        call cline%set('outfile',  outstk//'.txt')
        call cline%set('nptcls',   NHALF)
        call cline%set('nthr',     nthr)
        call cline%set('smpd',     SMPD)
        call cline%set('mskdiam',  MSKDIAM)
        call cline%set('pgrp',     'c1')
        call cline%set('ctf',      'yes')
        call cline%set('snr',      SIM_SNR)
        call cline%set('kv',       KV)
        call cline%set('cs',       CS)
        call cline%set('fraca',    FRACA)
        call cline%set('defocus',  DEFOCUS)
        call cline%set('dferr',    0.0)
        call cline%set('astigerr', 0.0)
        call cline%set('bfac',     0.0)
        call cline%set('bfacerr',  0.0)
        call cline%set('sherr',    0.0)
        call simulate%execute(cline)
        call cline%kill
        stack = simple_abspath(outstk)
        deallocate(simulate, cline)
    end function simulate_stack

    function interleave_stacks( stk_a, stk_b, outstk ) result( stack )
        class(string), intent(in) :: stk_a, stk_b, outstk
        type(string) :: stack
        type(image)  :: img
        integer :: i
        call img%new([BOX,BOX,1], SMPD, wthreads=.false.)
        do i = 1, NHALF
            call img%read(stk_a, i)
            call img%write(outstk, 2*i - 1)
            call img%read(stk_b, i)
            call img%write(outstk, 2*i)
        enddo
        call img%kill
        stack = simple_abspath(outstk)
    end function interleave_stacks

    !> Both truth states occur in both even/odd halfsets. Rows 2i-1 and 2i share their pose,
    !! so orientation coverage cannot become a state cue.
    subroutine build_flex_pca_phantom_project( projfile, stack, oritab )
        class(string), intent(in) :: projfile, stack, oritab
        class(sp_project), allocatable :: project
        class(oris),       allocatable :: poses
        type(ctfparams) :: ctf
        type(cmdline)   :: cline_compenv
        integer :: i, ipose, eo
        allocate(project, poses)
        call poses%new(NHALF, is_ptcl=.true.)
        call poses%read(oritab, [1,NHALF])
        ctf%smpd    = SMPD
        ctf%kv      = KV
        ctf%cs      = CS
        ctf%fraca   = FRACA
        ctf%dfx     = DEFOCUS
        ctf%dfy     = DEFOCUS
        ctf%angast  = 0.0
        ctf%ctfflag = CTFFLAG_YES
        call project%add_stk(stack, ctf)
        do i = 1, NPTCLS
            ipose = (i + 1) / 2
            eo     = mod(ipose, 2)
            call project%os_ptcl3D%set_euler(i, poses%get_euler(ipose))
            call project%os_ptcl3D%set_state(i, 1)
            call project%os_ptcl3D%set(i, 'eo',   real(eo))
            call project%os_ptcl3D%set(i, 'proj', 1.0)
            call project%os_ptcl2D%set_euler(i, poses%get_euler(ipose))
            call project%os_ptcl2D%set_state(i, 1)
            call project%os_ptcl2D%set(i, 'eo', real(eo))
        enddo
        call project%update_projinfo(projfile)
        ! the computing environment new_project writes; the distributed master's queue setup reads it
        call project%update_compenv(cline_compenv)
        call project%write(projfile)
        call poses%kill
        call project%kill
        deallocate(project, poses)
    end subroutine build_flex_pca_phantom_project

    !> The project's state weight set as a table: weights(nstates,nptcls) and the hard labels(nptcls)
    subroutine read_state_weight_table( project, field, nstates, weights, labels, status, message )
        class(sp_project),         intent(in)  :: project
        class(oris),               intent(in)  :: field
        integer,                   intent(out) :: nstates
        real(real32), allocatable, intent(out) :: weights(:,:)
        integer,      allocatable, intent(out) :: labels(:)
        integer,                   intent(out) :: status
        character(len=*),          intent(out) :: message
        type(state_weight_set)    :: weight_set
        real(real32), allocatable :: column(:)
        integer :: s
        nstates = 0
        call weight_set%new(project, field, status, message)
        if( status /= 0 ) return
        nstates = weight_set%get_nstates()
        allocate(weights(nstates, weight_set%get_nptcls()))
        do s = 1, nstates
            call weight_set%get_weights(s, column)
            weights(s,:) = column
        enddo
        call weight_set%get_labels(labels)
        call weight_set%kill
    end subroutine read_state_weight_table

    !> The common production command line. nparts=1 selects shared memory; nparts>1 selects
    !! the distributed master and its real worker processes through the local queue system.
    subroutine configure_flex_pca_phantom( cline, truth_mean, basis_backend, state_backend, nparts, nthr, nbins, box_crop )
        class(cmdline),    intent(inout) :: cline
        class(string),     intent(in)    :: truth_mean
        character(len=*),  intent(in)    :: basis_backend, state_backend
        integer,           intent(in)    :: nparts, nthr
        integer, optional, intent(in)    :: nbins, box_crop
        call cline%set('prg',               'flex_pca')
        call cline%set('projfile',          'phantom.simple')
        call cline%set('mkdir',             'no')
        call cline%set('oritype',           'ptcl3D')
        call cline%set('vol1',              truth_mean)
        call cline%set('pgrp',              'c1')
        call cline%set('mskdiam',           MSKDIAM)
        call cline%set('box_crop',          BOX)
        if( present(box_crop) ) call cline%set('box_crop', box_crop)
        call cline%set('lp',                8.0)
        call cline%set('nstates',           1)
        call cline%set('npreimages',        3)
        call cline%set('preimage_auto',     'yes')
        call cline%set('min_neff',          MIN_NEFF)
        call cline%set('neigs',             FLEX_PHANTOM_NCOMP)
        call cline%set('nkern',             1)
        call cline%set('state_axis',        1)
        call cline%set('n_probe_iters',     2)
        call cline%set('column_separation', 2)
        call cline%set('nbins',             1)
        if( present(nbins) ) call cline%set('nbins', nbins)
        call cline%set('rec_states',        'yes')
        call cline%set('rec_backend',       basis_backend)
        call cline%set('rec_states_backend', state_backend)
        call cline%set('maxits_pcg',        4)
        call cline%set('rtol',              0.0)
        call cline%set('objfun',            'euclid')
        call cline%set('sigma_est',         'global')
        call cline%set('umap',              'no')
        call cline%set('cache',             'no')
        call cline%set('qsys_name',         'local')
        call cline%set('nparts',            nparts)
        call cline%set('nthr',              nthr)
    end subroutine configure_flex_pca_phantom

    subroutine run_application_arm( basis_backend, state_backend, nbins, box_crop, truth_mean, result )
        character(len=*), intent(in)  :: basis_backend, state_backend
        integer,          intent(in)  :: nbins, box_crop
        class(string),    intent(in)  :: truth_mean
        type(arm_result), intent(out) :: result
        class(cmdline),            allocatable :: cline
        class(commander_flex_pca), allocatable :: flex_pca
        integer(timer_int_kind) :: started
        allocate(cline, flex_pca)
        result = arm_result()
        call configure_flex_pca_phantom(cline, truth_mean, basis_backend, state_backend, 1, NTHR, nbins=nbins, &
            &box_crop=box_crop)
        call set_fixed_seed(SEED_APP, propagate=.true.)
        started = tic()
        call flex_pca%execute(cline)
        result%seconds = real(toc(started))
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
        call assert_int(status, 0, trim(tag)//': delivered flex-weight set validates')
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
        do i = 1, NPTCLS
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
        result%label_accuracy = real(sum(maxval(truth_counts,dim=2))) / real(NPTCLS)
        call assert_true(result%label_accuracy >= LABEL_ACCURACY_MIN, &
            &trim(tag)//': published state labels separate the two truth modes from chance')
        call assert_int(nlabel_mismatch, 0, trim(tag)//': weight flags and project labels agree')
        call assert_int(nweight_mismatch, 0, trim(tag)//': hard labels select the maximum kernel weight')
        call assert_int(nselection_mismatch, 0, trim(tag)//': the weight set selects exactly the particles with state > 0')
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

end module simple_flex_pca_application_tester
