!@descr: unit tests for the Cartesian 3D strategy of a Cartesian pass (simple_strategy3D_cont)
! A one-state Cartesian calculator in a builder whose ptcl3D field holds the seeds, each slot the
! noise-free truth prediction centred on its stored shift; no files, no polar calculator. Pinned:
! the seed contract, the lifecycle, accepted and rejected commits (pose, corr_cart, improved flag,
! status, convergence fields), the euclid sigma2 residual at the committed pose, threaded batches
! equal to serial, particles outside the sample untouched, continuous pass fields, the athres_cont bound.
module simple_strategy3D_cont_tester
use simple_test_cftc_fixture, only: objective_gradient
use simple_core_module_api, only: dp, euler2m, string, rad2deg, del_file
!$ use omp_lib,             only: omp_get_max_threads, omp_get_num_procs, omp_set_num_threads
use simple_defs_ori,        only: N_PTCL_ORIPARAMS
use simple_builder,         only: builder
use simple_parameters,      only: parameters
use simple_ori,             only: ori
use simple_oris,            only: oris
use simple_cartft_calc,     only: cartft_calc
use simple_strategy3D,      only: strategy3D
use simple_strategy3D_alloc, only: s3D, clean_strategy3D, prep_strategy3D
use simple_eul_prob_tab_utils, only: cart_shift_prior_sigma2
use simple_strategy3D_srch, only: strategy3D_spec
use simple_strategy3D_cont, only: strategy3D_cont, cont_seed_is_valid, cont_reference_geometry_is_valid, check_cont_seeds
use simple_cartft_pose_opt, only: CARTFT_ACCEPTED
use simple_type_defs,       only: ctfparams, CTFFLAG_NO, OBJFUN_CC, OBJFUN_EUCLID
use simple_forked_process, only: forked_process, FORK_STATUS_FAILED, FORK_STATUS_RUNNING
use simple_cmdline, only: cmdline
use, intrinsic :: ieee_arithmetic, only: ieee_value, ieee_quiet_nan, ieee_positive_inf
use simple_test_utils
use simple_test_sigma2_fixture, only: write_test_sigma2_store
implicit none
private
public :: run_all_strategy3D_cont_tests

integer, parameter :: BOX = 16, NPTCLS = 6
integer, parameter :: KFROMTO(2) = [2, BOX/2 - 1]
! single-precision Euler angles of the stored poses limit the recovered rotation
real(dp), parameter :: ROTATION_TOL = 2.e-3_dp
real,     parameter :: SHIFT_TOL    = 5.e-3

type, extends(forked_process) :: invalid_prior_width_child
contains
    procedure :: execute => reject_invalid_prior_width
end type invalid_prior_width_child

type, extends(forked_process) :: prior_failure_master_child
contains
    procedure :: execute => reject_failed_prior_part
end type prior_failure_master_child

contains

    subroutine run_all_strategy3D_cont_tests()
        write(*,'(A)') '**** running all Cartesian 3D strategy tests ****'
        s3D%cont_shift_prior_sigma2 = 9._dp
        call clean_strategy3D()
        call assert_true(s3D%cont_shift_prior_sigma2 == 0._dp, 'Cartesian prior width must reset between passes')
        write(*,'(A)') 'test_shift_prior_width'
        call test_shift_prior_width()
        write(*,'(A)') 'test_shift_prior_angstrom_units'
        call test_shift_prior_angstrom_units()
        write(*,'(A)') 'test_invalid_prior_width (expected child rejections)'
        call test_invalid_prior_width()
        write(*,'(A)') 'test_prior_master_failure (expected worker and master rejections)'
        call test_prior_master_failure()
        call test_prior_strategy_commit()
        write(*,'(A)') 'test_seed_contract'
        call test_seed_contract()
        call test_reference_geometry()
        call test_selected_seeds()
        call test_particle_seed_contract()
        write(*,'(A)') 'test_lifecycle_and_rejection'
        call test_lifecycle_and_rejection()
        write(*,'(A)') 'test_commit_contract'
        call test_commit_contract()
        write(*,'(A)') 'test_threaded_batch'
        call test_threaded_batch()
        write(*,'(A)') 'test_outside_sample_untouched'
        call test_outside_sample_untouched()
        write(*,'(A)') 'test_followup_keeps_convergence_fields'
        call test_followup_keeps_convergence_fields()
        write(*,'(A)') 'test_rotation_bound'
        call test_rotation_bound()
    end subroutine run_all_strategy3D_cont_tests

    !> All active records, including zeros, contribute independently of sampling markers.
    subroutine test_shift_prior_width()
        type(oris) :: os
        call os%new(4, is_ptcl=.true.)
        call os%set_state(1, 1)
        call os%set_state(2, 1)
        call os%set_state(3, 1)
        call os%set_state(4, 0)
        call os%set(1, 'shincarg', 0.)
        call os%set(2, 'shincarg', 3.)
        call os%set(3, 'shincarg', 3.)
        call os%set(4, 'shincarg', 100.)
        call assert_true(cart_shift_prior_sigma2(os, 4._dp) == 3._dp, &
            &'Cartesian variance must retain zero increments and exclude inactive records')
        call os%set(2, 'shincarg', 0.)
        call os%set(3, 'shincarg', 0.)
        call assert_true(cart_shift_prior_sigma2(os, 8._dp) == 4._dp, &
            &'zero history must use the low-pass floor')
        call os%kill
    end subroutine test_shift_prior_width

    !> Three Angstroms spans 2, 1.5 and 1 native pixels at these samplings, independently of cropping.
    subroutine test_shift_prior_angstrom_units()
        real,     parameter :: SAMPLINGS(3) = [1.5, 2., 3.]
        real(dp), parameter :: VARIANCES(3) = [4._dp, 2.25_dp, 1._dp]
        integer,  parameter :: CROPS(2) = [BOX, BOX/2]
        class(builder),    allocatable, target :: b
        class(parameters), allocatable, target :: p
        integer :: i, j
        allocate(b, p)
        p%l_cart_refine = .true.
        p%cont_prior = 'shift'
        p%cont_prior_shift_sigma = 3.
        p%box = BOX
        p%lp = 12.
        do i = 1, size(SAMPLINGS)
            p%smpd = SAMPLINGS(i)
            do j = 1, size(CROPS)
                p%box_crop = CROPS(j)
                p%smpd_crop = p%smpd*real(BOX)/real(CROPS(j))
                call prep_strategy3D(p, b)
                call assert_true(s3D%cont_shift_prior_sigma2 == VARIANCES(i), &
                    &'fixed Angstrom width must convert once using native sampling, independently of crop')
                call assert_real(3., p%cont_prior_shift_sigma, 0., 'preparation must retain the command-line width in Angstroms')
            end do
        end do
        call b%spproj%os_ptcl3D%new(1, is_ptcl=.true.)
        b%spproj_field => b%spproj%os_ptcl3D
        call b%spproj_field%set_state(1, 1)
        call b%spproj_field%set(1, 'shincarg', 0.)
        p%smpd = 1.5
        p%box_crop = BOX
        p%smpd_crop = p%smpd
        p%cont_prior_shift_sigma = 0.
        call prep_strategy3D(p, b)
        call assert_true(s3D%cont_shift_prior_sigma2 == 4._dp, &
            &'automatic width keeps the native-pixel floor for a 12 A band at 1.5 A per pixel')
        p%cont_prior_shift_sigma = 1.5
        call prep_strategy3D(p, b)
        call assert_true(s3D%cont_shift_prior_sigma2 == 1._dp, &
            &'fixed 1.5 A is one native pixel here and must not inherit the automatic floor')
        p%cont_prior = 'none'
        call prep_strategy3D(p, b)
        call assert_true(s3D%cont_shift_prior_sigma2 == 0._dp, 'disabled prior clears the converted variance')
        call clean_strategy3D()
        call b%spproj%os_ptcl3D%kill
        nullify(b%spproj_field)
    end subroutine test_shift_prior_angstrom_units

    subroutine test_invalid_prior_width()
        type(invalid_prior_width_child) :: child
        class(cmdline), allocatable :: cline
        integer :: i
        allocate(cline)
        do i = 1, 3
            call cline%set('which_iter', i)
            call child%start(name=string('INVALID_PRIOR_WIDTH'), cline=cline)
            call child%await_final_status()
            call assert_int(FORK_STATUS_FAILED, child%status(), 'invalid active width increment must fail its worker')
            call assert_int(1, child%get_exit_code(), 'invalid width must fail through the reported error path')
            call child%destroy
        end do
        call cline%kill()
    end subroutine test_invalid_prior_width

    subroutine reject_invalid_prior_width( self, cline )
        class(invalid_prior_width_child), intent(inout) :: self
        class(cmdline), intent(inout) :: cline
        type(oris) :: os
        real(dp) :: variance
        real :: increment
        increment = -1.
        select case(cline%get_iarg('which_iter'))
            case(2)
                increment = ieee_value(1., ieee_quiet_nan)
            case(3)
                increment = ieee_value(1., ieee_positive_inf)
        end select
        call os%new(1, is_ptcl=.true.)
        call os%set_state(1, 1)
        call os%set(1, 'shincarg', increment)
        variance = cart_shift_prior_sigma2(os, 4._dp)
        write(*,'(a,es14.6)') 'Unexpected valid width: ', variance
        call os%kill
    end subroutine reject_invalid_prior_width

    !> Relay a real invalid-width worker exit through the production queue status boundary.
    !! The master must stop before merge; a watchdog makes a missed failure a bounded test.
    subroutine test_prior_master_failure()
        type(invalid_prior_width_child) :: worker
        type(prior_failure_master_child) :: master
        type(string) :: saved, root
        class(cmdline), allocatable :: cline
        integer :: nfail0, unit, i, state
        logical :: merged, timed_out
        nfail0 = tests_failed
        call enter_fixture('prior_master_failure', saved, root)
        allocate(cline)
        call cline%set('which_iter', 1)
        call worker%start(name=string('PRIOR_FAILED_PART'), cline=cline)
        call worker%await_final_status()
        call assert_int(FORK_STATUS_FAILED, worker%status(), 'invalid prior part must fail')
        open(newunit=unit, file='EXIT_CODE_JOB_1', status='replace', action='write')
        write(unit,'(i0)') worker%get_exit_code()
        close(unit)
        call worker%destroy
        call master%start(name=string('PRIOR_FAILED_MASTER'))
        timed_out = .true.
        do i = 1, 10
            state = master%status()
            if( state /= FORK_STATUS_RUNNING )then
                timed_out = .false.
                exit
            endif
            call sleep(1)
        end do
        if( timed_out ) call master%kill()
        call master%await_final_status()
        call assert_true(.not. timed_out, 'master must detect the failed prior part within ten seconds')
        call assert_int(FORK_STATUS_FAILED, master%status(), 'queue master must fail on the worker error')
        inquire(file='PRIOR_MERGE_REACHED', exist=merged)
        call assert_true(.not. merged, 'master must stop before merge after invalid prior input')
        call master%destroy
        call cline%kill()
        call leave_fixture(saved, root, nfail0)
    end subroutine test_prior_master_failure

    subroutine reject_failed_prior_part( self, cline )
        use simple_qsys_ctrl, only: qsys_ctrl
        use simple_qsys_local, only: qsys_local
        class(prior_failure_master_child), intent(inout) :: self
        class(cmdline), intent(inout) :: cline
        class(qsys_ctrl), allocatable :: ctrl
        class(qsys_local), allocatable, target :: queue
        integer, allocatable, target :: parts(:,:)
        logical, allocatable :: done(:), submitted(:)
        integer :: unit
        allocate(ctrl, queue, parts(1,2))
        parts(1,:) = [1,1]
        call queue%new()
        call ctrl%new(string('false'), queue, parts, [1,1], 1, .false.)
        call ctrl%get_jobs_status(done, submitted)
        done = .false.
        submitted = .true.
        call ctrl%set_jobs_status(done, submitted)
        call ctrl%schedule_jobs()
        open(newunit=unit, file='PRIOR_MERGE_REACHED', status='replace')
        close(unit)
        call ctrl%kill()
        call queue%kill()
    end subroutine reject_failed_prior_part

    !> Exercise the production strategy with prior on and independently recompute its data score.
    subroutine test_prior_strategy_commit()
        class(builder), allocatable, target :: b
        class(parameters), allocatable, target :: p
        type(ori) :: result
        real(dp) :: objective, gradient(5)
        integer :: slot
        allocate(b, p)
        call build_fixture(b, p, OBJFUN_EUCLID)
        p%cont_prior = 'pose'
        p%cont_prior_rot_sigma = 5.
        s3D%cont_shift_prior_sigma2 = 4._dp
        call run_particle(b, p, 3)
        call b%spproj_field%get_ori(3, result)
        slot = b%cftc%get_ptcl_slot(3)
        call objective_gradient(b%cftc, 1, .false., slot, real(result%get_mat(),dp), &
            &real(result%get_2Dshift()-[-0.08,0.06],dp), objective, gradient)
        call assert_true(abs(result%get('corr_cart')-real(b%cftc%score(slot,objective))) < 1.e-4, &
            &'prior strategy must publish data-only score at the returned pose')
        call result%kill
        call kill_fixture(b)
        call clean_strategy3D()
    end subroutine test_prior_strategy_commit

    ! FIXTURE

    subroutine build_volume( volume )
        real, allocatable, intent(out) :: volume(:,:,:)
        integer :: i, j, k
        allocate(volume(BOX,BOX,BOX))
        do k = 1, BOX
            do j = 1, BOX
                do i = 1, BOX
                    volume(i,j,k) = sin(0.11*real(i) + 0.17*real(j) - 0.07*real(k)) + 0.02*real(i*j - k) + 0.003*real(i*k)
                end do
            end do
        end do
    end subroutine build_volume

    real function truth_e3( i )
        integer, intent(in) :: i
        truth_e3 = 28. + 3.*real(i)
    end function truth_e3

    !> Parameters, a builder with NPTCLS seeds perturbed from their truth poses, and the
    !! Cartesian calculator with every particle's slot: particle i in slot NPTCLS+1-i, holding
    !! the truth prediction centred on the stored shift. Particle 2 has state zero.
    subroutine build_fixture( b, p, objfun )
        class(builder),    target, intent(inout) :: b
        class(parameters), target, intent(inout) :: p
        integer,                  intent(in)    :: objfun
        type(ctfparams)      :: no_ctf
        real, allocatable    :: volume(:,:,:)
        complex, allocatable :: observed(:,:)
        real    :: seed_shift(2), truth_shift(2)
        integer :: i, pinds(NPTCLS)
        p%oritype          = 'ptcl3D'
        p%cc_objfun        = objfun
        p%inpl_cont        = 'no'
        p%box              = BOX
        p%box_crop         = BOX
        p%trs              = 5.
        p%athres_cont      = 15.
        ! the polar searches' threshold, deliberately narrow: test_commit_contract recovers a 1.8
        ! degree seed error only if the Cartesian solve ignores it
        p%athres           = 0.5
        p%kfromto          = KFROMTO
        p%fromp            = 1
        p%top              = NPTCLS
        call b%spproj%os_ptcl3D%new(NPTCLS, is_ptcl=.true.)
        b%spproj_field => b%spproj%os_ptcl3D
        call b%pgrpsyms%new('c1')
        call build_volume(volume)
        call b%cftc%new(1, BOX, 1)
        call b%cftc%set_ref(1, .true.,  volume)
        call b%cftc%set_ref(1, .false., volume)
        call b%cftc%new_ptcls(NPTCLS)
        call b%cftc%set_sigma2(1, spread(spread(1., 1, BOX/2), 2, 2))
        allocate(observed(-BOX/2:BOX/2,-BOX/2:BOX/2))
        no_ctf%ctfflag = CTFFLAG_NO
        do i = 1, NPTCLS
            truth_shift = [0.31, -0.24] + 0.05*real(i)
            seed_shift  = [-0.08, 0.06]
            call b%spproj_field%set_euler(i, [20., 36.2, truth_e3(i) + 0.7])
            call b%spproj_field%set_shift(i, seed_shift)
            call b%spproj_field%set_state(i, merge(0, 1, i == 2))
            call b%spproj_field%set(i, 'eo',   real(mod(i,2)))
            call b%spproj_field%set(i, 'proj', 3.)
            call b%spproj_field%set(i, 'corr', 0.42)
            call b%cftc%predict(1, .true., real(euler2m([19., 37., truth_e3(i)]), dp), &
                &real(truth_shift - seed_shift, dp), observed)
            pinds(NPTCLS + 1 - i) = i
            if( objfun == OBJFUN_CC )then
                call b%cftc%set_ptcl(NPTCLS + 1 - i, observed, no_ctf, KFROMTO)
            else
                call b%cftc%set_ptcl(NPTCLS + 1 - i, observed, no_ctf, mod(i,2) + 1, KFROMTO)
            endif
        end do
        call b%cftc%set_ptcl_inds(pinds)
        if( objfun == OBJFUN_EUCLID )then
            ! the records the residual replaces, read from a unit store as a matcher reads them
            call write_test_sigma2_store('tmp_strategy3D_cont_tester_sigma2.bin', NPTCLS, BOX)
            call b%esig%new(p, string('tmp_strategy3D_cont_tester_sigma2.bin'), BOX)
            call b%esig%read_part
        endif
    end subroutine build_fixture

    subroutine kill_fixture( b )
        class(builder), intent(inout) :: b
        call b%cftc%kill
        call b%esig%kill
        call del_file('tmp_strategy3D_cont_tester_sigma2.bin')
        call b%pgrpsyms%kill
        call b%spproj%os_ptcl3D%kill
        nullify(b%spproj_field)
    end subroutine kill_fixture

    !> One particle through a class(strategy3D) pointer, as the matcher runs it.
    subroutine run_particle( b, p, iptcl )
        class(builder),    target, intent(inout) :: b
        class(parameters), target, intent(inout) :: p
        integer,                  intent(in)    :: iptcl
        class(strategy3D), pointer :: strat
        type(strategy3D_spec) :: spec
        spec%iptcl     = iptcl
        spec%iptcl_map = iptcl
        allocate(strategy3D_cont :: strat)
        call strat%new(p, spec, b)
        call strat%srch(b%spproj_field, 1)
        call strat%kill
        deallocate(strat)
    end subroutine run_particle

    pure function rotation_distance( a, c ) result( d )
        real(dp), intent(in) :: a(3,3), c(3,3)
        real(dp) :: d
        d = acos(max(-1._dp, min(1._dp, 0.5_dp*(sum(a*c) - 1._dp))))
    end function rotation_distance

    ! TESTS
    subroutine test_reference_geometry()
        use, intrinsic :: ieee_arithmetic, only: ieee_value, ieee_quiet_nan
        integer :: dims(3,3)
        real :: sampling(3)
        dims = 32
        sampling = 1.5
        call assert_true(cont_reference_geometry_is_valid(dims, sampling, 32, 1.5), 'valid reference geometry rejected')
        ! Both grids span 48 Angstrom: 32 * 1.5 = 16 * 3.
        dims = 16
        sampling = 3.
        call assert_true(cont_reference_geometry_is_valid(dims, sampling, 32, 1.5), 'compatible coarse reference rejected')
        dims(1,2) = 15
        call assert_false(cont_reference_geometry_is_valid(dims, sampling, 32, 1.5), 'mismatched half dimensions accepted')
        dims = 16
        sampling(3) = 2.
        call assert_false(cont_reference_geometry_is_valid(dims, sampling, 32, 1.5), 'mismatched half sampling accepted')
        sampling = 2.
        call assert_false(cont_reference_geometry_is_valid(dims, sampling, 32, 1.5), 'wrong physical extent accepted')
        sampling = 0.
        call assert_false(cont_reference_geometry_is_valid(dims, sampling, 32, 1.5), 'zero sampling accepted')
        sampling = 1.5
        sampling(2) = ieee_value(0., ieee_quiet_nan)
        call assert_false(cont_reference_geometry_is_valid(dims, sampling, 32, 1.5), 'NaN sampling accepted')
        sampling = 1.5
        dims = 0
        call assert_false(cont_reference_geometry_is_valid(dims, sampling, 32, 1.5), 'empty map geometry accepted')
    end subroutine test_reference_geometry

    ! An invalid row outside the selected state must not reject that state's continuation.
    subroutine test_selected_seeds()
        class(oris), allocatable :: particles
        type(ori) :: seed
        allocate(particles)
        call seed%new(is_ptcl=.true.)
        call particles%new(2, is_ptcl=.true.)
        call seed%set_euler([0., 0., 0.])
        call seed%set_shift([0., 0.])
        call seed%set('state', 1.)
        call seed%set('eo', 0.)
        call seed%set('proj', 1.)
        call particles%set_ori(1, seed)
        call seed%set('state', 2.)
        call seed%set('proj', 0.)
        call particles%set_ori(2, seed)
        call check_cont_seeds(particles, [.true., .false.])
        call assert_true(particles%get(2, 'proj') == 0., 'seed validation changed an excluded row')
        call seed%kill()
        call particles%kill()
    end subroutine test_selected_seeds

    ! Both binary labels are valid in packed particle records, including the zero-valued even half.
    subroutine test_particle_seed_contract()
        type(ori) :: seed
        call seed%new(is_ptcl=.true.)
        call seed%set_euler([0., 0., 0.])
        call seed%set_shift([0., 0.])
        call seed%set('state', 1.)
        call seed%set('proj', 1.)
        call seed%set('eo', 0.)
        call assert_true(cont_seed_is_valid(seed), 'a packed even-half seed was rejected')
        call seed%set('eo', 1.)
        call assert_true(cont_seed_is_valid(seed), 'a packed odd-half seed was rejected')
        call seed%set('eo', 0.5)
        call assert_false(cont_seed_is_valid(seed), 'a fractional packed half label was accepted')
        call seed%set('eo', -1.)
        call assert_false(cont_seed_is_valid(seed), 'an unassigned packed half label was accepted')
        call seed%kill()
    end subroutine test_particle_seed_contract

    ! E19 (seed validity, moved from the pose_cont strategy tester): the identity is a valid
    ! initialized pose; explicit state, half and projection, not nonzero Euler angles, define it
    subroutine test_seed_contract()
        type(ori) :: seed
        call seed%set_euler([0., 0., 0.])
        call seed%set_shift([1.25, -0.75])
        call seed%set('state', 1.)
        call seed%set('eo', 0.)
        call seed%set('proj', 1.)
        call assert_true(cont_seed_is_valid(seed), 'a valid identity seed was rejected')
        call seed%set_shift([0., 0.])
        call assert_true(cont_seed_is_valid(seed), 'zero shifts were rejected')
        call seed%delete_entry('e1')
        call assert_false(cont_seed_is_valid(seed), 'a missing Euler angle was accepted')
        call seed%set('e1', 0.)
        call seed%set('proj', 0.)
        call assert_false(cont_seed_is_valid(seed), 'a seed without projection direction was accepted')
        call seed%set('proj', 1.)
        call seed%delete_entry('eo')
        call assert_false(cont_seed_is_valid(seed), 'a missing half-set was accepted as even')
        call seed%set('eo', 0.5)
        call assert_false(cont_seed_is_valid(seed), 'a fractional half-set was rounded and accepted')
        call seed%set('eo', 2.)
        call assert_false(cont_seed_is_valid(seed), 'a seed without half-set was accepted')
        call seed%kill()
    end subroutine test_seed_contract

    ! N15: lifecycle through a class(strategy3D) pointer, repeated, with build%pftc never
    ! constructed; an identity seed is accepted as a seed; a state-zero particle is rejected
    ! and keeps state zero without a score
    subroutine test_lifecycle_and_rejection()
        class(builder),    allocatable, target :: b
        class(parameters), allocatable, target :: p
        integer :: icycle
        allocate(b, p)
        do icycle = 1, 2
            call build_fixture(b, p, OBJFUN_CC)
            call assert_false(b%pftc%exists(), 'the fixture constructed a polar calculator')
            call b%spproj_field%set_euler(1, [0., 0., 0.])
            call run_particle(b, p, 1)
            call assert_true(b%spproj_field%get(1, 'corr_cart') > 0., 'an identity seed was not refined')
            call run_particle(b, p, 2)
            call assert_true(b%spproj_field%get_state(2) == 0 .and. b%spproj_field%get(2, 'corr_cart') == 0., &
                &'a state-zero particle was not rejected')
            call assert_false(b%pftc%exists(), 'the Cartesian strategy constructed a polar calculator')
            call kill_fixture(b)
        end do
        deallocate(b, p)
    end subroutine test_lifecycle_and_rejection

    ! N16: an accepted solve commits the pose, corr_cart = cc at the result (the calculator's
    ! score of the objective evaluated in the test), the improved flag and the convergence
    ! fields of the seed-to-result pose update; a rejected solve (a particle whose observation is
    ! the noise-free prediction at its seed, so no pose improves on the seed) leaves the pose
    ! bit-identical and writes corr_cart at the seed; corr is unchanged in both. Under euclid
    ! the sigma owner receives the residual at the committed pose.
    subroutine test_commit_contract()
        class(builder),    allocatable, target :: b
        class(parameters), allocatable, target :: p
        type(ori) :: before, after
        real(dp)  :: objective, gradient(5), truth(3,3), expected_score
        real      :: native_before(2), sigma_shells_sum
        real, allocatable :: stored(:)
        complex, allocatable :: observed(:,:)
        type(ctfparams) :: no_ctf
        integer   :: islot, iobj
        integer, parameter :: objfuns(2) = [OBJFUN_CC, OBJFUN_EUCLID]
        allocate(b, p)
        do iobj = 1, 2
            call build_fixture(b, p, objfuns(iobj))
            ! accepted
            call b%spproj_field%get_ori(3, before)
            call run_particle(b, p, 3)
            call b%spproj_field%get_ori(3, after)
            truth = real(euler2m([19., 37., truth_e3(3)]), dp)
            call assert_true(after%get('pose_cont_improved') == 1., 'an improving solve was not flagged improved')
            call assert_true(after%get('pose_cont_status') == real(CARTFT_ACCEPTED), 'an accepted solve did not record its status')
            call assert_true(rotation_distance(real(after%get_mat(), dp), truth) < ROTATION_TOL .and. &
                &all(abs(after%get_2Dshift() - ([0.31, -0.24] + 0.15)) < SHIFT_TOL), 'an accepted solve did not commit the pose')
            islot = b%cftc%get_ptcl_slot(3)
            call objective_gradient(b%cftc, 1, .false., islot, real(after%get_mat(), dp), &
                &real(after%get_2Dshift() - before%get_2Dshift(), dp), objective, gradient)
            expected_score = b%cftc%score(islot, objective)
            call assert_true(abs(after%get('corr_cart') - real(expected_score)) < 1.e-4, &
                &'corr_cart is not the score at the committed pose')
            call assert_true(after%get('corr') == 0.42, 'a Cartesian pass wrote corr')
            call assert_true(after%get('dist') > 0. .and. after%get('frac') == 100. .and. &
                &after%get('mi_proj') == merge(1., 0., after%get('dist') <= p%angthres_mi_proj) .and. &
                &abs(after%get('shincarg') - norm2(after%get_2Dshift() - before%get_2Dshift())) < 1.e-5, &
                &'the convergence fields are not the seed-to-result pose update')
            if( objfuns(iobj) == OBJFUN_EUCLID )then
                stored = b%esig%get_sigma2_part(3)
                sigma_shells_sum = sum(stored)
                call assert_true(sigma_shells_sum > 0., 'the sigma owner received no residual under euclid')
            endif
            ! rejected: the seed is the optimum of its observation
            call b%spproj_field%get_ori(4, before)
            allocate(observed(-BOX/2:BOX/2,-BOX/2:BOX/2))
            no_ctf%ctfflag = CTFFLAG_NO
            islot = b%cftc%get_ptcl_slot(4)
            call b%cftc%predict(1, .true., real(before%get_mat(), dp), [0._dp, 0._dp], observed)
            if( objfuns(iobj) == OBJFUN_CC )then
                call b%cftc%set_ptcl(islot, observed, no_ctf, KFROMTO)
            else
                call b%cftc%set_ptcl(islot, observed, no_ctf, mod(4,2) + 1, KFROMTO)
            endif
            deallocate(observed)
            native_before = before%get_2Dshift()
            call run_particle(b, p, 4)
            call b%spproj_field%get_ori(4, after)
            call assert_true(after%get('pose_cont_improved') == 0. .and. all(after%get_euler() == before%get_euler()) .and. &
                &all(after%get_2Dshift() == native_before), 'a rejected solve changed the pose')
            ! the rejection outcomes are the codes above CARTFT_ACCEPTED
            call assert_true(after%get('pose_cont_status') > real(CARTFT_ACCEPTED), 'a rejected solve did not record its status')
            islot = b%cftc%get_ptcl_slot(4)
            call objective_gradient(b%cftc, 1, .true., islot, real(before%get_mat(), dp), [0._dp, 0._dp], objective, gradient)
            call assert_true(abs(after%get('corr_cart') - real(b%cftc%score(islot, objective))) < 1.e-4, &
                &'a rejected solve did not write corr_cart at the seed')
            call assert_true(after%get('corr') == 0.42, 'a rejected Cartesian solve wrote corr')
            call kill_fixture(b)
        end do
        call before%kill
        call after%kill
        deallocate(b, p)
    end subroutine test_commit_contract

    ! N27: the batch run through srch at every processor, 20 times, gives the poses, scores and
    ! flags of the run at one thread bit for bit. Expected value: the one-thread run.
    subroutine test_threaded_batch()
        integer, parameter :: NREPS = 20
        class(builder),    allocatable, target :: b
        class(parameters), allocatable, target :: p
        real    :: single(9,NPTCLS), team(9,NPTCLS)
        integer :: irep, ndiff, nthr_prev, nthr_full
        allocate(b, p)
        nthr_prev = 1
        nthr_full = 1
        !$ nthr_prev = omp_get_max_threads()
        !$ nthr_full = omp_get_num_procs()
        !$ call omp_set_num_threads(1)
        call run_batch(single)
        !$ call omp_set_num_threads(nthr_full)
        ndiff = 0
        do irep = 1, NREPS
            call run_batch(team)
            ndiff = ndiff + count(transfer(team, 0, size(team)) /= transfer(single, 0, size(single)))
        end do
        !$ call omp_set_num_threads(nthr_prev)
        call assert_int(0, ndiff, 'a batch at every processor differs from the run at one thread')
        deallocate(b, p)

      contains

        !> a fresh fixture's batch through srch, then its record
        subroutine run_batch( values )
            real, intent(out) :: values(9,NPTCLS)
            integer :: i
            call build_fixture(b, p, OBJFUN_CC)
            !$omp parallel do default(shared) private(i) schedule(static,1)
            do i = 1, NPTCLS
                call run_particle(b, p, i)
            end do
            !$omp end parallel do
            call record(values)
            call kill_fixture(b)
        end subroutine run_batch

        subroutine record( values )
            real, intent(out) :: values(9,NPTCLS)
            integer :: j
            do j = 1, NPTCLS
                values(1:3,j) = b%spproj_field%get_euler(j)
                values(4:5,j) = b%spproj_field%get_2Dshift(j)
                values(6,j)   = b%spproj_field%get(j, 'corr_cart')
                values(7,j)   = b%spproj_field%get(j, 'pose_cont_improved')
                values(8,j)   = b%spproj_field%get(j, 'dist')
                values(9,j)   = real(b%spproj_field%get_state(j))
            end do
        end subroutine record

    end subroutine test_threaded_batch

    ! N22: a pass over a sample (particles 3 and 5) leaves a particle outside it (4) bit for
    ! bit: pose, corr, corr_cart and flags, the whole record. Expected value: the input record.
    subroutine test_outside_sample_untouched()
        class(builder),    allocatable, target :: b
        class(parameters), allocatable, target :: p
        type(ori) :: o
        real      :: before(N_PTCL_ORIPARAMS), after(N_PTCL_ORIPARAMS)
        allocate(b, p)
        call build_fixture(b, p, OBJFUN_CC)
        call b%spproj_field%set(4, 'corr_cart', 0.123)
        call b%spproj_field%set(4, 'pose_cont_improved', 1.)
        call b%spproj_field%get_ori(4, o)
        call o%ori2prec(before)
        call run_particle(b, p, 3)
        call run_particle(b, p, 5)
        call b%spproj_field%get_ori(4, o)
        call o%ori2prec(after)
        call assert_true(all(after == before), 'a particle outside the pass sample changed')
        call assert_true(b%spproj_field%get(3, 'pose_cont_improved') == 1., 'the sampled particle was not refined')
        call o%kill
        call kill_fixture(b)
        deallocate(b, p)
    end subroutine test_outside_sample_untouched

    ! N23: with the convergence fields of a discrete search in the record, a continuous followup pass
    ! (l_cont_followup) commits pose, corr_cart and the flag and leaves dist, dist_inpl,
    ! shincarg, mi_proj, mi_state and frac as they were; a refine=cont pass writes them from
    ! the seed-to-result pose update. Expected values: the input record; the pose update computed here.
    subroutine test_followup_keeps_convergence_fields()
        character(len=9), parameter :: KEYS(6) = [character(len=9) :: 'dist', 'dist_inpl', 'shincarg', &
            &'mi_proj', 'mi_state', 'frac']
        real,             parameter :: DISCRETE(6) = [7.5, 3.25, 1.5, 0., 0., 55.]
        class(builder),    allocatable, target :: b
        class(parameters), allocatable, target :: p
        type(ori) :: before, after
        integer   :: i, k
        allocate(b, p)
        call build_fixture(b, p, OBJFUN_CC)
        do i = 3, 5, 2
            do k = 1, size(KEYS)
                call b%spproj_field%set(i, trim(KEYS(k)), DISCRETE(k))
            end do
        end do
        ! the continuous followup pass on particle 3
        p%pose_cont       = 'yes'
        p%l_cont_followup   = .true.
        call b%spproj_field%get_ori(3, before)
        call run_particle(b, p, 3)
        call b%spproj_field%get_ori(3, after)
        call assert_true(after%get('pose_cont_improved') == 1. .and. any(after%get_euler() /= before%get_euler()), &
            &'the continuous followup pass did not commit the refined pose')
        call assert_true(after%get('corr') == before%get('corr') .and. after%get('corr_cart') > 0., &
            &'the continuous followup pass did not write corr_cart alone')
        do k = 1, size(KEYS)
            call assert_true(after%get(trim(KEYS(k))) == DISCRETE(k), &
                &'the continuous followup pass changed the discrete '//trim(KEYS(k)))
        end do
        ! a refine=cont pass on particle 5
        p%pose_cont       = 'no'
        p%l_cont_followup   = .false.
        call b%spproj_field%get_ori(5, before)
        call run_particle(b, p, 5)
        call b%spproj_field%get_ori(5, after)
        call assert_true(abs(after%get('dist') - rad2deg(before.euldist.after)) < 1.e-3 .and. &
            &abs(after%get('shincarg') - norm2(after%get_2Dshift() - before%get_2Dshift())) < 1.e-5 .and. &
            &after%get('frac') == 100. .and. after%get('mi_state') == 1., &
            &'a refine=cont pass did not write the seed-to-result pose update')
        call before%kill
        call after%kill
        call kill_fixture(b)
        deallocate(b, p)
    end subroutine test_followup_keeps_convergence_fields

    ! The rotation bound of a Cartesian pass is athres_cont, not the polar athres (the fixture sets
    ! that to 0.5 degrees): with athres_cont = 0.5 the committed rotation stays within 0.5 degrees
    ! of the seed and the truth 1.8 degrees away is not reached.
    subroutine test_rotation_bound()
        real(dp), parameter :: BOUND = 0.5_dp*acos(-1._dp)/180._dp
        class(builder),    allocatable, target :: b
        class(parameters), allocatable, target :: p
        type(ori) :: before, after
        real(dp)  :: truth(3,3)
        allocate(b, p)
        call build_fixture(b, p, OBJFUN_CC)
        p%athres_cont = 0.5
        call b%spproj_field%get_ori(3, before)
        call run_particle(b, p, 3)
        call b%spproj_field%get_ori(3, after)
        truth = real(euler2m([19., 37., truth_e3(3)]), dp)
        call assert_true(rotation_distance(real(after%get_mat(), dp), real(before%get_mat(), dp)) <= BOUND + ROTATION_TOL, &
            &'the committed rotation exceeds athres_cont')
        call assert_true(rotation_distance(real(after%get_mat(), dp), truth) > BOUND, &
            &'a solve bounded by athres_cont reached a truth outside the bound')
        call kill_fixture(b)
        call before%kill
        call after%kill
        deallocate(b, p)
    end subroutine test_rotation_bound

end module simple_strategy3D_cont_tester
