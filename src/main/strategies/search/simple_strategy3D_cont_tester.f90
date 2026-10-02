!@descr: unit tests for the Cartesian 3D strategy of a Cartesian pass (simple_strategy3D_cont)
! A small fixture with no file and no polar calculator: a one-state Cartesian calculator in a
! builder whose ptcl3D field holds the seeds, each particle's slot holding the noise-free truth
! prediction centred on its stored shift (as the batch preparation centres it). Pinned: the seed
! validity contract (moved from the pose_cont strategy tester with E19), the lifecycle through a
! class(strategy3D) pointer with build%pftc never constructed, an identity seed accepted and a
! state-zero particle rejected (N15); an accepted solve committing pose, corr_cart, the improved
! flag and the convergence fields of the seed-to-result motion, and a rejected solve leaving the
! pose bit-identical with corr_cart written at the seed, corr unchanged in both (N16); under
! euclid the sigma owner receiving the residual at the committed pose; and a batch run by a team
! of three threads giving the serial poses, scores and flags (N27); a particle outside the
! pass's sample left bit-identical (N22); a polish pass (pose_cont=yes) leaving the convergence
! fields of the discrete search while a refine=cont pass writes them from the seed-to-result
! motion (N23); the rotation bound taken from athres_cont, not the polar athres.
module simple_strategy3D_cont_tester
use simple_core_module_api, only: dp, euler2m, string, rad2deg
use simple_defs_ori,        only: N_PTCL_ORIPARAMS
use simple_builder,         only: builder
use simple_parameters,      only: parameters
use simple_ori,             only: ori
use simple_cartft_calc,     only: cartft_calc
use simple_strategy3D,      only: strategy3D
use simple_strategy3D_srch, only: strategy3D_spec
use simple_strategy3D_cont, only: strategy3D_cont, cont_seed_is_valid
use simple_type_defs,       only: ctfparams, CTFFLAG_NO, OBJFUN_CC, OBJFUN_EUCLID
use simple_test_utils
implicit none
private
public :: run_all_strategy3D_cont_tests

integer, parameter :: BOX = 16, NPTCLS = 6
integer, parameter :: KFROMTO(2) = [2, BOX/2 - 1]
! single-precision Euler angles of the stored poses limit the recovered rotation
real(dp), parameter :: ROTATION_TOL = 2.e-3_dp
real,     parameter :: SHIFT_TOL    = 5.e-3

contains

    subroutine run_all_strategy3D_cont_tests()
        write(*,'(A)') '**** running all Cartesian 3D strategy tests ****'
        write(*,'(A)') 'test_seed_contract'
        call test_seed_contract()
        write(*,'(A)') 'test_lifecycle_and_rejection'
        call test_lifecycle_and_rejection()
        write(*,'(A)') 'test_commit_contract'
        call test_commit_contract()
        write(*,'(A)') 'test_threaded_batch'
        call test_threaded_batch()
        write(*,'(A)') 'test_outside_sample_untouched'
        call test_outside_sample_untouched()
        write(*,'(A)') 'test_polish_keeps_convergence_fields'
        call test_polish_keeps_convergence_fields()
        write(*,'(A)') 'test_rotation_bound'
        call test_rotation_bound()
    end subroutine run_all_strategy3D_cont_tests

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
        type(builder),    target, intent(inout) :: b
        type(parameters), target, intent(inout) :: p
        integer,                  intent(in)    :: objfun
        type(ctfparams)      :: no_ctf
        real, allocatable    :: volume(:,:,:), sigma2(:)
        complex, allocatable :: observed(:,:)
        real    :: seed_shift(2), truth_shift(2)
        integer :: i, pinds(NPTCLS)
        p%oritype          = 'ptcl3D'
        p%cc_objfun        = objfun
        p%inpl_cont        = 'no'
        p%projrec          = 'no'
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
        allocate(sigma2(0:BOX/2), source=1.)
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
                call b%cftc%set_ptcl(NPTCLS + 1 - i, observed, no_ctf, sigma2, KFROMTO)
            endif
        end do
        call b%cftc%set_ptcl_inds(pinds)
        if( objfun == OBJFUN_EUCLID )then
            call b%esig%new(p, string('tmp_strategy3D_cont_tester_sigma2.bin'), BOX)
            call b%esig%allocate_ptcls
        endif
    end subroutine build_fixture

    subroutine kill_fixture( b )
        type(builder), intent(inout) :: b
        call b%cftc%kill
        call b%esig%kill
        call b%pgrpsyms%kill
        call b%spproj%os_ptcl3D%kill
        nullify(b%spproj_field)
    end subroutine kill_fixture

    !> One particle through a class(strategy3D) pointer, as the matcher runs it.
    subroutine run_particle( b, p, iptcl )
        type(builder),    target, intent(inout) :: b
        type(parameters), target, intent(inout) :: p
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
        call seed%set('proj', 0.)
        call assert_false(cont_seed_is_valid(seed), 'a seed without projection direction was accepted')
        call seed%set('proj', 1.)
        call seed%set('eo', 2.)
        call assert_false(cont_seed_is_valid(seed), 'a seed without half-set was accepted')
        call seed%kill()
    end subroutine test_seed_contract

    ! N15: lifecycle through a class(strategy3D) pointer, repeated, with build%pftc never
    ! constructed; an identity seed is accepted as a seed; a state-zero particle is rejected
    ! and keeps state zero without a score
    subroutine test_lifecycle_and_rejection()
        type(builder),    target :: b
        type(parameters), target :: p
        integer :: icycle
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
    end subroutine test_lifecycle_and_rejection

    ! N16: an accepted solve commits the pose, corr_cart = cc at the result (the calculator's
    ! score of the objective evaluated in the test), the improved flag and the convergence
    ! fields of the seed-to-result motion; a rejected solve (a particle whose observation is
    ! the noise-free prediction at its seed, so no pose improves on the seed) leaves the pose
    ! bit-identical and writes corr_cart at the seed; corr is unchanged in both. Under euclid
    ! the sigma owner receives the residual at the committed pose.
    subroutine test_commit_contract()
        type(builder),    target :: b
        type(parameters), target :: p
        type(ori) :: before, after
        real(dp)  :: objective, gradient(5), truth(3,3), expected_score
        real      :: native_before(2), sigma_shells_sum
        real, allocatable :: stored(:), sigma2(:)
        complex, allocatable :: observed(:,:)
        type(ctfparams) :: no_ctf
        integer   :: islot, iobj
        integer, parameter :: objfuns(2) = [OBJFUN_CC, OBJFUN_EUCLID]
        do iobj = 1, 2
            call build_fixture(b, p, objfuns(iobj))
            ! accepted
            call b%spproj_field%get_ori(3, before)
            call run_particle(b, p, 3)
            call b%spproj_field%get_ori(3, after)
            truth = real(euler2m([19., 37., truth_e3(3)]), dp)
            call assert_true(after%get('pose_cont_improved') == 1., 'an improving solve was not flagged improved')
            call assert_true(rotation_distance(real(after%get_mat(), dp), truth) < ROTATION_TOL .and. &
                &all(abs(after%get_2Dshift() - ([0.31, -0.24] + 0.15)) < SHIFT_TOL), 'an accepted solve did not commit the pose')
            islot = b%cftc%get_ptcl_slot(3)
            call b%cftc%objective_gradient(1, .false., islot, real(after%get_mat(), dp), &
                &real(after%get_2Dshift() - before%get_2Dshift(), dp), objective, gradient)
            expected_score = b%cftc%score(islot, objective)
            call assert_true(abs(after%get('corr_cart') - real(expected_score)) < 1.e-4, &
                &'corr_cart is not the score at the committed pose')
            call assert_true(after%get('corr') == 0.42, 'a Cartesian pass wrote corr')
            call assert_true(after%get('dist') > 0. .and. after%get('frac') == 100. .and. &
                &after%get('mi_proj') == merge(1., 0., after%get('dist') <= p%angthres_mi_proj) .and. &
                &abs(after%get('shincarg') - norm2(after%get_2Dshift() - before%get_2Dshift())) < 1.e-5, &
                &'the convergence fields are not the seed-to-result motion')
            if( objfuns(iobj) == OBJFUN_EUCLID )then
                stored = b%esig%get_sigma2_part(3)
                sigma_shells_sum = sum(stored)
                call assert_true(sigma_shells_sum > 0., 'the sigma owner received no residual under euclid')
            endif
            ! rejected: the seed is the optimum of its observation
            call b%spproj_field%get_ori(4, before)
            allocate(observed(-BOX/2:BOX/2,-BOX/2:BOX/2))
            allocate(sigma2(0:BOX/2), source=1.)
            no_ctf%ctfflag = CTFFLAG_NO
            islot = b%cftc%get_ptcl_slot(4)
            call b%cftc%predict(1, .true., real(before%get_mat(), dp), [0._dp, 0._dp], observed)
            if( objfuns(iobj) == OBJFUN_CC )then
                call b%cftc%set_ptcl(islot, observed, no_ctf, KFROMTO)
            else
                call b%cftc%set_ptcl(islot, observed, no_ctf, sigma2, KFROMTO)
            endif
            deallocate(observed, sigma2)
            native_before = before%get_2Dshift()
            call run_particle(b, p, 4)
            call b%spproj_field%get_ori(4, after)
            call assert_true(after%get('pose_cont_improved') == 0. .and. all(after%get_euler() == before%get_euler()) .and. &
                &all(after%get_2Dshift() == native_before), 'a rejected solve changed the pose')
            islot = b%cftc%get_ptcl_slot(4)
            call b%cftc%objective_gradient(1, .true., islot, real(before%get_mat(), dp), [0._dp, 0._dp], objective, gradient)
            call assert_true(abs(after%get('corr_cart') - real(b%cftc%score(islot, objective))) < 1.e-4, &
                &'a rejected solve did not write corr_cart at the seed')
            call assert_true(after%get('corr') == 0.42, 'a rejected Cartesian solve wrote corr')
            call kill_fixture(b)
        end do
        call before%kill
        call after%kill
    end subroutine test_commit_contract

    ! N27: the batch run through srch by a team of three threads gives the poses, scores and
    ! flags of the serial run bit for bit. Expected value: the serial run.
    subroutine test_threaded_batch()
        type(builder),    target :: b
        type(parameters), target :: p
        real    :: serial(9,NPTCLS), team(9,NPTCLS)
        integer :: i
        call build_fixture(b, p, OBJFUN_CC)
        do i = 1, NPTCLS
            call run_particle(b, p, i)
        end do
        call record(serial)
        call kill_fixture(b)
        call build_fixture(b, p, OBJFUN_CC)
        !$omp parallel do num_threads(3) default(shared) private(i) schedule(static,1)
        do i = 1, NPTCLS
            call run_particle(b, p, i)
        end do
        !$omp end parallel do
        call record(team)
        call assert_true(all(team == serial), 'a threaded batch differs from the serial run')
        call kill_fixture(b)

      contains

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
        type(builder),    target :: b
        type(parameters), target :: p
        type(ori) :: o
        real      :: before(N_PTCL_ORIPARAMS), after(N_PTCL_ORIPARAMS)
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
    end subroutine test_outside_sample_untouched

    ! N23: with the convergence fields of a discrete search in the record, a polish pass
    ! (l_cont_polish) commits pose, corr_cart and the flag and leaves dist, dist_inpl,
    ! shincarg, mi_proj, mi_state and frac as they were; a refine=cont pass writes them from
    ! the seed-to-result motion. Expected values: the input record; the motion computed here.
    subroutine test_polish_keeps_convergence_fields()
        character(len=9), parameter :: KEYS(6) = [character(len=9) :: 'dist', 'dist_inpl', 'shincarg', &
            &'mi_proj', 'mi_state', 'frac']
        real,             parameter :: DISCRETE(6) = [7.5, 3.25, 1.5, 0., 0., 55.]
        type(builder),    target :: b
        type(parameters), target :: p
        type(ori) :: before, after
        integer   :: i, k
        call build_fixture(b, p, OBJFUN_CC)
        do i = 3, 5, 2
            do k = 1, size(KEYS)
                call b%spproj_field%set(i, trim(KEYS(k)), DISCRETE(k))
            end do
        end do
        ! the polish pass on particle 3
        p%pose_cont       = 'yes'
        p%l_cont_polish   = .true.
        call b%spproj_field%get_ori(3, before)
        call run_particle(b, p, 3)
        call b%spproj_field%get_ori(3, after)
        call assert_true(after%get('pose_cont_improved') == 1. .and. any(after%get_euler() /= before%get_euler()), &
            &'the polish pass did not commit the refined pose')
        call assert_true(after%get('corr') == before%get('corr') .and. after%get('corr_cart') > 0., &
            &'the polish pass did not write corr_cart alone')
        do k = 1, size(KEYS)
            call assert_true(after%get(trim(KEYS(k))) == DISCRETE(k), &
                &'the polish pass changed the discrete '//trim(KEYS(k)))
        end do
        ! a refine=cont pass on particle 5
        p%pose_cont       = 'no'
        p%l_cont_polish   = .false.
        call b%spproj_field%get_ori(5, before)
        call run_particle(b, p, 5)
        call b%spproj_field%get_ori(5, after)
        call assert_true(abs(after%get('dist') - rad2deg(before.euldist.after)) < 1.e-3 .and. &
            &abs(after%get('shincarg') - norm2(after%get_2Dshift() - before%get_2Dshift())) < 1.e-5 .and. &
            &after%get('frac') == 100. .and. after%get('mi_state') == 1., &
            &'a refine=cont pass did not write the seed-to-result motion')
        call before%kill
        call after%kill
        call kill_fixture(b)
    end subroutine test_polish_keeps_convergence_fields

    ! The rotation bound of a Cartesian pass is athres_cont, not the polar athres (the fixture sets
    ! that to 0.5 degrees): with athres_cont = 0.5 the committed rotation stays within 0.5 degrees
    ! of the seed and the truth 1.8 degrees away is not reached.
    subroutine test_rotation_bound()
        real(dp), parameter :: BOUND = 0.5_dp*acos(-1._dp)/180._dp
        type(builder),    target :: b
        type(parameters), target :: p
        type(ori) :: before, after
        real(dp)  :: truth(3,3)
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
    end subroutine test_rotation_bound

end module simple_strategy3D_cont_tester
