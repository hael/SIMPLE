!@descr: unit tests for the 3D probability tables: the seed shift of an assigned candidate, the dense part files, the in-plane samplers, the posterior draw of the dense and of the neighbourhood table, the n_eff control
! The seed table holds the particle-frame seed every candidate of a particle was scored with; an assigned, unrefined
! candidate records it as it is, whatever its in-plane angle. A refined candidate keeps the shift of its search.
! ref_assign draws every particle's reference from its tempered row posterior on synthetic dense tables of one and
! of several states, with ties, a flat and an invalid particle. A particle's state is the state of its closest
! reference; an exact tie goes to the lowest reference index.
module simple_eul_prob_tab_tester
use simple_core_module_api
use simple_builder,            only: builder
use simple_parameters,         only: parameters
use simple_eul_prob_tab,       only: eul_prob_tab
use simple_eul_prob_tab_neigh, only: eul_prob_tab_neigh
use simple_eul_prob_tab_utils, only: record_seed_shift, &
    &sample_bounded_dists, sample_power_dists, prob_candidate, prob_dense_entry, prob_refined_entry, &
    &posterior_temperature, posterior_row_stats, draw_posterior_row, draw_posterior_extras, prob_beta_at, prob_neff_readout, &
    &prob_insert_table, prob_insert_pose, prob_insert_poses, prob_insert_cone, &
    &PROB_NEFF_MIN_ROWS
use simple_rnd,                only: seed_rnd_fixed
use simple_sym,                only: sym
use simple_oris,               only: oris
use simple_type_defs,          only: ptcl_ref
use simple_test_utils
!$ use omp_lib, only: omp_get_max_threads, omp_get_num_procs, omp_set_num_threads
implicit none
private
public :: run_all_eul_prob_tab_tests

! eight in-plane rotations
integer, parameter :: NROTS = 8
real,    parameter :: SEED(2) = [1.5, -0.5]
! the synthetic dense tables of the assignment tests
integer, parameter :: NPTCLS_TAB = 150
integer, parameter :: SEED_TAB   = 20261016   !< the distances and the particles' angular distances
integer, parameter :: SEED_DRAW  = 20261017   !< the assignment's draws

contains

    subroutine run_all_eul_prob_tab_tests()
        write(*,'(A)') '**** running all 3D probability table tests ****'
        call test_record_seed_shift()
        write(*,'(A)') 'test_dense_part_round_trip'
        call test_dense_part_round_trip()
        write(*,'(A)') 'test_dense_read_threads'
        call test_dense_read_threads()
        write(*,'(A)') 'test_bounded_and_power_samplers'
        call test_bounded_and_power_samplers()
        write(*,'(A)') 'test_posterior_row_stats'
        call test_posterior_row_stats()
        write(*,'(A)') 'test_draw_posterior_row'
        call test_draw_posterior_row()
        write(*,'(A)') 'test_draw_posterior_extras'
        call test_draw_posterior_extras()
        write(*,'(A)') 'test_insert_cone'
        call test_insert_cone()
        write(*,'(A)') 'test_insert_table_round_trip'
        call test_insert_table_round_trip()
        write(*,'(A)') 'test_insert_extras_multistate'
        call test_insert_extras_multistate()
        write(*,'(A)') 'test_ref_assign_posterior'
        call test_ref_assign_posterior()
        write(*,'(A)') 'test_prob_beta_at'
        call test_prob_beta_at()
        write(*,'(A)') 'test_ref_assign_neigh_posterior'
        call test_ref_assign_neigh_posterior()
        write(*,'(A)') 'test_neff_readout'
        call test_neff_readout()
    end subroutine run_all_eul_prob_tab_tests

    !> The n_eff readout: 200 rows of 50 distances k s_r (k = 0..49, s_r = 0.01 r) at temperature 2, so row r's
    !! posterior at factor beta is geometric with ratio exp(-2 beta s_r). Reference values from neff_control_ref.py
    !! (numpy): one row at a = 0.5 has exp(H) = 5.493135; the mean over the rows is 49.9721 at beta 0.001, 6.023447
    !! at 0.5, 3.582882 at 1 and 2.243834 at 2, falling monotonically with the factor. Rows with fewer valid
    !! distances than nmin, a non-positive temperature or outside the stride are not kept, invalid distances are
    !! dropped
    subroutine test_neff_readout()
        integer, parameter :: NROWS_NC = 200, NCAND_NC = 50
        type(prob_neff_readout) :: neff_ro
        real    :: row(NCAND_NC), prev, neff
        integer :: r, k, nwanted
        logical :: l_monotone
        ! one row
        call neff_ro%new(1, 1, 2)
        row = [(0.25 * real(k - 1), k = 1,NCAND_NC)]   ! a = 2 * 0.25 = 0.5 at beta 1
        call neff_ro%offer(NCAND_NC, row, 2.)
        call assert_int(1, neff_ro%nkept(), 'n_eff readout: one row kept')
        call assert_real(5.493135, neff_ro%mean_neff(1.), 1.e-4, 'n_eff readout: exp(entropy) of a geometric row (numpy)')
        call assert_real(real(NCAND_NC), neff_ro%mean_neff(1.e-6), 1.e-2, 'n_eff readout: a cold factor of 1e-6 is the uniform row')
        call assert_real(1., neff_ro%mean_neff(1.e3), 1.e-4, 'n_eff readout: a huge factor is the argmin')
        call neff_ro%kill
        ! the 200 rows, all wanted (stride 1)
        call neff_ro%new(NROWS_NC, NROWS_NC, 6)
        do r = 1,NROWS_NC
            call assert_true(neff_ro%wanted(r), 'n_eff readout: every row wanted at stride 1')
            row = [(0.01 * real(r) * real(k - 1), k = 1,NCAND_NC)]
            call neff_ro%offer(NCAND_NC, row, 2.)
        enddo
        call assert_int(NROWS_NC, neff_ro%nkept(), 'n_eff readout: 200 rows kept')
        call assert_real(49.9721,  neff_ro%mean_neff(1.e-3), 1.e-2, 'n_eff readout: mean n_eff at 0.001 (numpy)')
        call assert_real(6.023447, neff_ro%mean_neff(0.5),   1.e-3, 'n_eff readout: mean n_eff at 0.5 (numpy)')
        call assert_real(3.582882, neff_ro%mean_neff(1.),    1.e-3, 'n_eff readout: mean n_eff at 1 (numpy)')
        call assert_real(2.243834, neff_ro%mean_neff(2.),    1.e-3, 'n_eff readout: mean n_eff at 2 (numpy)')
        l_monotone = .true.
        prev       = huge(prev)
        do k = -30,30
            neff = neff_ro%mean_neff(10.**(real(k)/10.))
            if( neff > prev ) l_monotone = .false.
            prev = neff
        enddo
        call assert_true(l_monotone, 'n_eff readout: the mean n_eff falls monotonically with the factor')
        call neff_ro%report(0.3) ! the log line, over 200 rows
        call neff_ro%kill
        ! doubling the temperatures halves the factor at which a mean is met: the mean at beta 0.5 and T 4 is
        ! the mean at beta 1 and T 2
        call neff_ro%new(NROWS_NC, NROWS_NC, 6)
        do r = 1,NROWS_NC
            row = [(0.01 * real(r) * real(k - 1), k = 1,NCAND_NC)]
            call neff_ro%offer(NCAND_NC, row, 4.)
        enddo
        call assert_real(3.582882, neff_ro%mean_neff(0.5), 1.e-3, 'n_eff readout: the mean scales with beta times the temperature')
        call neff_ro%kill
        ! the stride, the minimum row size, the temperature, invalid distances
        call neff_ro%new(1000, 100, 6)
        nwanted = 0
        do r = 1,1000
            if( neff_ro%wanted(r) ) nwanted = nwanted + 1
        enddo
        call assert_int(100, nwanted, 'n_eff readout: 100 of 1000 rows wanted at stride 10')
        call assert_true(neff_ro%wanted(1) .and. neff_ro%wanted(11) .and. .not. neff_ro%wanted(2), 'n_eff readout: the stride starts at row 1')
        row = [(0.25 * real(k - 1), k = 1,NCAND_NC)]
        call neff_ro%offer(5, row(1:5), 2.)
        call assert_int(0, neff_ro%nkept(), 'n_eff readout: a row of five candidates is not kept at nmin 6')
        call neff_ro%offer(NCAND_NC, row, 0.)
        call assert_int(0, neff_ro%nkept(), 'n_eff readout: a row without a temperature is not kept')
        row(7:NCAND_NC) = huge(1.)
        call neff_ro%offer(NCAND_NC, row, 2.)
        call assert_int(1, neff_ro%nkept(), 'n_eff readout: six valid distances among invalid ones are kept')
        call assert_real(6., neff_ro%mean_neff(1.e-6), 1.e-3, 'n_eff readout: the invalid distances are dropped')
        call assert_true(PROB_NEFF_MIN_ROWS > 1, 'n_eff readout: nothing is reported below PROB_NEFF_MIN_ROWS rows')
        call neff_ro%report(0.3) ! below the minimum: silent
        call neff_ro%kill
    end subroutine test_neff_readout

    !> The tempering ramp (prob_beta_at): beta0 at the first iteration, beta1 at the last, the geometric
    !! midpoint halfway (0.1 to 1 over 25 iterations: sqrt(0.1) at the 13th), the same factor every step;
    !! flat when the ends are equal; a one-iteration run at beta1; iterations outside the run clamped
    subroutine test_prob_beta_at()
        real    :: ratio, prev
        integer :: it
        logical :: l_ok
        call assert_real(0.1, prob_beta_at(1, 25, 0.1, 1.), 1.e-6, 'the ramp starts at beta0')
        call assert_real(1.0, prob_beta_at(25, 25, 0.1, 1.), 1.e-6, 'the ramp ends at beta1')
        call assert_real(sqrt(0.1), prob_beta_at(13, 25, 0.1, 1.), 1.e-6, 'the geometric midpoint halfway')
        l_ok  = .true.
        ratio = prob_beta_at(2, 25, 0.1, 1.) / prob_beta_at(1, 25, 0.1, 1.)
        prev  = prob_beta_at(1, 25, 0.1, 1.)
        do it = 2,25
            if( abs(prob_beta_at(it, 25, 0.1, 1.) / prev - ratio) > 1.e-5 ) l_ok = .false.
            prev = prob_beta_at(it, 25, 0.1, 1.)
        enddo
        call assert_true(l_ok, 'every iteration cools by the same factor')
        call assert_real(10.**(1./24.), ratio, 1.e-5, 'the factor is (beta1/beta0)^(1/(n-1))')
        call assert_real(0.3, prob_beta_at(7, 20, 0.3, 0.3), 1.e-6, 'equal ends give a flat run')
        call assert_real(1.0, prob_beta_at(1, 1, 0.1, 1.), 1.e-6, 'a one-iteration run is at beta1')
        call assert_real(0.1, prob_beta_at(0, 25, 0.1, 1.), 1.e-6, 'an iteration before the run is clamped to the start')
        call assert_real(1.0, prob_beta_at(30, 25, 0.1, 1.), 1.e-6, 'an iteration past the run is clamped to the end')
    end subroutine test_prob_beta_at

    !> the neighbourhood table's posterior assignment (ref_assign_neigh; the same draw on every
    !! neighbourhood mode): a synthetic candidate store of 150 particles over 3 states x 20 directions (1 to 12
    !! candidates each, some without an in-plane index; every particle has a valid one, since a particle without one
    !! is re-seeded through the calculator, which this test has not). At a huge temperature every particle takes the
    !! state of its closest valid candidate and the argmin within it, with its in-plane index and distance and the
    !! evaluated fraction. At T = 2 some particles leave the mode and none its state; the same seed gives the same
    !! assignment, another seed another.
    subroutine test_ref_assign_neigh_posterior()
        integer, parameter :: NSPACE = 20, NSTATES = 3, NP = 150
        class(parameters), allocatable, target :: params
        class(builder),    allocatable, target :: build
        type(eul_prob_tab_neigh) :: tab
        type(prob_candidate), allocatable :: cands(:)
        integer, allocatable :: counts(:), pinds(:), state_best(:), ref_best(:), run1(:), run2(:)
        real,    allocatable :: dist_best(:)
        integer :: i, k, n, nrefs, nwrong, nmoved, nstate_ok, iref, nvalid
        nrefs = NSPACE * NSTATES
        allocate(params, build)
        params%nspace          = NSPACE
        params%nstates         = NSTATES
        params%prob_athres     = 10.
        params%l_doshift       = .false.
        params%prob_neigh_mode = 'snhc'
        params%prob_beta       = 1.
        call seed_rnd_fixed(SEED_TAB)
        allocate(build%spproj_field)
        call build%spproj_field%new(NP, is_ptcl=.true.)
        do i = 1,NP
            call build%spproj_field%set_state(i, mod(i-1, NSTATES) + 1)
            call build%spproj_field%set(i, 'sampled', 1.)
            call build%spproj_field%set(i, 'dist',    20. * ran3())
        enddo
        pinds = [(i, i=1,NP)]
        allocate(counts(NP))
        do i = 1,NP
            counts(i) = 1 + int(ran3() * 11.)
        enddo
        allocate(cands(sum(counts)))
        n = 0
        do i = 1,NP
            do k = 1,counts(i)
                n = n + 1
                cands(n)%iref = 1 + int(ran3() * real(nrefs))       ! a full reference index (state, direction)
                cands(n)%inpl = 1 + int(ran3() * real(NROTS))
                cands(n)%dist = 4. * ran3()
                if( k > 1 .and. ran3() < 0.1 ) cands(n)%inpl = 0    ! not a candidate
            enddo
        enddo
        ! the expected argmin per particle: the closest valid candidate over all states
        allocate(state_best(NP), ref_best(NP), dist_best(NP), run1(NP), run2(NP))
        n = 0
        do i = 1,NP
            state_best(i) = 0
            ref_best(i)   = 0
            dist_best(i)  = huge(1.)
            do k = 1,counts(i)
                if( cands(n+k)%inpl <= 0 ) cycle
                if( cands(n+k)%dist < dist_best(i) )then
                    dist_best(i)  = cands(n+k)%dist
                    ref_best(i)   = cands(n+k)%iref
                    state_best(i) = (cands(n+k)%iref - 1) / NSPACE + 1
                endif
            enddo
            n = n + counts(i)
        enddo
        ! a huge temperature: the argmin within the argmin's state
        call set_store(1.e9)
        call seed_rnd_fixed(SEED_DRAW)
        call tab%ref_assign
        nwrong = 0
        n = 0
        do i = 1,NP
            if( state_best(i) > 0 )then
                iref = (tab%assgn_map(i)%istate - 1) * NSPACE + tab%assgn_map(i)%iproj
                if( iref /= ref_best(i) ) nwrong = nwrong + 1
                if( tab%assgn_map(i)%dist /= dist_best(i) ) nwrong = nwrong + 1
                if( tab%assgn_map(i)%inpl < 1 ) nwrong = nwrong + 1
                nvalid = 0
                do k = 1,counts(i)
                    if( cands(n+k)%inpl > 0 ) nvalid = nvalid + 1
                enddo
                if( abs(tab%assgn_map(i)%frac - 100. * real(nvalid) / real(nrefs)) > 1.e-4 ) nwrong = nwrong + 1
            endif
            n = n + counts(i)
        enddo
        call assert_int(0, nwrong, 'a huge temperature gives every particle its closest valid neighbourhood candidate')
        call tab%kill
        ! T = 2: some particles leave the mode, none its state; seeds
        call set_store(2.)
        call seed_rnd_fixed(SEED_DRAW)
        call tab%ref_assign
        nmoved    = 0
        nstate_ok = 0
        do i = 1,NP
            run1(i) = (tab%assgn_map(i)%istate - 1) * NSPACE + tab%assgn_map(i)%iproj
            if( state_best(i) == 0 ) cycle
            if( run1(i) /= ref_best(i) ) nmoved = nmoved + 1
            if( tab%assgn_map(i)%istate == state_best(i) ) nstate_ok = nstate_ok + 1
        enddo
        call assert_true(nmoved > 0 .and. nmoved < NP, 'at T = 2 some particles leave the mode')
        call assert_int(NP, nstate_ok, 'the state is the argmin''s state under the draw')
        call tab%kill
        call set_store(2.)
        call seed_rnd_fixed(SEED_DRAW)
        call tab%ref_assign
        run2 = (tab%assgn_map(:)%istate - 1) * NSPACE + tab%assgn_map(:)%iproj
        call assert_int(0, count(run2 /= run1), 'the same seed gives the same neighbourhood posterior assignment')
        call tab%kill
        call set_store(2.)
        call seed_rnd_fixed(SEED_DRAW + 1)
        call tab%ref_assign
        run2 = (tab%assgn_map(:)%istate - 1) * NSPACE + tab%assgn_map(:)%iproj
        call assert_true(any(run2 /= run1), 'another seed gives another neighbourhood posterior assignment')
        call tab%kill
        call build%spproj_field%kill
        deallocate(params, build, cands, counts, pinds, state_best, ref_best, dist_best, run1, run2)

    contains

        subroutine set_store( temp )
            real, intent(in) :: temp
            call tab%new_neigh_global(params, build, pinds)
            call tab%candidate_store%new_ragged(pinds, counts)
            tab%candidate_store%candidates = cands
            tab%temps = temp
        end subroutine set_store

    end subroutine test_ref_assign_neigh_posterior

    !> The draw from a row posterior (draw_posterior_row): a huge temperature gives a unique argmin every
    !! time and one of two equal minima with equal probability (both drawn over 50 draws, never another);
    !! temperature 0 draws uniformly over the valid distances and never an invalid one; at T = 1 the row
    !! [0, 1, 2] is drawn with p = [1, e^-1, e^-2]/Z = [0.665, 0.245, 0.090] to within 0.02 over 20000 draws;
    !! the same seed gives the same sequence; a row of invalid distances draws 0
    subroutine test_draw_posterior_row()
        integer, parameter :: NDRAWS = 20000
        real    :: dists(6), w(6), dist
        integer :: which, i, counts(6), seq1(50), seq2(50)
        logical :: l_ok
        ! huge temperature, a unique minimum: the argmin every time
        dists = [0.3, 0.1, 0.15, 0.5, huge(1.), 0.2]
        call seed_rnd_fixed(20261023)
        l_ok = .true.
        do i = 1,50
            call draw_posterior_row(6, dists, huge(1.), w, which, dist)
            l_ok = l_ok .and. which == 2 .and. dist == 0.1
        enddo
        call assert_true(l_ok, 'a huge temperature draws a unique argmin every time')
        ! huge temperature, two equal minima: either, with equal probability, never another
        dists = [0.3, 0.1, 0.1, 0.5, huge(1.), 0.2]
        counts = 0
        do i = 1,50
            call draw_posterior_row(6, dists, huge(1.), w, which, dist)
            counts(which) = counts(which) + 1
        enddo
        call assert_int(50, counts(2) + counts(3), 'a huge temperature draws only the equal minima')
        call assert_true(counts(2) > 0 .and. counts(3) > 0, 'a huge temperature draws both equal minima')
        ! temperature 0: uniform over the five valid distances, never the invalid one
        counts = 0
        do i = 1,NDRAWS
            call draw_posterior_row(6, dists, 0., w, which, dist)
            counts(which) = counts(which) + 1
        enddo
        call assert_int(0, counts(5), 'temperature 0 never draws an invalid distance')
        call assert_true(all(counts([1,2,3,4,6]) > NDRAWS/5 - 300) .and. all(counts([1,2,3,4,6]) < NDRAWS/5 + 300), &
            &'temperature 0 draws the valid distances uniformly')
        ! T = 1 on [0, 1, 2]
        counts = 0
        do i = 1,NDRAWS
            call draw_posterior_row(3, [0., 1., 2.], 1., w, which, dist)
            counts(which) = counts(which) + 1
        enddo
        call assert_true(abs(real(counts(1))/NDRAWS - 0.665241) < 0.02, 'the mode is drawn with its posterior mass')
        call assert_true(abs(real(counts(2))/NDRAWS - 0.244728) < 0.02, 'the second is drawn with e^-1 of the mode''s mass')
        call assert_true(abs(real(counts(3))/NDRAWS - 0.090031) < 0.02, 'the third is drawn with e^-2 of the mode''s mass')
        ! reproducibility from the seed
        call seed_rnd_fixed(20261024)
        do i = 1,50
            call draw_posterior_row(3, [0., 1., 2.], 1., w, seq1(i), dist)
        enddo
        call seed_rnd_fixed(20261024)
        do i = 1,50
            call draw_posterior_row(3, [0., 1., 2.], 1., w, seq2(i), dist)
        enddo
        call assert_true(all(seq1 == seq2), 'the same seed gives the same draws')
        call assert_true(any(seq1 /= seq1(1)), 'the draws at T = 1 are not all the same')
        ! no valid distance
        call draw_posterior_row(3, [huge(1.), huge(1.), huge(1.)], 1., w, which, dist)
        call assert_int(0, which, 'a row without a valid distance draws 0')
    end subroutine test_draw_posterior_row

    !> ref_assign on synthetic dense tables (one state of 200 directions, three states of 200, two states of one
    !! direction; random distances with the tie cases: particle 1 flat, particle 2 tied at its minimum between a
    !! reference of state 2 and one of state 3, particle 3 between state 1 and state 3, particle 4 with invalid
    !! distances only, particles 5 and 6 identical): with a huge
    !! temperature every particle whose row minimum is unique takes it; the flat particle 1 takes a uniformly
    !! drawn direction in state 1 (the state label is the lowest-index argmin's) and the invalid particle 4
    !! takes reference 1; the state-2/state-3 and state-1/state-3 ties of particles 2 and 3 go to the lower
    !! state, where their minimum is unique. At a temperature of 2 with distances in [0, 5] the draw leaves the
    !! mode for a fraction of the particles, and the state is always the argmin's; the draw is seeded, so two
    !! runs from the same seed agree and two seeds differ
    subroutine test_ref_assign_posterior()
        call check_posterior_assignment(nspace=200, nstates=1)
        call check_posterior_assignment(nspace=200, nstates=3)
        call check_posterior_assignment(nspace=1,  nstates=2)
    end subroutine test_ref_assign_posterior

    subroutine check_posterior_assignment( nspace, nstates )
        integer, intent(in) :: nspace, nstates
        class(parameters), allocatable, target :: params
        class(builder),    allocatable, target :: build
        type(eul_prob_tab) :: tab
        real,    allocatable :: dists(:,:)
        integer, allocatable :: pinds(:), iref_argmin(:), iref_run1(:), iref_run2(:), iref_run3(:)
        integer :: nrefs, i, iref, nsame, nstate_ok, nmoved
        character(len=64) :: label
        write(label,'(A,I0,A,I0)') 'nspace ', nspace, ', nstates ', nstates
        nrefs = nspace * nstates
        allocate(params, build)
        params%nspace      = nspace
        params%nstates     = nstates
        params%prob_athres = 10.
        params%l_doshift   = .false.
        params%prob_beta   = 1.
        call seed_rnd_fixed(SEED_TAB)
        allocate(build%spproj_field)
        call build%spproj_field%new(NPTCLS_TAB, is_ptcl=.true.)
        do i = 1,NPTCLS_TAB
            call build%spproj_field%set_state(i, mod(i-1, nstates) + 1)
            call build%spproj_field%set(i, 'sampled', 1.)
            call build%spproj_field%set(i, 'dist',    20. * ran3())
        enddo
        pinds = [(i, i=1,NPTCLS_TAB)]
        allocate(dists(nrefs,NPTCLS_TAB), iref_argmin(NPTCLS_TAB), iref_run1(NPTCLS_TAB), iref_run2(NPTCLS_TAB), iref_run3(NPTCLS_TAB))
        do i = 1,NPTCLS_TAB
            do iref = 1,nrefs
                dists(iref,i) = 5. * ran3()
            enddo
        enddo
        dists(:,1) = 1.
        if( nstates >= 3 )then
            dists(:,2) = 3. + dists(:,2)
            dists(nspace + min(5,nspace), 2)   = 0.01
            dists(2*nspace + min(7,nspace), 2) = 0.01
            dists(:,3) = 3. + dists(:,3)
            dists(min(10,nspace), 3)    = 0.01
            dists(2*nspace + 1, 3)      = 0.01
        endif
        dists(:,4) = huge(1.)
        dists(:,6) = dists(:,5)
        ! the argmin of every row, the lowest index among equal minima
        do i = 1,NPTCLS_TAB
            iref_argmin(i) = 1
            do iref = 2,nrefs
                if( dists(iref,i) < dists(iref_argmin(i),i) ) iref_argmin(i) = iref
            enddo
        enddo
        call tab%new(params, build, pinds)
        do i = 1,NPTCLS_TAB
            do iref = 1,nrefs
                call set_dist(tab, iref, i, dists(iref,i))
            enddo
        enddo
        ! a huge temperature: the argmin wherever the row minimum is unique (every row but the flat particle 1)
        tab%temps = huge(1.)
        call seed_rnd_fixed(SEED_DRAW)
        call tab%ref_assign
        nsame = 0
        do i = 2,NPTCLS_TAB
            if( assigned_ref(tab, i, nspace) == iref_argmin(i) ) nsame = nsame + 1
        enddo
        call assert_int(NPTCLS_TAB - 1, nsame, 'a huge temperature assigns every particle with a unique minimum its row argmin ('//trim(label)//')')
        call assert_int(1, tab%assgn_map(1)%istate, 'a flat particle takes the state of reference 1 ('//trim(label)//')')
        call assert_true(tab%assgn_map(1)%iproj >= 1 .and. tab%assgn_map(1)%iproj <= nspace, &
            &'a flat particle takes a direction of its state ('//trim(label)//')')
        call assert_int(1, tab%assgn_map(4)%istate, 'an invalid particle takes the state of reference 1 ('//trim(label)//')')
        call assert_int(1, tab%assgn_map(4)%iproj,  'an invalid particle takes reference 1 ('//trim(label)//')')
        if( nstates >= 3 )then
            call assert_int(2, tab%assgn_map(2)%istate, 'a state 2 / state 3 tie goes to state 2 ('//trim(label)//')')
            call assert_int(1, tab%assgn_map(3)%istate, 'a state 1 / state 3 tie goes to state 1 ('//trim(label)//')')
        endif
        ! the particles' temperatures: a draw that leaves the mode for some, never the argmin's state
        tab%temps = 2.
        call seed_rnd_fixed(SEED_DRAW)
        call tab%ref_assign
        nstate_ok = 0
        nmoved    = 0
        do i = 1,NPTCLS_TAB
            iref_run1(i) = assigned_ref(tab, i, nspace)
            if( tab%assgn_map(i)%istate == (iref_argmin(i) - 1)/nspace + 1 ) nstate_ok = nstate_ok + 1
            if( iref_run1(i) /= iref_argmin(i) ) nmoved = nmoved + 1
        enddo
        call assert_int(NPTCLS_TAB, nstate_ok, 'the state is the argmin''s state under the draw ('//trim(label)//')')
        ! the tempering factor: temps 4 at prob_beta 0.5 is temps 2 at 1, draw for draw
        params%prob_beta = 0.5
        tab%temps = 4.
        call seed_rnd_fixed(SEED_DRAW)
        call tab%ref_assign
        nsame = 0
        do i = 1,NPTCLS_TAB
            if( assigned_ref(tab, i, nspace) == iref_run1(i) ) nsame = nsame + 1
        enddo
        call assert_int(NPTCLS_TAB, nsame, 'prob_beta scales the temperature of the draw ('//trim(label)//')')
        params%prob_beta = 1.
        tab%temps = 2.
        if( nspace > 1 )then
            call assert_true(nmoved > 0, 'at temperature 2 on distances in [0,5] some particles leave the mode ('//trim(label)//')')
            call assert_true(nmoved < NPTCLS_TAB, 'at temperature 2 some particles stay at the mode ('//trim(label)//')')
        else
            call assert_int(0, nmoved, 'with one direction per state nothing moves ('//trim(label)//')')
        endif
        call seed_rnd_fixed(SEED_DRAW)
        call tab%ref_assign
        do i = 1,NPTCLS_TAB
            iref_run2(i) = assigned_ref(tab, i, nspace)
        enddo
        call assert_true(all(iref_run1 == iref_run2), 'the same seed gives the same assignment ('//trim(label)//')')
        if( nspace > 1 )then
            call seed_rnd_fixed(SEED_DRAW + 1)
            call tab%ref_assign
            do i = 1,NPTCLS_TAB
                iref_run3(i) = assigned_ref(tab, i, nspace)
            enddo
            call assert_true(any(iref_run1 /= iref_run3), 'another seed gives another assignment ('//trim(label)//')')
        endif
        call tab%kill
        call build%spproj_field%kill
        deallocate(build%spproj_field)
        deallocate(params, build)
    end subroutine check_posterior_assignment

    !> the reference index (state block, direction) a particle was assigned
    integer function assigned_ref( tab, i, nspace ) result( iref )
        type(eul_prob_tab), intent(in) :: tab
        integer,            intent(in) :: i, nspace
        iref = (tab%assgn_map(i)%istate - 1) * nspace + tab%assgn_map(i)%iproj
    end function assigned_ref

    !> The weighted insertion's table: a particle without extras keeps weight 1; add records the committed weight
    !! and the extras in order; the file round trip returns every field; the readout line prints; the part's
    !! subset keeps the part's particles in its order with their entries; the resolved-poses layout follows the
    !! subset (slots, counts, committed weights, zero weights until filled)
    subroutine test_insert_table_round_trip()
        type(prob_insert_table) :: tab, back, part
        type(prob_insert_poses) :: poses
        type(prob_insert_pose)  :: pose
        type(prob_candidate)    :: cands(2)
        type(string) :: fname
        integer :: i
        call tab%new(4, [11, 12, 13, 14], 3)
        cands(1) = prob_candidate(101, 7, 0.4, 0.5, 1.5, .true.)
        cands(2) = prob_candidate(102, 8, 0.6, -1.0, 2.0, .true.)
        call tab%add(2, 0.6, 2, cands, [0.3, 0.1])
        cands(1) = prob_candidate(103, 9, 0.2, 0.25, 0.0, .true.)
        call tab%add(4, 0.9, 1, cands(1:1), [0.1])
        call assert_int(4, tab%get_nptcls(), 'insert table: four particles')
        call assert_int(3, tab%get_kmax(),   'insert table: kmax 3')
        call assert_real(1., tab%get_committed_weight(1), 0., 'insert table: a particle without extras keeps weight 1')
        call assert_int(0, tab%get_nextra(1), 'insert table: and no extras')
        call assert_real(0.6, tab%get_committed_weight(2), 0., 'insert table: the committed weight recorded')
        call assert_int(2, tab%get_nextra(2), 'insert table: two extras recorded')
        pose = tab%get_extra(2, 2)
        call assert_true(pose%iref == 102 .and. pose%inpl == 8 .and. pose%x == -1.0 .and. pose%y == 2.0 .and. pose%w == 0.1, &
            &'insert table: the second extra of particle 2')
        call tab%report
        fname = string('prob_insert_table_test.dat')
        call tab%write(fname)
        call back%read(fname)
        call assert_int(4, back%get_nptcls(), 'insert table round trip: the particle count')
        call assert_int(3, back%get_kmax(),   'insert table round trip: kmax')
        do i = 1,4
            call assert_int(tab%get_pind(i), back%get_pind(i), 'insert table round trip: the particle index')
            call assert_real(tab%get_committed_weight(i), back%get_committed_weight(i), 0., 'insert table round trip: the committed weight')
            call assert_int(tab%get_nextra(i), back%get_nextra(i), 'insert table round trip: the extra count')
        enddo
        pose = back%get_extra(4, 1)
        call assert_true(pose%iref == 103 .and. pose%inpl == 9 .and. pose%x == 0.25 .and. pose%y == 0.0 .and. pose%w == 0.1, &
            &'insert table round trip: the extra of particle 4')
        pose = back%get_extra(2, 1)
        call assert_true(pose%iref == 101 .and. pose%inpl == 7 .and. pose%w == 0.3, 'insert table round trip: the first extra of particle 2')
        call del_file(fname)
        ! a part holding particles 14, 12 and 20 (the last unknown to the master's table), in that order
        call part%subset(tab, [14, 12, 20])
        call assert_int(3, part%get_nptcls(), 'insert table subset: the part''s particle count')
        call assert_int(14, part%get_pind(1), 'insert table subset: the part''s order')
        call assert_int(1, part%get_nextra(1), 'insert table subset: particle 14''s one extra')
        call assert_int(2, part%get_nextra(2), 'insert table subset: particle 12''s two extras')
        call assert_int(0, part%get_nextra(3), 'insert table subset: an unknown particle has no extra')
        call assert_real(0.9, part%get_committed_weight(1), 0., 'insert table subset: particle 14''s committed weight')
        call assert_real(1.0, part%get_committed_weight(3), 0., 'insert table subset: an unknown particle keeps weight 1')
        pose = part%get_extra(2, 2)
        call assert_true(pose%iref == 102 .and. pose%w == 0.1, 'insert table subset: particle 12''s second extra')
        ! the resolved-poses layout from the part's table
        call poses%setup(part, [14, 12, 20])
        call assert_int(3, poses%npose, 'insert poses: three slots')
        call assert_true(poses%has(14) .and. poses%has(12) .and. poses%has(20) .and. .not. poses%has(13), &
            &'insert poses: the part''s particles have entries, others do not')
        call assert_int(1, poses%first(1), 'insert poses: particle 14''s slot')
        call assert_int(2, poses%first(2), 'insert poses: particle 12''s slots follow')
        call assert_int(2, poses%loc(12),  'insert poses: the global index maps to the local one')
        call assert_real(0.6, poses%wc(2), 0., 'insert poses: the committed weight')
        call assert_true(all(poses%w == 0.), 'insert poses: the slots'' weights are zero until the search fills them')
        call poses%kill
        call part%kill
        call tab%kill
        call back%kill
    end subroutine test_insert_table_round_trip

    !> The weighted insertion's extras on a two-state dense table (the round-trip parts: 50 directions, 2 states,
    !! 40 particles, a tenth of the entries refined): every extra lies in the committed state's block, is a refined
    !! entry and is not the committed reference; a particle's weights sum to 1 and no extra is below
    !! PROB_INSERT_WMIN; extras do occur at a hot temperature; the master's table round-trips through the file the
    !! workers read; topk = 1 records no extra and weight 1 everywhere
    subroutine test_insert_extras_multistate()
        integer, parameter :: NSPACE_MS = 50, NSTATES_MS = 2, NPTCLS_MS = 40, NPARTS_MS = 2
        class(parameters), allocatable, target :: params
        class(builder),    allocatable, target :: build
        type(eul_prob_tab)      :: master
        type(prob_insert_table) :: back
        type(prob_insert_pose)  :: pose, pose_a
        type(prob_insert_cone)  :: cone
        type(prob_dense_entry), allocatable :: rows(:,:)
        real,                   allocatable :: shifts(:,:,:), seeds(:,:)
        logical,                allocatable :: refined(:,:)
        integer,                allocatable :: pinds(:)
        type(string) :: fname
        integer :: i, k, nextra_tot, nwrong_state, nwrong_refined, nwrong_committed, nwrong_weight, nwrong_dup, iref_committed, istate
        real    :: wsum
        call make_round_trip_parts(params, build, NSPACE_MS, NSTATES_MS, NPTCLS_MS, NPARTS_MS, rows, shifts, refined, seeds)
        pinds = [(i, i=NPTCLS_MS,1,-1)]
        params%topk = 3
        call master%new(params, build, pinds)
        call master%read_tabs_to_glob(string('prob_tab_tester_part'), NPARTS_MS, 1)
        master%temps = 1. ! hot rows: extras occur
        ! no cone set: no extra pose, the committed weight 1 everywhere
        call seed_rnd_fixed(SEED_DRAW)
        call master%ref_assign
        nextra_tot = 0
        do i = 1,NPTCLS_MS
            nextra_tot = nextra_tot + master%insert_tab%get_nextra(i)
            if( master%insert_tab%get_committed_weight(i) /= 1. ) nextra_tot = nextra_tot + 1
        enddo
        call assert_int(0, nextra_tot, 'no insertion cone: no extra pose and the committed weight 1 everywhere')
        ! a cone of a ten-thousandth of an Angstrom: no pose but the committed one passes, so no extra pose
        call cone%new([1.e-4, 1.e-4], 100., 1.)
        call master%set_insert_cone(cone)
        call seed_rnd_fixed(SEED_DRAW)
        call master%ref_assign
        nextra_tot = 0
        do i = 1,NPTCLS_MS
            nextra_tot = nextra_tot + master%insert_tab%get_nextra(i)
        enddo
        call assert_int(0, nextra_tot, 'a cone of the committed pose alone: no extra pose')
        ! a resolution past the rim: the whole sphere, every in-plane angle and shift pass, the extras are the state's
        ! refined entries
        call cone%new([1000., 1000.], 100., 1.)
        call master%set_insert_cone(cone)
        call seed_rnd_fixed(SEED_DRAW)
        call master%ref_assign
        nextra_tot = 0; nwrong_state = 0; nwrong_refined = 0; nwrong_committed = 0; nwrong_weight = 0; nwrong_dup = 0
        do i = 1,NPTCLS_MS
            istate         = master%assgn_map(i)%istate
            iref_committed = (istate - 1) * NSPACE_MS + master%assgn_map(i)%iproj
            wsum           = master%insert_tab%get_committed_weight(i)
            do k = 1,master%insert_tab%get_nextra(i)
                pose       = master%insert_tab%get_extra(i, k)
                nextra_tot = nextra_tot + 1
                wsum       = wsum + pose%w
                if( (pose%iref - 1) / NSPACE_MS + 1 /= istate )  nwrong_state     = nwrong_state + 1
                if( .not. refined(pose%iref, pinds(i)) )          nwrong_refined   = nwrong_refined + 1
                if( pose%iref == iref_committed )                 nwrong_committed = nwrong_committed + 1
                if( pose%w < 0.01 - 1.e-6 )                       nwrong_weight    = nwrong_weight + 1
                if( k == 2 )then
                    pose_a = master%insert_tab%get_extra(i, 1)
                    if( pose%iref == pose_a%iref ) nwrong_dup = nwrong_dup + 1
                endif
            enddo
            if( abs(wsum - 1.) > 1.e-5 ) nwrong_weight = nwrong_weight + 1
            if( master%insert_tab%get_nextra(i) > 2 ) nwrong_weight = nwrong_weight + 1
        enddo
        call assert_true(nextra_tot > 0, 'multi-state extras: hot rows give extra poses')
        call assert_int(0, nwrong_state,     'multi-state extras: every extra is in the committed state''s block')
        call assert_int(0, nwrong_refined,   'multi-state extras: every extra is a refined entry')
        call assert_int(0, nwrong_committed, 'multi-state extras: no extra is the committed reference')
        call assert_int(0, nwrong_dup,       'multi-state extras: the extras of a particle are distinct')
        call assert_int(0, nwrong_weight,    'multi-state extras: weights sum to 1, no extra below 1%, at most topk - 1 extras')
        ! the file the workers read
        fname = string('prob_insert_multistate_test.dat')
        call master%write_insert_table(fname)
        call back%read(fname)
        nwrong_weight = 0
        do i = 1,NPTCLS_MS
            if( back%get_pind(i) /= pinds(i) ) nwrong_weight = nwrong_weight + 1
            if( back%get_committed_weight(i) /= master%insert_tab%get_committed_weight(i) ) nwrong_weight = nwrong_weight + 1
            if( back%get_nextra(i) /= master%insert_tab%get_nextra(i) ) nwrong_weight = nwrong_weight + 1
            do k = 1,back%get_nextra(i)
                pose   = back%get_extra(i, k)
                pose_a = master%insert_tab%get_extra(i, k)
                if( pose%iref /= pose_a%iref .or. pose%inpl /= pose_a%inpl .or. pose%x /= pose_a%x .or. pose%y /= pose_a%y &
                    &.or. pose%w /= pose_a%w ) nwrong_weight = nwrong_weight + 1
            enddo
        enddo
        call assert_int(0, nwrong_weight, 'multi-state extras: the table round-trips through the file')
        call del_file(fname)
        call back%kill
        ! topk 1: the committed pose alone
        params%topk = 1
        call seed_rnd_fixed(SEED_DRAW)
        call master%ref_assign
        nextra_tot = 0
        do i = 1,NPTCLS_MS
            nextra_tot = nextra_tot + master%insert_tab%get_nextra(i)
            if( master%insert_tab%get_committed_weight(i) /= 1. ) nextra_tot = nextra_tot + 1
        enddo
        call assert_int(0, nextra_tot, 'topk 1: no extra pose and the committed weight 1 everywhere')
        call master%kill
        call cone%kill
        call kill_round_trip_parts(params, build, NPARTS_MS)
    end subroutine test_insert_extras_multistate

    !> The insertion cone: the half-angle whose arc at the rim of a particle of diameter diam is the resolution (10 A
    !! at 100 A: 0.2 rad, 11.46 degrees) and the shift bound in pixels; no resolution gives no cone; a resolution past
    !! the rim's half-circumference the whole sphere. On a 1000-direction C1 spiral with a 36-step rotation grid a
    !! neighbour within the half-angle at the same in-plane angle and shift is within, a far direction is not, the
    !! same direction a half-turn away in-plane is not (one step of 10 degrees is), a shift past the bound is not.
    !! Modulo C4, equatorial azimuths 1 and 89 degrees are 2 degrees apart and pass a 3-degree cone that C1 rejects
    subroutine test_insert_cone()
        type(prob_insert_cone) :: cone
        type(oris)             :: eul, eul2
        type(sym), target      :: se1, se4
        real    :: e1(3), e2(3), esym(3), euldist, inpldist
        integer :: i, inear, ifar
        call cone%new([10.], 100., 1.)
        call assert_real(11.459156, cone%get_angle(1), 1.e-3, 'insertion cone: 10 A at a 100 A particle is 11.46 degrees')
        call assert_true(cone%is_open(), 'insertion cone: open')
        call cone%new([0.], 100., 1.)
        call assert_true(.not. cone%is_open() .and. cone%get_angle(1) < 0., 'insertion cone: no resolution, no cone')
        call cone%new([1000.], 100., 1.)
        call assert_real(180., cone%get_angle(1), 1.e-3, 'insertion cone: a resolution past the rim is the whole sphere')
        ! the C1 spiral
        call se1%new('c1')
        call eul%new(1000, is_ptcl=.false.)
        call se1%build_refspiral(eul)
        call cone%new([10.], 100., 1.)
        call cone%set_space(eul, se1, 36)
        call assert_true(cone%has_space(), 'insertion cone: the space is set')
        inear = 0
        ifar  = 0
        e1    = eul%get_euler(1)
        do i = 2,1000
            e2 = eul%get_euler(i)
            call se1%sym_dists(e1, e2, esym, euldist, inpldist)
            if( inear == 0 .and. euldist > 3. .and. euldist < 11. ) inear = i
            if( ifar  == 0 .and. euldist > 20. ) ifar = i
            if( inear > 0 .and. ifar > 0 ) exit
        enddo
        call assert_true(inear > 0 .and. ifar > 0, 'insertion cone: a near and a far direction found on the spiral')
        call assert_true(cone%within(1, 1, 1, [0.,0.], inear, 1, [0.,0.]), 'insertion cone: a direction within the half-angle is within')
        call assert_true(.not. cone%within(1, 1, 1, [0.,0.], ifar, 1, [0.,0.]), 'insertion cone: a far direction is not')
        call assert_true(.not. cone%within(1, 1, 1, [0.,0.], 1, 19, [0.,0.]), 'insertion cone: the same direction a half-turn away in-plane is another pose')
        call assert_true(cone%within(1, 1, 1, [0.,0.], 1, 2, [0.,0.]), 'insertion cone: one rotation step (10 degrees) is within')
        call assert_true(.not. cone%within(1, 1, 1, [0.,0.], inear, 1, [11.,0.]), 'insertion cone: a shift past the bound is another pose')
        call assert_true(cone%within(1, 1, 1, [0.,0.], inear, 1, [6.,6.]), 'insertion cone: a shift within the bound is within')
        call assert_true(.not. cone%within(2, 1, 1, [0.,0.], inear, 1, [0.,0.]), 'insertion cone: a state without a cone admits nothing')
        ! modulo C4: azimuths 1 and 89 degrees on the equator
        call eul2%new(2, is_ptcl=.false.)
        call eul2%set_euler(1, [1., 90., 0.])
        call eul2%set_euler(2, [89., 90., 0.])
        call cone%new([2.618], 100., 1.)   ! 3 degrees
        call cone%set_space(eul2, se1, 36)
        call assert_true(.not. cone%within(1, 1, 1, [0.,0.], 2, 1, [0.,0.]), 'insertion cone: 88 degrees apart in C1')
        call se4%new('c4')
        call cone%set_space(eul2, se4, 36)
        call assert_true(cone%within(1, 1, 1, [0.,0.], 2, 1, [0.,0.]), 'insertion cone: 2 degrees apart modulo C4')
        call cone%kill
        call eul%kill
        call eul2%kill
        call se1%kill
        call se4%kill
    end subroutine test_insert_cone

    !> The extra poses of the weighted insertion (draw_posterior_extras): kmax = 1 keeps the committed entry
    !! alone at weight 1 and leaves the random stream untouched; the extras are drawn without replacement from
    !! the eligible entries only, never the committed one, with the weights the posterior at the given temperature over the kept
    !! set summing to 1; an extra below wmin is dropped; the committed entry is kept whatever its share; over many
    !! draws an extra is included with the probability of the remainder's renormalised posterior
    subroutine test_draw_posterior_extras()
        integer, parameter :: NDRAWS = 20000
        real    :: dists(6), work(6), w(3), dist, wrow(6)
        integer :: inds(3), nkept, i, k, nsecond, ndropped, seq1(20), seq2(20)
        logical :: mask(6), l_ok
        ! T = 1 on [0, 1, 2, 3, huge, 0.5], all eligible; the committed entry the mode
        dists = [0., 1., 2., 3., huge(1.), 0.5]
        mask  = .true.
        call seed_rnd_fixed(20261029)
        do i = 1,20
            call draw_posterior_row(3, [0., 1., 2.], 1., wrow, seq1(i), dist)
        enddo
        call seed_rnd_fixed(20261029)
        do i = 1,20
            call draw_posterior_extras(6, dists, 1., mask, 1, 1, 0.01, work, nkept, inds, w)
            call assert_int(1, nkept, 'kmax 1 keeps the committed entry alone')
            call assert_int(1, inds(1), 'kmax 1 keeps the committed entry')
            call assert_real(1., w(1), 0., 'kmax 1 gives the committed entry weight 1')
            call draw_posterior_row(3, [0., 1., 2.], 1., wrow, seq2(i), dist)
        enddo
        call assert_true(all(seq1 == seq2), 'kmax 1 draws no random number: the stream is the single draw''s')
        ! kmax 3: three distinct kept entries, the committed first, valid and eligible, weights summing to 1
        l_ok = .true.
        do i = 1,200
            call draw_posterior_extras(6, dists, 1., mask, 1, 3, 1.e-6, work, nkept, inds, w)
            l_ok = l_ok .and. nkept == 3 .and. inds(1) == 1 .and. inds(2) /= 1 .and. inds(3) /= 1 .and. inds(2) /= inds(3) &
                &.and. inds(2) /= 5 .and. inds(3) /= 5 .and. abs(sum(w(1:3)) - 1.) < 1.e-5
        enddo
        call assert_true(l_ok, 'kmax 3 keeps three distinct valid entries, the committed first, weights summing to 1')
        ! the weights are the posterior at the given temperature over the kept set: with entries {1 (d 0), 6 (d 0.5), 2 (d 1)} kept
        ! they are proportional to 1, e^-0.5, e^-1
        do i = 1,500
            call draw_posterior_extras(6, dists, 1., mask, 1, 3, 1.e-6, work, nkept, inds, w)
            if( nkept == 3 .and. ((inds(2) == 6 .and. inds(3) == 2) .or. (inds(2) == 2 .and. inds(3) == 6)) ) exit
        enddo
        call assert_true(i <= 500, 'the pair {6, 2} of extras comes up within 500 draws')
        do k = 2,3
            if( inds(k) == 6 ) call assert_real(exp(-0.5) / (1. + exp(-0.5) + exp(-1.)), w(k), 1.e-5, 'the weight of entry 6 is its share of the kept set')
            if( inds(k) == 2 ) call assert_real(exp(-1.0) / (1. + exp(-0.5) + exp(-1.)), w(k), 1.e-5, 'the weight of entry 2 is its share of the kept set')
        enddo
        call assert_real(1. / (1. + exp(-0.5) + exp(-1.)), w(1), 1.e-5, 'the committed entry''s weight is its share of the kept set')
        ! eligibility: only entry 6 eligible, so one extra at most, and it is entry 6
        mask = .false.
        mask(6) = .true.
        call draw_posterior_extras(6, dists, 1., mask, 1, 3, 1.e-6, work, nkept, inds, w)
        call assert_int(2, nkept, 'with one eligible entry there is one extra')
        call assert_int(6, inds(2), 'and it is the eligible entry')
        ! the committed entry itself ineligible changes nothing about its being kept first
        mask = .true.
        mask(1) = .false.
        call draw_posterior_extras(6, dists, 1., mask, 1, 2, 1.e-6, work, nkept, inds, w)
        call assert_true(nkept == 2 .and. inds(1) == 1 .and. inds(2) /= 1, 'the committed entry is kept first whatever its eligibility')
        ! the mass rule: a near-delta row at a cold temperature keeps the committed entry alone
        mask = .true.
        ndropped = 0
        do i = 1,100
            call draw_posterior_extras(6, dists, 10., mask, 1, 3, 0.01, work, nkept, inds, w)
            if( nkept == 1 ) ndropped = ndropped + 1
        enddo
        call assert_int(100, ndropped, 'at T = 10 every extra falls below 1% and is dropped: the committed entry alone')
        ! inclusion frequency: three entries at T = 1, p = (0.665, 0.245, 0.090); given the committed entry 1, the
        ! single extra is entry 2 with probability 0.245 / (0.245 + 0.090) = 0.731
        nsecond = 0
        do i = 1,NDRAWS
            call draw_posterior_extras(3, [0., 1., 2.], 1., [.true., .true., .true.], 1, 2, 1.e-6, work, nkept, inds, w)
            if( nkept == 2 .and. inds(2) == 2 ) nsecond = nsecond + 1
        enddo
        call assert_true(abs(real(nsecond) / NDRAWS - 0.731059) < 0.02, 'an extra is drawn from the remainder''s renormalised posterior')
        ! a committed entry 80 nats above the minimum: its weight underflows, the kept set is the committed entry alone
        call draw_posterior_extras(3, [0., 100., 2.], 1., [.true., .true., .true.], 2, 3, 0.01, work, nkept, inds, w)
        call assert_true(nkept >= 1 .and. inds(1) == 2 .and. abs(sum(w(1:nkept)) - 1.) < 1.e-5, 'a committed entry of no mass is still kept first, the weights still sum to 1')
        ! no valid distance: the committed entry alone
        call draw_posterior_extras(3, [huge(1.), huge(1.), huge(1.)], 1., [.true., .true., .true.], 1, 3, 0.01, work, nkept, inds, w)
        call assert_int(1, nkept, 'a row without a valid distance keeps the committed entry alone')
    end subroutine test_draw_posterior_extras

    !> The pose posterior over a table row: T = pi wsqsum / (2 pftsz) (15.707963 for 1000 and 100); at T = 1000
    !! the row [0, 0.001, 0.002, 0.5, invalid, 0.4] gives p = [e^0, e^-1, e^-2, ~0, 0, ~0]/Z = [0.665241,
    !! 0.244728, 0.090031, 0, 0, 0] (numpy), so the mode is reference 1 with that mass, exp(H) = 2.298819
    !! and three references hold 99% of the mass; at T = 0 the posterior is uniform over the five valid
    !! distances (exp(H) = 5, k99 = 5); a row of invalid distances gives nvalid 0 and zeros
    subroutine test_posterior_row_stats()
        real,    parameter :: ROW(6) = [0.0, 0.001, 0.002, 0.5, huge(1.), 0.4]
        real    :: w(6), pmax, neff, temp
        integer :: inds(6), imin, k99, nvalid, i
        logical :: l_sorted
        temp = posterior_temperature(1000.d0, 100)
        call assert_real(15.707963, temp, 1.e-4, 'posterior temperature pi wsqsum / (2 pftsz)')
        call posterior_row_stats(6, ROW, 1000., w, inds, imin, pmax, neff, k99, nvalid)
        call assert_int(5, nvalid, 'posterior: five valid distances')
        call assert_int(1, imin,   'posterior: the mode is the lowest distance')
        call assert_real(0.665241, pmax, 1.e-5, 'posterior: the mass of the mode')
        call assert_real(2.298819, neff, 1.e-4, 'posterior: exp(entropy)')
        call assert_int(3, k99, 'posterior: three references hold 99% of the mass')
        call assert_real(0.665241, w(6), 1.e-5, 'posterior: sorted ascending, the mode on top')
        call assert_int(1, inds(6), 'posterior: the mode''s reference on top')
        call assert_real(0.244728, w(5), 1.e-5, 'posterior: the second reference next')
        call assert_int(2, inds(5), 'posterior: its index')
        call assert_real(0., w(1), 1.e-7, 'posterior: the invalid distance carries no mass')
        l_sorted = .true.
        do i = 2,6
            if( w(i) < w(i-1) ) l_sorted = .false.
        end do
        call assert_true(l_sorted, 'posterior: the weights come back sorted ascending')
        call assert_real(1., sum(w), 1.e-5, 'posterior: the weights sum to one')
        call posterior_row_stats(6, ROW, 0., w, inds, imin, pmax, neff, k99, nvalid)
        call assert_real(0.2, pmax, 1.e-6, 'uniform posterior: the mode has 1/5')
        call assert_real(5.,  neff, 1.e-4, 'uniform posterior: exp(entropy) = 5')
        call assert_int(5, k99, 'uniform posterior: all five hold 99% of the mass')
        ! the tempering convention: a factor below one on the temperature heats (more entropy on the same row);
        ! at 500 the row's [0, 0.001, 0.002] weights are [1, e^-0.5, e^-1]/Z, pmax 0.506480 and exp(H) 2.773725 (numpy)
        call posterior_row_stats(6, ROW, 0.5 * 1000., w, inds, imin, pmax, neff, k99, nvalid)
        call assert_real(0.506480, pmax, 1.e-5, 'tempered posterior (factor 0.5): the mode''s mass falls')
        call assert_true(neff > 2.298819, 'a factor below one on the temperature raises the entropy')
        call assert_real(2.773725, neff, 1.e-4, 'tempered posterior (factor 0.5): exp(entropy)')
        call posterior_row_stats(3, [huge(1.), huge(1.), huge(1.)], 1000., w, inds, imin, pmax, neff, k99, nvalid)
        call assert_int(0, nvalid, 'no valid distance: nvalid 0')
        call assert_int(0, imin,   'no valid distance: no mode')
        call assert_real(0., neff, 1.e-7, 'no valid distance: exp(entropy) 0')
        call assert_int(0, k99,    'no valid distance: k99 0')
    end subroutine test_posterior_row_stats

    !> Two dense part files written by worker tables in batches, read by a master whose particle order is the reverse:
    !! every row, refined shift and seed comes back exactly, an unrefined entry without shift
    subroutine test_dense_part_round_trip()
        integer, parameter :: NSPACE_RT = 50, NSTATES_RT = 2, NPTCLS_RT = 40, NPARTS_RT = 2
        class(parameters), allocatable, target :: params
        class(builder),    allocatable, target :: build
        type(eul_prob_tab) :: master
        type(prob_dense_entry), allocatable :: rows(:,:)
        real,                   allocatable :: shifts(:,:,:), seeds(:,:)
        logical,                allocatable :: refined(:,:)
        integer,                allocatable :: pinds(:)
        type(prob_candidate) :: cand
        integer :: nrefs, i, ri, i_glob, nwrong
        nrefs = NSPACE_RT * NSTATES_RT
        call make_round_trip_parts(params, build, NSPACE_RT, NSTATES_RT, NPTCLS_RT, NPARTS_RT, rows, shifts, refined, seeds)
        pinds = [(i, i=NPTCLS_RT,1,-1)]
        call master%new(params, build, pinds)
        call master%read_tabs_to_glob(string('prob_tab_tester_part'), NPARTS_RT, 1)
        nwrong = 0
        do i = 1,NPTCLS_RT
            i_glob = NPTCLS_RT - i + 1
            if( any(master%seed_shifts(:,i_glob) /= seeds(:,i)) .or. .not. master%seed_has_sh(i_glob) ) nwrong = nwrong + 1
            if( master%temps(i_glob) /= 100. + real(i) ) nwrong = nwrong + 1
            do ri = 1,nrefs
                cand = master%dense_candidate(ri, i_glob)
                if( cand%iref /= ri .or. cand%dist /= rows(ri,i)%dist .or. cand%inpl /= rows(ri,i)%inpl )then
                    nwrong = nwrong + 1
                else if( refined(ri,i) )then
                    if( .not. cand%has_sh .or. cand%x /= shifts(1,ri,i) .or. cand%y /= shifts(2,ri,i) ) nwrong = nwrong + 1
                else if( cand%has_sh .or. cand%x /= 0. .or. cand%y /= 0. )then
                    nwrong = nwrong + 1
                endif
            enddo
        enddo
        call assert_int(0, nwrong, 'the dense part files read back exactly')
        call master%kill
        call kill_round_trip_parts(params, build, NPARTS_RT)
    end subroutine test_dense_part_round_trip

    !> The bounded and power samplers over an array of distances draw exactly as the ones over a distance function
    !! they replaced (kept here as reference_bounded and reference_power), from the same seed: 200 random sets of 1 to
    !! 300 distances with ties, a flat set, sentinel distances, every sample size and several angular bounds and powers
    subroutine test_bounded_and_power_samplers()
        real,    allocatable :: dists(:), pvec(:), pvec2(:), cdf(:)
        integer, allocatable :: inds(:), inds2(:)
        real    :: dist1, corr1, dist2, corr2, athres, power
        integer :: iset, n, nsample, which1, which2, i, nwrong_b, nwrong_p
        call seed_rnd_fixed(SEED_TAB)
        nwrong_b = 0
        nwrong_p = 0
        do iset = 1,200
            n = 1 + int(299. * ran3())
            allocate(dists(n), pvec(n), pvec2(n), cdf(n), inds(n), inds2(n))
            do i = 1,n
                dists(i) = 10. * ran3()
            enddo
            if( n > 4 ) dists(n/2:n/2+2) = dists(1)
            if( mod(iset,10) == 0 ) dists = 3.
            if( mod(iset,13) == 0 ) dists(1:(n+1)/2) = huge(1.)
            if( mod(iset,17) == 0 ) dists = 0.
            nsample = 1 + mod(iset-1, n)
            athres  = 1. + 179. * ran3()
            power   = 0.5 + 3. * ran3()
            call seed_rnd_fixed(SEED_DRAW + iset)
            call reference_bounded(n, athres, 10., dist1, which1, pvec, inds)
            call seed_rnd_fixed(SEED_DRAW + iset)
            call sample_bounded_dists(n, dists, athres, 10., dist2, which2, pvec2, inds2)
            if( which1 /= which2 .or. dist1 /= dist2 ) nwrong_b = nwrong_b + 1
            call seed_rnd_fixed(SEED_DRAW + iset)
            call reference_power(n, power, nsample, dist1, corr1, which1, pvec, inds)
            call seed_rnd_fixed(SEED_DRAW + iset)
            call sample_power_dists(n, dists, power, nsample, dist2, corr2, which2, pvec2, inds2, cdf)
            if( which1 /= which2 .or. dist1 /= dist2 .or. corr1 /= corr2 ) nwrong_p = nwrong_p + 1
            deallocate(dists, pvec, pvec2, cdf, inds, inds2)
        enddo
        call assert_int(0, nwrong_b, 'the array bounded sampler draws exactly as the function sampler')
        call assert_int(0, nwrong_p, 'the array power sampler draws exactly as the function sampler')

    contains

        real function dist_at( i )
            integer, intent(in) :: i
            dist_at = dists(i)
        end function dist_at

        !> the bounded sampler over a distance function, as it was
        subroutine reference_bounded( n_loc, athres_ub_in, prob_athres, dist, which, pvec_sorted, sorted_inds )
            integer, intent(in)    :: n_loc
            real,    intent(in)    :: athres_ub_in, prob_athres
            real,    intent(out)   :: dist
            integer, intent(out)   :: which
            real,    intent(inout) :: pvec_sorted(:)
            integer, intent(inout) :: sorted_inds(:)
            integer :: num_lb, num_ub, pick
            real    :: athres_ub, athres_lb, sum_pvec, rnd
            athres_ub = min(prob_athres, athres_ub_in)
            athres_lb = min(athres_ub / 10., 1.)
            num_ub    = min(n_loc,max(1,int(athres_ub * real(n_loc) / 180.)))
            num_lb    = 1 + floor(athres_lb / athres_ub * num_ub)
            call reference_select(n_loc, num_ub, pvec_sorted, sorted_inds)
            rnd      = ran3()
            sum_pvec = sum(pvec_sorted(1:num_ub))
            if( sum_pvec < TINY )then
                pick  = min(num_ub, 1 + floor(real(num_ub) * rnd))
                which = sorted_inds(pick)
            else
                pvec_sorted(1:num_ub) = pvec_sorted(1:num_ub) / sum_pvec
                if( rnd > sum(pvec_sorted(1:num_lb)) )then
                    which = sorted_inds(1)
                else
                    which = sorted_inds(num_ub)
                endif
            endif
            dist = dist_at(which)
        end subroutine reference_bounded

        !> the power sampler over a distance function, as it was
        subroutine reference_power( n_loc, pw, nsmpl, dist, corr, which, pvec_sorted, sorted_inds )
            integer, intent(in)    :: n_loc, nsmpl
            real,    intent(in)    :: pw
            real,    intent(out)   :: dist, corr
            integer, intent(out)   :: which
            real,    intent(inout) :: pvec_sorted(:)
            integer, intent(inout) :: sorted_inds(:)
            real    :: cdf_loc(max(1,nsmpl)), rnd, csum
            integer :: num_ub, j, rank, pick
            num_ub = min(n_loc,max(1,nsmpl))
            call reference_select(n_loc, num_ub, pvec_sorted, sorted_inds)
            do j = 1,num_ub
                cdf_loc(j) = exp(-pvec_sorted(num_ub-j+1))
            enddo
            if( all(cdf_loc(1:num_ub) < TINY) )then
                pick  = 1
                which = sorted_inds(pick)
                dist  = pvec_sorted(pick)
                corr  = cdf_loc(num_ub)
                return
            endif
            where( cdf_loc(1:num_ub) < TINY ) cdf_loc(1:num_ub) = 0.
            do j = 2,num_ub
                cdf_loc(j) = cdf_loc(j) + cdf_loc(j-1)
            enddo
            csum = cdf_loc(num_ub)
            cdf_loc(1:num_ub) = cdf_loc(1:num_ub) / csum
            rnd = ran3()
            rnd = min(1.,max(0.,1.-rnd**pw))
            rank = 0
            do j = 1,num_ub
                if( cdf_loc(j) > rnd )then
                    rank = j
                    exit
                endif
            enddo
            if( rank == 0 ) rank = num_ub
            pick  = num_ub - rank + 1
            which = sorted_inds(pick)
            dist  = pvec_sorted(pick)
            corr  = exp(-dist)
        end subroutine reference_power

        !> the num_ub lowest distances of the function, as the samplers selected them (bounded max-heap, sort)
        subroutine reference_select( n_loc, num_ub, pvec_sorted, sorted_inds )
            integer, intent(in)    :: n_loc, num_ub
            real,    intent(inout) :: pvec_sorted(:)
            integer, intent(inout) :: sorted_inds(:)
            real    :: dist_tmp
            integer :: j
            do j = 1,num_ub
                pvec_sorted(j) = dist_at(j)
                sorted_inds(j) = j
            enddo
            do j = num_ub / 2, 1, -1
                call heap_sift_down(pvec_sorted, sorted_inds, j, num_ub)
            enddo
            do j = num_ub + 1,n_loc
                dist_tmp = dist_at(j)
                if( dist_tmp < pvec_sorted(1) )then
                    pvec_sorted(1) = dist_tmp
                    sorted_inds(1) = j
                    call heap_sift_down(pvec_sorted, sorted_inds, 1, num_ub)
                endif
            enddo
            call hpsort(pvec_sorted(1:num_ub), sorted_inds(1:num_ub))
        end subroutine reference_select

        subroutine heap_sift_down( pv, si, root_in, heap_size )
            real,    intent(inout) :: pv(:)
            integer, intent(inout) :: si(:)
            integer, intent(in)    :: root_in, heap_size
            integer :: root, child, swap_i, tmp_ind
            real    :: tmp_val
            root = root_in
            do
                child = 2 * root
                if( child > heap_size ) exit
                swap_i = root
                if( pv(swap_i) < pv(child) ) swap_i = child
                if( child + 1 <= heap_size )then
                    if( pv(swap_i) < pv(child + 1) ) swap_i = child + 1
                endif
                if( swap_i == root ) exit
                tmp_val    = pv(root)
                pv(root)   = pv(swap_i)
                pv(swap_i) = tmp_val
                tmp_ind    = si(root)
                si(root)   = si(swap_i)
                si(swap_i) = tmp_ind
                root = swap_i
            enddo
        end subroutine heap_sift_down

    end subroutine test_bounded_and_power_samplers

    !> The parallel read of dense part files at one thread, at four and, for fourteen parts, 24 times at every processor:
    !! the same rows, refined lists and seeds, bit for bit, with changing reference counts
    subroutine test_dense_read_threads()
        integer, parameter :: NSPACE_TH = 40, NSTATES_TH = 3, NPTCLS_TH = 70
        integer, parameter :: NPARTS_TH(4) = [1, 3, 7, 14]
        class(parameters), allocatable, target :: params
        class(builder),    allocatable, target :: build
        type(eul_prob_tab) :: master1, master4
        type(prob_dense_entry), allocatable :: rows(:,:)
        real,                   allocatable :: shifts(:,:,:), seeds(:,:)
        logical,                allocatable :: refined(:,:)
        integer,                allocatable :: pinds(:)
        integer :: i, k, nthr_prev, nwrong, nthr_read, iread, nreads
        character(len=32) :: label
        nthr_prev = 1
        !$ nthr_prev = omp_get_max_threads()
        do k = 1,size(NPARTS_TH)
            write(label,'(I0,A)') NPARTS_TH(k), ' parts'
            call make_round_trip_parts(params, build, k*NSPACE_TH, NSTATES_TH, NPTCLS_TH, NPARTS_TH(k), rows, shifts,&
                &refined, seeds)
            ! the master's particles in a scrambled order
            pinds = [(mod(17*i, NPTCLS_TH) + 1, i=0,NPTCLS_TH-1)]
            !$ call omp_set_num_threads(1)
            call master1%new(params, build, pinds)
            call master1%read_tabs_to_glob(string('prob_tab_tester_part'), NPARTS_TH(k), 1)
            nthr_read = 4
            nreads    = 1
            if( NPARTS_TH(k) == 14 )then
                !$ nthr_read = omp_get_num_procs()
                nreads = 24
            endif
            !$ call omp_set_num_threads(nthr_read)
            call master4%new(params, build, pinds)
            nwrong = 0
            do iread = 1,nreads
                call master4%read_tabs_to_glob(string('prob_tab_tester_part'), NPARTS_TH(k), 1)
                nwrong = nwrong + count(master4%dists /= transpose(rows(:,pinds)%dist)) +&
                    &count(master4%inpls /= transpose(rows(:,pinds)%inpl)) +&
                    &count(master4%seed_shifts /= seeds(:,pinds)) + count(.not. master4%seed_has_sh) +&
                    &count_differences(master1, master4)
            enddo
            !$ call omp_set_num_threads(nthr_prev)
            call assert_int(0, nwrong, 'the parallel reads equal the generated rows and the read at one thread ('//&
                &trim(label)//')')
            call assert_int(36, master4%seed_nrots, 'the rotation grid of the parts is taken ('//trim(label)//')')
            call master1%kill
            call master4%kill
            call kill_round_trip_parts(params, build, NPARTS_TH(k))
        enddo
    end subroutine test_dense_read_threads

    !> the fields of two masters' read that differ: rows, refined lists, offsets and seeds
    integer function count_differences( a, b ) result( nwrong )
        type(eul_prob_tab), intent(in) :: a, b
        nwrong = count(a%dists /= b%dists) + count(a%inpls /= b%inpls) +&
            &count(a%refined_offsets /= b%refined_offsets) + count(a%seed_shifts /= b%seed_shifts) +&
            &count(a%seed_has_sh .neqv. b%seed_has_sh)
        if( size(a%refined) /= size(b%refined) )then
            nwrong = nwrong + 1
        else
            nwrong = nwrong + count(a%refined(:)%i /= b%refined(:)%i) + count(a%refined(:)%ri /= b%refined(:)%ri) +&
                &count(a%refined(:)%x /= b%refined(:)%x) + count(a%refined(:)%y /= b%refined(:)%y)
        endif
    end function count_differences

    !> NPARTS_RT dense part files prob_tab_tester_part<n>.dat of consecutive particles, written by worker tables in
    !! two batches each: random rows, a refined shift on about one entry in ten, a seed per particle
    subroutine make_round_trip_parts( params, build, nspace, nstates, nptcls, nparts, rows, shifts, refined, seeds )
        class(parameters), allocatable, target, intent(inout) :: params
        class(builder),    allocatable, target, intent(inout) :: build
        integer,                                intent(in)    :: nspace, nstates, nptcls, nparts
        type(prob_dense_entry), allocatable,    intent(out)   :: rows(:,:)
        real,                   allocatable,    intent(out)   :: shifts(:,:,:), seeds(:,:)
        logical,                allocatable,    intent(out)   :: refined(:,:)
        type(eul_prob_tab) :: worker
        type(prob_refined_entry), allocatable :: batch_refined(:)
        integer, allocatable :: parts(:,:)
        integer :: nrefs, i, ri, ipart, i_loc, nloc, ibatch, first, last, n
        nrefs = nspace * nstates
        allocate(params, build)
        params%nspace      = nspace
        params%nstates     = nstates
        params%prob_athres = 10.
        params%l_doshift   = .true.
        ! the point group and the direction space, as the assignment's masters build them (build_eulspace): the
        ! weighted insertion measures its candidates' poses in the direction space
        call build%pgrpsyms%new(trim(params%pgrp))
        call build%build_eulspace(params)
        call seed_rnd_fixed(SEED_TAB)
        allocate(build%spproj_field)
        call build%spproj_field%new(nptcls, is_ptcl=.true.)
        do i = 1,nptcls
            call build%spproj_field%set_state(i, mod(i-1, nstates) + 1)
        enddo
        allocate(rows(nrefs,nptcls), shifts(2,nrefs,nptcls), refined(nrefs,nptcls), seeds(2,nptcls))
        do i = 1,nptcls
            seeds(:,i) = [ran3() - 0.5, ran3() - 0.5]
            do ri = 1,nrefs
                rows(ri,i)     = prob_dense_entry(5. * ran3(), 1 + int(36. * ran3()))
                refined(ri,i) = ran3() < 0.1
                shifts(:,ri,i) = [ran3() - 0.5, ran3() - 0.5]
            enddo
        enddo
        parts = split_nobjs_even(nptcls, nparts)
        do ipart = 1,nparts
            nloc = parts(ipart,2) - parts(ipart,1) + 1
            call worker%new_worker(params, build, [(i, i=parts(ipart,1),parts(ipart,2))])
            worker%seed_nrots = 36
            worker%seed_shifts = seeds(:,parts(ipart,1):parts(ipart,2))
            worker%seed_has_sh = .true.
            worker%temps = [(100. + real(i), i=parts(ipart,1),parts(ipart,2))]
            call worker%begin_write(string('prob_tab_tester_part')//int2str(ipart)//'.dat')
            do ibatch = 1,2
                first = merge(1, nloc/2 + 1, ibatch == 1)
                last  = merge(nloc/2, nloc, ibatch == 1)
                allocate(batch_refined(count(refined(:,parts(ipart,1)+first-1:parts(ipart,1)+last-1))))
                n = 0
                do i_loc = first,last
                    i = parts(ipart,1) + i_loc - 1
                    do ri = 1,nrefs
                        if( .not. refined(ri,i) ) cycle
                        n = n + 1
                        batch_refined(n) = prob_refined_entry(i_loc, ri, shifts(1,ri,i), shifts(2,ri,i))
                    enddo
                enddo
                call worker%write_rows(first, rows(:,parts(ipart,1)+first-1:parts(ipart,1)+last-1), batch_refined)
                deallocate(batch_refined)
            enddo
            call worker%write_tab(string('prob_tab_tester_part')//int2str(ipart)//'.dat')
            call worker%kill
        enddo
    end subroutine make_round_trip_parts

    subroutine kill_round_trip_parts( params, build, nparts )
        class(parameters), allocatable, target, intent(inout) :: params
        class(builder),    allocatable, target, intent(inout) :: build
        integer,                                intent(in)    :: nparts
        integer :: ipart
        do ipart = 1,nparts
            call del_file(string('prob_tab_tester_part')//int2str(ipart)//'.dat')
        enddo
        call build%spproj_field%kill
        deallocate(build%spproj_field)
        call build%eulspace%kill
        call build%pgrpsyms%kill
        deallocate(params, build)
    end subroutine kill_round_trip_parts

    subroutine set_dist( tab, iref, i, dist )
        type(eul_prob_tab), intent(inout) :: tab
        integer,            intent(in)    :: iref, i
        real,               intent(in)    :: dist
        tab%dists(i,iref) = dist
        tab%inpls(i,iref) = 0
    end subroutine set_dist

    subroutine test_record_seed_shift()
        type(ptcl_ref) :: assgn
        integer :: inpl
        write(*,'(A)') 'test_record_seed_shift'
        ! the particle-frame seed is recorded whatever the angle
        do inpl = 1, NROTS
            assgn = unrefined(inpl)
            call record_seed_shift(assgn, SEED, .true., .true.)
            call assert_true(assgn%has_sh .and. assgn%x == SEED(1) .and. assgn%y == SEED(2), &
                &'the seed is recorded as it is, at every in-plane angle')
        enddo
        ! a refined candidate keeps its own shift, a particle without seed or without shift search is left alone
        assgn        = unrefined(3)
        assgn%has_sh = .true.
        assgn%x      = 0.25
        assgn%y      = -0.75
        call record_seed_shift(assgn, SEED, .true., .true.)
        call assert_true(assgn%x == 0.25 .and. assgn%y == -0.75, 'a refined candidate keeps the shift of its search')
        assgn = unrefined(3)
        call record_seed_shift(assgn, SEED, .false., .true.)
        call assert_false(assgn%has_sh, 'no seed, nothing recorded')
        call record_seed_shift(assgn, SEED, .true., .false.)
        call assert_false(assgn%has_sh, 'no shift search, nothing recorded')
    end subroutine test_record_seed_shift

    function unrefined( inpl ) result( assgn )
        integer, intent(in) :: inpl
        type(ptcl_ref) :: assgn
        assgn%inpl = inpl
    end function unrefined

end module simple_eul_prob_tab_tester
