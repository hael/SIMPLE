!@descr: unit tests for the 3D probability tables: the seed shift of an assigned candidate, the in-plane samplers, the posterior draw of the dense and of the neighbourhood table on the pipeline's batch, the weighted insertion (the top extras, the cone, the table), the n_eff readout
! The seed table holds the particle-frame seed every candidate of a particle was scored with; an assigned, unrefined
! candidate records it as it is, whatever its in-plane angle. A refined candidate keeps the shift of its search.
! The pipeline's batch assignment draws every particle's reference from its tempered row posterior on synthetic dense
! batches of one and of several states, with ties, a flat and an invalid particle. A particle's state is the state of its closest
! reference; an exact tie goes to the lowest reference index.
module simple_eul_prob_tab_tester
use simple_core_module_api
use simple_builder,            only: builder
use simple_parameters,         only: parameters
use simple_eul_prob_tab,       only: eul_prob_tab
use simple_eul_prob_tab_neigh, only: eul_prob_tab_neigh
use simple_eul_prob_tab_utils, only: record_seed_shift, &
    &sample_bounded_dists, sample_power_dists, prob_candidate, prob_dense_entry, prob_refined_entry, &
    &posterior_temperature, posterior_row_stats, draw_posterior_row, top_posterior_extras, prob_beta_at, prob_neff_readout, &
    &prob_insert_table, prob_insert_pose, prob_insert_cone, &
    &PROB_NEFF_MIN_ROWS, map_prior_weight, insert_pose_loss, update_screen_best
use simple_rnd,                only: seed_rnd_fixed, ran3, rnd_key, rnd_keyed
use, intrinsic :: iso_fortran_env, only: int64
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
        call test_screen_width_shortlist()
        write(*,'(A)') 'test_bounded_and_power_samplers'
        call test_bounded_and_power_samplers()
        write(*,'(A)') 'test_posterior_row_stats'
        call test_posterior_row_stats()
        write(*,'(A)') 'test_draw_posterior_row'
        call test_draw_posterior_row()
        write(*,'(A)') 'test_top_posterior_extras'
        call test_top_posterior_extras()
        write(*,'(A)') 'test_insert_pose_loss'
        call test_insert_pose_loss()
        write(*,'(A)') 'test_insert_cone'
        call test_insert_cone()
        write(*,'(A)') 'test_insert_table'
        call test_insert_table()
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
        call test_neff_report()
        call fix_test_seed()
    end subroutine run_all_eul_prob_tab_tests

    subroutine test_screen_width_shortlist()
        integer, parameter :: N = 12, K = 3
        real :: losses(N), best(K), work(N), cutoff
        logical :: measured(N)
        integer :: icase, i
        do icase = 1, 4
            select case(icase)
                case(1)
                    losses = [(real(i), i=1,N)]
                case(2)
                    losses = [(real(N+1-i), i=1,N)]
                case(3)
                    losses = 2.
                case(4)
                    losses = [5., 5., 5., 1., 7., 2., 2., 0., 8., 2., 0., 0.]
            end select
            best = huge(1.)
            do i = 1, N
                measured(i) = losses(i) <= maxval(best)
                call update_screen_best(losses(i), best)
            enddo
            work = losses
            do i = 1, K-1
                work(minloc(work,dim=1)) = huge(1.)
            enddo
            cutoff = minval(work)
            call assert_true(maxval(best) == cutoff, 'running cutoff equals the independent kth-smallest loss')
            call assert_true(all(measured .or. losses > cutoff), 'every candidate at or below the final cutoff has a width')
            if( icase == 1 ) call assert_int(K, count(measured), 'ascending losses only measure the initial shortlist')
            if( icase == 3 ) call assert_true(all(measured), 'ties are never pruned by the width shortcut')
        enddo
    end subroutine test_screen_width_shortlist

    subroutine test_neff_report()
        use, intrinsic :: iso_fortran_env, only: iostat_end
        type(prob_neff_readout) :: neff_ro
        real :: row(50)
        character(len=256) :: line
        integer :: r, saved_log, scratch_log, ios
        write(*,'(A)') 'test_neff_report'
        open(newunit=scratch_log, status='scratch', action='readwrite', iostat=ios)
        call assert_int(0, ios, 'n_eff report: open capture')
        if( ios /= 0 ) return
        saved_log = logfhandle
        row = 0.
        call neff_ro%new(PROB_NEFF_MIN_ROWS, PROB_NEFF_MIN_ROWS, 6)
        do r = 1,PROB_NEFF_MIN_ROWS
            call neff_ro%offer(size(row), row, 1.)
        enddo
        logfhandle = scratch_log
        call neff_ro%report(0.3)
        logfhandle = saved_log
        rewind(scratch_log)
        read(scratch_log,'(A)',iostat=ios) line
        call assert_int(0, ios, 'n_eff report: a line at the minimum row count')
        if( ios == 0 ) call assert_char('>>> POSTERIOR DRAW N_EFF: 50.00 at beta 0.30, 50.00 at beta 1.00 ('//&
            &int2str(PROB_NEFF_MIN_ROWS)//' rows)', trim(line), 'n_eff report: unpadded numbers with leading zeros')
        call neff_ro%kill
        close(scratch_log)
        open(newunit=scratch_log, status='scratch', action='readwrite', iostat=ios)
        call assert_int(0, ios, 'n_eff report: open below-threshold capture')
        if( ios /= 0 ) return
        call neff_ro%new(PROB_NEFF_MIN_ROWS - 1, PROB_NEFF_MIN_ROWS - 1, 6)
        do r = 1,PROB_NEFF_MIN_ROWS - 1
            call neff_ro%offer(size(row), row, 1.)
        enddo
        logfhandle = scratch_log
        call neff_ro%report(0.3)
        logfhandle = saved_log
        rewind(scratch_log)
        read(scratch_log,'(A)',iostat=ios) line
        call assert_int(iostat_end, ios, 'n_eff report: silent below the minimum row count')
        call neff_ro%kill
        close(scratch_log)
    end subroutine test_neff_report

    !> The n_eff readout: 200 rows of 50 distances k s_r (k = 0..49, s_r = 0.01 r) at temperature 2, so row r's
    !! posterior at factor beta is geometric with ratio exp(-2 beta s_r). Reference values from a numpy emulation:
    !! one row at a = 0.5 has exp(H) = 5.493135; the mean over the rows is 49.9721 at beta 0.001, 6.023447
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

    !> the neighbourhood table's posterior assignment (the pipeline's batch; the same draw on every
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
        call set_fixed_seed(SEED_DRAW, propagate=.true.)
        call assign_store
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
        call set_fixed_seed(SEED_DRAW, propagate=.true.)
        call assign_store
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
        call set_fixed_seed(SEED_DRAW, propagate=.true.)
        call assign_store
        run2 = (tab%assgn_map(:)%istate - 1) * NSPACE + tab%assgn_map(:)%iproj
        call assert_int(0, count(run2 /= run1), 'the same seed gives the same neighbourhood posterior assignment')
        call tab%kill
        call set_store(2.)
        call set_fixed_seed(SEED_DRAW + 1, propagate=.true.)
        call assign_store
        run2 = (tab%assgn_map(:)%istate - 1) * NSPACE + tab%assgn_map(:)%iproj
        call assert_true(any(run2 /= run1), 'another seed gives another neighbourhood posterior assignment')
        call tab%kill
        call build%spproj_field%kill
        deallocate(params, build, cands, counts, pinds, state_best, ref_best, dist_best, run1, run2)

    contains

        !> the pipeline's table over the particles, every particle's candidates appended to the first thread
        !! buffer in their order (as the fill appends them), the temperature set
        subroutine set_store( temp )
            real, intent(in) :: temp
            integer :: ip, kk, pos
            call tab%new_neigh_pipeline(params, build, pinds)
            pos = 0
            do ip = 1,NP
                do kk = 1,counts(ip)
                    pos = pos + 1
                    call tab%candidate_buffers(1)%append(ip, cands(pos))
                enddo
            enddo
            tab%temps = temp
        end subroutine set_store

        !> the one batch of every row assigned
        subroutine assign_store
            call tab%open_assignment
            call tab%assign_batch(1, NP)
            call tab%close_assignment
        end subroutine assign_store

    end subroutine test_ref_assign_neigh_posterior

    !> The draw from a row posterior (draw_posterior_row) with its variate given: against an independent
    !! cumulative-posterior oracle on hand-built rows (the mode's mass, ties at the minimum, an invalid entry,
    !! a row without a valid entry, the variate at the bin boundaries and at its ends), and over 20000 keyed
    !! variates the frequencies of [0, 1, 2] at T = 1 are [1, e^-1, e^-2]/Z = [0.665, 0.245, 0.090] within
    !! 0.02, temperature 0 draws uniformly over the valid distances, a huge temperature the argmin; the same
    !! variate gives the same draw
    subroutine test_draw_posterior_row()
        integer, parameter :: NDRAWS = 20000
        real     :: dists(6), w(6), dist
        real(dp) :: u, cdf(6)
        integer  :: which, i, counts(6), k
        logical  :: l_ok
        type(rnd_key) :: key
        ! the oracle: the cumulative posterior of a row at temperature T over its valid entries
        dists = [0.3, 0.1, 0.15, 0.5, huge(1.), 0.2]
        call oracle_cdf(6, dists, 1., cdf)
        l_ok = .true.
        do i = 1,6
            if( .not. valid(dists(i)) ) cycle
            ! a variate just below the entry's upper bound draws the entry; at the bound the next one
            u = cdf(i) - 1.d-9
            call draw_posterior_row(6, dists, 1., w, u, which, dist)
            l_ok = l_ok .and. which == i .and. dist == dists(i)
        enddo
        call assert_true(l_ok, 'a variate below an entry''s cumulative bound draws that entry')
        call draw_posterior_row(6, dists, 1., w, 0.d0, which, dist)
        call assert_int(1, which, 'the variate 0 draws the first valid entry')
        call draw_posterior_row(6, dists, 1., w, 1.d0 - 2.d0**(-53), which, dist)
        call assert_int(6, which, 'the largest variate draws the last valid entry')
        call draw_posterior_row(6, dists, 1., w, cdf(1), which, dist)
        call assert_int(1, which, 'a variate at the first entry''s bound draws the first (the bound belongs to its bin)')
        ! huge temperature, a unique minimum: the argmin at every variate
        l_ok = .true.
        do i = 0,50
            u = real(i,dp) / 51.d0
            call draw_posterior_row(6, dists, huge(1.), w, u, which, dist)
            l_ok = l_ok .and. which == 2 .and. dist == 0.1
        enddo
        call assert_true(l_ok, 'a huge temperature draws a unique argmin at every variate')
        ! huge temperature, two equal minima: either, by the variate's half, never another
        dists = [0.3, 0.1, 0.1, 0.5, huge(1.), 0.2]
        counts = 0
        do i = 0,49
            u = (real(i,dp) + 0.5d0) / 50.d0
            call draw_posterior_row(6, dists, huge(1.), w, u, which, dist)
            counts(which) = counts(which) + 1
        enddo
        call assert_int(50, counts(2) + counts(3), 'a huge temperature draws only the equal minima')
        call assert_int(25, counts(2), 'a huge temperature draws the first of two equal minima on the lower half of the variates')
        call assert_int(25, counts(3), 'and the second on the upper half')
        ! temperature 0 over keyed variates: uniform over the five valid distances, never the invalid one
        key%seed = 20261104_int64; key%family = 7
        counts = 0
        do i = 1,NDRAWS
            key%ordinal = i
            call draw_posterior_row(6, dists, 0., w, rnd_keyed(key), which, dist)
            counts(which) = counts(which) + 1
        enddo
        call assert_int(0, counts(5), 'temperature 0 never draws an invalid distance')
        call assert_true(all(counts([1,2,3,4,6]) > NDRAWS/5 - 300) .and. all(counts([1,2,3,4,6]) < NDRAWS/5 + 300), &
            &'temperature 0 draws the valid distances uniformly over keyed variates')
        ! T = 1 on [0, 1, 2] over keyed variates
        counts = 0
        do i = 1,NDRAWS
            key%ordinal = i
            call draw_posterior_row(3, [0., 1., 2.], 1., w, rnd_keyed(key), which, dist)
            counts(which) = counts(which) + 1
        enddo
        call assert_true(abs(real(counts(1))/NDRAWS - 0.665241) < 0.02, 'the mode is drawn with its posterior mass')
        call assert_true(abs(real(counts(2))/NDRAWS - 0.244728) < 0.02, 'the second is drawn with e^-1 of the mode''s mass')
        call assert_true(abs(real(counts(3))/NDRAWS - 0.090031) < 0.02, 'the third is drawn with e^-2 of the mode''s mass')
        ! the same variate gives the same draw, whatever was drawn before
        call draw_posterior_row(3, [0., 1., 2.], 1., w, 0.7d0, which, dist)
        k = which
        call seed_rnd_fixed(20261024)
        u = real(ran3(),dp)
        call draw_posterior_row(3, [0., 1., 2.], 1., w, 0.7d0, which, dist)
        call assert_int(k, which, 'the same variate gives the same draw')
        call assert_int(2, which, 'the variate 0.7 on [0, 1, 2] at T = 1 draws the second (the mode holds 0.665)')
        ! no valid distance
        call draw_posterior_row(3, [huge(1.), huge(1.), huge(1.)], 1., w, 0.5d0, which, dist)
        call assert_int(0, which, 'a row without a valid distance draws 0')

    contains

        logical function valid( d )
            real, intent(in) :: d
            valid = d < 0.1 * huge(d)
        end function valid

        !> the oracle: p_r = exp(-T (d_r - d_min)) over the valid entries, zero elsewhere, each weight a single
        !! precision number as the selector's, summed in double; cdf(r) = the cumulative mass up to and including
        !! r, over the total
        subroutine oracle_cdf( n, d, temp, c )
            integer,  intent(in)  :: n
            real,     intent(in)  :: d(n), temp
            real(dp), intent(out) :: c(n)
            real     :: p(n), dmin
            real(dp) :: z, arg
            integer  :: r
            dmin = huge(1.)
            do r = 1,n
                if( valid(d(r)) ) dmin = min(dmin, d(r))
            enddo
            p = 0.
            do r = 1,n
                if( .not. valid(d(r)) ) cycle
                arg = real(max(0.,temp),dp) * real(d(r) - dmin,dp)
                if( arg < 80.d0 ) p(r) = real(exp(-arg))
            enddo
            z = 0.d0
            do r = 1,n
                z = z + real(p(r),dp)
            enddo
            c(1) = real(p(1),dp) / z
            do r = 2,n
                c(r) = c(r-1) + real(p(r),dp) / z
            enddo
        end subroutine oracle_cdf

    end subroutine test_draw_posterior_row

    !> the pipeline's batch assignment on synthetic dense tables (one state of 200 directions, three states of 200, two states of one
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
        call tab%new_pipeline(params, build, pinds)
        ! a huge temperature: the argmin wherever the row minimum is unique (every row but the flat particle 1)
        tab%temps = huge(1.)
        call set_fixed_seed(SEED_DRAW, propagate=.true.)
        call assign_rows_all
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
        call set_fixed_seed(SEED_DRAW, propagate=.true.)
        call assign_rows_all
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
        call set_fixed_seed(SEED_DRAW, propagate=.true.)
        call assign_rows_all
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
        call set_fixed_seed(SEED_DRAW, propagate=.true.)
        call assign_rows_all
        do i = 1,NPTCLS_TAB
            iref_run2(i) = assigned_ref(tab, i, nspace)
        enddo
        call assert_true(all(iref_run1 == iref_run2), 'the same seed gives the same assignment ('//trim(label)//')')
        if( nspace > 1 )then
            call set_fixed_seed(SEED_DRAW + 1, propagate=.true.)
            call assign_rows_all
            do i = 1,NPTCLS_TAB
                iref_run3(i) = assigned_ref(tab, i, nspace)
            enddo
            call assert_true(any(iref_run1 /= iref_run3), 'another seed gives another assignment ('//trim(label)//')')
        endif
        call tab%kill
        call build%spproj_field%kill
        deallocate(build%spproj_field)
        deallocate(params, build)

    contains

        !> every row as one batch of the pipeline's table (the in-plane index 0, no refined entry), assigned
        subroutine assign_rows_all
            allocate(tab%batch_rows(nrefs,NPTCLS_TAB), tab%batch_refined(1,NPTCLS_TAB), tab%batch_nrefined(NPTCLS_TAB))
            do i = 1,NPTCLS_TAB
                do iref = 1,nrefs
                    tab%batch_rows(iref,i) = prob_dense_entry(dists(iref,i), 0)
                enddo
            enddo
            tab%batch_nrefined = 0
            call tab%open_assignment
            call tab%assign_batch(1, NPTCLS_TAB)
            call tab%close_assignment
        end subroutine assign_rows_all

    end subroutine check_posterior_assignment

    !> the reference index (state block, direction) a particle was assigned
    integer function assigned_ref( tab, i, nspace ) result( iref )
        type(eul_prob_tab), intent(in) :: tab
        integer,            intent(in) :: i, nspace
        iref = (tab%assgn_map(i)%istate - 1) * nspace + tab%assgn_map(i)%iproj
    end function assigned_ref

    !> The weighted insertion's table: a particle without extras keeps weight 1; add records the committed weight
    !! and the extras in order, and counts the cone's occupancy (particle 2 had 5 eligible entries, above the cap
    !! of 2; particle 4 had 1; the mean over the four particles added is 1.5 and a quarter of them were capped);
    !! the readout line prints
    subroutine test_insert_table()
        type(prob_insert_table) :: tab
        type(prob_insert_pose)  :: pose
        type(prob_candidate)    :: cands(2)
        call tab%new(4, [11, 12, 13, 14], 3)
        cands(1) = prob_candidate(101, 7, 0.4, 0.5, 1.5, .true.)
        cands(2) = prob_candidate(102, 8, 0.6, -1.0, 2.0, .true.)
        call tab%add(2, 0.6, 2, cands, [0.3, 0.1], 5)
        cands(1) = prob_candidate(103, 9, 0.2, 0.25, 0.0, .true.)
        call tab%add(4, 0.9, 1, cands(1:1), [0.1], 1)
        call tab%add(1, 1.0, 0, cands(1:0), [real ::], 0)
        call tab%add(3, 1.0, 0, cands(1:0), [real ::], 0)
        call assert_int(4, tab%get_nptcls(), 'insert table: four particles')
        call assert_int(3, tab%get_kmax(),   'insert table: kmax 3')
        call assert_real(1.5,  tab%get_mean_eligible(), 1.e-6, 'insert table: the eligible entries per particle')
        call assert_real(0.25, tab%get_frac_capped(),   1.e-6, 'insert table: the fraction of particles above the cap')
        call assert_real(1., tab%get_committed_weight(1), 0., 'insert table: a particle without extras keeps weight 1')
        call assert_int(0, tab%get_nextra(1), 'insert table: and no extras')
        call assert_real(0.6, tab%get_committed_weight(2), 0., 'insert table: the committed weight recorded')
        call assert_int(2, tab%get_nextra(2), 'insert table: two extras recorded')
        pose = tab%get_extra(2, 2)
        call assert_true(pose%iref == 102 .and. pose%inpl == 8 .and. pose%x == -1.0 .and. pose%y == 2.0 .and. pose%w == 0.1, &
            &'insert table: the second extra of particle 2')
        call tab%report
        call tab%kill
    end subroutine test_insert_table

    !> The weighted insertion's extras on a two-state dense batch (50 directions, 2 states, 40 particles in reverse
    !! order, a tenth of the entries refined, their distances lowered by 50 so that the committed draw lands on a
    !! refined entry, the condition for extras): every extra lies in the committed state's block, is a refined entry
    !! and is not the committed reference; a particle's weights sum to 1, its extras come in descending weight and
    !! none is below PROB_INSERT_WMIN; extras do occur at a hot temperature; a second assignment gives the same
    !! extras (the selection is deterministic); topk = 1 records no extra and weight 1 everywhere
    subroutine test_insert_extras_multistate()
        integer, parameter :: NSPACE_MS = 50, NSTATES_MS = 2, NPTCLS_MS = 40
        class(parameters), allocatable, target :: params
        class(builder),    allocatable, target :: build
        type(eul_prob_tab)      :: master
        type(prob_insert_pose)  :: pose, pose_a
        type(prob_insert_cone)  :: cone
        type(prob_dense_entry), allocatable :: rows(:,:)
        real,                   allocatable :: shifts(:,:,:), seeds(:,:)
        logical,                allocatable :: refined(:,:)
        integer,                allocatable :: pinds(:), nextra_first(:)
        integer :: i, k, nextra_tot, nwrong_state, nwrong_refined, nwrong_committed, nwrong_weight, nwrong_dup, iref_committed
        integer :: istate, nwrong_order, nwrong_again
        real    :: wsum
        call make_dense_rows(params, build, NSPACE_MS, NSTATES_MS, NPTCLS_MS, rows, shifts, refined, seeds)
        pinds = [(i, i=NPTCLS_MS,1,-1)]
        params%topk = 3
        call master%new_pipeline(params, build, pinds)
        master%seed_nrots = 36
        master%temps = 1. ! hot rows: extras occur
        ! no cone set: no extra pose, the committed weight 1 everywhere
        call set_fixed_seed(SEED_DRAW, propagate=.true.)
        call assign_master
        nextra_tot = 0
        do i = 1,NPTCLS_MS
            nextra_tot = nextra_tot + master%insert_tab%get_nextra(i)
            if( master%insert_tab%get_committed_weight(i) /= 1. ) nextra_tot = nextra_tot + 1
        enddo
        call assert_int(0, nextra_tot, 'no insertion cone: no extra pose and the committed weight 1 everywhere')
        ! a cone of a ten-thousandth of an Angstrom: no pose but the committed one passes, so no extra pose
        call cone%new([1.e-4, 1.e-4], 100., 1.)
        call master%set_insert_cone(cone)
        call set_fixed_seed(SEED_DRAW, propagate=.true.)
        call assign_master
        nextra_tot = 0
        do i = 1,NPTCLS_MS
            nextra_tot = nextra_tot + master%insert_tab%get_nextra(i)
        enddo
        call assert_int(0, nextra_tot, 'a cone of the committed pose alone: no extra pose')
        ! a resolution past the rim: the whole sphere, every in-plane angle and shift pass, the extras are the state's
        ! refined entries
        call cone%new([1000., 1000.], 100., 1.)
        call master%set_insert_cone(cone)
        call set_fixed_seed(SEED_DRAW, propagate=.true.)
        call assign_master
        nextra_tot = 0; nwrong_state = 0; nwrong_refined = 0; nwrong_committed = 0; nwrong_weight = 0; nwrong_dup = 0
        nwrong_order = 0
        allocate(nextra_first(NPTCLS_MS), source=0)
        do i = 1,NPTCLS_MS
            istate          = master%assgn_map(i)%istate
            iref_committed  = (istate - 1) * NSPACE_MS + master%assgn_map(i)%iproj
            wsum            = master%insert_tab%get_committed_weight(i)
            nextra_first(i) = master%insert_tab%get_nextra(i)
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
                    if( pose%iref == pose_a%iref ) nwrong_dup   = nwrong_dup + 1
                    if( pose%w    >  pose_a%w    ) nwrong_order = nwrong_order + 1
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
        call assert_int(0, nwrong_order,     'multi-state extras: the extras of a particle come in descending weight')
        call assert_int(0, nwrong_weight,    'multi-state extras: weights sum to 1, no extra below 1%, at most topk - 1 extras')
        call assert_true(master%insert_tab%get_mean_eligible() > 0., 'multi-state extras: the cone''s occupancy is counted')
        ! the same assignment again: the same extras (the selection is deterministic, the draws reseeded)
        call set_fixed_seed(SEED_DRAW, propagate=.true.)
        call assign_master
        nwrong_again = 0
        do i = 1,NPTCLS_MS
            if( master%insert_tab%get_nextra(i) /= nextra_first(i) ) nwrong_again = nwrong_again + 1
        enddo
        call assert_int(0, nwrong_again, 'multi-state extras: a second assignment selects the same extras')
        deallocate(nextra_first)
        ! topk 1: the committed pose alone
        params%topk = 1
        call set_fixed_seed(SEED_DRAW, propagate=.true.)
        call assign_master
        nextra_tot = 0
        do i = 1,NPTCLS_MS
            nextra_tot = nextra_tot + master%insert_tab%get_nextra(i)
            if( master%insert_tab%get_committed_weight(i) /= 1. ) nextra_tot = nextra_tot + 1
        enddo
        call assert_int(0, nextra_tot, 'topk 1: no extra pose and the committed weight 1 everywhere')
        call master%kill
        call cone%kill
        call build%spproj_field%kill
        deallocate(build%spproj_field)
        call build%eulspace%kill
        call build%pgrpsyms%kill
        deallocate(params, build, rows, shifts, refined, seeds, pinds)

    contains

        !> the batch of every row from the synthetic rows (in the table's particle order), the refined entries with
        !! their shifts and their distances lowered by 50, the seeds; assigned
        subroutine assign_master
            integer :: ii, ri, n, nrefs
            nrefs = NSPACE_MS * NSTATES_MS
            allocate(master%batch_rows(nrefs,NPTCLS_MS), master%batch_refined(nrefs,NPTCLS_MS), master%batch_nrefined(NPTCLS_MS))
            master%batch_nrefined = 0
            do ii = 1,NPTCLS_MS
                n = 0
                do ri = 1,nrefs
                    master%batch_rows(ri,ii) = rows(ri,pinds(ii))
                    if( .not. refined(ri,pinds(ii)) ) cycle
                    master%batch_rows(ri,ii)%dist = rows(ri,pinds(ii))%dist - 50.
                    n = n + 1
                    master%batch_refined(n,ii) = prob_refined_entry(ii, ri, shifts(1,ri,pinds(ii)), shifts(2,ri,pinds(ii)), &
                        &master%batch_rows(ri,ii)%dist, -1.)
                enddo
                master%batch_nrefined(ii)  = n
                master%seed_shifts(:,ii)   = seeds(:,pinds(ii))
                master%seed_has_sh(ii)     = .true.
            enddo
            call master%open_assignment
            call master%assign_batch(1, NPTCLS_MS)
            call master%close_assignment
        end subroutine assign_master
    end subroutine test_insert_extras_multistate

    !> The insertion cone: a rotation by alpha and a shift ds between two poses move the rim of a particle of
    !! diameter diam by at most diam alpha / 2 + ds, which the cone bounds by the resolution (10 A at 100 A: the
    !! half-angle at equal shifts is 0.2 rad, 11.46 degrees; the shift bound at equal rotations 10 pixels at 1 A);
    !! no resolution gives no cone; a resolution past the rim's half-circumference the whole sphere. On a
    !! 1000-direction C1 spiral with a 36-step rotation grid, from an equatorial direction: a direction whose
    !! rotation at the same in-plane angle and shift is within the half-angle is within, a far one is not, the same
    !! direction a half-turn away in-plane is not (one step of 10 degrees, 8.7 A at the rim, is), a shift past the
    !! bound is not, one within it is, and the rotation step plus a 2-pixel shift (10.7 A) is not while plus 1 pixel
    !! (9.7 A) it is. The rotation is the geodesic angle, not the direction and psi angles apart: at the pole, the
    !! directions [1,1,0] and [181,1,0] at the same in-plane index are 2 degrees apart by their normals and the same
    !! in psi, yet half a turn apart as rotations, so not within; with the second's in-plane angle a half-turn away
    !! they are 2 degrees apart as rotations, so within. Modulo C4, equatorial azimuths 1 and 89 degrees are 2
    !! degrees apart and pass a 3-degree cone that C1 rejects
    subroutine test_insert_cone()
        type(prob_insert_cone) :: cone
        type(oris)             :: eul, eul2
        type(sym)              :: se1, se4
        real    :: e1(3), e2(3), geo
        integer :: i, ic, inear, ifar
        call cone%new([10.], 100., 1.)
        call assert_real(11.459156, cone%get_angle(1), 1.e-3, 'insertion cone: 10 A at a 100 A particle is 11.46 degrees')
        call assert_true(cone%is_open(), 'insertion cone: open')
        call cone%new([0.], 100., 1.)
        call assert_true(.not. cone%is_open() .and. cone%get_angle(1) < 0., 'insertion cone: no resolution, no cone')
        call cone%new([1000.], 100., 1.)
        call assert_real(180., cone%get_angle(1), 1.e-3, 'insertion cone: a resolution past the rim is the whole sphere')
        ! the C1 spiral, from the direction nearest the equator
        call se1%new('c1')
        call eul%new(1000, is_ptcl=.false.)
        call se1%build_refspiral(eul)
        call cone%new([10.], 100., 1.)
        call cone%set_space(eul, se1, 36)
        call assert_true(cone%has_space(), 'insertion cone: the space is set')
        ic = 1
        e1 = eul%get_euler(ic)
        do i = 2,1000
            e2 = eul%get_euler(i)
            if( abs(e2(2) - 90.) < abs(e1(2) - 90.) )then
                ic = i
                e1 = e2
            endif
        enddo
        e1(3) = 360.
        inear = 0
        ifar  = 0
        do i = 1,1000
            if( i == ic ) cycle
            e2    = eul%get_euler(i)
            e2(3) = 360.
            geo   = geodesic_deg(e1, e2)
            if( inear == 0 .and. geo > 3.  .and. geo < 11. ) inear = i
            if( ifar  == 0 .and. geo > 20. ) ifar = i
            if( inear > 0 .and. ifar > 0 ) exit
        enddo
        call assert_true(inear > 0 .and. ifar > 0, 'insertion cone: a near and a far direction found on the spiral')
        call assert_true(cone%within(1, ic, 1, [0.,0.], inear, 1, [0.,0.]), 'insertion cone: a rotation within the half-angle is within')
        call assert_true(.not. cone%within(1, ic, 1, [0.,0.], ifar, 1, [0.,0.]), 'insertion cone: a far direction is not')
        call assert_true(.not. cone%within(1, ic, 1, [0.,0.], ic, 19, [0.,0.]), 'insertion cone: the same direction a half-turn away in-plane is another pose')
        call assert_true(cone%within(1, ic, 1, [0.,0.], ic, 2, [0.,0.]), 'insertion cone: one rotation step (10 degrees, 8.7 A at the rim) is within')
        call assert_true(.not. cone%within(1, ic, 1, [0.,0.], ic, 1, [11.,0.]), 'insertion cone: a shift past the bound is another pose')
        call assert_true(cone%within(1, ic, 1, [0.,0.], ic, 1, [6.,6.]), 'insertion cone: a shift within the bound is within')
        call assert_true(.not. cone%within(1, ic, 1, [0.,0.], ic, 2, [2.,0.]), 'insertion cone: the rotation step plus a 2-pixel shift (10.7 A) is another pose')
        call assert_true(cone%within(1, ic, 1, [0.,0.], ic, 2, [1.,0.]), 'insertion cone: the rotation step plus a 1-pixel shift (9.7 A) is within')
        call assert_true(.not. cone%within(2, ic, 1, [0.,0.], inear, 1, [0.,0.]), 'insertion cone: a state without a cone admits nothing')
        ! the polar pair: 2 degrees apart by their normals, half a turn apart as rotations
        call eul2%new(2, is_ptcl=.false.)
        call eul2%set_euler(1, [1.,   1., 0.])
        call eul2%set_euler(2, [181., 1., 0.])
        call cone%new([10.], 100., 1.)
        call cone%set_space(eul2, se1, 36)
        call assert_true(.not. cone%within(1, 1, 1, [0.,0.], 2, 1, [0.,0.]), 'insertion cone: the polar pair at the same in-plane index is half a turn apart')
        call assert_true(cone%within(1, 1, 1, [0.,0.], 2, 19, [0.,0.]), 'insertion cone: the polar pair with the in-plane angle a half-turn away is 2 degrees apart')
        ! modulo C4: azimuths 1 and 89 degrees on the equator
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

    contains

        !> the angle between two poses as rotations, degrees
        real function geodesic_deg( ea, eb )
            real, intent(in) :: ea(3), eb(3)
            real :: m(3,3)
            m = matmul(euler2m(ea), transpose(euler2m(eb)))
            geodesic_deg = rad2deg(acos(max(-1., min(1., 0.5 * (m(1,1) + m(2,2) + m(3,3) - 1.)))))
        end function geodesic_deg

    end subroutine test_insert_cone

    !> The extra poses of the weighted insertion (top_posterior_extras): kmax = 1 keeps the committed entry alone at
    !! weight 1; no random number is drawn at any kmax, so the stream is the single draw's; the extras are the
    !! eligible entries of highest posterior at the given temperature, in descending order, never the committed
    !! entry nor an invalid one, with the weights the posterior over the kept set summing to 1 (the committed
    !! entry's its own share, zero when it has no mass); an extra below wmin is dropped with every lighter one; on a
    !! tie the earlier entry stays; nelig counts the eligible valid entries other than the committed one
    !> the review's T2: the insertion weighs poses by the data term with the map's prior on the pose's increment, in
    !! the dense and the neighbourhood tables alike (insert_pose_loss); two poses of equal data term 20 and 22 px from
    !! the stored shift under sigma = 10 px share the kept set 0.603483 / 0.396517 at any temperature
    subroutine test_insert_pose_loss()
        real     :: d(2), w(2), temps(2)
        real(dp) :: prior_w
        integer  :: inds(2), nkept, nelig, it
        prior_w = map_prior_weight(10., 1.)
        call assert_real(0.005, real(prior_w), 1.e-7, 'the map''s prior''s weight is 1 / (2 sigma^2)')
        temps = [0.5, 2.]
        do it = 1, 2
            d(1) = insert_pose_loss(0.3, [20.d0, 0.d0], prior_w, temps(it))
            d(2) = insert_pose_loss(0.3, [0.d0, -22.d0], prior_w, temps(it))
            call top_posterior_extras(2, d, temps(it), [.true., .true.], 1, 2, 1.e-6, nkept, inds, w, nelig)
            call assert_int(2, nkept, 'both poses kept')
            call assert_real(0.603483, w(1), 1.e-5, 'the pose moved 20 px weighs 0.603483')
            call assert_real(0.396517, w(2), 1.e-5, 'the pose moved 22 px weighs 0.396517')
        enddo
        call assert_real(0.3, insert_pose_loss(0.3, [20.d0, 0.d0], prior_w, 0.), 0., &
            &'without a temperature the loss is the data term')
    end subroutine test_insert_pose_loss

    subroutine test_top_posterior_extras()
        real    :: dists(6), w(5), z
        integer :: inds(5), nkept, nelig
        logical :: mask(6)
        ! T = 1 on [0, 1, 2, 3, huge, 0.5], all eligible; the committed entry the mode
        dists = [0., 1., 2., 3., huge(1.), 0.5]
        mask  = .true.
        call top_posterior_extras(6, dists, 1., mask, 1, 1, 0.01, nkept, inds, w, nelig)
        call assert_int(1, nkept, 'kmax 1 keeps the committed entry alone')
        call assert_int(1, inds(1), 'kmax 1 keeps the committed entry')
        call assert_real(1., w(1), 0., 'kmax 1 gives the committed entry weight 1')
        call assert_int(4, nelig, 'kmax 1 still counts the eligible valid entries other than the committed one')
        ! kmax 3: the two eligible entries of highest weight, 6 (d 0.5) then 2 (d 1), the committed first, the weights
        ! their shares of the kept set, the invalid entry never selected
        call top_posterior_extras(6, dists, 1., mask, 1, 3, 1.e-6, nkept, inds, w, nelig)
        z = 1. + exp(-0.5) + exp(-1.)
        call assert_int(3, nkept, 'kmax 3 keeps the committed entry and two extras')
        call assert_true(inds(1) == 1 .and. inds(2) == 6 .and. inds(3) == 2, 'the extras are the entries of highest weight in descending order')
        call assert_real(1. / z,        w(1), 1.e-5, 'the committed entry''s weight is its share of the kept set')
        call assert_real(exp(-0.5) / z, w(2), 1.e-5, 'the weight of entry 6 is its share of the kept set')
        call assert_real(exp(-1.0) / z, w(3), 1.e-5, 'the weight of entry 2 is its share of the kept set')
        call assert_int(4, nelig, 'four eligible valid entries other than the committed one')
        ! a cap above the eligible count takes them all
        call top_posterior_extras(3, [0., 1., 2.], 1., [.true., .true., .true.], 1, 5, 1.e-6, nkept, inds, w, nelig)
        call assert_true(nkept == 3 .and. nelig == 2 .and. inds(2) == 2 .and. inds(3) == 3, 'a cap above the eligible count takes every eligible entry')
        ! the order is by weight, not by index
        call top_posterior_extras(5, [0., 3., 1., 2., 0.5], 1., [.true., .true., .true., .true., .true.], 1, 4, 1.e-6, nkept, inds, w, nelig)
        call assert_true(nkept == 4 .and. inds(2) == 5 .and. inds(3) == 3 .and. inds(4) == 4, 'the extras come in descending weight whatever their row order')
        ! the review's M4: an excluded row far better than the kept set changes nothing (the normaliser is local)
        call top_posterior_extras(3, [0., 0.9, 0.901], 100., [.false., .true., .true.], 2, 2, 1.e-6, nkept, inds, w, nelig)
        call assert_true(nkept == 2 .and. inds(1) == 2 .and. inds(2) == 3, 'an excluded better row keeps the eligible extra')
        call assert_real(0.524980, w(1), 1.e-5, 'the committed entry''s weight is its conditional share')
        call assert_real(0.475020, w(2), 1.e-5, 'the extra''s weight is its conditional share')
        call top_posterior_extras(3, [5., 0.9, 0.901], 100., [.false., .true., .true.], 2, 2, 1.e-6, nkept, inds, w, nelig)
        call assert_real(0.524980, w(1), 1.e-5, 'changing an excluded row leaves the kept set''s weights')
        ! a tie keeps the earlier entries
        call top_posterior_extras(4, [0., 1., 1., 1.], 1., [.true., .true., .true., .true.], 1, 3, 1.e-6, nkept, inds, w, nelig)
        call assert_true(nkept == 3 .and. inds(2) == 2 .and. inds(3) == 3, 'equal weights keep the earlier entries')
        ! eligibility: only entry 6 eligible, so one extra at most, and it is entry 6
        mask = .false.
        mask(6) = .true.
        call top_posterior_extras(6, dists, 1., mask, 1, 3, 1.e-6, nkept, inds, w, nelig)
        call assert_int(2, nkept, 'with one eligible entry there is one extra')
        call assert_int(6, inds(2), 'and it is the eligible entry')
        call assert_int(1, nelig, 'and the occupancy is one')
        ! the committed entry itself ineligible changes nothing about its being kept first
        mask = .true.
        mask(1) = .false.
        call top_posterior_extras(6, dists, 1., mask, 1, 2, 1.e-6, nkept, inds, w, nelig)
        call assert_true(nkept == 2 .and. inds(1) == 1 .and. inds(2) == 6, 'the committed entry is kept first whatever its eligibility')
        ! the mass rule: a near-delta row at a cold temperature keeps the committed entry alone
        mask = .true.
        call top_posterior_extras(6, dists, 10., mask, 1, 3, 0.01, nkept, inds, w, nelig)
        call assert_true(nkept == 1 .and. w(1) == 1., 'at T = 10 every extra falls below 1% and is dropped: the committed entry alone')
        ! a committed entry 100 nats above the minimum: its weight underflows to zero, it is still kept first, and the
        ! particle's mass is on its extras, the weights summing to 1
        call top_posterior_extras(3, [0., 100., 2.], 1., [.true., .true., .true.], 2, 3, 0.01, nkept, inds, w, nelig)
        call assert_true(nkept == 3 .and. inds(1) == 2 .and. inds(2) == 1 .and. inds(3) == 3, 'a committed entry of no mass is kept first, its extras after it')
        call assert_true(w(1) == 0. .and. abs(sum(w(1:3)) - 1.) < 1.e-5, 'the committed entry of no mass has weight zero and the weights sum to 1')
        ! no valid distance: the committed entry alone
        call top_posterior_extras(3, [huge(1.), huge(1.), huge(1.)], 1., [.true., .true., .true.], 1, 3, 0.01, nkept, inds, w, nelig)
        call assert_true(nkept == 1 .and. nelig == 0, 'a row without a valid distance keeps the committed entry alone')
    end subroutine test_top_posterior_extras

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

    !> synthetic dense rows of nptcls particles over nspace directions and nstates states: random rows, a refined
    !! shift on about one entry in ten, a seed per particle; the point group and the direction space as the
    !! assignment needs them (the weighted insertion measures its candidates' poses in the direction space)
    subroutine make_dense_rows( params, build, nspace, nstates, nptcls, rows, shifts, refined, seeds )
        class(parameters), allocatable, target, intent(inout) :: params
        class(builder),    allocatable, target, intent(inout) :: build
        integer,                                intent(in)    :: nspace, nstates, nptcls
        type(prob_dense_entry), allocatable,    intent(out)   :: rows(:,:)
        real,                   allocatable,    intent(out)   :: shifts(:,:,:), seeds(:,:)
        logical,                allocatable,    intent(out)   :: refined(:,:)
        integer :: nrefs, i, ri
        nrefs = nspace * nstates
        allocate(params, build)
        params%nspace      = nspace
        params%nstates     = nstates
        params%prob_athres = 10.
        params%l_doshift   = .true.
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
    end subroutine make_dense_rows

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
