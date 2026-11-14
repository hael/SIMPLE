!@descr: unit tests for the 2D probability tables: the pipeline's dense and sparse selection (simple_eul_prob_tab2D)
! A batch of candidates for three particles and two classes goes through the pipeline's table (the thread buffers
! the fill appends to, the batch's rows or candidate store, the per-particle draw from the tempered posterior at a
! sharp temperature); the assignment must give each particle its best class, distance, in-plane index, peak count
! and fraction: dense (refine=prob) and sparse (refine=prob_neigh). The draw on synthetic tables with ties, flat
! and invalid rows at several temperatures and seeds. The posterior rows a student's teacher reads from a kept dense
! batch, and the assignment from those rows.
module simple_eul_prob_tab2D_tester
use, intrinsic :: iso_fortran_env, only: int64
use simple_core_module_api
use simple_builder,            only: builder
use simple_parameters,         only: parameters
use simple_eul_prob_tab2D,     only: eul_prob_tab2D
use simple_eul_prob_tab_utils, only: prob_candidate, prob_row_record, inpl_quadrature, prob_insert_inpl, &
    &temperature_or_argmin
use simple_type_defs,          only: ptcl_ref
use simple_rnd,                only: seed_rnd_fixed, rnd_scope_key, rnd_draw_key, rnd_keyed, rnd_key, &
    &RND_PROG_REFINE2D, RND_PASS_ALIGN, RND_FAM_INPLQ, RND_FAM_ASSIGN
use simple_test_utils
implicit none
private
public :: run_all_eul_prob_tab2D_tests
#include "simple_local_flags.inc"

type(parameters), target :: params
type(builder),    target :: build
type(eul_prob_tab2D)     :: table
integer, parameter :: NPTCLS = 3, NCLASSES = 2
integer :: pinds(NPTCLS) = [10,20,30]   ! project indices of the three particles

contains

    subroutine run_all_eul_prob_tab2D_tests()
        write(*,'(A)') '**** running all 2D probability table tests ****'
        params%ncls         = NCLASSES
        params%npeaks_inpl  = NCLASSES
        params%l_doshift    = .false.
        params%prob_beta    = 1.
        call test_dense_batch()
        call test_dense_posterior_rows()
        call test_sparse_batch()
        call test_sparse_uneven_batches()
        call test_sparse_shared_best_class()
        call test_dense_posterior_assignment()
        call test_sparse_posterior_assignment()
        call test_inpl_quadrature()
        call test_inpl_phase()
        call test_inpl_insert_channel()
        call table%kill
        call fix_test_seed()
    end subroutine run_all_eul_prob_tab2D_tests

    !> the dense posterior draw on a synthetic table of 150 particles and 20 classes with ties (a flat particle,
    !! identical particles, a class tie at the minimum, invalid distances): at a huge temperature every particle with a unique finite minimum takes its
    !! argmin (its distance and in-plane index), a class tie one of the tied classes, the flat particle any class
    !! and the invalid one class 1; a particle without a temperature takes its argmin too; at T = 2 on distances
    !! in [0, 4] some particles leave the mode; the same seed gives the same assignment and another seed another
    subroutine test_dense_posterior_assignment()
        integer, parameter :: NP = 150, NC = 20
        type(eul_prob_tab2D) :: tab
        real,    allocatable :: dists(:,:)
        integer, allocatable :: pinds_t(:), icls_a(:), icls_b(:)
        integer :: i, icls, nwrong, loc(1), ncls_prev, npeaks_prev
        logical :: l_unique(NP)
        character(len=len(params%refine)) :: refine_prev
        write(*,'(A)') 'test_dense_posterior_assignment'
        ncls_prev   = params%ncls
        npeaks_prev = params%npeaks_inpl
        refine_prev = params%refine
        params%ncls        = NC
        params%npeaks_inpl = NC
        params%refine      = 'prob'
        pinds_t = [(i, i=1,NP)]
        allocate(dists(NC,NP))
        call seed_rnd_fixed(20261023)
        do i = 1, NP
            do icls = 1, NC
                dists(icls,i) = 4. * ran3()
            enddo
        enddo
        dists(:,1)  = 1.5                  ! flat
        dists(:,8)  = dists(:,7)           ! identical particles
        dists(:,9)  = max(dists(:,9), 0.5) ! the tie below is the only minimum, whatever the random fill drew
        dists(3,9)  = 0.01                 ! a class tie at the minimum
        dists(11,9) = 0.01
        dists(5,12) = huge(1.)             ! invalid distances
        dists(:,13) = huge(1.)
        l_unique = .true.
        l_unique([1,9,13]) = .false.
        ! a huge temperature: the argmin
        call set_fixed_seed(5001, propagate=.true.)
        call assign_dense(tab, 1.e9)
        nwrong = 0
        do i = 1, NP
            if( l_unique(i) )then
                loc = minloc(dists(:,i))
                if( tab%assgn_map(i)%icls /= loc(1) ) nwrong = nwrong + 1
                if( tab%assgn_map(i)%dist /= dists(loc(1),i) ) nwrong = nwrong + 1
                if( tab%assgn_map(i)%inpl /= loc(1) ) nwrong = nwrong + 1
            endif
        enddo
        call assert_int(0, nwrong, 'a huge temperature gives every particle with a unique minimum its argmin class')
        call assert_true(tab%assgn_map(9)%icls == 3 .or. tab%assgn_map(9)%icls == 11, 'a class tie draws one of the tied classes')
        call assert_true(tab%assgn_map(1)%icls >= 1 .and. tab%assgn_map(1)%icls <= NC, 'the flat particle takes a class')
        call assert_int(1, tab%assgn_map(13)%icls, 'the particle without a valid distance takes class 1')
        call assert_true(all(tab%assgn_map(:)%frac == 100.) .and. all(tab%assgn_map(:)%npeaks == 0), 'frac 100 and no peak count')
        icls_a = tab%assgn_map(:)%icls
        ! no temperature from the fill: the argmin, as at the huge temperature
        call set_fixed_seed(5001, propagate=.true.)
        call assign_dense(tab, 0.)
        call assert_int(0, count(pack(tab%assgn_map(:)%icls, l_unique) /= pack(icls_a, l_unique)),&
            &'a particle without a temperature takes its argmin')
        ! T = 2: some particles leave the mode, at the same seed the same ones
        call set_fixed_seed(5002, propagate=.true.)
        call assign_dense(tab, 2.)
        icls_a = tab%assgn_map(:)%icls
        nwrong = 0
        do i = 1, NP
            if( .not. l_unique(i) ) cycle
            loc = minloc(dists(:,i))
            if( icls_a(i) /= loc(1) ) nwrong = nwrong + 1
            if( tab%assgn_map(i)%dist /= dists(icls_a(i),i) ) nwrong = nwrong + NP
        enddo
        call assert_true(nwrong > 0 .and. nwrong < NP, 'at T = 2 some particles leave the mode, with the drawn class''s distance')
        call set_fixed_seed(5002, propagate=.true.)
        call assign_dense(tab, 2.)
        call assert_int(0, count(tab%assgn_map(:)%icls /= icls_a), 'the same seed gives the same dense posterior assignment')
        call set_fixed_seed(5003, propagate=.true.)
        call assign_dense(tab, 2.)
        icls_b = tab%assgn_map(:)%icls
        call assert_true(any(icls_b /= icls_a), 'another seed gives another dense posterior assignment')
        call tab%kill
        params%ncls        = ncls_prev
        params%npeaks_inpl = npeaks_prev
        params%refine      = refine_prev

    contains

        !> the pipeline's table over the particles, every row's candidates (one per class, the in-plane index the
        !! class) appended to the first thread buffer as the fill appends them, the temperature set, the batch of
        !! every row assigned
        subroutine assign_dense( t, temp )
            type(eul_prob_tab2D), intent(inout) :: t
            real,                 intent(in)    :: temp
            type(prob_candidate) :: candidate
            integer :: ip, ic
            call t%new_pipeline(params, build, pinds_t)
            do ip = 1, NP
                do ic = 1, NC
                    candidate      = prob_candidate()
                    candidate%iref = ic
                    candidate%inpl = ic
                    candidate%dist = dists(ic,ip)
                    call t%candidate_buffers(1)%append(ip, candidate)
                enddo
            enddo
            t%temps = temp
            call t%open_assignment
            call t%assign_batch(1, NP)
            call t%close_assignment
        end subroutine assign_dense

    end subroutine test_dense_posterior_assignment

    !> the sparse posterior draw on a synthetic candidate batch of 300 particles and 20 classes (1 to 20 candidates
    !! per particle, ties, candidates without an in-plane index): at a huge temperature every particle takes a
    !! candidate at its smallest valid distance, with the peak count of the candidates with an in-plane index and
    !! their fraction of the classes; at T = 2 some particles leave the mode; the same seed gives the same
    !! assignment and another seed another
    subroutine test_sparse_posterior_assignment()
        integer, parameter :: NP = 300, NC = 20
        class(eul_prob_tab2D), allocatable :: tab
        type(prob_candidate),  allocatable :: cands(:)
        integer, allocatable :: counts(:), pinds_t(:), icls_a(:), icls_b(:)
        integer :: i, k, n, nwrong, npeaks, ncls_prev, npeaks_prev
        real    :: best
        character(len=len(params%refine)) :: refine_prev
        write(*,'(A)') 'test_sparse_posterior_assignment'
        ncls_prev   = params%ncls
        npeaks_prev = params%npeaks_inpl
        refine_prev = params%refine
        params%ncls        = NC
        params%npeaks_inpl = NC
        params%refine      = 'prob_neigh'
        pinds_t = [(i, i=1,NP)]
        call seed_rnd_fixed(20261023)
        allocate(counts(NP))
        do i = 1, NP
            counts(i) = 1 + int(ran3() * real(NC-1))
        enddo
        counts(7:8) = NC
        allocate(cands(sum(counts)))
        n = 0
        do i = 1, NP
            do k = 1, counts(i)
                n = n + 1
                cands(n)%iref = 1 + mod(k * 7 + i, NC)
                cands(n)%inpl = 1 + int(ran3() * 35.)
                cands(n)%dist = 4. * ran3()
                if( k > 1 .and. ran3() < 0.05 ) cands(n)%inpl = 0  ! not a candidate
            enddo
        enddo
        n = sum(counts(1:6))
        cands(n+1:n+NC)%dist = 0.7                                 ! a flat particle
        cands(n+1:n+NC)%inpl = 3
        allocate(tab)
        ! a huge temperature: a candidate at the smallest valid distance
        call set_fixed_seed(6001, propagate=.true.)
        call assign_sparse(tab, 1.e9)
        nwrong = 0
        n = 0
        do i = 1, NP
            best   = huge(1.)
            npeaks = 0
            do k = 1, counts(i)
                if( cands(n+k)%inpl > 0 )then
                    npeaks = npeaks + 1
                    best   = min(best, cands(n+k)%dist)
                endif
            enddo
            if( tab%assgn_map(i)%dist /= best ) nwrong = nwrong + 1
            if( tab%assgn_map(i)%npeaks /= npeaks ) nwrong = nwrong + 1
            if( abs(tab%assgn_map(i)%frac - 100. * real(npeaks) / real(NC)) > 1.e-4 ) nwrong = nwrong + 1
            if( tab%assgn_map(i)%inpl < 1 ) nwrong = nwrong + 1
            n = n + counts(i)
        enddo
        call assert_int(0, nwrong, 'a huge temperature gives every particle a candidate at its smallest valid distance')
        call tab%kill
        ! T = 2: some particles leave the mode, at the same seed the same ones
        call set_fixed_seed(6002, propagate=.true.)
        call assign_sparse(tab, 2.)
        icls_a = tab%assgn_map(:)%icls
        nwrong = 0
        n = 0
        do i = 1, NP
            best = huge(1.)
            do k = 1, counts(i)
                if( cands(n+k)%inpl > 0 ) best = min(best, cands(n+k)%dist)
            enddo
            if( tab%assgn_map(i)%dist > best ) nwrong = nwrong + 1
            n = n + counts(i)
        enddo
        call assert_true(nwrong > 0 .and. nwrong < NP, 'at T = 2 some particles leave the mode')
        call tab%kill
        call set_fixed_seed(6002, propagate=.true.)
        call assign_sparse(tab, 2.)
        call assert_int(0, count(tab%assgn_map(:)%icls /= icls_a), 'the same seed gives the same sparse posterior assignment')
        call tab%kill
        call set_fixed_seed(6003, propagate=.true.)
        call assign_sparse(tab, 2.)
        icls_b = tab%assgn_map(:)%icls
        call assert_true(any(icls_b /= icls_a), 'another seed gives another sparse posterior assignment')
        call tab%kill
        deallocate(tab)
        params%ncls        = ncls_prev
        params%npeaks_inpl = npeaks_prev
        params%refine      = refine_prev

    contains

        !> the pipeline's table over the particles, every particle's candidates appended to the first thread
        !! buffer in their order (the order the fill gives the batch's store), the temperature set, the batch of
        !! every row assigned
        subroutine assign_sparse( t, temp )
            class(eul_prob_tab2D), intent(inout) :: t
            real,                  intent(in)    :: temp
            integer :: ip, k, pos
            call t%new_pipeline(params, build, pinds_t)
            pos = 0
            do ip = 1, NP
                do k = 1, counts(ip)
                    pos = pos + 1
                    call t%candidate_buffers(1)%append(ip, cands(pos))
                enddo
            enddo
            t%temps = temp
            call t%open_assignment
            call t%assign_batch(1, NP)
            call t%close_assignment
        end subroutine assign_sparse

    end subroutine test_sparse_posterior_assignment

    !> refine=prob: every particle carries a candidate per class; at the sharp temperature the draw is the smaller
    !! distance, and a refined candidate keeps its shift through the batch's rows
    subroutine test_dense_batch()
        type(prob_candidate) :: candidates(6)
        integer :: particle_indices(6)
        write(*,'(A)') 'test_dense_batch'
        params%refine = 'prob'
        call set_candidate(candidates(1),1,0.1)
        call set_candidate(candidates(2),2,0.9)
        call set_candidate(candidates(3),1,0.8)
        call set_candidate(candidates(4),2,0.2)
        call set_candidate(candidates(5),1,0.3)
        call set_candidate(candidates(6),2,0.4)
        candidates(4)%x      =  1.5
        candidates(4)%y      = -2.5
        candidates(4)%has_sh = .true.
        particle_indices = [1,1,2,2,3,3]
        call assign_batch_of(particle_indices, candidates)
        call assert_assignment('dense batch',[1,2,1],[0.1,0.2,0.3],[0,0,0],[100.,100.,100.])
        call assert_shifts('dense batch',[0.,1.5,0.],[0.,-2.5,0.],[.false.,.true.,.false.])
        call table%kill
    end subroutine test_dense_batch

    !> refine=prob_neigh: the batch carries only the evaluated candidates, one or two per particle
    subroutine test_sparse_batch()
        type(prob_candidate) :: candidates(4)
        integer :: particle_indices(4)
        write(*,'(A)') 'test_sparse_batch'
        params%refine = 'prob_neigh'
        call set_candidate(candidates(1),1,0.1)
        call set_candidate(candidates(2),2,0.2)
        call set_candidate(candidates(3),1,0.3)
        call set_candidate(candidates(4),2,0.4)
        particle_indices = [1,2,3,3]
        call assign_batch_of(particle_indices, candidates)
        call assert_assignment('sparse batch',[1,2,1],[0.1,0.2,0.3],[1,1,2],[50.,50.,100.])
        call table%kill
    end subroutine test_sparse_batch

    !> Uneven batches preserve global IDs, seed/refined shifts and untouched rows outside the current batch.
    subroutine test_sparse_uneven_batches()
        integer, parameter :: NP = 11, BS = 4
        class(eul_prob_tab2D), allocatable :: tab
        type(prob_candidate) :: candidate
        integer :: ids(NP), i, first, last, best_class
        real    :: expected_x, expected_y
        write(*,'(A)') 'test_sparse_uneven_batches'
        allocate(tab)
        params%refine    = 'prob_neigh'
        params%l_doshift = .true.
        ids = [(100 + 7*i, i=1,NP)]
        call tab%new_pipeline(params, build, ids)
        tab%temps = 1.e9
        do i = 1,NP
            tab%seed_shifts(:,i) = [real(i), -real(i)]
            tab%seed_has_sh(i) = .true.
        enddo
        call tab%open_assignment
        do first = 1,NP,BS
            last = min(NP, first + BS - 1)
            do i = first,last
                best_class = 1 + mod(i, NCLASSES)
                candidate = prob_candidate()
                candidate%iref = best_class
                candidate%inpl = best_class
                candidate%dist = 0.1
                if( mod(i,2) == 0 )then
                    candidate%has_sh = .true.
                    candidate%x = 0.5 * real(i)
                    candidate%y = 0.25 * real(i)
                endif
                call tab%candidate_buffers(1)%append(i, candidate)
                candidate = prob_candidate()
                candidate%iref = 3 - best_class
                candidate%inpl = 3 - best_class
                candidate%dist = 0.9
                call tab%candidate_buffers(1)%append(i, candidate)
            enddo
            call tab%assign_batch(first, last)
            if( last < NP ) call assert_true(all(tab%assgn_map(last+1:)%icls == 0), 'sparse batch leaves later rows untouched')
        enddo
        call tab%close_assignment
        do i = 1,NP
            call assert_int(ids(i), tab%assgn_map(i)%pind, 'uneven sparse batches preserve global particle IDs')
            call assert_int(1 + mod(i,NCLASSES), tab%assgn_map(i)%icls, 'uneven sparse batches draw the unique minimum')
            call assert_int(2, tab%assgn_map(i)%npeaks, 'uneven sparse batches preserve candidate counts')
            expected_x = real(i)
            expected_y = -real(i)
            if( mod(i,2) == 0 )then
                expected_x = 0.5 * real(i)
                expected_y = 0.25 * real(i)
            endif
            call assert_true(tab%assgn_map(i)%has_sh, 'uneven sparse batches retain seeded and refined shifts')
            call assert_true(tab%assgn_map(i)%x == expected_x .and. tab%assgn_map(i)%y == expected_y, &
                &'uneven sparse batches map shifts to the right particle')
        enddo
        call tab%kill
        deallocate(tab)
        params%l_doshift = .false.
    end subroutine test_sparse_uneven_batches

    !> two particles share their best class: at the sharp temperature both draw it, no class forced to be populated
    subroutine test_sparse_shared_best_class()
        type(prob_candidate) :: candidates(6)
        integer :: particle_indices(6)
        write(*,'(A)') 'test_sparse_shared_best_class'
        params%refine = 'prob_neigh'
        call set_candidate(candidates(1),1,0.5)
        call set_candidate(candidates(2),2,0.1)
        call set_candidate(candidates(3),1,0.4)
        call set_candidate(candidates(4),2,0.2)
        call set_candidate(candidates(5),1,0.3)
        call set_candidate(candidates(6),2,0.6)
        particle_indices = [1,1,2,2,3,3]
        call assign_batch_of(particle_indices, candidates)
        call assert_assignment('sparse shared best class',[2,2,1],[0.1,0.2,0.3],[2,2,2],[100.,100.,100.])
        call table%kill
    end subroutine test_sparse_shared_best_class

    !> The in-plane quadrature of one row's posterior (the 2D insertion's selector) against closed forms, on rows of
    !! n = 360 rotations at temperature 1: a flat row gives K points at the cumulative ranks (j - 1 + v) / K, so with
    !! K = 8 and v = 0.3 the rotations 14 + 45 (j - 1), weights 1/8, and an angular response that vanishes on every
    !! harmonic but the multiples of 8 (the rule averages over a cyclic group of order K, not the circle: the harmonics
    !! 8 and 16 survive at full amplitude); a row with all its mass in one rotation gives that rotation at weight 1, at
    !! a huge temperature too; two equal rotations split the points; a row across the cut of the index (rotations 359,
    !! 360, 1, 2) keeps its circular mean; masses 4/7, 2/7, 1/7 in three rotations are recovered exactly by K = 7 at v =
    !! 0.5 (the points fall at the cumulative ranks 1/14 .. 13/14) with every harmonic equal to the full row's; a row
    !! without a valid distance gives no entry and one valid distance that one. A Gaussian row of 10 rotations' width
    !! (n_eff 41.3): at K = 64 the response on the harmonics up to 8 is within 0.02 of the full row's at two phases
    !! (an emulation of the rule in Python, with the same row, gave a largest difference of 0.0134); K = 32 on a flat
    !! row leaves the multiples of 32 at 0.988 and every other harmonic up to 20 below 0.05 (emulation: 0.035). The
    !! points are snapped to the rotation grid, so the cancellation is exact only where K divides n (K = 36 at the
    !! phase 0.37: rotations 4 + 10 (j - 1), every harmonic up to 50 but the multiples of 36 below 1e-5) and
    !! otherwise leaves a residual that grows with the harmonic: K = 64 keeps 12.5 % of harmonic 40 at any phase
    !! (emulation: 0.1250 at the phases 0.3 and 0.7), and no harmonic up to 50 above that
    subroutine test_inpl_quadrature()
        integer,  parameter :: N = 360, KMAX = 64
        real,     parameter :: TOL_HARMONIC = 1.e-5, TOL_GAUSS = 0.02, TOL_QUANT = 0.05
        real     :: dists(N), w(KMAX), neff, full(N)
        integer  :: inds(KMAX), nent, m, j, k
        real(dp) :: hq, hf
        logical  :: l_ok
        write(*,'(A)') 'test_inpl_quadrature'
        ! a flat row, K = 8 at v = 0.3
        dists = 1.
        call inpl_quadrature(N, dists, 1., 8, 0.3_dp, nent, inds, w, neff)
        call assert_int(8, nent, 'a flat row gives K entries')
        l_ok = .true.
        do j = 1,nent
            if( inds(j) /= 14 + 45*(j-1) ) l_ok = .false.
            if( abs(w(j) - 0.125) > 1.e-7 ) l_ok = .false.
        enddo
        call assert_true(l_ok, 'a flat row at K = 8, v = 0.3 gives the rotations 14 + 45 (j - 1) at weight 1/8')
        call assert_real(1., sum(w(1:nent)), 1.e-6, 'the weights of the entries sum to one')
        call assert_real(360., neff, 1.e-2, 'the effective number of rotations of a flat row is its length')
        l_ok = .true.
        do m = 1,20
            hq = abs(harmonic(m, N, nent, inds, w))
            if( mod(m,8) /= 0 )then
                if( hq > TOL_HARMONIC ) l_ok = .false.
            else
                if( abs(hq - 1.d0) > TOL_HARMONIC ) l_ok = .false.
            endif
        enddo
        call assert_true(l_ok, 'a flat row at K = 8 annihilates every harmonic but the multiples of 8 (they survive at full amplitude)')
        ! a flat row, K = 32 at v = 0.7 (the points quantised to the rotation grid)
        call inpl_quadrature(N, dists, 1., 32, 0.7_dp, nent, inds, w, neff)
        l_ok = .true.
        do m = 1,20
            if( abs(harmonic(m, N, nent, inds, w)) > TOL_QUANT ) l_ok = .false.
        enddo
        call assert_true(l_ok, 'a flat row at K = 32 leaves no harmonic up to 20 above the quantisation of the points')
        call assert_real(0.988, real(abs(harmonic(32, N, nent, inds, w))), 5.e-3, 'a flat row at K = 32 leaves harmonic 32 at 0.988')
        ! K dividing n at a phase that is no tie: exactly equally spaced
        dists = 1.
        call inpl_quadrature(N, dists, 1., 36, 0.37_dp, nent, inds, w, neff)
        l_ok = nent == 36
        do j = 1,nent
            if( inds(j) /= 4 + 10*(j-1) ) l_ok = .false.
        enddo
        call assert_true(l_ok, 'a flat row at K = 36 (dividing 360) at a phase that is no tie gives the rotations 4 + 10 (j - 1)')
        l_ok = .true.
        do m = 1,50
            if( mod(m,36) /= 0 .and. abs(harmonic(m, N, nent, inds, w)) > TOL_HARMONIC ) l_ok = .false.
        enddo
        call assert_true(l_ok, 'a flat row at K = 36 annihilates every harmonic up to 50 but the multiples of 36')
        ! K = 64 does not divide n: the snapping of the points to the rotation grid leaves 12.5 % of harmonic 40
        call inpl_quadrature(N, dists, 1., 64, 0.7_dp, nent, inds, w, neff)
        call assert_real(0.125, real(abs(harmonic(40, N, nent, inds, w))), 2.e-3, &
            &'a flat row at K = 64 keeps 12.5 percent of harmonic 40: the points are snapped to the rotation grid')
        l_ok = .true.
        do m = 1,50
            if( mod(m,64) /= 0 .and. abs(harmonic(m, N, nent, inds, w)) > 0.13 ) l_ok = .false.
        enddo
        call assert_true(l_ok, 'a flat row at K = 64 leaves no harmonic up to 50 above 13 percent')
        ! all the mass in one rotation: at the temperature 1 and at the argmin temperature
        dists = 50.
        dists(37) = 0.
        call inpl_quadrature(N, dists, 1., 8, 0.3_dp, nent, inds, w, neff)
        call assert_true(nent == 1 .and. inds(1) == 37 .and. abs(w(1) - 1.) < 1.e-6, 'a row with all its mass in one rotation gives it at weight 1')
        call inpl_quadrature(N, dists, temperature_or_argmin(0.), 8, 0.3_dp, nent, inds, w, neff)
        call assert_true(nent == 1 .and. inds(1) == 37 .and. abs(w(1) - 1.) < 1.e-6, 'a huge temperature gives the argmin at weight 1')
        ! two equal rotations
        dists = 100.
        dists(10)  = 0.
        dists(100) = 0.
        call inpl_quadrature(N, dists, 1., 4, 0.5_dp, nent, inds, w, neff)
        call assert_true(nent == 2 .and. inds(1) == 10 .and. inds(2) == 100 .and. abs(w(1) - 0.5) < 1.e-6 .and. abs(w(2) - 0.5) < 1.e-6, &
            &'two equal rotations split the points evenly, merged into two entries')
        ! a row across the cut of the index
        dists = 100.
        dists([359,360,1,2]) = 0.
        call inpl_quadrature(N, dists, 1., 4, 0.5_dp, nent, inds, w, neff)
        call assert_true(nent == 4 .and. all(inds(1:4) == [1,2,359,360]), 'a row across the cut of the index gives its four rotations')
        full = 0.
        full([359,360,1,2]) = 0.25
        hq = abs(harmonic(1, N, nent, inds, w) - harmonic_row(1, N, full))
        call assert_real(0., real(hq), TOL_HARMONIC, 'a row across the cut keeps its circular mean (harmonic 1 equals the full row)')
        ! masses 4/7, 2/7, 1/7: recovered exactly by K = 7 at v = 0.5
        dists = 100.
        dists(100) = 0.
        dists(101) = log(2.)
        dists(102) = 2. * log(2.)
        call inpl_quadrature(N, dists, 1., 7, 0.5_dp, nent, inds, w, neff)
        call assert_true(nent == 3 .and. all(inds(1:3) == [100,101,102]), 'masses 4/7, 2/7, 1/7 give three entries')
        call assert_true(abs(w(1) - 4./7.) < 1.e-6 .and. abs(w(2) - 2./7.) < 1.e-6 .and. abs(w(3) - 1./7.) < 1.e-6, &
            &'masses 4/7, 2/7, 1/7 are recovered exactly at K = 7')
        call assert_real(2.60049, neff, 1.e-3, 'the effective number of rotations is exp of the entropy (4/7, 2/7, 1/7)')
        full = 0.
        full(100) = 4./7.
        full(101) = 2./7.
        full(102) = 1./7.
        l_ok = .true.
        do m = 1,20
            if( abs(harmonic(m, N, nent, inds, w) - harmonic_row(m, N, full)) > TOL_HARMONIC ) l_ok = .false.
        enddo
        call assert_true(l_ok, 'a row of masses 4/7, 2/7, 1/7 has every harmonic up to 20 equal to the full row''s')
        ! invalid and partly invalid rows
        dists = huge(1.)
        call inpl_quadrature(N, dists, 1., 8, 0.3_dp, nent, inds, w, neff)
        call assert_int(0, nent, 'a row without a valid distance gives no entry')
        dists(21) = 3.
        call inpl_quadrature(N, dists, 1., 8, 0.3_dp, nent, inds, w, neff)
        call assert_true(nent == 1 .and. inds(1) == 21 .and. abs(w(1) - 1.) < 1.e-6, 'a row with one valid distance gives that rotation at weight 1')
        ! a Gaussian row of width 10 rotations, K = 64 at two phases
        do k = 1,N
            j = min(abs(k - 100), N - abs(k - 100))
            dists(k) = 0.5 * (real(j) / 10.)**2
        enddo
        call posterior_of(N, dists, 1., full)
        do j = 1,2
            call inpl_quadrature(N, dists, 1., 64, merge(0.3_dp, 0.9_dp, j == 1), nent, inds, w, neff)
            l_ok = .true.
            do m = 1,8
                if( abs(harmonic(m, N, nent, inds, w) - harmonic_row(m, N, full)) > TOL_GAUSS ) l_ok = .false.
            enddo
            call assert_true(l_ok, 'a Gaussian row at K = 64 has its harmonics up to 8 within the emulated error of the full row''s')
        enddo
        call assert_real(41.327, neff, 1.e-2, 'the effective number of rotations of the Gaussian row')
        hf = sum(real(full,dp))
        call assert_real(1., real(hf), 1.e-5, 'the reference posterior of the test sums to one')
    end subroutine test_inpl_quadrature

    !> the response of an entry set on harmonic m: sum_j w_j exp(i m theta_j), theta_j = 2 pi (ind_j - 1) / n
    complex(dp) function harmonic( m, n, nent, inds, w ) result( h )
        integer, intent(in) :: m, n, nent, inds(:)
        real,    intent(in) :: w(:)
        integer :: j
        h = cmplx(0.d0,0.d0,dp)
        do j = 1,nent
            h = h + real(w(j),dp) * exp(cmplx(0.d0, real(m,dp) * 2.d0 * DPI * real(inds(j) - 1,dp) / real(n,dp), dp))
        enddo
    end function harmonic

    !> the same response of a full posterior over the n rotations
    complex(dp) function harmonic_row( m, n, p ) result( h )
        integer, intent(in) :: m, n
        real,    intent(in) :: p(n)
        integer :: r
        h = cmplx(0.d0,0.d0,dp)
        do r = 1,n
            h = h + real(p(r),dp) * exp(cmplx(0.d0, real(m,dp) * 2.d0 * DPI * real(r - 1,dp) / real(n,dp), dp))
        enddo
    end function harmonic_row

    !> p_r proportional to exp(-temp (dist_r - dist_min)), normalised
    subroutine posterior_of( n, dists, temp, p )
        integer, intent(in)  :: n
        real,    intent(in)  :: dists(n), temp
        real,    intent(out) :: p(n)
        p = exp(-temp * (dists - minval(dists)))
        p = p / sum(p)
    end subroutine posterior_of

    !> The quadrature's phase is one keyed variate per particle in [0,1): the same key gives the same phase, another
    !! particle another, and the family is its own (the assignment draw of the same particle differs)
    subroutine test_inpl_phase()
        type(rnd_key) :: scope
        real(dp)      :: v1, v1b, v2, va
        integer       :: ip
        logical       :: l_range
        write(*,'(A)') 'test_inpl_phase'
        scope = rnd_scope_key(20261113_int64, RND_PROG_REFINE2D, RND_PASS_ALIGN, 3)
        v1  = rnd_keyed(rnd_draw_key(scope, 1, RND_FAM_INPLQ))
        v1b = rnd_keyed(rnd_draw_key(scope, 1, RND_FAM_INPLQ))
        v2  = rnd_keyed(rnd_draw_key(scope, 2, RND_FAM_INPLQ))
        va  = rnd_keyed(rnd_draw_key(scope, 1, RND_FAM_ASSIGN))
        call assert_true(v1 == v1b, 'the same key gives the same phase')
        call assert_true(v1 /= v2, 'another particle gets another phase')
        call assert_true(v1 /= va, 'the phase family is not the assignment family')
        l_range = .true.
        do ip = 1,1000
            va = rnd_keyed(rnd_draw_key(scope, ip, RND_FAM_INPLQ))
            if( va < 0.d0 .or. va >= 1.d0 ) l_range = .false.
        enddo
        call assert_true(l_range, 'the phases of 1000 particles lie in [0,1)')
    end subroutine test_inpl_phase

    !> The batch's side channel: entries and weights per slot as set, a slot with no row or a collapsed row has no
    !! entries (the committed pose alone), the counters add the batch's rows, entries (a collapsed row counts one)
    !! and collapsed rows and persist over batches while begin_batch clears the slots
    subroutine test_inpl_insert_channel()
        type(prob_insert_inpl) :: ins
        integer,        pointer :: nent(:)
        real,           pointer :: e3(:,:), w(:,:)
        integer(int64)          :: nrows, npose, nsingle
        real(dp)                :: neff_sum
        write(*,'(A)') 'test_inpl_insert_channel'
        call ins%new(4, 5)
        call ins%begin_batch
        call ins%set_spread(2, 3, [10., 20., 30.], [0.5, 0.25, 0.25], 7.)
        call ins%set_single(3, 1.5)
        call ins%end_batch(4)
        call ins%get_view(nent, e3, w)
        call assert_true(nent(1) == 0 .and. nent(2) == 3 .and. nent(3) == 0 .and. nent(4) == 0, 'a spread slot has its entries, the others none')
        call assert_true(all(abs(e3(1:3,2) - [10.,20.,30.]) < 1.e-6) .and. all(abs(w(1:3,2) - [0.5,0.25,0.25]) < 1.e-7), 'the entries are the angles and weights set')
        call ins%get_counts(nrows, npose, nsingle, neff_sum)
        call assert_true(nrows == 2_int64 .and. npose == 4_int64 .and. nsingle == 1_int64, 'two rows, four entries, one collapsed')
        call assert_real(8.5, real(neff_sum), 1.e-5, 'the effective numbers of the rows add')
        call ins%begin_batch
        call ins%get_view(nent, e3, w)
        call assert_true(all(nent == 0), 'a new batch starts every slot at the committed pose')
        call ins%set_spread(1, 2, [1., 2.], [0.5, 0.5], 2.)
        call ins%end_batch(2)
        call ins%get_counts(nrows, npose, nsingle, neff_sum)
        call assert_true(nrows == 3_int64 .and. npose == 6_int64 .and. nsingle == 1_int64, 'the counters persist over batches')
        call ins%kill
    end subroutine test_inpl_insert_channel

    !> The teacher's rows (keep_batch_rows, get_posterior_row): each row's posterior over the classes at beta times the
    !! particle's temperature, the masses in descending order summing to one with no tail when every class is kept,
    !! the ratio exp(-T (L_2 - L_1)), the mode's class, distance, in-plane index and refined shift; a particle without
    !! a temperature read at its argmin; and the batch assigned from the kept rows as from the buffers
    subroutine test_dense_posterior_rows()
        real, parameter :: TEMP = 2., BETA = 0.5
        type(prob_candidate)  :: candidates(6)
        type(prob_row_record) :: rec
        integer :: particle_indices(6), j
        real    :: dists(NCLASSES), t_eff
        write(*,'(A)') 'test_dense_posterior_rows'
        params%refine = 'prob'
        call set_candidate(candidates(1),1,0.1)
        call set_candidate(candidates(2),2,0.9)
        call set_candidate(candidates(3),1,0.8)
        call set_candidate(candidates(4),2,0.2)
        call set_candidate(candidates(5),1,0.3)
        call set_candidate(candidates(6),2,0.4)
        candidates(4)%x      =  1.5
        candidates(4)%y      = -2.5
        candidates(4)%has_sh = .true.
        particle_indices = [1,1,2,2,3,3]
        call table%new_pipeline(params, build, pinds)
        table%temps = TEMP
        table%temps(3) = 0.
        do j = 1,size(candidates)
            call table%candidate_buffers(1)%append(particle_indices(j), candidates(j))
        enddo
        call table%keep_batch_rows(1, NPTCLS)
        t_eff = BETA * TEMP
        call table%get_posterior_row(1, BETA, 4, rec, dists)
        call assert_int(pinds(1), rec%iptcl, 'the row names its particle')
        call assert_real(t_eff, rec%temp, 1.e-6, 'the row is read at beta times the particle''s temperature')
        call assert_int(2, rec%nkept, 'every valid class is kept when kmax allows')
        call assert_true(rec%inds(1) == 1 .and. rec%inds(2) == 2, 'the classes in descending mass')
        call assert_real(1., rec%w(1) + rec%w(2), 1.e-5, 'the masses sum to one')
        call assert_real(0., rec%wtail, 1.e-5, 'and leave no tail')
        call assert_real(exp(-t_eff * 0.8), rec%w(2) / rec%w(1), 1.e-4, 'the masses follow exp(-T (L - L_min))')
        call assert_true(rec%imin == 1 .and. abs(rec%lmin - 0.1) < 1.e-6 .and. rec%inpl_mode == 1, 'the mode''s class, loss and rotation')
        call assert_true(.not. rec%has_sh, 'an unrefined mode carries no shift')
        call assert_true(all(abs(dists - [0.1, 0.9]) < 1.e-6), 'the whole row''s losses come back')
        call table%get_posterior_row(2, BETA, 1, rec)
        call assert_true(rec%nkept == 1 .and. rec%inds(1) == 2, 'kmax caps the kept classes')
        call assert_real(1. - rec%w(1), rec%wtail, 1.e-5, 'the tail holds what the kept class leaves')
        call assert_true(rec%has_sh .and. rec%x_mode == 1.5 .and. rec%y_mode == -2.5, 'a refined mode carries its shift')
        call table%get_posterior_row(3, BETA, 4, rec)
        call assert_true(rec%imin == 1 .and. abs(rec%w(1) - 1.) < 1.e-6, 'a particle without a temperature is read at its argmin')
        table%temps = 1.e9
        call table%open_assignment
        call table%assign_batch(1, NPTCLS)
        call table%close_assignment
        call assert_assignment('dense batch from kept rows',[1,2,1],[0.1,0.2,0.3],[0,0,0],[100.,100.,100.])
        call assert_shifts('dense batch from kept rows',[0.,1.5,0.],[0.,-2.5,0.],[.false.,.true.,.false.])
        call rec%kill
        call table%kill
    end subroutine test_dense_posterior_rows

    subroutine set_candidate( candidate, icls, dist )
        type(prob_candidate), intent(out) :: candidate
        integer,              intent(in)  :: icls
        real,                 intent(in)  :: dist
        candidate%iref = icls
        candidate%inpl = icls
        candidate%dist = dist
    end subroutine set_candidate

    !> the module's table as the pipeline's over the three particles: the candidates appended to the first thread
    !! buffer in the given order (as the fill appends them), the temperature sharp enough that the draw is the argmin
    !! (the distances differ by 0.1 or more), the one batch assigned
    subroutine assign_batch_of( particle_indices, candidates )
        integer,              intent(in) :: particle_indices(:)
        type(prob_candidate), intent(in) :: candidates(:)
        integer :: j
        if( size(particle_indices) /= size(candidates) ) THROW_HARD('candidate batch test size mismatch')
        call table%new_pipeline(params, build, pinds)
        table%temps = 1.e9
        do j = 1,size(candidates)
            call table%candidate_buffers(1)%append(particle_indices(j), candidates(j))
        enddo
        call table%open_assignment
        call table%assign_batch(1, NPTCLS)
        call table%close_assignment
    end subroutine assign_batch_of

    subroutine assert_assignment( label, expected_classes, expected_dists, expected_npeaks, expected_fracs )
        character(len=*), intent(in) :: label
        integer,          intent(in) :: expected_classes(NPTCLS)
        real,             intent(in) :: expected_dists(NPTCLS)
        integer,          intent(in) :: expected_npeaks(NPTCLS)
        real,             intent(in) :: expected_fracs(NPTCLS)
        call assert_true(all(table%assgn_map(:)%pind == pinds), label//': particle indices')
        call assert_true(all(table%assgn_map(:)%icls == expected_classes), label//': best class per particle')
        call assert_true(all(abs(table%assgn_map(:)%dist - expected_dists) <= 1.e-6), label//': distance of the best class')
        call assert_true(all(table%assgn_map(:)%inpl == expected_classes), label//': in-plane index of the best class')
        call assert_true(all(table%assgn_map(:)%npeaks == expected_npeaks), label//': peak count')
        call assert_true(all(abs(table%assgn_map(:)%frac - expected_fracs) <= 1.e-6), label//': evaluated fraction')
    end subroutine assert_assignment

    subroutine assert_shifts( label, expected_x, expected_y, expected_has_sh )
        character(len=*), intent(in) :: label
        real,             intent(in) :: expected_x(NPTCLS), expected_y(NPTCLS)
        logical,          intent(in) :: expected_has_sh(NPTCLS)
        call assert_true(all(table%assgn_map(:)%x == expected_x) .and. all(table%assgn_map(:)%y == expected_y),&
            &label//': shifts of the assigned candidates')
        call assert_true(all(table%assgn_map(:)%has_sh .eqv. expected_has_sh), label//': shift flags')
    end subroutine assert_shifts

end module simple_eul_prob_tab2D_tester
