!@descr: unit tests for the 2D probability tables: the pipeline's dense and sparse selection (simple_eul_prob_tab2D)
! A batch of candidates for three particles and two classes goes through the pipeline's table (the thread buffers
! the fill appends to, the batch's rows or candidate store, the per-particle draw from the tempered posterior at a
! sharp temperature); the assignment must give each particle its best class, distance, in-plane index, peak count
! and fraction: dense (refine=prob) and sparse (refine=prob_neigh). The draw on synthetic tables with ties, flat
! and invalid rows at several temperatures and seeds.
module simple_eul_prob_tab2D_tester
use simple_core_module_api
use simple_builder,            only: builder
use simple_parameters,         only: parameters
use simple_eul_prob_tab2D,     only: eul_prob_tab2D
use simple_eul_prob_tab_utils, only: prob_candidate
use simple_type_defs,          only: ptcl_ref
use simple_rnd,                only: seed_rnd_fixed
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
        call test_sparse_batch()
        call test_sparse_shared_best_class()
        call test_dense_posterior_assignment()
        call test_sparse_posterior_assignment()
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
