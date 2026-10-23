!@descr: unit tests for the streamed 2D probability tables: dense and sparse merge and assignment (simple_eul_prob_tab2D)
! A worker's candidate stream for three particles and two classes is written through the table's
! writer, merged into the global table, assigned and written; the assignment file must give each
! particle its best class, distance, in-plane index, peak count and fraction: dense (refine=prob)
! and sparse (refine=prob_snhc). The per-particle posterior draw (prob_assign=posterior) on the
! synthetic tables of the reference tests, and the particle temperatures through the part files.
module simple_eul_prob_tab2D_tester
use simple_core_module_api
use simple_builder,            only: builder
use simple_parameters,         only: parameters
use simple_eul_prob_tab2D,     only: eul_prob_tab2D
use simple_eul_prob_tab_utils, only: prob_candidate, prob_candidate_store, sample_likelihood_index
use simple_type_defs,          only: ptcl_ref
use simple_rnd,                only: seed_rnd_fixed
use simple_test_utils
!$ use omp_lib,                only: omp_get_max_threads, omp_get_num_procs, omp_set_num_threads
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
        call test_dense_stream()
        call test_sparse_stream()
        call test_sparse_single_neighborhood()
        call test_dense_assign_against_reference()
        call test_sparse_assign_against_reference()
        call test_dense_posterior_assignment()
        call test_sparse_posterior_assignment()
        call test_parallel_read('prob')
        call test_parallel_read('prob_snhc')
        call table%kill
    end subroutine run_all_eul_prob_tab2D_tests

    !> prob_assign=posterior on the dense synthetic table of test_dense_assign_against_reference (ties, a flat
    !! particle, invalid distances): at a huge temperature every particle with a unique finite minimum takes its
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
        params%l_prob_assign_posterior = .true.
        params%prob_beta = 1.
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
        dists(3,9)  = 0.01                 ! a class tie at the minimum
        dists(11,9) = 0.01
        dists(5,12) = huge(1.)             ! invalid distances
        dists(:,13) = huge(1.)
        l_unique = .true.
        l_unique([1,9,13]) = .false.
        ! a huge temperature: the argmin
        call tab%new(params, build, pinds_t)
        call set_table(tab, 1.e9)
        call seed_rnd_fixed(5001)
        call tab%ref_assign
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
        call set_table(tab, 0.)
        call seed_rnd_fixed(5001)
        call tab%ref_assign
        call assert_int(0, count(pack(tab%assgn_map(:)%icls, l_unique) /= pack(icls_a, l_unique)),&
            &'a particle without a temperature takes its argmin')
        ! T = 2: some particles leave the mode, at the same seed the same ones
        call set_table(tab, 2.)
        call seed_rnd_fixed(5002)
        call tab%ref_assign
        icls_a = tab%assgn_map(:)%icls
        nwrong = 0
        do i = 1, NP
            if( .not. l_unique(i) ) cycle
            loc = minloc(dists(:,i))
            if( icls_a(i) /= loc(1) ) nwrong = nwrong + 1
            if( tab%assgn_map(i)%dist /= dists(icls_a(i),i) ) nwrong = nwrong + NP
        enddo
        call assert_true(nwrong > 0 .and. nwrong < NP, 'at T = 2 some particles leave the mode, with the drawn class''s distance')
        call set_table(tab, 2.)
        call seed_rnd_fixed(5002)
        call tab%ref_assign
        call assert_int(0, count(tab%assgn_map(:)%icls /= icls_a), 'the same seed gives the same dense posterior assignment')
        call set_table(tab, 2.)
        call seed_rnd_fixed(5003)
        call tab%ref_assign
        icls_b = tab%assgn_map(:)%icls
        call assert_true(any(icls_b /= icls_a), 'another seed gives another dense posterior assignment')
        call tab%kill
        params%l_prob_assign_posterior = .false.
        params%ncls        = ncls_prev
        params%npeaks_inpl = npeaks_prev
        params%refine      = refine_prev

    contains

        subroutine set_table( t, temp )
            type(eul_prob_tab2D), intent(inout) :: t
            real,                 intent(in)    :: temp
            integer :: ip, ic
            do ip = 1, NP
                do ic = 1, NC
                    t%rows(ic,ip)%inpl = ic
                    t%rows(ic,ip)%dist = dists(ic,ip)
                enddo
            enddo
            t%temps = temp
        end subroutine set_table

    end subroutine test_dense_posterior_assignment

    !> prob_assign=posterior on the sparse synthetic store of test_sparse_assign_against_reference (1 to 20
    !! candidates per particle, ties, candidates without an in-plane index): at a huge temperature every particle
    !! takes a candidate at its smallest valid distance, with the peak count and fraction of the frontier's record;
    !! at T = 2 some particles leave the mode; the same seed gives the same assignment and another seed another
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
        params%refine      = 'prob_snhc'
        params%l_prob_assign_posterior = .true.
        params%prob_beta = 1.
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
        call set_store(tab, 1.e9)
        call seed_rnd_fixed(6001)
        call tab%ref_assign
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
        call set_store(tab, 2.)
        call seed_rnd_fixed(6002)
        call tab%ref_assign
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
        call set_store(tab, 2.)
        call seed_rnd_fixed(6002)
        call tab%ref_assign
        call assert_int(0, count(tab%assgn_map(:)%icls /= icls_a), 'the same seed gives the same sparse posterior assignment')
        call tab%kill
        call set_store(tab, 2.)
        call seed_rnd_fixed(6003)
        call tab%ref_assign
        icls_b = tab%assgn_map(:)%icls
        call assert_true(any(icls_b /= icls_a), 'another seed gives another sparse posterior assignment')
        call tab%kill
        deallocate(tab)
        params%l_prob_assign_posterior = .false.
        params%ncls        = ncls_prev
        params%npeaks_inpl = npeaks_prev
        params%refine      = refine_prev

    contains

        subroutine set_store( t, temp )
            class(eul_prob_tab2D), intent(inout) :: t
            real,                  intent(in)    :: temp
            call t%new(params, build, pinds_t)
            call t%candidate_store%new_ragged(pinds_t, counts)
            t%candidate_store%candidates = cands
            t%temps = temp
        end subroutine set_store

    end subroutine test_sparse_posterior_assignment

    !> refine=prob: every particle carries a candidate per class; the smaller distance wins, and a refined
    !! candidate keeps its shift through the dense part file
    subroutine test_dense_stream()
        type(prob_candidate) :: candidates(6)
        integer :: particle_indices(6)
        write(*,'(A)') 'test_dense_stream'
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
        call write_candidate_stream('prob2d_dense_part1.dat',particle_indices,candidates)
        call table%new(params,build,pinds)
        call table%read_tabs_to_glob(string('prob2d_dense_part'),1,1)
        call table%ref_assign
        call table%write_assignment(string('prob2d_dense_assignment.dat'))
        call assert_assignment('prob2d_dense_assignment.dat',[1,2,1],[0.1,0.2,0.3],&
            &[0,0,0],[100.,100.,100.])
        call assert_shifts('prob2d_dense_assignment.dat',[0.,1.5,0.],[0.,-2.5,0.],[.false.,.true.,.false.])
        call table%kill
        call del_file('prob2d_dense_part1.dat')
        call del_file('prob2d_dense_assignment.dat')
    end subroutine test_dense_stream

    !> refine=prob_snhc: the stream carries only the evaluated candidates, one or two per particle
    subroutine test_sparse_stream()
        type(prob_candidate) :: candidates(4)
        integer :: particle_indices(4)
        write(*,'(A)') 'test_sparse_stream'
        params%refine = 'prob_snhc'
        call set_candidate(candidates(1),1,0.1)
        call set_candidate(candidates(2),2,0.2)
        call set_candidate(candidates(3),1,0.3)
        call set_candidate(candidates(4),2,0.4)
        particle_indices = [1,2,3,3]
        call write_candidate_stream('prob2d_sparse_part1.dat',particle_indices,candidates)
        call table%new(params,build,pinds)
        call table%read_tabs_to_glob(string('prob2d_sparse_part'),1,1)
        call table%ref_assign
        call table%write_assignment(string('prob2d_sparse_assignment.dat'))
        call assert_assignment('prob2d_sparse_assignment.dat',[1,2,1],[0.1,0.2,0.3],&
            &[1,1,2],[50.,50.,100.])
        call table%kill
        call del_file('prob2d_sparse_part1.dat')
        call del_file('prob2d_sparse_assignment.dat')
    end subroutine test_sparse_stream

    !> A one-class support still selects the globally best class frontier.
    subroutine test_sparse_single_neighborhood()
        type(prob_candidate) :: candidates(6)
        integer :: particle_indices(6)
        write(*,'(A)') 'test_sparse_single_neighborhood'
        params%refine = 'prob_snhc'
        call set_candidate(candidates(1),1,0.5)
        call set_candidate(candidates(2),2,0.1)
        call set_candidate(candidates(3),1,0.4)
        call set_candidate(candidates(4),2,0.2)
        call set_candidate(candidates(5),1,0.3)
        call set_candidate(candidates(6),2,0.6)
        particle_indices = [1,1,2,2,3,3]
        call write_candidate_stream('prob2d_sparse_single_part1.dat',particle_indices,candidates)
        call table%new(params,build,pinds)
        call table%read_tabs_to_glob(string('prob2d_sparse_single_part'),1,1)
        call table%ref_assign
        call table%write_assignment(string('prob2d_sparse_single_assignment.dat'))
        call assert_assignment('prob2d_sparse_single_assignment.dat',[2,2,1],[0.1,0.2,0.3],&
            &[2,2,2],[100.,100.,100.])
        call table%kill
        call del_file('prob2d_sparse_single_part1.dat')
        call del_file('prob2d_sparse_single_assignment.dat')
    end subroutine test_sparse_single_neighborhood

    !> The dense assignment against the one it replaced (reference_assign2D, which rescans the frontier every step):
    !! synthetic tables of 150 particles and 20 classes with ties (a flat particle, identical particles, a class
    !! tie at the minimum, invalid distances), several draw seeds, the same class for every particle; and the
    !! assignment at one thread and 20 times at every processor, the same classes and distances
    subroutine test_dense_assign_against_reference()
        integer, parameter :: NP = 150, NC = 20, NSEEDS = 4
        type(eul_prob_tab2D) :: tab
        real,    allocatable :: dists(:,:)
        integer, allocatable :: icls_ref(:), pinds_t(:), icls1(:)
        real,    allocatable :: dist1(:)
        integer :: i, icls, iseed, nwrong, irep, nthr_prev, nthr_full, ncls_prev, npeaks_prev
        character(len=len(params%refine)) :: refine_prev
        write(*,'(A)') 'test_dense_assign_against_reference'
        ncls_prev   = params%ncls
        npeaks_prev = params%npeaks_inpl
        refine_prev = params%refine
        params%ncls        = NC
        params%npeaks_inpl = NC
        params%refine      = 'prob'
        nthr_prev = 1
        nthr_full = 1
        !$ nthr_prev = omp_get_max_threads()
        !$ nthr_full = omp_get_num_procs()
        pinds_t = [(i, i=1,NP)]
        allocate(dists(NC,NP))
        call seed_rnd_fixed(20261019)
        do i = 1, NP
            do icls = 1, NC
                dists(icls,i) = 4. * ran3()
            enddo
        enddo
        dists(:,1)  = 1.5                  ! flat
        dists(:,8)  = dists(:,7)           ! identical particles
        dists(3,9)  = 0.01                 ! a class tie at the minimum
        dists(11,9) = 0.01
        dists(5,12) = huge(1.)             ! invalid distances
        dists(:,13) = huge(1.)
        nwrong = 0
        do iseed = 1, NSEEDS
            call tab%new(params, build, pinds_t)
            call set_table(tab)
            call seed_rnd_fixed(1000 + iseed)
            call tab%ref_assign
            call seed_rnd_fixed(1000 + iseed)
            call reference_assign2D(dists, tab%nhood_sz, icls_ref)
            nwrong = nwrong + count(tab%assgn_map(:)%icls /= icls_ref)
            call tab%kill
        enddo
        call assert_int(0, nwrong, 'the dense 2D assignment equals the rescanning assignment it replaced')
        call seed_rnd_fixed(1001)
        call reference_assign2D(dists, NC, icls1)
        call seed_rnd_fixed(1002)
        call reference_assign2D(dists, NC, icls_ref)
        call assert_true(any(icls1 /= icls_ref), 'the draws change the dense 2D assignment')
        ! one thread against every processor
        call tab%new(params, build, pinds_t)
        call set_table(tab)
        !$ call omp_set_num_threads(1)
        call seed_rnd_fixed(2000)
        call tab%ref_assign
        icls1 = tab%assgn_map(:)%icls
        dist1 = tab%assgn_map(:)%dist
        !$ call omp_set_num_threads(nthr_full)
        nwrong = 0
        do irep = 1, 20
            call set_table(tab)
            call seed_rnd_fixed(2000)
            call tab%ref_assign
            nwrong = nwrong + count(tab%assgn_map(:)%icls /= icls1) + count(tab%assgn_map(:)%dist /= dist1)
        enddo
        !$ call omp_set_num_threads(nthr_prev)
        call assert_int(0, nwrong, 'the dense 2D assignment at every processor equals the one at one thread')
        call tab%kill
        params%ncls        = ncls_prev
        params%npeaks_inpl = npeaks_prev
        params%refine      = refine_prev

    contains

        subroutine set_table( t )
            type(eul_prob_tab2D), intent(inout) :: t
            integer :: ip, ic
            do ip = 1, NP
                do ic = 1, NC
                    t%rows(ic,ip)%inpl = ic
                    t%rows(ic,ip)%dist = dists(ic,ip)
                enddo
            enddo
        end subroutine set_table

    end subroutine test_dense_assign_against_reference

    !> The dense 2D assignment as it was (every class active): each class sorts the particles; each step scans every
    !! class's head to its first available particle and draws among the heads. icls_out(i) is particle i's class.
    subroutine reference_assign2D( dists, nhood, icls_out )
        real,                 intent(in)  :: dists(:,:)
        integer,              intent(in)  :: nhood
        integer, allocatable, intent(out) :: icls_out(:)
        integer, allocatable :: stab(:,:), heads(:), eligible(:), order(:)
        real,    allocatable :: work(:), cdists(:), sorted(:)
        logical, allocatable :: avail(:)
        integer :: nc, np, ic, ip, neligible, chosen, iptcl
        nc = size(dists,1)
        np = size(dists,2)
        allocate(stab(np,nc), heads(nc), eligible(nc), order(nc), work(np), cdists(nc), sorted(nc), avail(np))
        allocate(icls_out(np), source=0)
        do ic = 1, nc
            stab(:,ic) = [(ip, ip=1,np)]
            work = dists(ic,:)
            call hpsort(work, stab(:,ic))
        enddo
        heads = 1
        avail = .true.
        do while( any(avail) )
            neligible = 0
            do ic = 1, nc
                do while( heads(ic) <= np )
                    if( avail(stab(heads(ic),ic)) ) exit
                    heads(ic) = heads(ic) + 1
                enddo
                if( heads(ic) <= np )then
                    neligible = neligible + 1
                    eligible(neligible) = ic
                    cdists(neligible)   = dists(ic,stab(heads(ic),ic))
                endif
            enddo
            if( neligible == 0 ) exit
            call sample_likelihood_index(neligible, cdists, min(nhood, neligible), huge(1.0), sorted, order, chosen)
            ic    = eligible(chosen)
            iptcl = stab(heads(ic),ic)
            avail(iptcl)    = .false.
            icls_out(iptcl) = ic
        enddo
    end subroutine reference_assign2D

    !> The sparse 2D assignment against the one it replaced (reference_sparse_assign2D: class lists gathered and
    !! sorted serially, heads read through the store): synthetic stores of 300 particles and 20 classes with 1 to 20
    !! candidates each, ties (distances on a coarse grid, a flat particle, identical particles), candidates without an
    !! in-plane index, several draw seeds; and the assignment at one thread and 20 times at every processor
    subroutine test_sparse_assign_against_reference()
        integer, parameter :: NP = 300, NC = 20, NSEEDS = 4
        class(eul_prob_tab2D), allocatable :: tab
        type(prob_candidate),  allocatable :: cands(:)
        integer, allocatable :: counts(:), pinds_t(:), pos_ref(:), npeaks_ref(:), icls1(:), inpl1(:), npeaks1(:)
        real,    allocatable :: frac_ref(:), dist1(:), frac1(:)
        integer :: i, k, n, iseed, irep, nwrong, nthr_prev, nthr_full, ncls_prev, npeaks_prev
        character(len=len(params%refine)) :: refine_prev
        write(*,'(A)') 'test_sparse_assign_against_reference'
        ncls_prev   = params%ncls
        npeaks_prev = params%npeaks_inpl
        refine_prev = params%refine
        params%ncls        = NC
        params%npeaks_inpl = NC
        params%refine      = 'prob_snhc'
        nthr_prev = 1
        nthr_full = 1
        !$ nthr_prev = omp_get_max_threads()
        !$ nthr_full = omp_get_num_procs()
        pinds_t = [(i, i=1,NP)]
        call seed_rnd_fixed(20261019)
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
                cands(n)%dist = real(nint(10. * ran3())) / 10.    ! ties
                if( k > 1 .and. ran3() < 0.05 ) cands(n)%inpl = 0  ! not a candidate
            enddo
        enddo
        n = sum(counts(1:6))
        cands(n+1:n+NC)%dist = 0.7                                 ! a flat particle
        cands(n+NC+1:n+2*NC) = cands(n+1:n+NC)                     ! and an identical one
        cands(n+1:n+2*NC)%inpl = 3
        allocate(tab)
        nwrong = 0
        do iseed = 1, NSEEDS
            call set_store(tab)
            call seed_rnd_fixed(3000 + iseed)
            call tab%ref_assign
            call seed_rnd_fixed(3000 + iseed)
            call reference_sparse_assign2D(tab%candidate_store, NC, tab%nhood_sz, pos_ref, frac_ref, npeaks_ref)
            nwrong = nwrong + count(tab%assgn_map(:)%icls /= cands(pos_ref)%iref) +&
                &count(tab%assgn_map(:)%inpl /= cands(pos_ref)%inpl) + count(tab%assgn_map(:)%dist /= cands(pos_ref)%dist) +&
                &count(tab%assgn_map(:)%npeaks /= npeaks_ref) + count(tab%assgn_map(:)%frac /= frac_ref)
            call tab%kill
        enddo
        call assert_int(0, nwrong, 'the sparse 2D assignment equals the serial assignment it replaced')
        ! one thread against every processor
        !$ call omp_set_num_threads(1)
        call set_store(tab)
        call seed_rnd_fixed(4000)
        call tab%ref_assign
        icls1   = tab%assgn_map(:)%icls
        inpl1   = tab%assgn_map(:)%inpl
        dist1   = tab%assgn_map(:)%dist
        frac1   = tab%assgn_map(:)%frac
        npeaks1 = tab%assgn_map(:)%npeaks
        call tab%kill
        !$ call omp_set_num_threads(nthr_full)
        nwrong = 0
        do irep = 1, 20
            call set_store(tab)
            call seed_rnd_fixed(4000)
            call tab%ref_assign
            nwrong = nwrong + count(tab%assgn_map(:)%icls /= icls1) + count(tab%assgn_map(:)%inpl /= inpl1) +&
                &count(tab%assgn_map(:)%dist /= dist1) + count(tab%assgn_map(:)%frac /= frac1) +&
                &count(tab%assgn_map(:)%npeaks /= npeaks1)
            call tab%kill
        enddo
        !$ call omp_set_num_threads(nthr_prev)
        call assert_int(0, nwrong, 'the sparse 2D assignment at every processor equals the one at one thread')
        deallocate(tab)
        params%ncls        = ncls_prev
        params%npeaks_inpl = npeaks_prev
        params%refine      = refine_prev

    contains

        subroutine set_store( t )
            class(eul_prob_tab2D), intent(inout) :: t
            call t%new(params, build, pinds_t)
            call t%candidate_store%new_ragged(pinds_t, counts)
            t%candidate_store%candidates = cands
        end subroutine set_store

    end subroutine test_sparse_assign_against_reference

    !> the sparse 2D assignment as it was before its class lists were sorted in parallel: the reference of
    !! test_sparse_assign_against_reference. Out: the store position of each particle's candidate, its fraction and
    !! peak count. Every class is active.
    subroutine reference_sparse_assign2D( store, nclasses, nhood_sz, pos_out, frac_out, npeaks_out )
        type(prob_candidate_store), intent(in)  :: store
        integer,                    intent(in)  :: nclasses, nhood_sz
        integer, allocatable,       intent(out) :: pos_out(:), npeaks_out(:)
        real,    allocatable,       intent(out) :: frac_out(:)
        integer, allocatable :: class_counts(:), class_offsets(:), class_fill(:), class_pos(:), active_classes(:)
        integer, allocatable :: class_candidates(:), candidate_particles(:), order(:), work_candidates(:)
        integer, allocatable :: sel_classes(:), sel_pos(:), likelihood_order(:)
        real,    allocatable :: class_dist(:), work_d(:), sel_dists(:), likelihood_dists(:)
        logical, allocatable :: ptcl_avail(:)
        integer :: nptcls, i, k, pos, idx, icls, assigned_ptcl, assigned_idx, nleft, nsel, total, nactive, maxclass
        integer :: last_class, assigned_candidate
        real    :: best_dist, huge_val
        nptcls = size(store%offsets) - 1
        allocate(pos_out(nptcls), npeaks_out(nptcls), source=0)
        allocate(frac_out(nptcls), source=0.)
        huge_val = huge(1.0) / 2.0
        allocate(class_counts(nclasses), class_offsets(nclasses+1), class_fill(nclasses), class_pos(nclasses))
        class_counts = 0
        do i = 1,nptcls
            do pos = int(store%offsets(i)), int(store%offsets(i+1))-1
                if( .not. valid(pos) ) cycle
                icls = store%candidates(pos)%iref
                class_counts(icls) = class_counts(icls) + 1
            enddo
        enddo
        nactive = count(class_counts > 0)
        allocate(active_classes(nactive))
        active_classes = pack((/(icls, icls=1,nclasses)/), class_counts > 0)
        class_offsets(1) = 1
        do icls = 1,nclasses
            class_offsets(icls+1) = class_offsets(icls) + class_counts(icls)
        enddo
        total = class_offsets(nclasses+1) - 1
        allocate(class_candidates(total), candidate_particles(size(store%candidates)), source=0)
        class_fill = class_offsets(1:nclasses)
        do i = 1,nptcls
            do pos = int(store%offsets(i)), int(store%offsets(i+1))-1
                candidate_particles(pos) = i
                if( .not. valid(pos) ) cycle
                icls = store%candidates(pos)%iref
                idx = class_fill(icls)
                class_candidates(idx) = pos
                class_fill(icls) = idx + 1
            enddo
        enddo
        maxclass = max(1, maxval(class_counts))
        allocate(work_d(maxclass), order(maxclass), work_candidates(maxclass))
        do idx = 1,nactive
            icls = active_classes(idx)
            k    = class_counts(icls)
            if( k <= 1 ) cycle
            pos = class_offsets(icls)
            work_candidates(1:k) = class_candidates(pos:pos+k-1)
            do i = 1,k
                work_d(i) = store%candidates(work_candidates(i))%dist
            enddo
            order(1:k) = (/(i,i=1,k)/)
            call hpsort(work_d(1:k), order(1:k))
            do i = 1,k
                class_candidates(pos+i-1) = work_candidates(order(i))
            enddo
        enddo
        allocate(class_dist(nclasses), sel_pos(nclasses), sel_classes(nactive), sel_dists(nactive),&
            &likelihood_dists(nactive), likelihood_order(nactive), ptcl_avail(nptcls))
        ptcl_avail = .true.
        nleft      = nptcls
        class_pos  = 1
        class_dist = huge(1.0)
        sel_pos    = 0
        nsel       = 0
        do idx = 1,nactive
            icls = active_classes(idx)
            call advance_head(icls)
            call sync_class(icls)
        enddo
        do while( nleft > 0 )
            if( nsel == 0 ) exit
            call sample_likelihood_index(nsel, sel_dists, min(nhood_sz, nsel), huge_val, likelihood_dists,&
                &likelihood_order, assigned_idx)
            icls = sel_classes(assigned_idx)
            assigned_candidate = class_candidates(class_offsets(icls) + class_pos(icls) - 1)
            assigned_ptcl = candidate_particles(assigned_candidate)
            ptcl_avail(assigned_ptcl) = .false.
            nleft = nleft - 1
            call record(assigned_ptcl, assigned_candidate)
            do idx = int(store%offsets(assigned_ptcl)), int(store%offsets(assigned_ptcl+1))-1
                if( .not. valid(idx) ) cycle
                icls = store%candidates(idx)%iref
                if( class_pos(icls) > class_counts(icls) ) cycle
                pos = class_offsets(icls)
                assigned_candidate = class_candidates(pos + class_pos(icls) - 1)
                if( candidate_particles(assigned_candidate) /= assigned_ptcl ) cycle
                call advance_head(icls)
                call sync_class(icls)
            enddo
        enddo
        do i = 1,nptcls
            if( .not. ptcl_avail(i) ) cycle
            assigned_candidate = 0
            best_dist = huge(1.0)
            do pos = int(store%offsets(i)), int(store%offsets(i+1))-1
                if( .not. valid(pos) ) cycle
                if( store%candidates(pos)%dist < best_dist )then
                    best_dist          = store%candidates(pos)%dist
                    assigned_candidate = pos
                endif
            enddo
            call record(i, assigned_candidate)
        enddo

    contains

        logical function valid( p )
            integer, intent(in) :: p
            valid = store%candidates(p)%iref >= 1 .and. store%candidates(p)%iref <= nclasses
            if( valid ) valid = store%candidates(p)%inpl > 0
        end function valid

        subroutine record( iptcl, p )
            integer, intent(in) :: iptcl, p
            integer :: q, npeaks
            npeaks = 0
            do q = int(store%offsets(iptcl)), int(store%offsets(iptcl+1))-1
                if( store%candidates(q)%iref < 1 .or. store%candidates(q)%iref > nclasses ) cycle
                if( store%candidates(q)%inpl > 0 ) npeaks = npeaks + 1
            enddo
            pos_out(iptcl)    = p
            npeaks_out(iptcl) = npeaks
            frac_out(iptcl)   = 100. * real(npeaks) / real(nclasses)
        end subroutine record

        subroutine advance_head( icls_loc )
            integer, intent(in) :: icls_loc
            integer :: cand
            do while( class_pos(icls_loc) <= class_counts(icls_loc) )
                cand = class_candidates(class_offsets(icls_loc) + class_pos(icls_loc) - 1)
                if( ptcl_avail(candidate_particles(cand)) )then
                    class_dist(icls_loc) = store%candidates(cand)%dist
                    return
                endif
                class_pos(icls_loc) = class_pos(icls_loc) + 1
            enddo
            class_dist(icls_loc) = huge(1.0)
        end subroutine advance_head

        subroutine sync_class( icls_loc )
            integer, intent(in) :: icls_loc
            if( sel_pos(icls_loc) > 0 )then
                if( class_dist(icls_loc) >= huge(1.0) )then
                    pos = sel_pos(icls_loc)
                    if( pos < nsel )then
                        last_class = sel_classes(nsel)
                        sel_classes(pos)    = last_class
                        sel_dists(pos)      = sel_dists(nsel)
                        sel_pos(last_class) = pos
                    endif
                    sel_pos(icls_loc) = 0
                    nsel = nsel - 1
                else
                    sel_dists(sel_pos(icls_loc)) = class_dist(icls_loc)
                endif
            else
                if( class_dist(icls_loc) < huge(1.0) )then
                    nsel = nsel + 1
                    sel_classes(nsel)  = icls_loc
                    sel_dists(nsel)    = class_dist(icls_loc)
                    sel_pos(icls_loc)  = nsel
                endif
            endif
        end subroutine sync_class

    end subroutine reference_sparse_assign2D

    !> The part reads run in parallel over the parts: synthetic part files (7 parts of interleaved particles, 20
    !! classes, refined entries with shifts in the dense layout, 1 to 20 candidates per particle in the stream layout,
    !! seeds) read at one thread and 20 times at every processor; the master's table must agree bit for bit
    subroutine test_parallel_read( refine )
        character(len=*), intent(in) :: refine
        integer, parameter :: NP = 700, NC = 20, NPARTS = 7, NREP = 20
        class(eul_prob_tab2D), allocatable :: worker, ref, tab
        type(prob_candidate) :: candidate
        integer, allocatable :: pinds_t(:), pinds_part(:)
        integer :: i, i_loc, ipart, k, ncand, irep, nwrong, nthr_prev, nthr_full, ncls_prev, npeaks_prev
        character(len=len(params%refine)) :: refine_prev
        character(len=32) :: fname
        write(*,'(A)') 'test_parallel_read '//refine
        ncls_prev   = params%ncls
        npeaks_prev = params%npeaks_inpl
        refine_prev = params%refine
        params%ncls        = NC
        params%npeaks_inpl = NC
        params%refine      = refine
        nthr_prev = 1
        nthr_full = 1
        !$ nthr_prev = omp_get_max_threads()
        !$ nthr_full = omp_get_num_procs()
        pinds_t = [(100 + 3*i, i=1,NP)]
        call seed_rnd_fixed(20261019)
        allocate(worker, ref, tab)
        do ipart = 1, NPARTS
            pinds_part = pack(pinds_t, [(mod(i-1,NPARTS) == ipart-1, i=1,NP)])
            call worker%new_worker(params, build, pinds_part)
            worker%seed_nrots = 36
            do i_loc = 1, size(pinds_part)
                worker%seed_shifts(:,i_loc) = [ran3() - 0.5, ran3() - 0.5]
                worker%seed_has_sh(i_loc)   = ran3() > 0.3
                worker%temps(i_loc)         = real(pinds_part(i_loc))
            enddo
            write(fname,'(A,I1,A)') 'prob2d_par_part', ipart, '.dat'
            call worker%begin_write(string(trim(fname)))
            do i_loc = 1, size(pinds_part)
                ncand = NC
                if( refine == 'prob_snhc' ) ncand = 1 + int(ran3() * real(NC-1))
                do k = 1, ncand
                    candidate        = prob_candidate()
                    candidate%iref   = k
                    if( refine == 'prob_snhc' ) candidate%iref = 1 + mod(k * 7 + i_loc, NC)
                    candidate%inpl   = 1 + int(ran3() * 35.)
                    candidate%dist   = 4. * ran3()
                    if( k <= 2 )then
                        candidate%x      = ran3() - 0.5
                        candidate%y      = ran3() - 0.5
                        candidate%has_sh = .true.
                    endif
                    call worker%candidate_buffers(1)%append(i_loc, candidate)
                enddo
            enddo
            call worker%write_tab(string(trim(fname)))
            call worker%kill
        enddo
        !$ call omp_set_num_threads(1)
        call ref%new(params, build, pinds_t)
        call ref%read_tabs_to_glob(string('prob2d_par_part'), NPARTS, 1)
        !$ call omp_set_num_threads(nthr_full)
        nwrong = 0
        do irep = 1, NREP
            call tab%new(params, build, pinds_t)
            call tab%read_tabs_to_glob(string('prob2d_par_part'), NPARTS, 1)
            nwrong = nwrong + count(tab%seed_shifts /= ref%seed_shifts) + count(tab%seed_has_sh .neqv. ref%seed_has_sh)
            nwrong = nwrong + count(tab%temps /= ref%temps)
            if( tab%seed_nrots /= ref%seed_nrots ) nwrong = nwrong + 1
            if( refine == 'prob' )then
                nwrong = nwrong + count(tab%rows%dist /= ref%rows%dist) + count(tab%rows%inpl /= ref%rows%inpl)
                nwrong = nwrong + count(tab%refined_offsets /= ref%refined_offsets)
                if( size(tab%refined) /= size(ref%refined) )then
                    nwrong = nwrong + 1
                else
                    nwrong = nwrong + count(tab%refined%i /= ref%refined%i) + count(tab%refined%ri /= ref%refined%ri) +&
                        &count(tab%refined%x /= ref%refined%x) + count(tab%refined%y /= ref%refined%y)
                endif
            else
                nwrong = nwrong + count(tab%candidate_store%offsets /= ref%candidate_store%offsets)
                if( size(tab%candidate_store%candidates) /= size(ref%candidate_store%candidates) )then
                    nwrong = nwrong + 1
                else
                    associate( a => tab%candidate_store%candidates, b => ref%candidate_store%candidates )
                    nwrong = nwrong + count(a%iref /= b%iref) + count(a%inpl /= b%inpl) + count(a%dist /= b%dist) +&
                        &count(a%x /= b%x) + count(a%y /= b%y) + count(a%has_sh .neqv. b%has_sh)
                    end associate
                endif
            endif
            call tab%kill
        enddo
        !$ call omp_set_num_threads(nthr_prev)
        call assert_int(0, nwrong, 'the 2D part read ('//refine//') at every processor equals the read at one thread')
        call assert_int(0, count(ref%temps /= real(pinds_t)), 'the particle temperatures land at their master particles')
        if( refine == 'prob' )then
            call assert_int(2*NP, size(ref%refined), 'the dense master holds every refined entry')
        else
            call assert_true(all(ref%candidate_store%offsets(2:) > ref%candidate_store%offsets(:NP)),&
                &'every particle of the sparse master has candidates')
        endif
        call ref%kill
        do ipart = 1, NPARTS
            write(fname,'(A,I1,A)') 'prob2d_par_part', ipart, '.dat'
            call del_file(trim(fname))
        enddo
        deallocate(worker, ref, tab)
        params%ncls        = ncls_prev
        params%npeaks_inpl = npeaks_prev
        params%refine      = refine_prev
    end subroutine test_parallel_read

    subroutine set_candidate( candidate, icls, dist )
        type(prob_candidate), intent(out) :: candidate
        integer,              intent(in)  :: icls
        real,                 intent(in)  :: dist
        candidate%iref = icls
        candidate%inpl = icls
        candidate%dist = dist
    end subroutine set_candidate

    !> a worker's part file through the table's own writer: the candidates in one chunk, in the given order
    subroutine write_candidate_stream( fname, particle_indices, candidates )
        character(len=*),     intent(in) :: fname
        integer,              intent(in) :: particle_indices(:)
        type(prob_candidate), intent(in) :: candidates(:)
        class(eul_prob_tab2D), allocatable :: worker
        integer :: j
        if( size(particle_indices) /= size(candidates) ) THROW_HARD('candidate stream test size mismatch')
        allocate(worker)
        call worker%new_worker(params,build,pinds)
        call worker%begin_write(string(fname))
        do j = 1,size(candidates)
            call worker%candidate_buffers(1)%append(particle_indices(j), candidates(j))
        enddo
        call worker%write_tab(string(fname))
        call worker%kill
        deallocate(worker)
    end subroutine write_candidate_stream

    subroutine assert_assignment( fname, expected_classes, expected_dists, expected_npeaks, expected_fracs )
        character(len=*), intent(in) :: fname
        integer,          intent(in) :: expected_classes(NPTCLS)
        real,             intent(in) :: expected_dists(NPTCLS)
        integer,          intent(in) :: expected_npeaks(NPTCLS)
        real,             intent(in) :: expected_fracs(NPTCLS)
        type(ptcl_ref) :: assignments(NPTCLS)
        integer :: funit, io_stat, nptcls_file
        call fopen(funit,string(fname),access='STREAM',action='READ',status='OLD',iostat=io_stat)
        call fileiochk('simple_eul_prob_tab2D_tester; read assignment '//fname,io_stat)
        read(funit,pos=1) nptcls_file
        call assert_int(NPTCLS, nptcls_file, fname//': particle count')
        if( nptcls_file /= NPTCLS )then
            call fclose(funit)
            return
        endif
        read(funit,pos=sizeof(nptcls_file)+1) assignments
        call fclose(funit)
        call assert_true(all(assignments(:)%pind == pinds), fname//': particle indices')
        call assert_true(all(assignments(:)%icls == expected_classes), fname//': best class per particle')
        call assert_true(all(abs(assignments(:)%dist - expected_dists) <= 1.e-6), fname//': distance of the best class')
        call assert_true(all(assignments(:)%inpl == expected_classes), fname//': in-plane index of the best class')
        call assert_true(all(assignments(:)%npeaks == expected_npeaks), fname//': peak count')
        call assert_true(all(abs(assignments(:)%frac - expected_fracs) <= 1.e-6), fname//': evaluated fraction')
    end subroutine assert_assignment

    subroutine assert_shifts( fname, expected_x, expected_y, expected_has_sh )
        character(len=*), intent(in) :: fname
        real,             intent(in) :: expected_x(NPTCLS), expected_y(NPTCLS)
        logical,          intent(in) :: expected_has_sh(NPTCLS)
        type(ptcl_ref) :: assignments(NPTCLS)
        integer :: funit, io_stat, nptcls_file
        call fopen(funit,string(fname),access='STREAM',action='READ',status='OLD',iostat=io_stat)
        call fileiochk('simple_eul_prob_tab2D_tester; read assignment '//fname,io_stat)
        read(funit,pos=1) nptcls_file
        read(funit,pos=sizeof(nptcls_file)+1) assignments
        call fclose(funit)
        call assert_true(all(assignments(:)%x == expected_x) .and. all(assignments(:)%y == expected_y),&
            &fname//': shifts of the assigned candidates')
        call assert_true(all(assignments(:)%has_sh .eqv. expected_has_sh), fname//': shift flags')
    end subroutine assert_shifts

end module simple_eul_prob_tab2D_tester
