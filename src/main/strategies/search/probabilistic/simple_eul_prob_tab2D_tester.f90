!@descr: unit tests for the streamed 2D probability tables: dense and sparse merge and assignment (simple_eul_prob_tab2D)
! A worker's candidate stream (header, particle indices, seed-shift table, candidates) is written
! by hand for three particles and two classes, merged into the global table, assigned and
! written; the assignment file must give each particle its best class, distance, in-plane
! index, peak count and fraction: dense (refine=prob) and sparse (refine=prob_snhc).
module simple_eul_prob_tab2D_tester
use, intrinsic :: iso_fortran_env, only: int64
use simple_core_module_api
use simple_builder,            only: builder
use simple_parameters,         only: parameters
use simple_eul_prob_tab2D,     only: eul_prob_tab2D
use simple_eul_prob_tab_utils, only: prob_candidate, write_seed_shift_table, sample_likelihood_index
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
        call table%kill
    end subroutine run_all_eul_prob_tab2D_tests

    !> refine=prob: every particle carries a candidate per class; the smaller distance wins
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
        particle_indices = [1,1,2,2,3,3]
        call write_candidate_stream('prob2d_dense_part1.dat',particle_indices,candidates)
        call table%new(params,build,pinds)
        call table%read_tabs_to_glob(string('prob2d_dense_part'),1,1)
        call table%ref_assign
        call table%write_assignment(string('prob2d_dense_assignment.dat'))
        call assert_assignment('prob2d_dense_assignment.dat',[1,2,1],[0.1,0.2,0.3],&
            &[0,0,0],[100.,100.,100.])
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
                    t%loc_tab(ic,ip)%iref = ic
                    t%loc_tab(ic,ip)%inpl = ic
                    t%loc_tab(ic,ip)%dist = dists(ic,ip)
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

    subroutine set_candidate( candidate, icls, dist )
        type(prob_candidate), intent(out) :: candidate
        integer,              intent(in)  :: icls
        real,                 intent(in)  :: dist
        candidate%iref = icls
        candidate%inpl = icls
        candidate%dist = dist
    end subroutine set_candidate

    subroutine write_candidate_stream( fname, particle_indices, candidates )
        character(len=*),     intent(in) :: fname
        integer,              intent(in) :: particle_indices(:)
        type(prob_candidate), intent(in) :: candidates(:)
        integer :: funit, io_stat, chunk_n
        integer(int64) :: header(4), addr
        real    :: seed_shifts(2,NPTCLS)
        logical :: seed_has_sh(NPTCLS)
        if( size(particle_indices) /= size(candidates) ) THROW_HARD('candidate stream test size mismatch')
        seed_shifts = 0.
        seed_has_sh = .false.
        chunk_n = size(candidates)
        header = [int(NCLASSES,int64),int(NPTCLS,int64),int(chunk_n,int64),1_int64]
        call fopen(funit,string(fname),access='STREAM',action='WRITE',status='REPLACE',iostat=io_stat)
        call fileiochk('simple_eul_prob_tab2D_tester; write stream '//fname,io_stat)
        write(funit,pos=1) header
        addr = sizeof(header) + 1
        write(funit,pos=addr) pinds
        addr = addr + sizeof(pinds)
        call write_seed_shift_table(funit,addr,8,seed_shifts,seed_has_sh)
        write(funit,pos=addr) chunk_n
        addr = addr + sizeof(chunk_n)
        write(funit,pos=addr) particle_indices
        addr = addr + sizeof(particle_indices)
        write(funit,pos=addr) candidates
        call fclose(funit)
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

end module simple_eul_prob_tab2D_tester
