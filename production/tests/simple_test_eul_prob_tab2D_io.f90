!@descr: validates streamed dense/sparse 2D probability-table merge and assignment
program simple_test_eul_prob_tab2D_io
use, intrinsic :: iso_fortran_env, only: int64
use simple_core_module_api
use simple_builder,            only: builder
use simple_parameters,         only: parameters
use simple_eul_prob_tab2D,     only: eul_prob_tab2D
use simple_eul_prob_tab_utils, only: prob_candidate, write_seed_shift_table
use simple_type_defs,          only: ptcl_ref
implicit none
#include "simple_local_flags.inc"

type(parameters), target :: params
type(builder),    target :: build
type(eul_prob_tab2D)     :: table
integer, parameter :: NPTCLS = 3, NCLASSES = 2
integer :: pinds(NPTCLS)

pinds = [10,20,30]
params%ncls         = NCLASSES
params%npeaks_inpl  = NCLASSES
params%l_doshift    = .false.

call test_dense_stream
call test_sparse_stream
call table%kill
call simple_end('**** SIMPLE_EUL_PROB_TAB2D_IO TEST NORMAL STOP ****')

contains

    subroutine test_dense_stream
        type(prob_candidate) :: candidates(6)
        integer :: particle_indices(6)
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

    subroutine test_sparse_stream
        type(prob_candidate) :: candidates(4)
        integer :: particle_indices(4)
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
        call fileiochk('simple_test_eul_prob_tab2D_io; write stream '//fname,io_stat)
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
        integer :: funit, io_stat, nptcls_file, i
        call fopen(funit,string(fname),access='STREAM',action='READ',status='OLD',iostat=io_stat)
        call fileiochk('simple_test_eul_prob_tab2D_io; read assignment '//fname,io_stat)
        read(funit,pos=1) nptcls_file
        if( nptcls_file /= NPTCLS ) THROW_HARD('2D probability assignment particle-count mismatch')
        read(funit,pos=sizeof(nptcls_file)+1) assignments
        call fclose(funit)
        do i = 1,NPTCLS
            if( assignments(i)%pind /= pinds(i) ) THROW_HARD('2D probability assignment particle mismatch')
            if( assignments(i)%icls /= expected_classes(i) ) THROW_HARD('2D probability assignment class mismatch')
            if( abs(assignments(i)%dist-expected_dists(i)) > 1.e-6 ) THROW_HARD('2D probability assignment distance mismatch')
            if( assignments(i)%inpl /= expected_classes(i) ) THROW_HARD('2D probability assignment in-plane mismatch')
            if( assignments(i)%npeaks /= expected_npeaks(i) ) THROW_HARD('2D probability assignment peak-count mismatch')
            if( abs(assignments(i)%frac-expected_fracs(i)) > 1.e-6 ) THROW_HARD('2D probability assignment fraction mismatch')
        enddo
    end subroutine assert_assignment

end program simple_test_eul_prob_tab2D_io
