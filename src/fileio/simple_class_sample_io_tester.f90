!@descr: unit test routines for the class-sampling checkpoint file (simple_class_sample_io)
! The ragged class_sample array written by abinitio/refine3D and read back by the matcher and the split
! checkpoint: field-by-field round trip, empty classes, replacement of a previously allocated array.
module simple_class_sample_io_tester
use simple_test_utils          ! assertions etc.
use simple_type_defs,          only: class_sample
use simple_string,             only: string
use simple_string_utils,       only: int2str
use simple_syslib,             only: del_file, file_exists
use simple_class_sample_io,    only: write_class_samples, read_class_samples, deallocate_class_samples
implicit none
private
public :: run_all_class_sample_io_tests

character(len=*), parameter :: CS_FILE = 'class_sample_io_tester.bin'
real,             parameter :: TOL     = 1.0e-6

contains

    subroutine run_all_class_sample_io_tests()
        write(*,'(A)') '**** running all class sample I/O tests ****'
        call test_ragged_roundtrip()
        call test_empty_class_roundtrip()
        call test_read_replaces_previous_array()
        call del_file(string(CS_FILE))
    end subroutine run_all_class_sample_io_tests

    ! three classes of different population; every field survives the write/read cycle, in order
    subroutine make_ragged( cs )
        type(class_sample), allocatable, intent(inout) :: cs(:)
        call deallocate_class_samples(cs)
        allocate(cs(3))
        cs(1)%clsind  = 7
        cs(2)%clsind  = 2
        cs(3)%clsind  = 11
        cs(1)%pop     = 1
        cs(2)%pop     = 2
        cs(3)%pop     = 4
        cs(1)%nsample = 1
        cs(2)%nsample = 1
        cs(3)%nsample = 3
        allocate(cs(1)%pinds(1), source=[42])
        allocate(cs(2)%pinds(2), source=[5, 9])
        allocate(cs(3)%pinds(4), source=[100, 3, 77, 12])
        allocate(cs(1)%ccs(1),   source=[0.75])
        allocate(cs(2)%ccs(2),   source=[0.5, -0.25])
        allocate(cs(3)%ccs(4),   source=[0.9, 0.8, 0.1, -0.6])
    end subroutine make_ragged

    subroutine assert_entry_equal( expected, actual, label )
        type(class_sample), intent(in) :: expected, actual
        character(len=*),   intent(in) :: label
        integer :: j
        call assert_int(expected%clsind,  actual%clsind,  label//': clsind')
        call assert_int(expected%pop,     actual%pop,     label//': pop')
        call assert_int(expected%nsample, actual%nsample, label//': nsample')
        call assert_true(allocated(actual%pinds) .and. allocated(actual%ccs), label//': pinds and ccs allocated')
        if( .not. (allocated(actual%pinds) .and. allocated(actual%ccs)) ) return
        call assert_int(size(expected%pinds), size(actual%pinds), label//': pinds size')
        call assert_int(size(expected%ccs),   size(actual%ccs),   label//': ccs size')
        if( size(actual%pinds) /= size(expected%pinds) .or. size(actual%ccs) /= size(expected%ccs) ) return
        do j = 1,size(expected%pinds)
            call assert_int(expected%pinds(j), actual%pinds(j), label//': pinds('//int2str(j)//')')
            call assert_real(expected%ccs(j),  actual%ccs(j), TOL, label//': ccs('//int2str(j)//')')
        end do
    end subroutine assert_entry_equal

    subroutine test_ragged_roundtrip()
        type(class_sample), allocatable :: cs(:), cs_read(:)
        integer :: i
        write(*,'(A)') 'test_ragged_roundtrip'
        call make_ragged(cs)
        call del_file(string(CS_FILE))
        call write_class_samples(cs, string(CS_FILE))
        call assert_true(file_exists(string(CS_FILE)), 'write_class_samples creates the file')
        call read_class_samples(cs_read, string(CS_FILE))
        call assert_true(allocated(cs_read), 'read_class_samples allocates the array')
        if( .not. allocated(cs_read) ) return
        call assert_int(3, size(cs_read), 'one entry per class comes back')
        if( size(cs_read) /= 3 ) return
        do i = 1,3
            call assert_entry_equal(cs(i), cs_read(i), 'class '//int2str(i))
        end do
        call deallocate_class_samples(cs)
        call deallocate_class_samples(cs_read)
    end subroutine test_ragged_roundtrip

    ! a class without members (pinds never allocated, as get_class_sample_stats leaves it) is stored as its
    ! three scalars; it comes back with pop 0 and zero-sized index/correlation arrays, and does not disturb
    ! the populated classes around it
    subroutine test_empty_class_roundtrip()
        type(class_sample), allocatable :: cs(:), cs_read(:)
        write(*,'(A)') 'test_empty_class_roundtrip'
        allocate(cs(3))
        cs(1)%clsind = 1
        cs(1)%pop    = 2
        allocate(cs(1)%pinds(2), source=[8, 6])
        allocate(cs(1)%ccs(2),   source=[0.3, 0.2])
        cs(2)%clsind = 4   ! empty: no pinds/ccs
        cs(3)%clsind = 5
        cs(3)%pop    = 1
        cs(3)%nsample = 1
        allocate(cs(3)%pinds(1), source=[13])
        allocate(cs(3)%ccs(1),   source=[0.95])
        call del_file(string(CS_FILE))
        call write_class_samples(cs, string(CS_FILE))
        call read_class_samples(cs_read, string(CS_FILE))
        call assert_int(3, size(cs_read), 'the empty class keeps its slot')
        if( size(cs_read) /= 3 ) return
        call assert_entry_equal(cs(1), cs_read(1), 'class before the empty one')
        call assert_int(4, cs_read(2)%clsind,  'empty class: clsind survives')
        call assert_int(0, cs_read(2)%pop,     'empty class: pop is 0')
        call assert_int(0, cs_read(2)%nsample, 'empty class: nsample is 0')
        call assert_true(allocated(cs_read(2)%pinds), 'empty class: pinds allocated on read')
        if( allocated(cs_read(2)%pinds) ) call assert_int(0, size(cs_read(2)%pinds), 'empty class: pinds has zero size')
        call assert_true(allocated(cs_read(2)%ccs), 'empty class: ccs allocated on read')
        if( allocated(cs_read(2)%ccs) ) call assert_int(0, size(cs_read(2)%ccs), 'empty class: ccs has zero size')
        call assert_entry_equal(cs(3), cs_read(3), 'class after the empty one')
        call deallocate_class_samples(cs)
        call deallocate_class_samples(cs_read)
    end subroutine test_empty_class_roundtrip

    ! reading into an already allocated array replaces it wholesale (no stale entries, no leak of old sizes)
    subroutine test_read_replaces_previous_array()
        type(class_sample), allocatable :: cs(:), cs_read(:)
        write(*,'(A)') 'test_read_replaces_previous_array'
        allocate(cs_read(5))
        cs_read(1)%clsind = 99
        allocate(cs_read(1)%pinds(3), source=[1, 2, 3])
        allocate(cs_read(1)%ccs(3),   source=[1., 1., 1.])
        call make_ragged(cs)
        call del_file(string(CS_FILE))
        call write_class_samples(cs, string(CS_FILE))
        call read_class_samples(cs_read, string(CS_FILE))
        call assert_int(3, size(cs_read), 'the previous 5-entry array is replaced by the 3 on disk')
        if( size(cs_read) /= 3 ) return
        call assert_entry_equal(cs(1), cs_read(1), 'first class after replacement')
        call assert_entry_equal(cs(3), cs_read(3), 'last class after replacement')
        call deallocate_class_samples(cs)
        call deallocate_class_samples(cs_read)
        call assert_false(allocated(cs_read), 'deallocate_class_samples releases the array')
    end subroutine test_read_replaces_previous_array

end module simple_class_sample_io_tester
