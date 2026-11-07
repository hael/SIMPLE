!@descr: closed-form capacity and geometry checks for matcher memory planning
module simple_matcher_memory_tester
use, intrinsic :: iso_fortran_env, only: int64
use simple_matcher_memory, only: matcher_batch_capacity, matcher_memory_terms
use simple_parameters,    only: parameters
use simple_test_utils,    only: assert_int, assert_true
implicit none
private
public :: run_all_matcher_memory_tests

contains

    subroutine run_all_matcher_memory_tests()
        integer(int64), parameter :: GIB = 1024_int64**3, SLOT = 256_int64*1024_int64
        class(parameters), allocatable :: params
        integer(int64) :: fixed, slot_bytes, larger_fixed, larger_slot
        write(*,'(A)') '**** running matcher memory tests ****'
        call assert_int(7, matcher_batch_capacity(1000000,1000,3*GIB+7*SLOT+17,3*GIB,SLOT), &
            &'capacity uses 64-bit byte counts and discards a partial slot')
        call assert_int(5, matcher_batch_capacity(5,1000,10*GIB,0_int64,SLOT), 'capacity is bounded by the sample')
        call assert_int(12, matcher_batch_capacity(1000,12,10*GIB,0_int64,SLOT), 'capacity preserves the preferred batch')
        call assert_int(1, matcher_batch_capacity(0,12,10*GIB,0_int64,SLOT), 'an empty part keeps one calculator slot')
        call assert_int(1, matcher_batch_capacity(1000,12,GIB,2*GIB,SLOT), 'fixed storage over budget leaves one slot')
        call assert_int(12, matcher_batch_capacity(1000,12,0_int64,GIB,SLOT), 'unknown memory keeps the preferred batch')
        allocate(params)
        params%box = 128
        params%box_crop = 64
        params%box_croppd = 128
        params%pftsz = 48
        params%kfromto = [2,20]
        params%fromp = 1
        params%top = 1000
        call matcher_memory_terms(params,50,1000,4,0,fixed,slot_bytes)
        call matcher_memory_terms(params,100,1000,4,0,larger_fixed,larger_slot)
        call assert_true(larger_fixed > fixed, 'more references increase fixed storage')
        call assert_true(larger_slot > slot_bytes, 'more references increase row storage per particle')
        call matcher_memory_terms(params,50,1000,8,0,larger_fixed,larger_slot)
        call assert_true(larger_fixed > fixed, 'larger teams reserve more thread scratch')
        call assert_true(larger_slot == slot_bytes, 'thread count does not change bytes per particle')
        call matcher_memory_terms(params,50,1000,4,params%box_crop,larger_fixed,larger_slot)
        call assert_true(larger_fixed > fixed, '2D reserves the class sums')
        call assert_true(larger_slot > slot_bytes, '2D reserves batch images and candidate staging')
        call matcher_memory_terms(params,50,1000,4,params%box,fixed,slot_bytes)
        call assert_true(fixed > larger_fixed .and. slot_bytes > larger_slot, 'full-size class sums reserve more memory than cropped sums')
        deallocate(params)
    end subroutine run_all_matcher_memory_tests

end module simple_matcher_memory_tester
