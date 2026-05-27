!@descr: coarray matrix-sum parallelization strategy test
module simple_coarray_matrix_sum_strategy
use simple_core_module_api
use simple_cmdline, only: cmdline
implicit none

public :: coarray_matrix_sum_strategy
public :: coarray_matrix_sum_coarray_strategy
public :: create_coarray_matrix_sum_strategy
private
#include "simple_local_flags.inc"

type, abstract :: coarray_matrix_sum_strategy
contains
    procedure(init_interface),     deferred :: initialize
    procedure(exec_interface),     deferred :: execute
    procedure(finalize_interface), deferred :: finalize_run
    procedure(cleanup_interface),  deferred :: cleanup
    procedure(endmsg_interface),   deferred :: end_message
end type coarray_matrix_sum_strategy

type, extends(coarray_matrix_sum_strategy) :: coarray_matrix_sum_coarray_strategy
    integer :: n = 4
contains
    procedure :: initialize   => coarray_initialize
    procedure :: execute      => coarray_execute
    procedure :: finalize_run => coarray_finalize_run
    procedure :: cleanup      => coarray_cleanup
    procedure :: end_message  => coarray_end_message
end type coarray_matrix_sum_coarray_strategy

abstract interface
    subroutine init_interface( self, cline )
        import :: coarray_matrix_sum_strategy, cmdline
        class(coarray_matrix_sum_strategy), intent(inout) :: self
        class(cmdline),                     intent(inout) :: cline
    end subroutine init_interface

    subroutine exec_interface( self )
        import :: coarray_matrix_sum_strategy
        class(coarray_matrix_sum_strategy), intent(inout) :: self
    end subroutine exec_interface

    subroutine finalize_interface( self )
        import :: coarray_matrix_sum_strategy
        class(coarray_matrix_sum_strategy), intent(inout) :: self
    end subroutine finalize_interface

    subroutine cleanup_interface( self )
        import :: coarray_matrix_sum_strategy
        class(coarray_matrix_sum_strategy), intent(inout) :: self
    end subroutine cleanup_interface

    function endmsg_interface( self ) result( msg )
        import :: coarray_matrix_sum_strategy
        class(coarray_matrix_sum_strategy), intent(in) :: self
        character(len=:), allocatable :: msg
    end function endmsg_interface
end interface

contains

    function create_coarray_matrix_sum_strategy( cline ) result( strategy )
        class(cmdline), intent(in) :: cline
        class(coarray_matrix_sum_strategy), allocatable :: strategy
        allocate(coarray_matrix_sum_coarray_strategy :: strategy)
    end function create_coarray_matrix_sum_strategy

    subroutine coarray_initialize( self, cline )
        class(coarray_matrix_sum_coarray_strategy), intent(inout) :: self
        class(cmdline),                             intent(inout) :: cline
        self%n = 4
        if( cline%defined('box') ) self%n = cline%get_iarg('box')
        if( self%n < 1 ) THROW_HARD('box must be > 0 for coarray matrix sum test')
    end subroutine coarray_initialize

    subroutine coarray_execute( self )
        class(coarray_matrix_sum_coarray_strategy), intent(inout) :: self
#ifdef USE_COARRAYS
        integer, allocatable :: matrix(:,:)[:], matrix_sum(:,:)
        integer :: img, expected_element, expected_total
        allocate(matrix(self%n,self%n)[*], matrix_sum(self%n,self%n))
        matrix = this_image()
        sync all
        if( this_image() == 1 )then
            matrix_sum = 0
            do img = 1,num_images()
                matrix_sum = matrix_sum + matrix(:,:)[img]
            end do
            expected_element = (num_images() * (num_images() + 1)) / 2
            expected_total   = self%n * self%n * expected_element
            if( any(matrix_sum /= expected_element) )then
                THROW_HARD('coarray matrix sum contains unexpected values')
            endif
            write(logfhandle,'(a,i0)') '>>> COARRAY MATRIX SIZE N: ', self%n
            write(logfhandle,'(a,i0)') '>>> COARRAY IMAGES:        ', num_images()
            write(logfhandle,'(a,i0)') '>>> MATRIX ELEMENT SUM:    ', expected_element
            write(logfhandle,'(a,i0)') '>>> TOTAL MATRIX SUM:      ', sum(matrix_sum)
            if( sum(matrix_sum) /= expected_total )then
                THROW_HARD('coarray total matrix sum mismatch')
            endif
        endif
        sync all
        deallocate(matrix, matrix_sum)
#else
        THROW_HARD('coarray_matrix_sum strategy requires building SIMPLE with USE_COARRAYS=ON')
#endif
    end subroutine coarray_execute

    subroutine coarray_finalize_run( self )
        class(coarray_matrix_sum_coarray_strategy), intent(inout) :: self
    end subroutine coarray_finalize_run

    subroutine coarray_cleanup( self )
        class(coarray_matrix_sum_coarray_strategy), intent(inout) :: self
    end subroutine coarray_cleanup

    function coarray_end_message( self ) result( msg )
        class(coarray_matrix_sum_coarray_strategy), intent(in) :: self
        character(len=:), allocatable :: msg
#ifdef USE_COARRAYS
        if( this_image() == 1 )then
            msg = '**** SIMPLE_TEST_COARRAY_MATRIX_SUM NORMAL STOP ****'
        else
            msg = ''
        endif
#else
        msg = '**** SIMPLE_TEST_COARRAY_MATRIX_SUM NORMAL STOP ****'
#endif
    end function coarray_end_message

end module simple_coarray_matrix_sum_strategy
