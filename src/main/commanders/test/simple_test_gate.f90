!@descr: pass/fail checks and metrics against floors of a workflow gate, logged and tabulated in a TSV file
! Every check and metric is logged as PASS or FAIL and written as one row of the
! table (name, value, floor, pass); a metric with a non-finite value fails. A
! metric that is reported without a bound has the floor NO_FLOOR (report).
! The gate has passed when no check or metric failed.
module simple_test_gate
use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
use simple_core_module_api
implicit none

public :: test_gate, NO_FLOOR
private
#include "simple_local_flags.inc"

real, parameter :: NO_FLOOR = -1. !< the floor recorded for a metric that is reported only

type :: test_gate
    private
    integer :: funit   = 0
    logical :: l_open  = .false.
    integer :: nfailed = 0
  contains
    procedure :: new
    procedure :: check
    procedure :: metric
    procedure :: report
    procedure :: passed
    procedure :: kill
end type test_gate

contains

    !> a gate writing its table to fname, which is replaced
    subroutine new( self, fname )
        class(test_gate), intent(inout) :: self
        class(string),    intent(in)    :: fname
        integer :: io_stat
        call self%kill
        call fopen(self%funit, file=fname, status='REPLACE', action='WRITE', iostat=io_stat)
        if( io_stat /= 0 ) THROW_HARD('cannot open the gate table '//fname%to_char())
        self%l_open = .true.
        write(self%funit,'(A)') 'name'//achar(9)//'value'//achar(9)//'floor'//achar(9)//'pass'
    end subroutine new

    !> a pass/fail check, tabulated as value 1 or 0 against floor 1
    subroutine check( self, name, ok )
        class(test_gate), intent(inout) :: self
        character(len=*), intent(in)    :: name
        logical,          intent(in)    :: ok
        if( self%l_open ) write(self%funit,'(A,A,I0,A,A,A,A)') name, achar(9), merge(1,0,ok), achar(9), '1', &
            &achar(9), trim(merge('yes','no ',ok))
        write(logfhandle,'(a)') '    '//trim(merge('PASS:','FAIL:',ok))//' '//name
        if( .not. ok ) self%nfailed = self%nfailed + 1
    end subroutine check

    !> a metric and the floor it is held to; ok is the caller's verdict
    subroutine metric( self, name, val, floor, ok )
        class(test_gate), intent(inout) :: self
        character(len=*), intent(in)    :: name
        real,             intent(in)    :: val, floor
        logical,          intent(in)    :: ok
        logical :: l_ok
        l_ok = ok .and. ieee_is_finite(val)
        if( self%l_open ) write(self%funit,'(A,A,F12.5,A,F12.5,A,A)') name, achar(9), val, achar(9), floor, &
            &achar(9), trim(merge('yes','no ',l_ok))
        write(logfhandle,'(a,f12.5,a,f12.5)') '    '//trim(merge('PASS:','FAIL:',l_ok))//' '//name//' = ', val, &
            &' floor ', floor
        if( .not. l_ok ) self%nfailed = self%nfailed + 1
    end subroutine metric

    !> a metric reported without a bound (fails only when not finite)
    subroutine report( self, name, val )
        class(test_gate), intent(inout) :: self
        character(len=*), intent(in)    :: name
        real,             intent(in)    :: val
        call self%metric(name, val, NO_FLOOR, .true.)
    end subroutine report

    logical function passed( self )
        class(test_gate), intent(in) :: self
        passed = self%nfailed == 0
    end function passed

    subroutine kill( self )
        class(test_gate), intent(inout) :: self
        if( self%l_open ) call fclose(self%funit)
        self%funit   = 0
        self%l_open  = .false.
        self%nfailed = 0
    end subroutine kill

end module simple_test_gate
