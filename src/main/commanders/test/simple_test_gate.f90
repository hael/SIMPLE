!@descr: pass/fail checks and metrics against floors of a workflow gate, logged and tabulated in a TSV file
! Every check and metric is logged as PASS or FAIL and written as one row of the
! table (name, value, floor, pass, sense); a metric with a non-finite value fails.
! The sense says how the floor bounds the value: min (value >= floor), max (value
! <= floor, a ceiling), check (1 or 0 against 1) or report (no bound, floor NO_FLOOR).
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
    procedure, private :: write_row
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
        write(self%funit,'(A)') 'name'//achar(9)//'value'//achar(9)//'floor'//achar(9)//'pass'//achar(9)//'sense'
    end subroutine new

    !> a pass/fail check, tabulated as value 1 or 0 against floor 1
    subroutine check( self, name, ok )
        class(test_gate), intent(inout) :: self
        character(len=*), intent(in)    :: name
        logical,          intent(in)    :: ok
        if( self%l_open )then
            write(self%funit,'(A,A,I0,A,A,A,A,A,A)') name, achar(9), merge(1,0,ok), achar(9), '1', &
                &achar(9), trim(merge('yes','no ',ok)), achar(9), 'check'
            ! row by row, so a killed run keeps its table
            flush(self%funit)
        endif
        write(logfhandle,'(a)') '    '//trim(merge('PASS:','FAIL:',ok))//' '//name
        flush(logfhandle)
        if( .not. ok ) self%nfailed = self%nfailed + 1
    end subroutine check

    !> a metric and the floor it is held to; ok is the caller's verdict. upper: the floor is a
    !! ceiling (smaller is better)
    subroutine metric( self, name, val, floor, ok, upper )
        class(test_gate),  intent(inout) :: self
        character(len=*),  intent(in)    :: name
        real,              intent(in)    :: val, floor
        logical,           intent(in)    :: ok
        logical, optional, intent(in)    :: upper
        character(len=3) :: sense
        logical :: l_ok
        l_ok  = ok .and. ieee_is_finite(val)
        sense = 'min'
        if( present(upper) )then
            if( upper ) sense = 'max'
        endif
        call self%write_row(name, val, floor, l_ok, sense)
    end subroutine metric

    !> a metric reported without a bound (fails only when not finite)
    subroutine report( self, name, val )
        class(test_gate), intent(inout) :: self
        character(len=*), intent(in)    :: name
        real,             intent(in)    :: val
        call self%write_row(name, val, NO_FLOOR, ieee_is_finite(val), 'report')
    end subroutine report

    subroutine write_row( self, name, val, floor, ok, sense )
        class(test_gate), intent(inout) :: self
        character(len=*), intent(in)    :: name, sense
        real,             intent(in)    :: val, floor
        logical,          intent(in)    :: ok
        if( self%l_open )then
            write(self%funit,'(A,A,ES13.6,A,ES13.6,A,A,A,A)') name, achar(9), val, achar(9), floor, &
                &achar(9), trim(merge('yes','no ',ok)), achar(9), sense
            ! row by row, so a killed run keeps its table
            flush(self%funit)
        endif
        write(logfhandle,'(a,g0.6,a,g0.6)') '    '//trim(merge('PASS:','FAIL:',ok))//' '//name//' = ', val, &
            &' floor ', floor
        flush(logfhandle)
        if( .not. ok ) self%nfailed = self%nfailed + 1
    end subroutine write_row

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
