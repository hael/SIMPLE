!@descr: for all geometry tests
module simple_commanders_test_geometry
use simple_commanders_api
implicit none
#include "simple_local_flags.inc"

type, extends(commander_base) :: commander_test_angres
  contains
    procedure :: execute      => exec_test_angres
end type commander_test_angres

contains

!> angular resolution of the projection-direction spiral as a function of its size.
!! find_angres is the largest third-nearest-neighbour distance (degrees) over the
!! spiral. The reference values are from the original sweep (500..20000 in steps
!! of 500), truncated to a ladder that keeps the O(n^2) cost modest.
subroutine exec_test_angres( self, cline )
    use simple_test_utils, only: begin_test_suite, end_test_suite, assert_true, assert_real, report_summary
    class(commander_test_angres), intent(inout) :: self
    class(cmdline),               intent(inout) :: cline
    integer, parameter :: NDIRS(4)      = [500, 1000, 2000, 4000]
    real,    parameter :: ANGRES_REF(4) = [9.95311069, 7.01425362, 4.97815180, 3.54364777]
    real,    parameter :: REL_TOL       = 0.02
    type(oris) :: os
    real       :: angres(size(NDIRS))
    integer    :: i
    logical    :: test_failed
    call begin_test_suite('angular resolution of the spiral')
    do i = 1,size(NDIRS)
        call os%new(NDIRS(i), is_ptcl=.false.)
        call os%spiral
        angres(i) = os%find_angres()
        write(logfhandle,'(A,I6,A,F10.5)') 'ndirs=', NDIRS(i), '  angres=', angres(i)
        call assert_real(ANGRES_REF(i), angres(i), REL_TOL*ANGRES_REF(i), 'angres matches the recorded value')
        call os%kill
    end do
    ! more directions -> finer angular sampling
    do i = 2,size(NDIRS)
        call assert_true(angres(i) < angres(i-1), 'angres decreases with the number of directions')
    end do
    ! the sampling interval scales as ndirs^(-1/2): four times the directions halves it
    call assert_real(2.0, angres(1)/angres(3), 0.2, 'angres(500)/angres(2000) ~ 2')
    call assert_real(2.0, angres(2)/angres(4), 0.2, 'angres(1000)/angres(4000) ~ 2')
    call end_test_suite
    call report_summary(failed=test_failed)
    if( test_failed ) error stop 1
    call simple_end('**** SIMPLE_TEST_ANGRES NORMAL STOP ****')
end subroutine exec_test_angres

end module simple_commanders_test_geometry
