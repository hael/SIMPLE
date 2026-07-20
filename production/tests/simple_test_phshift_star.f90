program simple_test_phshift_star
use simple_core_module_api
use simple_starproject_utils, only: RELION_PHASE_MULT
use simple_test_utils,       only: assert_real, report_summary, tests_failed
implicit none

real, parameter :: PHASE_DEG = 45.
real            :: phase_internal, phase_exported

! STAR stores rlnPhaseShift in degrees; SIMPLE stores the CTF phase in radians.
! The same multiplier is used by micrograph, particle-2D, and particle-3D maps.
phase_internal = PHASE_DEG * RELION_PHASE_MULT
call assert_real(PI/4., phase_internal, 1.e-6, &
    &'RELION phase-shift import converts degrees to SIMPLE radians')

phase_exported = phase_internal / RELION_PHASE_MULT
call assert_real(PHASE_DEG, phase_exported, 1.e-6, &
    &'RELION phase-shift export converts SIMPLE radians to degrees')

call report_summary()
if( tests_failed > 0 ) error stop 1
end program simple_test_phshift_star
