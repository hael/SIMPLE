program simple_test_phshift_star
use simple_core_module_api
use simple_starproject_utils, only: RELION_PHASE_DEG2RAD
use simple_test_utils,       only: assert_real, report_summary, tests_failed
implicit none

real, parameter :: PHASE_DEG = 45., BEYOND_PI_DEG = 250.
real            :: phase_internal, phase_exported

! STAR stores rlnPhaseShift in degrees; SIMPLE stores the CTF phase in radians.
! The same conversion factor is used by micrograph, particle-2D, and
! particle-3D maps.
phase_internal = PHASE_DEG * RELION_PHASE_DEG2RAD
call assert_real(PI/4., phase_internal, 1.e-6, &
    &'RELION phase-shift import converts degrees to SIMPLE radians')

phase_exported = phase_internal / RELION_PHASE_DEG2RAD
call assert_real(PHASE_DEG, phase_exported, 1.e-6, &
    &'RELION phase-shift export converts SIMPLE radians to degrees')

! RELION does not constrain rlnPhaseShift to [0,180): its aberration fit can move a
! constant gamma offset into that column. Such a value must survive import and
! export unchanged, because folding it at pi would negate the transfer function
! relative to what RELION computes from the same STAR file.
phase_internal = canonical_phshift(BEYOND_PI_DEG * RELION_PHASE_DEG2RAD)
call assert_real(BEYOND_PI_DEG * RELION_PHASE_DEG2RAD, phase_internal, 1.e-6, &
    &'a RELION phase shift beyond 180 degrees is not folded on import')

phase_exported = phase_internal / RELION_PHASE_DEG2RAD
call assert_real(BEYOND_PI_DEG, phase_exported, 1.e-4, &
    &'a phase shift beyond 180 degrees round-trips back to RELION degrees')

call report_summary()
if( tests_failed > 0 ) error stop 1
end program simple_test_phshift_star
