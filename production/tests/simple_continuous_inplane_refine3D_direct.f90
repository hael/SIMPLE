module continuous_inplane_refine3D_direct_test
implicit none
private
public :: run_direct_route_contract

contains

subroutine run_direct_route_contract()
use simple_strategy3D_alloc, only: clean_strategy3D, s3D
use simple_strategy3D_srch, only: strategy3D_srch
implicit none

type(strategy3D_srch) :: search
real, parameter :: expected_shift(2) = [1.25, -0.75]

! polish-only policy: legacy post-selection refinement is NEVER bypassed;
! the bypass predicate no longer exists (removing it from the type is the
! compile-time guarantee), and the joint solve only polishes the committed
! pose afterwards
search%continuous_active = .true.
if( search%joint_evaluation_invalid(.true.) ) &
    &error stop 'valid joint no-improvement result was marked invalid'
if( .not. search%joint_evaluation_invalid(.false.) ) &
    &error stop 'invalid joint result was not identified for candidate retention'

allocate(s3D%proj_space_shift(2,2,1), s3D%proj_space_corrs(2,1), &
    &s3D%proj_space_inplcoords(2,1), s3D%proj_space_inplinds(2,1), &
    &s3D%proj_space_inplvalid(2,1))
s3D%proj_space_shift      = 0.
s3D%proj_space_corrs      = -huge(1.)
s3D%proj_space_inplcoords = 0.
s3D%proj_space_inplinds   = 0
s3D%proj_space_inplvalid  = .false.
s3D%proj_space_corrs(2,1)      = 0.625
s3D%proj_space_inplcoords(2,1) = 17.
s3D%proj_space_inplinds(2,1)   = 17
search%ithr = 1

call search%store_discrete_seed_solution(2, 23, 0.75, expected_shift)
if( any(abs(s3D%proj_space_shift(:,2,1) - expected_shift) > 1.e-6) ) &
    &error stop 'no-improvement route did not store the selected seed shift'
if( abs(s3D%proj_space_corrs(2,1) - 0.75) > 1.e-6 ) &
    &error stop 'no-improvement route did not store the selected seed score'
if( s3D%proj_space_inplinds(2,1) /= 23 .or. &
    &abs(s3D%proj_space_inplcoords(2,1) - 23.) > 1.e-6 ) &
    &error stop 'no-improvement route did not replace the incoming grid angle'
if( s3D%proj_space_inplvalid(2,1) ) &
    &error stop 'finite no-improvement was incorrectly marked continuous-valid'
if( any(s3D%proj_space_shift(:,1,1) /= 0.) ) &
    &error stop 'retaining the selected shift changed a different reference'

call clean_strategy3D

write(*,'(a)') &
    &'REFINE3D_INPL_CONT_DIRECT BYPASS/NO_IMPROVEMENT/SCORE_SHIFT/INVALID: PASS'
end subroutine run_direct_route_contract

end module continuous_inplane_refine3D_direct_test
