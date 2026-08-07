module continuous_inplane_refine3D_joint_test
implicit none
private
public :: run_joint_state_contract

contains

subroutine run_joint_state_contract()
use simple_core_module_api, only: dp
use simple_strategy3D_alloc, only: clean_strategy3D, s3D
use simple_strategy3D_srch, only: strategy3D_srch
implicit none

type(strategy3D_srch) :: search

allocate(s3D%proj_space_shift(2,2,1), s3D%proj_space_corrs(2,1), &
    &s3D%proj_space_inplcoords(2,1), s3D%proj_space_inplinds(2,1), &
    &s3D%proj_space_inplvalid(2,1))
s3D%proj_space_shift      = 0.
s3D%proj_space_corrs      = -huge(1.)
s3D%proj_space_inplcoords = 0.
s3D%proj_space_inplinds   = 0
s3D%proj_space_inplvalid  = .false.
search%ithr = 1

call search%store_solution(2, 17, 0.5, sh=[1.25, -0.75])
call search%store_continuous_solution(2, 18, 17.625_dp, 0.6, [1.1, -0.7])

if( s3D%proj_space_inplinds(2,1) /= 18 ) &
    &error stop 'joint result did not retain its nearest legacy index'
if( abs(s3D%proj_space_inplcoords(2,1) - 17.625) > 1.e-6 ) &
    &error stop 'joint result did not retain its continuous coordinate'
if( .not. s3D%proj_space_inplvalid(2,1) ) &
    &error stop 'joint result was not marked valid'
if( abs(s3D%proj_space_corrs(2,1) - 0.6) > 1.e-6 ) &
    &error stop 'joint result score was not committed'
if( any(abs(s3D%proj_space_shift(:,2,1) - [1.1, -0.7]) > 1.e-6) ) &
    &error stop 'joint result shift was not committed'
if( s3D%proj_space_corrs(1,1) > -huge(1.)/2. ) &
    &error stop 'joint result changed a different reference candidate'

call clean_strategy3D

write(*,'(a)') 'REFINE3D_INPL_CONT_JOINT SELECTED_REFERENCE/STORE: PASS'
end subroutine run_joint_state_contract

end module continuous_inplane_refine3D_joint_test
