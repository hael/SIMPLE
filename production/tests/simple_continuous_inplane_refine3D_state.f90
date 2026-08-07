module continuous_inplane_refine3D_state_test
implicit none
private
public :: run_search_state_contract

contains

subroutine run_search_state_contract()
use simple_strategy3D_alloc, only: clean_strategy3D, s3D, &
    &seed_continuous_inplane_candidate
use simple_strategy3D_srch, only: strategy3D_srch
implicit none

real :: coordinate
logical :: valid
type(strategy3D_srch) :: search

call seed_continuous_inplane_candidate(1, coordinate, valid)
if( coordinate /= 1. ) error stop 'first grid index was not preserved as the continuous seed'
if( valid ) error stop 'grid seed was incorrectly marked as a valid continuous result'

call seed_continuous_inplane_candidate(288, coordinate, valid)
if( coordinate /= 288. ) error stop 'upper grid index was not preserved as the continuous seed'
if( valid ) error stop 'upper grid seed was incorrectly marked as a valid continuous result'

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
if( search%nsolns /= 1 ) error stop 'first stored candidate did not increment solution count'
if( s3D%proj_space_inplinds(2,1) /= 17 ) error stop 'stored integer in-plane index changed'
if( s3D%proj_space_inplcoords(2,1) /= 17. ) &
    &error stop 'stored continuous coordinate did not retain the grid seed'
if( s3D%proj_space_inplvalid(2,1) ) &
    &error stop 'legacy stored candidate was marked continuous-valid'
if( any(s3D%proj_space_shift(:,2,1) /= [1.25, -0.75]) ) &
    &error stop 'stored candidate shift changed'

! A rejected candidate must not disturb an accepted continuous result.
s3D%proj_space_inplcoords(2,1) = 17.25
s3D%proj_space_inplvalid(2,1)  = .true.
call search%store_solution(2, 18, 0.4, sh=[2., 2.])
if( s3D%proj_space_inplcoords(2,1) /= 17.25 .or. &
    &.not. s3D%proj_space_inplvalid(2,1) ) &
    &error stop 'rejected candidate disturbed continuous state'

! An improving legacy candidate must replace and invalidate stale continuous state.
call search%store_solution(2, 19, 0.6, sh=[0.5, -0.5])
if( s3D%proj_space_inplinds(2,1) /= 19 .or. &
    &s3D%proj_space_inplcoords(2,1) /= 19. .or. &
    &s3D%proj_space_inplvalid(2,1) ) &
    &error stop 'improving grid candidate did not reset continuous state'

call clean_strategy3D

write(*,'(a)') 'REFINE3D_INPL_CONT_STATE GRID_SEED/STORE/REPLACE: PASS'
end subroutine run_search_state_contract

end module continuous_inplane_refine3D_state_test
