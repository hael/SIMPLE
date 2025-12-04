! concrete strategy3D: refinement
module simple_strategy3D_greedy_smpl
include 'simple_lib.f08'
use simple_strategy3D_alloc
use simple_strategy3D_utils
use simple_parameters,       only: params_glob
use simple_builder,          only: build_glob
use simple_strategy3D,       only: strategy3D
use simple_strategy3D_srch,  only: strategy3D_spec
use simple_polarft_calc, only: pftc_glob
implicit none

public :: strategy3D_greedy_smpl
private
#include "simple_local_flags.inc"

type, extends(strategy3D) :: strategy3D_greedy_smpl
contains
    procedure :: new         => new_greedy_smpl
    procedure :: srch        => srch_greedy_smpl
    procedure :: oris_assign => oris_assign_greedy_smpl
    procedure :: kill        => kill_greedy_smpl
end type strategy3D_greedy_smpl

contains

    subroutine new_greedy_smpl( self, spec )
        class(strategy3D_greedy_smpl), intent(inout) :: self
        class(strategy3D_spec),     intent(inout) :: spec
        call self%s%new(spec)
        self%spec = spec
    end subroutine new_greedy_smpl

    subroutine srch_greedy_smpl( self, ithr )
        use simple_eul_prob_tab, only: angle_sampling, eulprob_dist_switch
        class(strategy3D_greedy_smpl), intent(inout) :: self
        integer,                       intent(in)    :: ithr
        integer :: iref, isample, loc(1), inds(self%s%nrots)
        real    :: inpl_corrs(self%s%nrots), sorted_corrs(self%s%nrots)
        if( build_glob%spproj_field%get_state(self%s%iptcl) > 0 )then
            ! set thread index
            self%s%ithr = ithr
            ! prep
            call self%s%prep4srch
            ! shift search on previous best reference
            call self%s%inpl_srch_first
            ! search
            do isample=1,self%s%nrefs
                iref = s3D%srch_order(isample,self%s%ithr)  ! set the stochastic reference index
                if( s3D%state_exists(s3D%proj_space_state(iref)) )then
                    ! identify the top scoring in-plane angle
                    if( params_glob%l_sh_first )then
                        call pftc_glob%gen_objfun_vals(iref, self%s%iptcl, self%s%xy_first, inpl_corrs)
                    else
                        call pftc_glob%gen_objfun_vals(iref, self%s%iptcl, [0.,0.],         inpl_corrs)
                    endif
                    loc = angle_sampling(eulprob_dist_switch(inpl_corrs), sorted_corrs, inds, s3D%smpl_inpl_athres(s3D%proj_space_state(iref)))
                    call self%s%store_solution(iref, loc(1), inpl_corrs(loc(1)))
                end if
            end do
            ! in greedy mode, we evaluate all refs
            self%s%nrefs_eval = self%s%nrefs
            ! take care of the in-planes
            call self%s%inpl_srch
            ! prepare weights and orientations
            call self%oris_assign
        else
            call build_glob%spproj_field%reject(self%s%iptcl)
        endif
    end subroutine srch_greedy_smpl

    subroutine oris_assign_greedy_smpl( self )
        class(strategy3D_greedy_smpl), intent(inout) :: self
        call extract_peak_ori(self%s)
    end subroutine oris_assign_greedy_smpl

    subroutine kill_greedy_smpl( self )
        class(strategy3D_greedy_smpl), intent(inout) :: self
        call self%s%kill
    end subroutine kill_greedy_smpl

end module simple_strategy3D_greedy_smpl
