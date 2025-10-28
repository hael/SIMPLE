! concrete strategy3D: greedy refinement
module simple_strategy3D_greedy_sub
include 'simple_lib.f08'
use simple_strategy3D_alloc
use simple_strategy3D_utils
use simple_parameters,       only: params_glob
use simple_builder,          only: build_glob
use simple_strategy3D,       only: strategy3D
use simple_strategy3D_srch,  only: strategy3D_srch, strategy3D_spec
use simple_polarft_corrcalc, only: pftcc_glob
implicit none

public :: strategy3D_greedy_sub
private
#include "simple_local_flags.inc"

type, extends(strategy3D) :: strategy3D_greedy_sub
contains
    procedure :: new         => new_greedy_sub
    procedure :: srch        => srch_greedy_sub
    procedure :: kill        => kill_greedy_sub
    procedure :: oris_assign => oris_assign_greedy_sub
end type strategy3D_greedy_sub

contains

    subroutine new_greedy_sub( self, spec )
        class(strategy3D_greedy_sub), intent(inout) :: self
        class(strategy3D_spec),       intent(inout) :: spec
        call self%s%new(spec)
        self%spec = spec
    end subroutine new_greedy_sub

    subroutine srch_greedy_sub( self, ithr )
        use simple_eul_prob_tab, only: angle_sampling, eulprob_dist_switch
        class(strategy3D_greedy_sub), intent(inout) :: self
        integer,                      intent(in)    :: ithr
        integer   :: iref, isample, loc(1), iproj, ipeak, inds(self%s%nrots)
        real      :: inpl_corrs(self%s%nrots), sorted_corrs(self%s%nrots)
        logical   :: lnns(params_glob%nspace)
        type(ori) :: o
        if( build_glob%spproj_field%get_state(self%s%iptcl) > 0 )then
            ! set thread index
            self%s%ithr = ithr
            ! prep
            call self%s%prep4srch
            ! shift search on previous best reference
            call self%s%inpl_srch_first
            ! search
            do isample=1,self%s%nrefs_sub
                iref = s3D%srch_order_sub(isample,self%s%ithr) ! set the reference index
                if( s3D%state_exists(s3D%proj_space_state(iref)) )then
                    ! identify the top scoring in-plane angle
                    if( params_glob%l_sh_first )then
                        call pftcc_glob%gencorrs(iref, self%s%iptcl, self%s%xy_first, inpl_corrs)
                    else
                        call pftcc_glob%gencorrs(iref, self%s%iptcl, inpl_corrs)
                    endif
                    loc = maxloc(inpl_corrs)
                    call self%s%store_solution(iref, loc(1), inpl_corrs(loc(1)))
                endif
            end do
            ! prepare peak orientations
            call extract_peak_oris(self%s)
            ! construct multi-neighborhood search space from subspace peaks
            lnns = .false.
            do ipeak = 1, self%s%npeaks
                call self%s%opeaks%get_ori(ipeak, o)
                call build_glob%pgrpsyms%nearest_proj_neighbors(build_glob%eulspace, o, params_glob%athres, lnns)
            end do
            ! include the previous best ori in the multi-neighborhood search
            call build_glob%pgrpsyms%nearest_proj_neighbors(build_glob%eulspace, self%s%o_prev, params_glob%athres, lnns)
            ! count the number of nearest neighbors
            self%s%nnn = count(lnns)
            ! search
            do iproj=1,params_glob%nspace
                if( .not. lnns(iproj) ) cycle
                iref = (self%s%prev_state - 1) * params_glob%nspace + iproj
                if( s3D%state_exists(s3D%proj_space_state(iref)) )then
                    ! identify the top scoring in-plane angle
                    if( params_glob%l_sh_first )then
                        call pftcc_glob%gencorrs(iref, self%s%iptcl, self%s%xy_first, inpl_corrs)
                    else
                        call pftcc_glob%gencorrs(iref, self%s%iptcl, inpl_corrs)
                    endif
                    if( params_glob%l_prob_inpl )then
                        loc = angle_sampling(eulprob_dist_switch(inpl_corrs), sorted_corrs, inds, s3D%smpl_inpl_athres(s3D%proj_space_state(iref)))
                    else
                        loc = maxloc(inpl_corrs)
                    endif
                    call self%s%store_solution(iref, loc(1), inpl_corrs(loc(1)))
                endif
            end do
            ! in greedy mode, we evaluate all refs
            self%s%nrefs_eval = self%s%nnn
            ! take care of the in-planes
            call self%s%inpl_srch_peaks
            ! prepare orientation
            call self%oris_assign
            ! cleanup
            call o%kill
        else
            call build_glob%spproj_field%reject(self%s%iptcl)
        endif
    end subroutine srch_greedy_sub

    subroutine oris_assign_greedy_sub( self )
        class(strategy3D_greedy_sub), intent(inout) :: self
        call extract_peak_ori(self%s)
    end subroutine oris_assign_greedy_sub

    subroutine kill_greedy_sub( self )
        class(strategy3D_greedy_sub), intent(inout) :: self
        call self%s%kill
    end subroutine kill_greedy_sub

end module simple_strategy3D_greedy_sub
