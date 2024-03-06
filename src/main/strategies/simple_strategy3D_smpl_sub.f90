! concrete strategy3D: stochastic top sampling
module simple_strategy3D_smpl_sub
include 'simple_lib.f08'
use simple_strategy3D_alloc
use simple_strategy3D_utils
use simple_parameters,       only: params_glob
use simple_builder,          only: build_glob
use simple_strategy3D,       only: strategy3D
use simple_strategy3D_srch,  only: strategy3D_srch, strategy3D_spec
use simple_polarft_corrcalc, only: pftcc_glob
use simple_eul_prob_tab,     only: eulprob_dist_switch
implicit none

public :: strategy3D_smpl_sub
private
#include "simple_local_flags.inc"

type, extends(strategy3D) :: strategy3D_smpl_sub
    type(strategy3D_srch) :: s
    type(strategy3D_spec) :: spec
contains
    procedure :: new         => new_smpl_sub
    procedure :: srch        => srch_smpl_sub
    procedure :: oris_assign => oris_assign_smpl_sub
    procedure :: kill        => kill_smpl_sub
end type strategy3D_smpl_sub

contains

    subroutine new_smpl_sub( self, spec )
        class(strategy3D_smpl_sub), intent(inout) :: self
        class(strategy3D_spec), intent(inout) :: spec
        call self%s%new(spec)
        self%spec = spec
    end subroutine new_smpl_sub

    subroutine srch_smpl_sub( self, ithr )
        class(strategy3D_smpl_sub), intent(inout) :: self
        integer,                intent(in)    :: ithr
        integer   :: iref, isample, loc(1), iproj, irot, ipeak
        integer   :: inds(self%s%nrots)
        logical   :: lnns(params_glob%nspace)
        real      :: inpl_corrs(self%s%nrots), sorted_corrs(self%s%nrots), corr
        type(ori) :: o
        ! execute search
        if( build_glob%spproj_field%get_state(self%s%iptcl) > 0 )then
            ! set thread index
            self%s%ithr = ithr
            ! prep
            call self%s%prep4srch
            ! search
            do isample=1,self%s%nrefs_sub
                iref = s3D%srch_order_sub(isample,self%s%ithr) ! set the reference index
                if( s3D%state_exists(s3D%proj_space_state(iref)) )then
                    ! identify the top scoring in-plane angle
                    call pftcc_glob%gencorrs(iref, self%s%iptcl, inpl_corrs)
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
            ! count the number of nearest neighbors
            self%s%nnn = count(lnns)
            ! search
            do iref=1,self%s%nrefs
                iproj = iref - (self%s%prev_state - 1) * params_glob%nspace
                if( .not. lnns(iproj) ) cycle
                if( s3D%state_exists( s3D%proj_space_state(iref) ) )then
                    call pftcc_glob%gencorrs(iref, self%s%iptcl, inpl_corrs)
                    irot = greedy_sampling(eulprob_dist_switch(inpl_corrs), sorted_corrs, inds, s3D%smpl_inpl_ns)
                    call self%s%store_solution(iref, irot, inpl_corrs(irot))
                endif
            enddo
            ! we evaluate all refs
            self%s%nrefs_eval = self%s%nnn
            ! take care of the in-planes
            call self%s%inpl_srch_peaks
            ! prepare orientation
            call self%oris_assign
        else
            call build_glob%spproj_field%reject(self%s%iptcl)
        endif
    end subroutine srch_smpl_sub

    subroutine oris_assign_smpl_sub( self )
        class(strategy3D_smpl_sub), intent(inout) :: self
        call extract_peak_ori(self%s)
    end subroutine oris_assign_smpl_sub

    subroutine kill_smpl_sub( self )
        class(strategy3D_smpl_sub), intent(inout) :: self
        call self%s%kill
    end subroutine kill_smpl_sub

end module simple_strategy3D_smpl_sub
    