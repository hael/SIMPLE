! concrete strategy3D: refinement
module simple_strategy3D_shc_smpl
include 'simple_lib.f08'
use simple_strategy3D_alloc
use simple_strategy3D_utils
use simple_parameters,       only: params_glob
use simple_builder,          only: build_glob
use simple_strategy3D,       only: strategy3D
use simple_strategy3D_srch,  only: strategy3D_srch, strategy3D_spec
use simple_polarft_corrcalc, only: pftcc_glob
implicit none

public :: strategy3D_shc_smpl
private
#include "simple_local_flags.inc"

type, extends(strategy3D) :: strategy3D_shc_smpl
    type(strategy3D_srch) :: s
    type(strategy3D_spec) :: spec
contains
    procedure          :: new         => new_shc_smpl
    procedure          :: srch        => srch_shc_smpl
    procedure          :: oris_assign => oris_assign_shc_smpl
    procedure          :: kill        => kill_shc_smpl
end type strategy3D_shc_smpl

contains

    subroutine new_shc_smpl( self, spec )
        class(strategy3D_shc_smpl), intent(inout) :: self
        class(strategy3D_spec),     intent(inout) :: spec
        call self%s%new(spec)
        self%spec = spec
    end subroutine new_shc_smpl

    subroutine srch_shc_smpl( self, ithr )
        use simple_eul_prob_tab, only: eulprob_dist_switch
        class(strategy3D_shc_smpl), intent(inout) :: self
        integer,                    intent(in)    :: ithr
        integer :: iref, isample, loc(1), inds(self%s%nrots)
        real    :: inpl_corrs(self%s%nrots), sorted_corrs(self%s%nrots)
        ! execute search
        if( build_glob%spproj_field%get_state(self%s%iptcl) > 0 )then
            ! set thread index
            self%s%ithr = ithr
            ! prep
            call self%s%prep4srch
            ! initialize, ctd
            self%s%nbetter    = 0
            self%s%nrefs_eval = 0
            do isample=1,self%s%nrefs
                iref = s3D%srch_order(self%s%ithr,isample)  ! set the stochastic reference index
                if( s3D%state_exists( s3D%proj_space_state(iref) ) )then
                    ! identify the top scoring in-plane angle
                    call pftcc_glob%gencorrs(iref, self%s%iptcl, inpl_corrs)
                    loc = greedy_sampling(eulprob_dist_switch(inpl_corrs), sorted_corrs, inds, s3D%smpl_inpl_ns)
                    call self%s%store_solution(iref, loc(1), inpl_corrs(loc(1)))
                    ! update nbetter to keep track of how many improving solutions we have identified
                    if( inpl_corrs(loc(1)) > self%s%prev_corr ) self%s%nbetter = self%s%nbetter + 1
                    ! keep track of how many references we are evaluating
                    self%s%nrefs_eval = self%s%nrefs_eval + 1
                end if
                ! exit condition
                if( self%s%nbetter > 0 ) exit
            end do
            call self%s%inpl_srch ! search shifts
            ! prepare weights and orientations
            call self%oris_assign
        else
            call build_glob%spproj_field%reject(self%s%iptcl)
        endif
    end subroutine srch_shc_smpl

    subroutine oris_assign_shc_smpl( self )
        class(strategy3D_shc_smpl), intent(inout) :: self
        call extract_peak_ori(self%s)
    end subroutine oris_assign_shc_smpl

    subroutine kill_shc_smpl( self )
        class(strategy3D_shc_smpl), intent(inout) :: self
        call self%s%kill
    end subroutine kill_shc_smpl

end module simple_strategy3D_shc_smpl
