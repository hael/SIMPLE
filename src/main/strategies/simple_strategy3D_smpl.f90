! concrete strategy3D: stochastic top sampling
module simple_strategy3D_smpl
include 'simple_lib.f08'
use simple_strategy3D_alloc
use simple_strategy3D_utils
use simple_parameters,       only: params_glob
use simple_builder,          only: build_glob
use simple_strategy3D,       only: strategy3D
use simple_strategy3D_srch,  only: strategy3D_srch, strategy3D_spec
use simple_polarft_corrcalc, only: pftcc_glob
use simple_eul_prob_tab,     only: eulprob_dist_switch, angle_sampling
implicit none

public :: strategy3D_smpl
private
#include "simple_local_flags.inc"

type, extends(strategy3D) :: strategy3D_smpl
    type(strategy3D_srch) :: s
    type(strategy3D_spec) :: spec
contains
    procedure :: new         => new_smpl
    procedure :: srch        => srch_smpl
    procedure :: oris_assign => oris_assign_smpl
    procedure :: kill        => kill_smpl
end type strategy3D_smpl

contains

    subroutine new_smpl( self, spec )
        class(strategy3D_smpl), intent(inout) :: self
        class(strategy3D_spec), intent(inout) :: spec
        call self%s%new(spec)
        self%spec = spec
    end subroutine new_smpl

    subroutine srch_smpl( self, ithr )
        class(strategy3D_smpl), intent(inout) :: self
        integer,                intent(in)    :: ithr
        integer :: iref, locs(self%s%nrefs), inds(self%s%nrots), inds_ref(self%s%nrefs), irot
        real    :: inpl_corrs(self%s%nrots), ref_corrs(self%s%nrefs), sorted_corrs(self%s%nrots), sorted_ref(self%s%nrefs), corr
        ! execute search
        if( build_glob%spproj_field%get_state(self%s%iptcl) > 0 )then
            ! init threaded search arrays
            self%s%ithr = ithr
            call self%s%prep4srch
            ! search
            ref_corrs = TINY
            do iref=1,self%s%nprojs
                if( s3D%state_exists( s3D%proj_space_state(iref) ) )then
                    call pftcc_glob%gencorrs(iref, self%s%iptcl, inpl_corrs)
                    irot = angle_sampling(eulprob_dist_switch(inpl_corrs), sorted_corrs, inds,&
                          &s3D%smpl_inpl_athres(s3D%proj_space_state(iref)))
                    locs(iref)      = irot
                    ref_corrs(iref) = inpl_corrs(irot)
                endif
            enddo
            self%s%nrefs_eval = self%s%nrefs
            iref = angle_sampling(eulprob_dist_switch(ref_corrs), sorted_ref, inds_ref,&
                   &s3D%smpl_refs_athres(build_glob%spproj_field%get_state(self%s%iptcl)))
            irot = locs(iref)
            corr = ref_corrs(iref)
            if( self%s%doshift )then
                call self%s%inpl_srch(ref=iref)
                ! checking if shift search is good
                if( s3D%proj_space_inplinds(iref, self%s%ithr) < 1 )then
                    call assign_ori(self%s, iref, irot, corr, [0.,0.])
                else
                    irot = s3D%proj_space_inplinds(iref, self%s%ithr)
                    corr = s3D%proj_space_corrs(   iref, self%s%ithr)
                    call assign_ori(self%s, iref, irot, corr, s3D%proj_space_shift(:,iref,self%s%ithr))
                endif
            else
                call assign_ori(self%s, iref, irot, corr, [0.,0.])
            endif
        else
            call build_glob%spproj_field%reject(self%s%iptcl)
        endif
    end subroutine srch_smpl

    subroutine oris_assign_smpl( self )
        class(strategy3D_smpl), intent(inout) :: self
        call extract_peak_ori(self%s)
    end subroutine oris_assign_smpl

    subroutine kill_smpl( self )
        class(strategy3D_smpl), intent(inout) :: self
        call self%s%kill
    end subroutine kill_smpl

end module simple_strategy3D_smpl
    