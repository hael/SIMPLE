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
        integer :: iref, locs(self%s%nrefs), inds(self%s%nrots)
        real    :: inpl_corrs(self%s%nrots), ref_corrs(self%s%nrefs), sorted_corrs(self%s%nrots)
        ! execute search
        if( build_glob%spproj_field%get_state(self%s%iptcl) > 0 )then
            ! init threaded search arrays
            self%s%ithr = ithr
            call prep_strategy3D_thread(ithr)
            do iref=1,self%s%nrefs
                if( s3D%state_exists( s3D%proj_space_state(iref) ) )then
                    ! identify the top scoring in-plane angle
                    call pftcc_glob%gencorrs(iref, self%s%iptcl, inpl_corrs)
                    locs(iref)      = reverse_multinomal(inpl_corrs, sorted_corrs, inds, s3D%smpl_inpl_ns, params_glob%l_reg_uni)
                    ref_corrs(iref) = inpl_corrs(locs(iref))
                endif
            enddo
            iref = reverse_multinomal(ref_corrs, s3D%smpl_refs_ns, params_glob%l_reg_uni)
            call assign_ori(self%s, iref, locs(iref), ref_corrs(iref))
            self%s%nrefs_eval = self%s%nrefs
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
    