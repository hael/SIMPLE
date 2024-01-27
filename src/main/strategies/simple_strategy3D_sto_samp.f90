! concrete strategy3D: stochastic top sampling
module simple_strategy3D_sto_samp
include 'simple_lib.f08'
use simple_strategy3D_alloc
use simple_strategy3D_utils
use simple_parameters,       only: params_glob
use simple_builder,          only: build_glob
use simple_strategy3D,       only: strategy3D
use simple_strategy3D_srch,  only: strategy3D_srch, strategy3D_spec
use simple_polarft_corrcalc, only: pftcc_glob
implicit none

public :: strategy3D_sto_samp
private
#include "simple_local_flags.inc"

type, extends(strategy3D) :: strategy3D_sto_samp
    type(strategy3D_srch) :: s
    type(strategy3D_spec) :: spec
contains
    procedure          :: new         => new_sto_samp
    procedure          :: srch        => srch_sto_samp
    procedure          :: oris_assign => oris_assign_sto_samp
    procedure          :: kill        => kill_sto_samp
end type strategy3D_sto_samp

contains

    subroutine new_sto_samp( self, spec )
        class(strategy3D_sto_samp), intent(inout) :: self
        class(strategy3D_spec),     intent(inout) :: spec
        call self%s%new(spec)
        self%spec = spec
    end subroutine new_sto_samp

    subroutine srch_sto_samp( self, ithr )
        class(strategy3D_sto_samp), intent(inout) :: self
        integer,                    intent(in)    :: ithr
        integer :: iref, locs(self%s%nrefs)
        real    :: inpl_corrs(self%s%nrots), ref_corrs(self%s%nrefs)
        ! execute search
        if( build_glob%spproj_field%get_state(self%s%iptcl) > 0 )then
            ! set thread index
            self%s%ithr = ithr
            ! prep
            call self%s%prep4srch
            do iref=1,self%s%nrefs
                if( s3D%state_exists( s3D%proj_space_state(iref) ) )then
                    ! identify the top scoring in-plane angle
                    call pftcc_glob%gencorrs(iref, self%s%iptcl, inpl_corrs)
                    locs(iref)      = reverse_multinomal(inpl_corrs, int(params_glob%reg_athres * self%s%nrots / 180.))
                    ref_corrs(iref) = inpl_corrs(locs(iref))
                endif
            enddo
            iref              = reverse_multinomal(ref_corrs, int(params_glob%reg_athres * self%s%nrefs / 180.))
            self%s%nrefs_eval = self%s%nrefs
            call assign_ori(self%s, iref, locs(iref), ref_corrs(iref))
        else
            call build_glob%spproj_field%reject(self%s%iptcl)
        endif
    end subroutine srch_sto_samp

    subroutine oris_assign_sto_samp( self )
        class(strategy3D_sto_samp), intent(inout) :: self
        call extract_peak_ori(self%s)
    end subroutine oris_assign_sto_samp

    subroutine kill_sto_samp( self )
        class(strategy3D_sto_samp),   intent(inout) :: self
        call self%s%kill
    end subroutine kill_sto_samp

end module simple_strategy3D_sto_samp
