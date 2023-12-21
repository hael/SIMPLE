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
        integer :: iref, loc, inpl_inds(self%s%nrots), irnd, ref_inds(self%s%nrefs), j
        real    :: inpl_corrs(self%s%nrots), rnd_num, ref_corrs(self%s%nrefs)
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
                    inpl_inds = (/(j,j=1,self%s%nrots)/)
                    call hpsort(inpl_corrs, inpl_inds)
                    call random_number(rnd_num)
                    irnd = 1 + floor(params_glob%reg_nrots * rnd_num)
                    loc  = inpl_inds(irnd)
                    ref_corrs(iref) = inpl_corrs(irnd)
                endif
            enddo
            ref_inds = (/(iref,iref=1,self%s%nrefs)/)
            call hpsort(ref_corrs, ref_inds)
            call random_number(rnd_num)
            irnd = 1 + floor(params_glob%reg_nrots * rnd_num)
            iref = ref_inds(irnd)
            self%s%nrefs_eval = self%s%nrefs
            call assign_ori(self%s, iref, loc, ref_corrs(irnd))
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
