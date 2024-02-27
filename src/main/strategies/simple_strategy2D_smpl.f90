! concrete strategy3D: stochastic top sampling
module simple_strategy2D_smpl
include 'simple_lib.f08'
use simple_strategy2D_alloc  ! singleton
use simple_strategy2D,       only: strategy2D
use simple_strategy2D_srch,  only: strategy2D_spec
use simple_builder,          only: build_glob
use simple_polarft_corrcalc, only: pftcc_glob
use simple_parameters,       only: params_glob
implicit none

public :: strategy2D_smpl
private
#include "simple_local_flags.inc"

type, extends(strategy2D) :: strategy2D_smpl
contains
    procedure :: new  => new_smpl
    procedure :: srch => srch_smpl
    procedure :: kill => kill_smpl
end type strategy2D_smpl

contains

    subroutine new_smpl( self, spec )
        class(strategy2D_smpl), intent(inout) :: self
        class(strategy2D_spec), intent(inout) :: spec
        call self%s%new(spec)
        self%spec = spec
    end subroutine new_smpl

    subroutine srch_smpl( self )
        class(strategy2D_smpl), intent(inout) :: self
        integer :: iref, locs(self%s%nrefs), inds(self%s%nrots), irot
        real    :: inpl_corrs(self%s%nrots), ref_corrs(self%s%nrefs), sorted_corrs(self%s%nrots), corr
        if( build_glob%spproj_field%get_state(self%s%iptcl) > 0 )then
            call self%s%prep4srch
            ! search

            ! ref_corrs HOW TO INIT???
            ! s2D%smpl_inpl_ns HOW TO SET???
            ! s2D%smpl_refs_ns HOW TO SET???

            do iref=1,self%s%nrefs
                if( s2D%cls_pops(iref) == 0 )cycle      
                call pftcc_glob%gencorrs(iref, self%s%iptcl, inpl_corrs)
                irot = reverse_multinomal(inpl_corrs, sorted_corrs, inds, s2D%smpl_inpl_ns)
                locs(iref)      = irot
                ref_corrs(iref) = inpl_corrs(irot)
            enddo
            self%s%nrefs_eval = self%s%nrefs
            iref = reverse_multinomal(ref_corrs, s2D%smpl_refs_ns)
            irot = locs(iref)
            corr = ref_corrs(iref)
            self%s%best_class = iref
            self%s%best_rot   = irot
            self%s%best_corr  = corr
            call self%s%inpl_srch
            call self%s%store_solution
        else
            call build_glob%spproj_field%reject(self%s%iptcl)
        endif
    end subroutine srch_smpl

    subroutine kill_smpl( self )
        class(strategy2D_smpl), intent(inout) :: self
        call self%s%kill
    end subroutine kill_smpl

end module simple_strategy2D_smpl