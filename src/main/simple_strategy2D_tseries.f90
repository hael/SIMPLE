module simple_strategy2D_tseries
include 'simple_lib.f08'
use simple_strategy2D_alloc
use simple_strategy2D,       only: strategy2D
use simple_strategy2D_srch,  only: strategy2D_srch, strategy2D_spec
use simple_builder,          only: build_glob
use simple_polarft_corrcalc, only: pftcc_glob
implicit none

public :: strategy2D_tseries
private

#include "simple_local_flags.inc"

real, parameter :: hwidth_ang = 15. ! degrees, < 180.

type, extends(strategy2D) :: strategy2D_tseries
    type(strategy2D_srch) :: s
    type(strategy2D_spec) :: spec
contains
    procedure :: new  => new_tseries
    procedure :: srch => srch_tseries
    procedure :: kill => kill_tseries
end type strategy2D_tseries

contains

    subroutine new_tseries( self, spec )
        class(strategy2D_tseries), intent(inout) :: self
        class(strategy2D_spec),   intent(inout) :: spec
        call self%s%new( spec )
        self%spec = spec
    end subroutine new_tseries

    subroutine srch_tseries( self )
        class(strategy2D_tseries), intent(inout) :: self
        logical :: rotmask(self%s%nrots)
        integer :: iref,inpl_ind,hwidth_rot
        real    :: corrs(self%s%nrots),inpl_corr,corr
        if( build_glob%spproj_field%get_state(self%s%iptcl) > 0 )then
            call self%s%prep4srch
            corr    = -huge(corr)
            hwidth_rot = max(1,floor(real(self%s%nrots)*hwidth_ang/360.))
            rotmask    = .true.
            rotmask(hwidth_rot+1:self%s%nrots-hwidth_rot) = .false.
            ! previous best
            iref = self%s%prev_class
            call per_ref_srch
            ! first backwards
            do iref = self%s%prev_class-1,1,-1
                if( s2D%cls_pops(iref) > 0 )then
                    call per_ref_srch
                    exit
                endif
            enddo
            ! first forward
            do iref = self%s%prev_class+1,self%s%nrefs,1
                if( s2D%cls_pops(iref) > 0 )then
                    call per_ref_srch
                    exit
                endif
            enddo
            self%s%nrefs_eval = self%s%nrefs
            call self%s%inpl_srch ! shift only search
            call self%s%store_solution
        else
            call build_glob%spproj_field%reject(self%s%iptcl)
        endif
        contains

            subroutine per_ref_srch
                if( s2D%cls_pops(iref) == 0 )return
                call pftcc_glob%gencorrs(iref, self%s%iptcl, corrs)
                inpl_ind  = maxloc(corrs, dim=1, mask=rotmask)
                inpl_corr = corrs(inpl_ind)
                if( inpl_corr >= corr )then
                    corr              = inpl_corr
                    self%s%best_class = iref
                    self%s%best_corr  = inpl_corr
                    self%s%best_rot   = inpl_ind
                endif
            end subroutine per_ref_srch

    end subroutine srch_tseries

    subroutine kill_tseries( self )
        class(strategy2D_tseries), intent(inout) :: self
        call self%s%kill
    end subroutine kill_tseries

end module simple_strategy2D_tseries
