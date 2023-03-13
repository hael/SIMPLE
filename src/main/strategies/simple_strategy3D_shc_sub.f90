! concrete strategy3D: shc refinement
module simple_strategy3D_shc_sub
include 'simple_lib.f08'
use simple_strategy3D_alloc  ! use all in there
use simple_strategy3D_utils  ! use all in there
use simple_strategy3D,       only: strategy3D
use simple_strategy3D_srch,  only: strategy3D_srch, strategy3D_spec
use simple_parameters,       only: params_glob
use simple_builder,          only: build_glob
use simple_polarft_corrcalc, only: pftcc_glob
implicit none

public :: strategy3D_shc_sub
private

#include "simple_local_flags.inc"

type, extends(strategy3D) :: strategy3D_shc_sub
    type(strategy3D_srch) :: s
    type(strategy3D_spec) :: spec
contains
    procedure :: new         => new_shc_sub
    procedure :: srch        => srch_shc_sub
    procedure :: kill        => kill_shc_sub
    procedure :: oris_assign => oris_assign_shc_sub
end type strategy3D_shc_sub

contains

    subroutine new_shc_sub( self, spec )
        class(strategy3D_shc_sub), intent(inout) :: self
        class(strategy3D_spec),    intent(inout) :: spec
        call self%s%new(spec)
        self%spec = spec
    end subroutine new_shc_sub

    subroutine srch_shc_sub( self, ithr )
        class(strategy3D_shc_sub), intent(inout) :: self
        integer,                   intent(in)    :: ithr
        integer   :: iref, isample, loc(1), iproj, ipeak
        real      :: inpl_corrs(self%s%nrots)
        logical   :: lnns(params_glob%nspace)
        type(ori) :: o
        if( build_glob%spproj_field%get_state(self%s%iptcl) > 0 )then
            ! set thread index
            self%s%ithr = ithr
            ! prep
            call self%s%prep4srch
            ! initialize, ctd
            self%s%nbetter    = 0
            self%s%nrefs_eval = 0
            ! search
            do isample=1,self%s%nrefs_sub
                iref = s3D%srch_order_sub(self%s%ithr,isample) ! set the reference index
                if( s3D%state_exists(s3D%proj_space_state(iref)) )then
                    ! identify the top scoring in-plane angle
                    call pftcc_glob%gencorrs(iref, self%s%iptcl, inpl_corrs)
                    loc = maxloc(inpl_corrs)
                    call self%s%store_solution(iref, loc(1), inpl_corrs(loc(1)))
                    ! update nbetter to keep track of how many improving solutions we have identified
                    if( inpl_corrs(loc(1)) > self%s%prev_corr ) self%s%nbetter = self%s%nbetter + 1
                    ! keep track of how many references we are evaluating
                    self%s%nrefs_eval = self%s%nrefs_eval + 1
                    ! exit condition
                    if( self%s%nbetter > 0 .and. self%s%nrefs_eval >= self%s%npeaks ) exit
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
            do iproj=1,params_glob%nspace
                if( .not. lnns(iproj) ) cycle
                iref = (self%s%prev_state - 1) * params_glob%nspace + iproj
                if( s3D%state_exists(s3D%proj_space_state(iref)) )then
                    ! identify the top scoring in-plane angle
                    call pftcc_glob%gencorrs(iref, self%s%iptcl, inpl_corrs)
                    loc = maxloc(inpl_corrs)
                    call self%s%store_solution(iref, loc(1), inpl_corrs(loc(1)))
                endif
            end do
            ! we evaluate all refs
            self%s%nrefs_eval = self%s%nnn
            ! take care of the in-planes
            call self%s%inpl_srch_peaks
            ! prepare orientation
            call self%oris_assign
        else
            call build_glob%spproj_field%reject(self%s%iptcl)
        endif
    end subroutine srch_shc_sub

    subroutine oris_assign_shc_sub( self )
        class(strategy3D_shc_sub), intent(inout) :: self
        call extract_peak_ori(self%s)
    end subroutine oris_assign_shc_sub

    subroutine kill_shc_sub( self )
        class(strategy3D_shc_sub), intent(inout) :: self
        call self%s%kill
    end subroutine kill_shc_sub

end module simple_strategy3D_shc_sub
