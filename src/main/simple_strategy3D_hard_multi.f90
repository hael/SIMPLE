! concrete strategy3D: probabilistic multi-state refinement with hard orientation assignment
module simple_strategy3D_hard_multi
include 'simple_lib.f08'
use simple_strategy3D_alloc
use simple_strategy3D_utils
use simple_parameters,       only: params_glob
use simple_builder,          only: build_glob
use simple_strategy3D,       only: strategy3D
use simple_strategy3D_srch,  only: strategy3D_srch, strategy3D_spec
use simple_polarft_corrcalc, only: pftcc_glob
implicit none

public :: strategy3D_hard_multi
private

#include "simple_local_flags.inc"

type, extends(strategy3D) :: strategy3D_hard_multi
    type(strategy3D_srch) :: s
    type(strategy3D_spec) :: spec
contains
    procedure          :: new         => new_hard_multi
    procedure          :: srch        => srch_hard_multi
    procedure          :: oris_assign => oris_assign_hard_multi
    procedure          :: kill        => kill_hard_multi
end type strategy3D_hard_multi

contains

    subroutine new_hard_multi( self, spec, npeaks )
        class(strategy3D_hard_multi), intent(inout) :: self
        class(strategy3D_spec),  intent(inout) :: spec
        integer,                 intent(in)    :: npeaks
        call self%s%new(spec, npeaks)
        self%spec = spec
    end subroutine new_hard_multi

    subroutine srch_hard_multi( self, ithr )
        class(strategy3D_hard_multi), intent(inout) :: self
        integer,                      intent(in)    :: ithr
        integer :: iref,isample,nrefs
        real    :: inpl_corrs(self%s%nrots)
        real    :: inpl_corrs_mir(self%s%nrots)
        ! execute search
        if( build_glob%spproj_field%get_state(self%s%iptcl) > 0 )then
            ! set thread index
            self%s%ithr = ithr
            if( self%s%neigh )then
                ! initialize
                call self%s%prep4srch(build_glob%nnmat)
                nrefs = self%s%nnnrefs
            else
                ! initialize
                call self%s%prep4srch()
                nrefs = self%s%nrefs
            endif
            ! initialize, ctd
            self%s%nrefs_eval = 0
            do isample=1,nrefs
                iref = s3D%srch_order(self%s%ithr,isample) ! set the stochastic reference index
                call per_ref_srch                               ! actual search
            end do
            ! in this mode mode, we evaluate all refs
            self%s%nrefs_eval = nrefs
            call sort_corrs(self%s)  ! sort in correlation projection direction space
            call self%s%inpl_srch ! search shifts
            ! prepare weights and orientations
            call self%oris_assign
        else
            call build_glob%spproj_field%reject(self%s%iptcl)
        endif
        DebugPrint  '>>> strategy3D_hard_multi :: FINISHED STOCHASTIC SEARCH'

        contains

            subroutine per_ref_srch
                integer :: loc(MAXNINPLPEAKS), loc_mir(MAXNINPLPEAKS)
                integer :: iref_mir
                if( s3D%state_exists( s3D%proj_space_state(iref) ) )then
                    if (mir_projns) then
                        if (s3D%proj_space_corrs_calcd(self%s%ithr,iref)) then
                            s3D%proj_space_corrs_srchd(self%s%ithr,iref) = .true.
                        else
                            call pftcc_glob%gencorrs_mir(iref, self%s%iptcl, inpl_corrs, inpl_corrs_mir)
                            ! identify the MAXNINPLPEAKS top scoring in-planes
                            loc      = maxnloc(inpl_corrs,     MAXNINPLPEAKS)
                            loc_mir  = maxnloc(inpl_corrs_mir, MAXNINPLPEAKS)
                            iref_mir = s3D%proj_mirror_idx(iref)
                            call self%s%store_solution(iref,     loc,     inpl_corrs(loc),         .true. )
                            call self%s%store_solution(iref_mir, loc_mir, inpl_corrs_mir(loc_mir), .false.)
                        end if
                    else
                        ! calculate in-plane correlations
                        call pftcc_glob%gencorrs(iref, self%s%iptcl, inpl_corrs)
                        ! identify the MAXNINPLPEAKS top scoring in-planes
                        loc = maxnloc(inpl_corrs, MAXNINPLPEAKS)
                        ! stash
                        call self%s%store_solution(iref, loc, inpl_corrs(loc), .true.)
                    end if
                endif
            end subroutine per_ref_srch

    end subroutine srch_hard_multi

    subroutine oris_assign_hard_multi( self )
        use simple_ori,  only: ori
        class(strategy3D_hard_multi), intent(inout) :: self
        type(ori) :: osym
        real      :: dist_inpl, euldist, updatecnt
        integer   :: best_loc(1)
        ! extract peak info
        updatecnt = build_glob%spproj_field%get(self%s%iptcl, 'updatecnt')
        call prob_select_peak( self%s, params_glob%tau, updatecnt )
        best_loc(1) = 1 ! by definition
        ! angular distances
        call build_glob%pgrpsyms%sym_dists( build_glob%spproj_field%get_ori(self%s%iptcl),&
            &s3D%o_peaks(self%s%iptcl)%get_ori(best_loc(1)), osym, euldist, dist_inpl )
        ! generate convergence stats
        call convergence_stats_multi(self%s, best_loc, euldist)
        ! set the distances before we update the orientation
        if( build_glob%spproj_field%isthere(self%s%iptcl,'dist') )then
            call build_glob%spproj_field%set(self%s%iptcl, 'dist', 0.5*euldist + 0.5*build_glob%spproj_field%get(self%s%iptcl,'dist'))
        else
            call build_glob%spproj_field%set(self%s%iptcl, 'dist', euldist)
        endif
        call build_glob%spproj_field%set(self%s%iptcl, 'dist_inpl', dist_inpl)
        ! all the other stuff
        call build_glob%spproj_field%set_euler(self%s%iptcl, s3D%o_peaks(self%s%iptcl)%get_euler(best_loc(1)))
        call build_glob%spproj_field%set_shift(self%s%iptcl, s3D%o_peaks(self%s%iptcl)%get_2Dshift(best_loc(1)))
        call build_glob%spproj_field%set(self%s%iptcl, 'state',     s3D%o_peaks(self%s%iptcl)%get(best_loc(1),'state'))
        call build_glob%spproj_field%set(self%s%iptcl, 'frac',      100.)
        call build_glob%spproj_field%set(self%s%iptcl, 'corr',      s3D%o_peaks(self%s%iptcl)%get(best_loc(1),'corr'))
        call build_glob%spproj_field%set(self%s%iptcl, 'specscore', self%s%specscore)
        call build_glob%spproj_field%set(self%s%iptcl, 'ow',        1.0)
        call build_glob%spproj_field%set(self%s%iptcl, 'proj',      s3D%o_peaks(self%s%iptcl)%get(best_loc(1),'proj'))
        call build_glob%spproj_field%set(self%s%iptcl, 'spread',    0.0)
        call build_glob%spproj_field%set(self%s%iptcl, 'npeaks',    1.0)
        DebugPrint   '>>> strategy3D_hard_multi :: EXECUTED oris_assign_hard_multi'
    end subroutine oris_assign_hard_multi

    subroutine kill_hard_multi( self )
        class(strategy3D_hard_multi),   intent(inout) :: self
        call self%s%kill
    end subroutine kill_hard_multi

end module simple_strategy3D_hard_multi
