! concrete strategy3D: probabilistic multi-state refinement
module simple_strategy3D_multi
include 'simple_lib.f08'
use simple_strategy3D_alloc
use simple_strategy3D_utils
use simple_parameters,       only: params_glob
use simple_builder,          only: build_glob
use simple_strategy3D,       only: strategy3D
use simple_strategy3D_srch,  only: strategy3D_srch, strategy3D_spec
use simple_polarft_corrcalc, only: pftcc_glob
implicit none

public :: strategy3D_multi
private

#include "simple_local_flags.inc"

type, extends(strategy3D) :: strategy3D_multi
    type(strategy3D_srch) :: s
    type(strategy3D_spec) :: spec
contains
    procedure          :: new         => new_multi
    procedure          :: srch        => srch_multi
    procedure          :: oris_assign => oris_assign_multi
    procedure          :: kill        => kill_multi
end type strategy3D_multi

contains

    subroutine new_multi( self, spec, npeaks )
        class(strategy3D_multi), intent(inout) :: self
        class(strategy3D_spec),  intent(inout) :: spec
        integer,                 intent(in)    :: npeaks
        call self%s%new( spec, npeaks )
        self%spec = spec
    end subroutine new_multi

    subroutine srch_multi( self, ithr )
        class(strategy3D_multi), intent(inout) :: self
        integer,                 intent(in)    :: ithr
        integer :: iref,isample,nrefs
        real    :: inpl_corrs(self%s%nrots)
        real    :: areal
        real    :: inpl_corrs_mir(self%s%nrots)
        ! execute search
        if( build_glob%spproj_field%get_state(self%s%iptcl) > 0 )then
            ! set thread index
            self%s%ithr = ithr
            if( self%s%neigh )then
                call self%s%prep4srch(build_glob%nnmat)
                nrefs = self%s%nnnrefs
            else
                ! initialize
                call self%s%prep4srch()
                nrefs = self%s%nrefs
            endif
            ! initialize, ctd
            self%s%nbetter    = 0
            self%s%nrefs_eval = 0
            do isample=1,nrefs
                iref = s3D%srch_order(self%s%ithr,isample)  ! set the stochastic reference index
                call per_ref_srch                           ! actual search
                if( self%s%nbetter >= self%s%npeaks ) exit  ! exit condition
            end do
            call sort_corrs(self%s)  ! sort in correlation projection direction space
            call self%s%inpl_srch    ! search shifts
            ! prepare weights and orientations
            call self%oris_assign
        else
            call build_glob%spproj_field%reject(self%s%iptcl)
        endif
        DebugPrint  '>>> STRATEGY3D_MULTI :: FINISHED STOCHASTIC SEARCH'

    contains

        subroutine per_ref_srch
            real    :: best_inpl_corr
            integer :: loc(MAXNINPLPEAKS)
            integer :: loc_mir(MAXNINPLPEAKS)
            integer :: iref_mir
            real    :: best_inpl_corr_mir
            if( s3D%state_exists( s3D%proj_space_state(iref) ) )then
                if (mir_projns) then
                    if (s3D%proj_space_corrs_calcd(self%s%ithr,iref)) then
                        best_inpl_corr = s3D%proj_space_corrs(self%s%ithr,iref,1)
                        s3D%proj_space_corrs_srchd(self%s%ithr,iref) = .true.
                    else
                        call pftcc_glob%gencorrs_mir(iref, self%s%iptcl, inpl_corrs, inpl_corrs_mir)
                        ! identify the MAXNINPLPEAKS top scoring in-planes
                        loc                = maxnloc(inpl_corrs,     MAXNINPLPEAKS)
                        loc_mir            = maxnloc(inpl_corrs_mir, MAXNINPLPEAKS)
                        best_inpl_corr     = inpl_corrs(loc(1))
                        best_inpl_corr_mir = inpl_corrs_mir(loc_mir(1))
                        iref_mir           = s3D%proj_mirror_idx(iref)
                        call self%s%store_solution(iref,     loc,     inpl_corrs(loc),         .true. )
                        call self%s%store_solution(iref_mir, loc_mir, inpl_corrs_mir(loc_mir), .false.)
                    end if
                else
                    ! identify the MAXNINPLPEAKS top scoring in-planes
                    call pftcc_glob%gencorrs(iref, self%s%iptcl, inpl_corrs)
                    loc                = maxnloc(inpl_corrs, MAXNINPLPEAKS)
                    best_inpl_corr     = inpl_corrs(loc(1))
                    call self%s%store_solution(iref, loc, inpl_corrs(loc), .true.)
                end if
                ! update nbetter to keep track of how many improving solutions we have identified
                if( self%s%npeaks == 1 )then
                    if( best_inpl_corr > self%s%prev_corr ) self%s%nbetter = self%s%nbetter + 1
                else
                    if( best_inpl_corr >= self%s%prev_corr ) self%s%nbetter = self%s%nbetter + 1
                endif
                ! keep track of how many references we are evaluating
                self%s%nrefs_eval = self%s%nrefs_eval + 1
            end if
        end subroutine per_ref_srch

    end subroutine srch_multi

    !>  \brief retrieves and preps npeaks orientations for reconstruction
    subroutine oris_assign_multi( self )
        use simple_ori,  only: ori
        class(strategy3D_multi), intent(inout) :: self
        type(ori) :: osym
        real      :: corrs(self%s%npeaks), ws(self%s%npeaks)
        real      :: wcorr, frac, ang_spread, dist_inpl, euldist
        integer   :: best_loc(1), neff_states, state
        logical   :: included(self%s%npeaks)
        ! extract peak info
        call extract_peaks(self%s, corrs, multistates=.true.)
        ! stochastic weights
        call corrs2softmax_weights(self%s, self%s%npeaks, corrs, params_glob%tau, ws, included, best_loc, wcorr )
        ! state reweighting
        call states_reweight(self%s, self%s%npeaks, ws, included, state, best_loc, wcorr)
        ! angular standard deviation
        ang_spread = estimate_ang_spread(self%s)
        ! angular distances
        call build_glob%pgrpsyms%sym_dists( build_glob%spproj_field%get_ori(self%s%iptcl),&
            &s3D%o_peaks(self%s%iptcl)%get_ori(best_loc(1)), osym, euldist, dist_inpl )
        ! generate convergence stats
        call convergence_stats_multi( self%s, best_loc, euldist )
        ! fraction of search space scanned
        neff_states = count(s3D%state_exists)
        if( self%s%neigh )then
            frac = 100.*real(self%s%nrefs_eval) / real(self%s%nnn * neff_states)
        else
            frac = 100.*real(self%s%nrefs_eval) / real(self%s%nprojs * neff_states)
        endif
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
        call build_glob%spproj_field%set(self%s%iptcl, 'state',     real(state))
        call build_glob%spproj_field%set(self%s%iptcl, 'frac',      frac)
        call build_glob%spproj_field%set(self%s%iptcl, 'corr',      wcorr)
        call build_glob%spproj_field%set(self%s%iptcl, 'specscore', self%s%specscore)
        call build_glob%spproj_field%set(self%s%iptcl, 'ow',        s3D%o_peaks(self%s%iptcl)%get(best_loc(1),'ow'))
        call build_glob%spproj_field%set(self%s%iptcl, 'proj',      s3D%o_peaks(self%s%iptcl)%get(best_loc(1),'proj'))
        call build_glob%spproj_field%set(self%s%iptcl, 'spread',    ang_spread)
        call build_glob%spproj_field%set(self%s%iptcl, 'npeaks',    real(self%s%npeaks_eff))
        DebugPrint   '>>> STRATEGY3D_MULTI :: EXECUTED ORIS_ASSIGN_MULTI'
    end subroutine oris_assign_multi

    subroutine kill_multi( self )
        class(strategy3D_multi),   intent(inout) :: self
        call self%s%kill
    end subroutine kill_multi

end module simple_strategy3D_multi
