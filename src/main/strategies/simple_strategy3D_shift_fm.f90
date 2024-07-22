! concrete strategy3D: exhaustive top sampling
module simple_strategy3D_shift_fm
include 'simple_lib.f08'
use simple_strategy3D_alloc
use simple_strategy3D_utils
use simple_parameters,       only: params_glob
use simple_builder,          only: build_glob
use simple_strategy3D,       only: strategy3D
use simple_strategy3D_srch,  only: strategy3D_srch, strategy3D_spec
use simple_polarft_corrcalc, only: pftcc_glob
implicit none

public :: strategy3D_shift_fm
private
#include "simple_local_flags.inc"

type, extends(strategy3D) :: strategy3D_shift_fm
    type(strategy3D_srch) :: s
    type(strategy3D_spec) :: spec
contains
    procedure :: new         => new_shift_fm
    procedure :: srch        => srch_shift_fm
    procedure :: oris_assign => oris_assign_shift_fm
    procedure :: kill        => kill_shift_fm
end type strategy3D_shift_fm

contains

    subroutine new_shift_fm( self, spec )
        class(strategy3D_shift_fm), intent(inout) :: self
        class(strategy3D_spec), intent(inout) :: spec
        call self%s%new(spec)
        self%spec = spec
    end subroutine new_shift_fm

    subroutine srch_shift_fm( self, ithr )
        class(strategy3D_shift_fm), intent(inout) :: self
        integer,                    intent(in)    :: ithr
        integer   :: iref, isample, irot, nfound
        integer   :: inds(self%s%nrefs), rot_inds(self%s%nrefs)
        real      :: inpl_abscorrs(pftcc_glob%get_pftsz()), abscorrs(self%s%nrefs)
        real      :: offset(2), score
        logical   :: found
        ! execute search
        if( build_glob%spproj_field%get_state(self%s%iptcl) > 0 )then
            ! set thread index
            self%s%ithr = ithr
            ! prep
            call self%s%prep4srch
            ! search
            do isample=1,self%s%nrefs
                iref = s3D%srch_order(isample,self%s%ithr)
                if( s3D%state_exists(s3D%proj_space_state(iref)) )then
                    call pftcc_glob%gencorrs_mag(iref, self%s%iptcl, inpl_abscorrs, kweight=.true.)
                    irot           = maxloc(inpl_abscorrs,dim=1)
                    rot_inds(iref) = irot
                    abscorrs(iref) = inpl_abscorrs(irot)
                    inds(iref)     = iref
                else
                    abscorrs(iref) = 0.
                    inds(iref)     = 0
                    rot_inds(iref) = 0
                endif
            end do
            call hpsort(abscorrs, inds)
            nfound = 0
            do isample = floor(0.95*real(self%s%nrefs)),self%s%nrefs
                iref = inds(isample)
                irot = rot_inds(iref)
                call self%s%fm_shsrch_obj%minimize(iref, self%s%iptcl, found, irot, score, offset )
                if( found ) nfound = nfound + 1
                call self%s%store_solution(iref, irot, score, sh=offset)
            enddo
            self%s%nrefs_eval = self%s%nrefs
            call self%oris_assign
            if( params_glob%part == 1 )print *,self%s%iptcl, nfound
        else
            call build_glob%spproj_field%reject(self%s%iptcl)
        endif
    end subroutine srch_shift_fm

    ! in-plane search only version
    ! subroutine srch_shift_fm( self, ithr )
    !     class(strategy3D_shift_fm), intent(inout) :: self
    !     integer,                    intent(in)    :: ithr
    !     integer   :: iref, isample, irot
    !     integer   :: inds(self%s%nrefs), rot_inds(self%s%nrefs)
    !     real      :: inpl_abscorrs(pftcc_glob%get_pftsz()), abscorrs(self%s%nrefs)
    !     real      :: offset(2), score
    !     ! execute search
    !     if( build_glob%spproj_field%get_state(self%s%iptcl) > 0 )then
    !         ! set thread index
    !         self%s%ithr = ithr
    !         ! prep
    !         call self%s%prep4srch
    !         ! search
    !         call pftcc_glob%gencorrs_abs(self%s%prev_ref, self%s%iptcl, inpl_abscorrs, kweight=.true.)
    !         irot = maxloc(inpl_abscorrs,dim=1)
    !         call self%s%fm_shsrch_obj%minimize(self%s%prev_ref, self%s%iptcl, irot, score, offset )
    !         call self%s%store_solution(self%s%prev_ref, irot, score, sh=offset)
    !         self%s%nrefs_eval = self%s%nrefs
    !         call self%oris_assign
    !     else
    !         call build_glob%spproj_field%reject(self%s%iptcl)
    !     endif
    ! end subroutine srch_shift_fm

    subroutine oris_assign_shift_fm( self )
        class(strategy3D_shift_fm), intent(inout) :: self
        call extract_peak_ori(self%s)
    end subroutine oris_assign_shift_fm

    subroutine kill_shift_fm( self )
        class(strategy3D_shift_fm), intent(inout) :: self
        call self%s%kill
    end subroutine kill_shift_fm

end module simple_strategy3D_shift_fm
        