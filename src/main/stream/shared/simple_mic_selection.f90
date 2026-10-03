!@descr: micrograph rejection on CTF resolution, ice fraction and astigmatism thresholds
!==============================================================================
! MODULE: simple_mic_selection
!
! PURPOSE:
!   The one micrograph threshold rule. Until now it was written out in the
!   `selection` commander (simple_commanders_project_core), twice in stream
!   p01 and twice in stream p04, with different comparisons: p04 rejected at
!   threshold-0.001, `selection` ignored astigmatism, p01 tested whether the
!   key existed anywhere in the segment. The rule here:
!
!   - a micrograph with state > 0 is rejected (state = 0) when any supplied
!     threshold is strictly exceeded by its own value for that key;
!   - a threshold that is not supplied is not applied;
!   - a micrograph without the key is not rejected on it;
!   - already rejected micrographs are left alone and not counted.
!
!   Rejection is one-way: raising a threshold later does not restore
!   micrographs, as in the code it replaces.
!
!   reject_mics_without_particles is the rule applied after picking and
!   extraction (no particles, or no readable box file), written out twice in
!   stream p03 until now.
!
! HOME:
!   In src/main/stream/shared for now; it belongs beside the micrograph
!   segment (oris or sp_project) once the `selection` commander calls it.
!
! TESTS:
!   simple_mic_selection_tester
!==============================================================================
module simple_mic_selection
use simple_string, only: string
use simple_fileio, only: file_exists
use simple_oris,   only: oris
implicit none

public :: reject_mics_by_thresholds, reject_mics_without_particles
private

contains

    !> Sets state 0 on every micrograph of @p os_mic with state > 0 that exceeds a supplied threshold;
    !! @p nrejected is the number newly rejected.
    subroutine reject_mics_by_thresholds( os_mic, nrejected, ctfres, icefrac, astig )
        class(oris),    intent(inout) :: os_mic
        integer,        intent(out)   :: nrejected
        real, optional, intent(in)    :: ctfres, icefrac, astig
        integer :: imic
        nrejected = 0
        do imic = 1,os_mic%get_noris()
            if( os_mic%get_state(imic) <= 0 ) cycle
            if( exceeds('ctfres', ctfres) .or. exceeds('icefrac', icefrac) .or. exceeds('astig', astig) )then
                call os_mic%set_state(imic, 0)
                nrejected = nrejected + 1
            endif
        enddo

    contains

        logical function exceeds( key, threshold )
            character(len=*), intent(in) :: key
            real, optional,   intent(in) :: threshold
            exceeds = .false.
            if( .not. present(threshold) ) return
            if( .not. os_mic%isthere(imic, key) ) return
            exceeds = os_mic%get(imic, key) > threshold
        end function exceeds

    end subroutine reject_mics_by_thresholds

    !> After picking and extraction: sets state 0 on every micrograph of @p os_mic with state > 0
    !! that has no particles or no readable box file; rows are kept so micrograph/stack indices hold.
    !! @p nrejected is the number newly rejected.
    subroutine reject_mics_without_particles( os_mic, nrejected )
        class(oris), intent(inout) :: os_mic
        integer,     intent(out)   :: nrejected
        type(string) :: boxfile
        integer      :: imic
        logical      :: l_reject
        nrejected = 0
        do imic = 1,os_mic%get_noris()
            if( os_mic%get_state(imic) <= 0 ) cycle
            l_reject = .not. os_mic%isthere(imic, 'nptcls')
            if( .not. l_reject ) l_reject = os_mic%get(imic, 'nptcls') <= 0.
            if( .not. l_reject ) l_reject = .not. os_mic%isthere(imic, 'boxfile')
            if( .not. l_reject )then
                boxfile  = os_mic%get_str(imic, 'boxfile')
                l_reject = boxfile%strlen() == 0
                if( .not. l_reject ) l_reject = .not. file_exists(boxfile)
            endif
            if( l_reject )then
                call os_mic%set_state(imic, 0)
                nrejected = nrejected + 1
            endif
        enddo
    end subroutine reject_mics_without_particles

end module simple_mic_selection
