!@descr: retired program names that still run, mapped to the programs that replaced them
module simple_ui_legacy_names
use simple_defs, only: logfhandle
implicit none

public :: NLEGACY_PRG_NAMES, LEGACY_PRG_NAMES, CURRENT_PRG_NAMES, LEGACY_PRG_PRIVATE
public :: canonical_prg_name, report_legacy_prg_name
private

! Renamed programs keep their old names as retired names, so that scripts,
! logged command lines and job descriptions written before a rename keep
! running. Retired names are not registered in the UI: prg=list and the JSON
! interface show current names only. nice_lite migration 0006 carries the
! same table for stored NICE jobs.
! 2026-10-03: abinitio2D/abinitio3D became solve2D/solve3D (their variants with
! them, abinitio2D_stream became pool2D); cluster2D became refine2D, the 2D
! counterpart of refine3D. LEGACY_PRG_PRIVATE marks simple_private_exec programs,
! which have no UI registration.
integer,          parameter :: NLEGACY_PRG_NAMES = 12
character(len=*), parameter :: LEGACY_PRG_NAMES(NLEGACY_PRG_NAMES) = [character(len=26) :: &
    &'abinitio2D',        'abinitio2D_chunks', 'abinitio2D_stream',          &
    &'abinitio3D',        'abinitio3D_cavgs',  'abinitio3D_addon',           &
    &'abinitio3D_addon_snapshots',             'abinitio3D_nano',            &
    &'abinitio3D_stream', 'cluster2D',         'cluster2D_distr',            &
    &'cluster2D_nano']
character(len=*), parameter :: CURRENT_PRG_NAMES(NLEGACY_PRG_NAMES) = [character(len=23) :: &
    &'solve2D',           'solve2D_chunks',    'pool2D',                     &
    &'solve3D',           'solve3D_cavgs',     'solve3D_addon',              &
    &'solve3D_addon_snapshots',                'solve3D_nano',               &
    &'solve3D_stream',    'refine2D',          'refine2D_distr',             &
    &'refine2D_nano']
logical,          parameter :: LEGACY_PRG_PRIVATE(NLEGACY_PRG_NAMES) = [ &
    &.false.,             .false.,             .false.,                      &
    &.false.,             .false.,             .false.,                      &
    &.false.,                                  .false.,                      &
    &.false.,             .true.,              .true.,                       &
    &.false.]

contains

    !> current name of a retired program; any other name is returned unchanged
    function canonical_prg_name( name ) result( canon )
        character(len=*), intent(in)  :: name
        character(len=:), allocatable :: canon
        integer :: i
        canon = trim(adjustl(name))
        do i = 1, NLEGACY_PRG_NAMES
            if( canon == trim(LEGACY_PRG_NAMES(i)) )then
                canon = trim(CURRENT_PRG_NAMES(i))
                return
            endif
        end do
    end function canonical_prg_name

    !> one line on the log when a retired program name was used
    subroutine report_legacy_prg_name( name )
        character(len=*), intent(in)  :: name
        character(len=:), allocatable :: canon
        canon = canonical_prg_name(name)
        if( canon /= trim(adjustl(name)) )then
            write(logfhandle,'(a)') '>>> '//trim(adjustl(name))//' was renamed '//canon//'; running '//canon
        endif
    end subroutine report_legacy_prg_name

end module simple_ui_legacy_names
