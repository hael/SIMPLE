!@descr: abinitio3D_addon validation against the base solution per state: FSC verdict, map correlation at base resolution, cohort cross-check and stage limits
! Per state: union vs base FSC (verdict on the FSC=0.143 shell move beyond SHELL_TOL), correlation in
! the mskdiam soft sphere up to the base resolution, docking below DOCK_CORR_FLOOR, optional cohort-only
! check (addon_diag), stage limits; kind=stage|state|cohort text records that read restores.
module simple_abinitio3D_addon_report
use simple_core_module_api
use simple_volpair_metrics, only: compare_volpair
use simple_dock_vols,       only: dock_vols
use simple_ori,             only: ori
use simple_oris,            only: oris
implicit none

public :: abinitio3D_addon_report, ADDON_REPORT_FNAME
private
#include "simple_local_flags.inc"

character(len=*), parameter :: ADDON_REPORT_FNAME = 'abinitio3D_addon_report.txt'
integer,          parameter :: SHELL_TOL       = 1     !< FSC=0.143 shells within which a state is UNCHANGED
real,             parameter :: FSC_CRIT        = 0.143 !< the resolution criterion the verdict uses
real,             parameter :: DOCK_CORR_FLOOR = 0.9   !< in-frame correlation below which the union map is docked
real,             parameter :: DOCK_HP         = 100.  !< docking band, high-pass (A)
real,             parameter :: DOCK_LP_MIN     = 8.    !< docking band, low-pass: the base resolution, no finer (A)
integer,          parameter :: VERDICT_LEN     = 12
character(len=*), parameter :: V_NONE = 'NOT_COMPARED', V_UP = 'IMPROVED', V_SAME = 'UNCHANGED', V_DOWN = 'REGRESSED'

type :: state_record
    integer :: pop_union = 0, pop_frozen = 0, dshell = 0
    logical :: l_fsc = .false., l_map = .false., l_docked = .false., l_cohort = .false.
    real    :: res05_base = 0., res0143_base = 0., res05_union = 0., res0143_union = 0., fsc_gain = 0.
    real    :: lp_corr = 0., corr = 0., dock_angle = 0., dock_shift = 0., dock_corr = 0.
    real    :: res05_cohort = 0., res0143_cohort = 0., corr_cohort = 0.
    character(len=VERDICT_LEN) :: verdict = V_NONE
end type state_record

type :: stage_record
    logical :: l_set = .false., l_base = .false.
    real    :: lp_addon = 0., lpstop_addon = 0., lp_base = 0., lpstop_base = 0.
end type stage_record

type :: abinitio3D_addon_report
    private
    integer :: nstates = 0, nstages = 0
    real    :: smpd = 0., mskdiam = 0.
    type(state_record), allocatable :: states(:)
    type(stage_record), allocatable :: stages(:)
  contains
    procedure :: new
    procedure :: set_stage
    procedure :: set_populations
    procedure :: compare_fsc
    procedure :: compare_maps
    procedure :: compare_cohort
    procedure :: get_verdict
    procedure :: get_dshell
    procedure :: get_corr
    procedure :: get_dock_angle
    procedure :: get_cohort_res0143
    procedure :: any_regressed
    procedure :: print
    procedure :: write
    procedure :: read
    procedure :: kill
    procedure, private :: check_state
    procedure, private :: corr_lp
    procedure, private :: dock
end type abinitio3D_addon_report

contains

    !> nstates states, stages 1..nstages, maps sampled at smpd (A) and masked
    !! at the base run's mskdiam (A)
    subroutine new( self, nstates, nstages, smpd, mskdiam )
        class(abinitio3D_addon_report), intent(inout) :: self
        integer,                        intent(in)    :: nstates, nstages
        real,                           intent(in)    :: smpd, mskdiam
        call self%kill
        self%nstates = nstates
        self%nstages = nstages
        self%smpd    = smpd
        self%mskdiam = mskdiam
        allocate(self%states(nstates), self%stages(nstages))
    end subroutine new

    !> the limits the add-on emitted at a stage and the base run's (a negative
    !! base limit means the base run did not run the stage)
    subroutine set_stage( self, istage, lp_addon, lpstop_addon, lp_base, lpstop_base )
        class(abinitio3D_addon_report), intent(inout) :: self
        integer,                        intent(in)    :: istage
        real,                           intent(in)    :: lp_addon, lpstop_addon, lp_base, lpstop_base
        if( istage < 1 .or. istage > self%nstages ) THROW_HARD('stage out of range; set_stage')
        self%stages(istage)%l_set        = .true.
        self%stages(istage)%l_base       = lp_base >= 0.
        self%stages(istage)%lp_addon     = lp_addon
        self%stages(istage)%lpstop_addon = lpstop_addon
        self%stages(istage)%lp_base      = lp_base
        self%stages(istage)%lpstop_base  = lpstop_base
    end subroutine set_stage

    subroutine set_populations( self, s, pop_union, pop_frozen )
        class(abinitio3D_addon_report), intent(inout) :: self
        integer,                        intent(in)    :: s, pop_union, pop_frozen
        call self%check_state(s)
        self%states(s)%pop_union  = pop_union
        self%states(s)%pop_frozen = pop_frozen
    end subroutine set_populations

    !> The union FSC against the base run's, both on the resolution axis res:
    !! resolutions, the move of the FSC=0.143 shell, the mean FSC gain up to the
    !! base run's FSC=0.143 shell and the verdict (curves of another length are
    !! not compared)
    subroutine compare_fsc( self, s, fsc_base, fsc_union, res )
        class(abinitio3D_addon_report), intent(inout) :: self
        integer,                        intent(in)    :: s
        real,                           intent(in)    :: fsc_base(:), fsc_union(:), res(:)
        integer :: kbase, kunion
        call self%check_state(s)
        if( size(fsc_base) /= size(res) .or. size(fsc_union) /= size(res) ) return
        associate( st => self%states(s) )
            call get_resolution(fsc_base,  res, st%res05_base,  st%res0143_base)
            call get_resolution(fsc_union, res, st%res05_union, st%res0143_union)
            kbase     = get_find_at_crit(fsc_base,  FSC_CRIT)
            kunion    = get_find_at_crit(fsc_union, FSC_CRIT)
            st%dshell   = kunion - kbase
            st%fsc_gain = sum(fsc_union(1:kbase) - fsc_base(1:kbase)) / real(kbase)
            if( st%dshell > SHELL_TOL )then
                st%verdict = V_UP
            else if( st%dshell < -SHELL_TOL )then
                st%verdict = V_DOWN
            else
                st%verdict = V_SAME
            endif
            st%l_fsc = .true.
        end associate
    end subroutine compare_fsc

    !> The union map against the base map, in the base frame, up to the base
    !! run's FSC=0.143 resolution (Nyquist before compare_fsc); docked when the
    !! correlation falls below DOCK_CORR_FLOOR
    subroutine compare_maps( self, s, base_fname, union_fname )
        class(abinitio3D_addon_report), intent(inout) :: self
        integer,                        intent(in)    :: s
        class(string),                  intent(in)    :: base_fname, union_fname
        real, allocatable :: fsc(:), res(:)
        logical :: ok
        call self%check_state(s)
        associate( st => self%states(s) )
            st%lp_corr = self%corr_lp(s)
            call compare_volpair(base_fname, union_fname, self%mskdiam, st%lp_corr, st%corr, fsc, res, ok)
            st%l_map = ok
        end associate
        if( self%states(s)%l_map .and. self%states(s)%corr < DOCK_CORR_FLOOR ) call self%dock(s, base_fname, union_fname)
    end subroutine compare_maps

    !> The cohort-only map against the base map: their FSC resolutions (the
    !! two come from disjoint particles) and their correlation up to the base
    !! run's FSC=0.143 resolution
    subroutine compare_cohort( self, s, base_fname, cohort_fname )
        class(abinitio3D_addon_report), intent(inout) :: self
        integer,                        intent(in)    :: s
        class(string),                  intent(in)    :: base_fname, cohort_fname
        real, allocatable :: fsc(:), res(:)
        logical :: ok
        call self%check_state(s)
        associate( st => self%states(s) )
            call compare_volpair(base_fname, cohort_fname, self%mskdiam, self%corr_lp(s), st%corr_cohort, fsc, res, ok)
            if( .not. ok ) return
            call get_resolution(fsc, res, st%res05_cohort, st%res0143_cohort)
            st%l_cohort = .true.
        end associate
    end subroutine compare_cohort

    function get_verdict( self, s ) result( verdict )
        class(abinitio3D_addon_report), intent(in) :: self
        integer,                        intent(in) :: s
        character(len=VERDICT_LEN) :: verdict
        call self%check_state(s)
        verdict = self%states(s)%verdict
    end function get_verdict

    integer function get_dshell( self, s )
        class(abinitio3D_addon_report), intent(in) :: self
        integer,                        intent(in) :: s
        call self%check_state(s)
        get_dshell = self%states(s)%dshell
    end function get_dshell

    !> the union map's correlation with the base map (0 when not compared)
    real function get_corr( self, s )
        class(abinitio3D_addon_report), intent(in) :: self
        integer,                        intent(in) :: s
        call self%check_state(s)
        get_corr = self%states(s)%corr
    end function get_corr

    !> the rotation (deg) that docked the union map onto the base map, -1 when
    !! the maps were not docked
    real function get_dock_angle( self, s )
        class(abinitio3D_addon_report), intent(in) :: self
        integer,                        intent(in) :: s
        call self%check_state(s)
        get_dock_angle = -1.
        if( self%states(s)%l_docked ) get_dock_angle = self%states(s)%dock_angle
    end function get_dock_angle

    !> the FSC=0.143 resolution of the cohort-only map against the base map
    !! (A; 0 when not compared)
    real function get_cohort_res0143( self, s )
        class(abinitio3D_addon_report), intent(in) :: self
        integer,                        intent(in) :: s
        call self%check_state(s)
        get_cohort_res0143 = self%states(s)%res0143_cohort
    end function get_cohort_res0143

    logical function any_regressed( self )
        class(abinitio3D_addon_report), intent(in) :: self
        integer :: s
        any_regressed = .false.
        do s = 1, self%nstates
            if( self%states(s)%verdict == V_DOWN ) any_regressed = .true.
        enddo
    end function any_regressed

    subroutine print( self )
        class(abinitio3D_addon_report), intent(in) :: self
        character(len=*), parameter :: PFX = '>>> ABINITIO3D_ADDON REPORT '
        integer :: i, s
        do i = 1, self%nstages
            associate( sg => self%stages(i) )
                if( .not. sg%l_set ) cycle
                if( sg%l_base )then
                    write(logfhandle,'(A,A,I0,4(A,F7.2),A)') PFX, 'STAGE ', i, ': LP ', sg%lp_addon, ' LPSTOP ', &
                        &sg%lpstop_addon, ' (BASE RUN LP ', sg%lp_base, ' LPSTOP ', sg%lpstop_base, ')'
                else
                    write(logfhandle,'(A,A,I0,2(A,F7.2),A)') PFX, 'STAGE ', i, ': LP ', sg%lp_addon, ' LPSTOP ', &
                        &sg%lpstop_addon, ' (NOT RUN BY THE BASE RUN)'
                endif
            end associate
        enddo
        do s = 1, self%nstates
            associate( st => self%states(s) )
                write(logfhandle,'(A,A,I0,A,I0,A,I0,A)') PFX, 'STATE ', s, ': UNION POPULATION ', st%pop_union, &
                    &' (FROZEN ', st%pop_frozen, ')'
                if( st%l_fsc )then
                    write(logfhandle,'(A,A,I0,2(A,F7.2),A,SP,I0,SS,A,2(A,F7.2),A,F7.3,A,A)') PFX, 'STATE ', s, &
                        &': FSC=0.143 ', st%res0143_union, ' A (BASE ', st%res0143_base, ' A, ', st%dshell, ' SHELLS)', &
                        &', FSC=0.5 ', st%res05_union, ' A (BASE ', st%res05_base, ' A), MEAN FSC GAIN ', st%fsc_gain, &
                        &': ', trim(st%verdict)
                else
                    write(logfhandle,'(A,A,I0,A)') PFX, 'STATE ', s, ': FSC NOT COMPARED WITH THE BASE RUN'
                endif
                if( st%l_map )then
                    write(logfhandle,'(A,A,I0,A,F7.2,A,F7.4)') PFX, 'STATE ', s, ': CORRELATION WITH THE BASE MAP TO ', &
                        &st%lp_corr, ' A ', st%corr
                    if( st%l_docked ) write(logfhandle,'(A,A,I0,3(A,F7.2))') PFX, 'STATE ', s, &
                        &': DOCKED ONTO THE BASE MAP: ROTATION ', st%dock_angle, ' DEG, SHIFT ', st%dock_shift, &
                        &' A, CORRELATION ', st%dock_corr
                else
                    write(logfhandle,'(A,A,I0,A)') PFX, 'STATE ', s, ': MAP NOT COMPARED WITH THE BASE MAP'
                endif
                if( st%l_cohort ) write(logfhandle,'(A,A,I0,3(A,F7.2),A,F7.4)') PFX, 'STATE ', s, &
                    &': COHORT-ONLY MAP VS BASE MAP FSC=0.5 ', st%res05_cohort, ' A, FSC=0.143 ', st%res0143_cohort, &
                    &' A, CORRELATION TO ', st%lp_corr, ' A ', st%corr_cohort
            end associate
        enddo
    end subroutine print

    !> one text record per listed stage, per state and per compared cohort
    subroutine write( self, fname )
        class(abinitio3D_addon_report), intent(in) :: self
        class(string),                  intent(in) :: fname
        type(oris) :: os
        integer    :: i, s, n
        n = count(self%stages(:)%l_set) + self%nstates + count(self%states(:)%l_cohort)
        call os%new(n, is_ptcl=.false.)
        n = 0
        do i = 1, self%nstages
            associate( sg => self%stages(i) )
                if( .not. sg%l_set ) cycle
                n = n + 1
                call os%set(n, 'kind',         'stage')
                call os%set(n, 'stage',        i)
                call os%set(n, 'base_run',     merge(1, 0, sg%l_base))
                call os%set(n, 'lp_addon',     sg%lp_addon)
                call os%set(n, 'lpstop_addon', sg%lpstop_addon)
                call os%set(n, 'lp_base',      sg%lp_base)
                call os%set(n, 'lpstop_base',  sg%lpstop_base)
            end associate
        enddo
        do s = 1, self%nstates
            associate( st => self%states(s) )
                n = n + 1
                call os%set(n, 'kind',          'state')
                call os%set(n, 'state',         s)
                call os%set(n, 'pop_union',     st%pop_union)
                call os%set(n, 'pop_frozen',    st%pop_frozen)
                call os%set(n, 'verdict',       trim(st%verdict))
                call os%set(n, 'fsc_compared',  merge(1, 0, st%l_fsc))
                call os%set(n, 'res0143_union', st%res0143_union)
                call os%set(n, 'res0143_base',  st%res0143_base)
                call os%set(n, 'res05_union',   st%res05_union)
                call os%set(n, 'res05_base',    st%res05_base)
                call os%set(n, 'dshell0143',    st%dshell)
                call os%set(n, 'fsc_gain',      st%fsc_gain)
                call os%set(n, 'map_compared',  merge(1, 0, st%l_map))
                call os%set(n, 'corr_lp',       st%lp_corr)
                call os%set(n, 'corr',          st%corr)
                call os%set(n, 'docked',        merge(1, 0, st%l_docked))
                call os%set(n, 'dock_angle',    st%dock_angle)
                call os%set(n, 'dock_shift',    st%dock_shift)
                call os%set(n, 'dock_corr',     st%dock_corr)
            end associate
        enddo
        do s = 1, self%nstates
            associate( st => self%states(s) )
                if( .not. st%l_cohort ) cycle
                n = n + 1
                call os%set(n, 'kind',           'cohort')
                call os%set(n, 'state',          s)
                call os%set(n, 'res0143_cohort', st%res0143_cohort)
                call os%set(n, 'res05_cohort',   st%res05_cohort)
                call os%set(n, 'corr_cohort',    st%corr_cohort)
            end associate
        enddo
        call os%write(fname)
        call os%kill
    end subroutine write

    !> the report as write left it (sampling and mask diameter are not
    !! recorded: a read report compares nothing more)
    subroutine read( self, fname )
        class(abinitio3D_addon_report), intent(inout) :: self
        class(string),                  intent(in)    :: fname
        type(oris)   :: os
        type(string) :: rec_kind, verdict
        integer      :: i, n, s, nstates, nstages
        if( .not. file_exists(fname) ) THROW_HARD('no abinitio3D_addon report: '//fname%to_char())
        n = nlines(fname)
        call os%new(n, is_ptcl=.false.)
        call os%read(fname)
        nstates = 0
        nstages = 0
        do i = 1, n
            rec_kind = os%get_str(i, 'kind')
            if( rec_kind%to_char() == 'state' ) nstates = max(nstates, os%get_int(i, 'state'))
            if( rec_kind%to_char() == 'stage' ) nstages = max(nstages, os%get_int(i, 'stage'))
        enddo
        call self%new(nstates, nstages, 0., 0.)
        do i = 1, n
            rec_kind = os%get_str(i, 'kind')
            select case(rec_kind%to_char())
                case('stage')
                    call self%set_stage(os%get_int(i, 'stage'), os%get(i, 'lp_addon'), os%get(i, 'lpstop_addon'), &
                        &os%get(i, 'lp_base'), os%get(i, 'lpstop_base'))
                case('state')
                    s = os%get_int(i, 'state')
                    associate( st => self%states(s) )
                        st%pop_union     = os%get_int(i, 'pop_union')
                        st%pop_frozen    = os%get_int(i, 'pop_frozen')
                        verdict          = os%get_str(i, 'verdict')
                        st%verdict       = verdict%to_char()
                        st%l_fsc         = os%get_int(i, 'fsc_compared') == 1
                        st%res0143_union = os%get(i, 'res0143_union')
                        st%res0143_base  = os%get(i, 'res0143_base')
                        st%res05_union   = os%get(i, 'res05_union')
                        st%res05_base    = os%get(i, 'res05_base')
                        st%dshell        = os%get_int(i, 'dshell0143')
                        st%fsc_gain      = os%get(i, 'fsc_gain')
                        st%l_map         = os%get_int(i, 'map_compared') == 1
                        st%lp_corr       = os%get(i, 'corr_lp')
                        st%corr          = os%get(i, 'corr')
                        st%l_docked      = os%get_int(i, 'docked') == 1
                        st%dock_angle    = os%get(i, 'dock_angle')
                        st%dock_shift    = os%get(i, 'dock_shift')
                        st%dock_corr     = os%get(i, 'dock_corr')
                    end associate
                case('cohort')
                    s = os%get_int(i, 'state')
                    associate( st => self%states(s) )
                        st%l_cohort       = .true.
                        st%res0143_cohort = os%get(i, 'res0143_cohort')
                        st%res05_cohort   = os%get(i, 'res05_cohort')
                        st%corr_cohort    = os%get(i, 'corr_cohort')
                    end associate
            end select
        enddo
        call os%kill
        call rec_kind%kill
        call verdict%kill
    end subroutine read

    subroutine kill( self )
        class(abinitio3D_addon_report), intent(inout) :: self
        if( allocated(self%states) ) deallocate(self%states)
        if( allocated(self%stages) ) deallocate(self%stages)
        self%nstates = 0
        self%nstages = 0
        self%smpd    = 0.
        self%mskdiam = 0.
    end subroutine kill

    ! private

    subroutine check_state( self, s )
        class(abinitio3D_addon_report), intent(in) :: self
        integer,                        intent(in) :: s
        if( s < 1 .or. s > self%nstates ) THROW_HARD('state out of range; abinitio3D_addon_report')
    end subroutine check_state

    !> the correlation band's low-pass limit: the base run's FSC=0.143
    !! resolution once compare_fsc found one, else Nyquist (0)
    real function corr_lp( self, s )
        class(abinitio3D_addon_report), intent(in) :: self
        integer,                        intent(in) :: s
        corr_lp = 0.
        if( self%states(s)%l_fsc .and. self%states(s)%res0143_base > 0. ) corr_lp = self%states(s)%res0143_base
    end function corr_lp

    !> the union map docked onto the base map: rotation (deg), shift (A) and
    !! docked correlation
    subroutine dock( self, s, base_fname, union_fname )
        class(abinitio3D_addon_report), intent(inout) :: self
        integer,                        intent(in)    :: s
        class(string),                  intent(in)    :: base_fname, union_fname
        type(dock_vols) :: docker
        type(ori)       :: o_dock, o_ident
        real :: eul(3), shift(3), cc
        call docker%new(base_fname, union_fname, self%smpd, DOCK_HP, max(DOCK_LP_MIN, self%corr_lp(s)), self%mskdiam)
        call docker%srch()
        call docker%get_dock_info(eul, shift, cc)
        call docker%kill()
        call o_dock%new_ori(.false.)
        call o_ident%new_ori(.false.)
        call o_dock%set_euler(eul)
        call o_ident%set_euler([0., 0., 0.])
        self%states(s)%l_docked   = .true.
        self%states(s)%dock_angle = rad2deg(o_dock%geodesic_dist_trace(o_ident))
        self%states(s)%dock_shift = norm2(shift) * self%smpd
        self%states(s)%dock_corr  = cc
        call o_dock%kill
        call o_ident%kill
    end subroutine dock

end module simple_abinitio3D_addon_report
