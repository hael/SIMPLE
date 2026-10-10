!@descr: partitioning, step projects and report of emulate_solve3D_stream, the offline replay of the stream's solve3D and solve3D_addon cycle
! The selected particles (ptcl2D state > 0, in row order) split into the base and consecutive add-on chunks;
! a step project is the prepared source with the rows of later chunks deselected. Contract:
! doc/implementation_notes/planned/emulate_solve3D_stream.md sec. 3, 4.2, 4.3 and 6.
module simple_solve3D_stream_emulation
use simple_core_module_api
use simple_sp_project, only: sp_project
implicit none

public :: partition_selected, prepare_emulation_source, select_emulation_step
public :: emulation_step, emulation_report, EMULATION_REPORT_FNAME
private
#include "simple_local_flags.inc"

character(len=*), parameter :: EMULATION_REPORT_FNAME = 'emulate_solve3D_stream_report.txt'
integer,          parameter :: MAX_STATE_ARTIFACTS    = 20  !< states whose map and FSC entries the source preparation drops
integer,          parameter :: VERDICT_LEN            = 12

!> one executed step: the base run (kind 'base') or an add-on (kind 'addon')
type :: emulation_step
    character(len=5)      :: kind       = 'base'
    integer               :: ipart      = 0     !< chunk index: 0 for the base
    integer               :: nfrozen    = 0     !< selected particles of the adopted results the step started from (0 for the base)
    integer               :: nadded     = 0     !< selected particles the step searches
    integer               :: last_stage = 0     !< the stage the step ran last (from its manifest)
    integer               :: box_crop   = 0
    real                  :: lp         = 0.    !< that stage's low-pass limit (A)
    real                  :: seconds    = 0.    !< wall time of the step
    logical               :: l_adopted  = .true.
    logical               :: l_regressed = .false. !< the add-on report holds a REGRESSED state
    character(len=STDLEN) :: projfile   = ''
    real,              allocatable :: res0143(:), res05(:)   !< per state, A (0: not available)
    real,              allocatable :: corr(:), res_cohort(:) !< per state, add-on steps (0: not compared)
    integer,           allocatable :: dshell(:)
    character(len=VERDICT_LEN), allocatable :: verdict(:)
end type emulation_step

!> the steps run so far, written as text after every step
type :: emulation_report
    private
    character(len=:), allocatable :: command, source
    integer :: nselected = 0, nbase = 0, naddon = 0, nchunks = 0, nstates = 0
    logical :: l_rollback = .true.
    integer :: nsteps = 0
    type(emulation_step), allocatable :: steps(:)
  contains
    procedure :: new
    procedure :: add_step
    procedure :: get_nsteps
    procedure :: write
    procedure :: kill
end type emulation_report

contains

    !> Split the selected rows into the base (chunk 0) and consecutive add-on chunks.
    !! The first nbase selected particles are the base; the rest go in chunks of
    !! naddon particles, the last taking the remainder when it is at least naddon
    !! (so it lies in [naddon, 2*naddon)), and being the only chunk, however
    !! small, when the rest is below naddon. part(i) is the chunk of row i, -1
    !! for a row that is not selected. A chunk of fewer than nmin particles, and a
    !! base that leaves nothing to add, are refused (status /= 0, msg names the
    !! defect).
    subroutine partition_selected( selected, nbase, naddon, nmin, part, nchunks, status, msg )
        logical,              intent(in)  :: selected(:)
        integer,              intent(in)  :: nbase, naddon, nmin
        integer, allocatable, intent(out) :: part(:)
        integer,              intent(out) :: nchunks, status
        character(len=*),     intent(out) :: msg
        integer :: n, nsel, nrem, i, rank
        status  = 1
        msg     = ''
        nchunks = 0
        n       = size(selected)
        nsel    = count(selected)
        allocate(part(n), source=-1)
        if( nbase < nmin )then
            write(msg,'(A,I0,A,I0,A)') 'nptcls_base (', nbase, ') is below the minimum of ', nmin, ' particles (5 per state)'
            return
        endif
        if( naddon < nmin )then
            write(msg,'(A,I0,A,I0,A)') 'nptcls_addon (', naddon, ') is below the minimum of ', nmin, ' particles (5 per state)'
            return
        endif
        if( nbase >= nsel )then
            write(msg,'(A,I0,A,I0,A)') 'nptcls_base (', nbase, ') leaves nothing to add: the project has ', nsel, &
                &' selected particles'
            return
        endif
        nrem = nsel - nbase
        if( nrem < nmin )then
            write(msg,'(A,I0,A,I0,A)') 'only ', nrem, ' selected particles remain after the base run, below the minimum of ', &
                &nmin, ' particles (5 per state)'
            return
        endif
        nchunks = max(1, nrem / naddon)
        rank    = 0
        do i = 1, n
            if( .not. selected(i) ) cycle
            rank = rank + 1
            if( rank <= nbase )then
                part(i) = 0
            else
                part(i) = min(nchunks, (rank - nbase - 1) / naddon + 1)
            endif
        enddo
        status = 0
    end subroutine partition_selected

    !> The source copy every step project is cut from: no 3D solution
    !! (alignment, update counts, class 3D, state maps and FSCs, state weights)
    !! and no run registration, so a step project is a plain project with the
    !! chosen particles selected. The 2D solution stays. ptcl3D states follow
    !! ptcl2D's selection (1 selected, 0 not).
    subroutine prepare_emulation_source( spproj )
        class(sp_project), intent(inout) :: spproj
        integer :: i, s
        call spproj%os_cls3D%kill
        if( spproj%os_ptcl3D%get_noris() > 0 )then
            call spproj%os_ptcl3D%delete_3Dalignment(keepshifts=.true.)
            call spproj%os_ptcl3D%clean_entry('updatecnt', 'sampled')
            call spproj%os_ptcl3D%delete_entry('res')
            call spproj%os_ptcl3D%delete_entry('res05')
            do i = 1, spproj%os_ptcl3D%get_noris()
                if( spproj%os_ptcl2D%get_state(i) > 0 )then
                    call spproj%os_ptcl3D%set_state(i, 1)
                else
                    call spproj%os_ptcl3D%set_state(i, 0)
                endif
            enddo
        endif
        do s = 1, MAX_STATE_ARTIFACTS
            call spproj%remove_state_artifacts_from_osout(s)
        enddo
        call spproj%remove_state_weights_from_osout
        if( spproj%projinfo%isthere(1, 'sigma2_state')     ) call spproj%projinfo%delete_entry('sigma2_state')
        if( spproj%projinfo%isthere(1, 'solve3D_manifest') ) call spproj%projinfo%delete_entry('solve3D_manifest')
        if( spproj%projinfo%isthere(1, 'solve3D_run_id')   ) call spproj%projinfo%delete_entry('solve3D_run_id')
    end subroutine prepare_emulation_source

    !> The project of the step that has chunks 0..kmax: the rows of later
    !! chunks, and of unselected particles, are deselected in ptcl2D and ptcl3D.
    subroutine select_emulation_step( spproj, part, kmax )
        class(sp_project), intent(inout) :: spproj
        integer,           intent(in)    :: part(:), kmax
        integer :: i
        if( size(part) /= spproj%os_ptcl2D%get_noris() .or. size(part) /= spproj%os_ptcl3D%get_noris() )then
            THROW_HARD('the partition does not cover the rows of the project; select_emulation_step')
        endif
        do i = 1, size(part)
            if( part(i) >= 0 .and. part(i) <= kmax ) cycle
            call spproj%os_ptcl2D%set_state(i, 0)
            call spproj%os_ptcl3D%set_state(i, 0)
        enddo
    end subroutine select_emulation_step

    ! ---- report ----------------------------------------------------------------

    subroutine new( self, command, source, nselected, nbase, naddon, nchunks, nstates, l_rollback )
        class(emulation_report), intent(inout) :: self
        character(len=*),        intent(in)    :: command, source
        integer,                 intent(in)    :: nselected, nbase, naddon, nchunks, nstates
        logical,                 intent(in)    :: l_rollback
        call self%kill
        self%command    = trim(command)
        self%source     = trim(source)
        self%nselected  = nselected
        self%nbase      = nbase
        self%naddon     = naddon
        self%nchunks    = nchunks
        self%nstates    = nstates
        self%l_rollback = l_rollback
        allocate(self%steps(0))
    end subroutine new

    subroutine add_step( self, step )
        class(emulation_report), intent(inout) :: self
        type(emulation_step),    intent(in)    :: step
        type(emulation_step), allocatable :: tmp(:)
        if( .not. allocated(self%steps) ) allocate(self%steps(0))
        allocate(tmp(self%nsteps + 1))
        if( self%nsteps > 0 ) tmp(1:self%nsteps) = self%steps(1:self%nsteps)
        tmp(self%nsteps + 1) = step
        call move_alloc(tmp, self%steps)
        self%nsteps = self%nsteps + 1
    end subroutine add_step

    integer function get_nsteps( self )
        class(emulation_report), intent(in) :: self
        get_nsteps = self%nsteps
    end function get_nsteps

    !> The text report, rewritten whole: what was run, a table of the steps, then
    !! the per-state resolutions and verdicts of each step.
    subroutine write( self, fname )
        class(emulation_report), intent(in) :: self
        class(string),           intent(in) :: fname
        character(len=STDLEN) :: line
        character(len=16)     :: sec1000, kindres
        integer :: funit, io_stat, i, s
        call fopen(funit, file=fname, status='replace', action='write', iostat=io_stat)
        if( io_stat /= 0 )then
            THROW_WARN('could not write the emulation report: '//fname%to_char())
            return
        endif
        write(funit,'(A)') 'EMULATE_SOLVE3D_STREAM REPORT'
        write(funit,'(A)') 'source project: '//trim(self%source)
        write(funit,'(A)') 'command: '//trim(self%command)
        write(funit,'(A,I0,A,I0,A,I0,A,I0,A,I0,A,A)') 'selected particles: ', self%nselected, '  base: ', self%nbase, &
            &'  add-on chunk: ', self%naddon, '  add-on steps: ', self%nchunks, '  states: ', self%nstates, &
            &'  rollback: ', merge('yes', 'no ', self%l_rollback)
        write(funit,'(A)') ''
        write(funit,'(A)') 'step kind   chunk   frozen   added   union  stage  lp(A)  box   wall(s)  s/1000added  result'
        do i = 1, self%nsteps
            associate( st => self%steps(i) )
                sec1000 = '-'
                if( st%kind == 'addon' .and. st%nadded > 0 ) write(sec1000,'(F12.1)') st%seconds / (real(st%nadded) / 1000.)
                kindres = 'adopted'
                if( .not. st%l_adopted ) kindres = 'ROLLED BACK'
                write(line,'(I4,1X,A5,1X,I6,1X,I8,1X,I7,1X,I7,1X,I5,1X,F6.2,1X,I4,1X,F9.1,1X,A12,2X,A)') i - 1, st%kind, &
                    &st%ipart, st%nfrozen, st%nadded, st%nfrozen + st%nadded, st%last_stage, st%lp, st%box_crop, st%seconds, &
                    &adjustr(sec1000(1:12)), trim(kindres)
                write(funit,'(A)') trim(line)
            end associate
        enddo
        write(funit,'(A)') ''
        write(funit,'(A)') 'per state (step, state: FSC=0.143 / FSC=0.5 resolution of the result in A; add-on steps: verdict against the'
        write(funit,'(A)') 'frozen solution, FSC=0.143 shell shift, map correlation, cohort-only resolution)'
        do i = 1, self%nsteps
            associate( st => self%steps(i) )
                if( .not. allocated(st%res0143) ) cycle
                do s = 1, size(st%res0143)
                    if( st%kind == 'addon' .and. allocated(st%verdict) )then
                        write(line,'(A,I0,A,I0,A,F7.2,A,F7.2,A,A,A,I0,A,F7.4,A,F7.2)') 'step ', i - 1, ' state ', s, &
                            &': FSC0.143 ', st%res0143(s), '  FSC0.5 ', st%res05(s), '  ', trim(st%verdict(s)), &
                            &'  dshell ', st%dshell(s), '  corr ', st%corr(s), '  cohort-only ', st%res_cohort(s)
                    else
                        write(line,'(A,I0,A,I0,A,F7.2,A,F7.2)') 'step ', i - 1, ' state ', s, &
                            &': FSC0.143 ', st%res0143(s), '  FSC0.5 ', st%res05(s)
                    endif
                    write(funit,'(A)') trim(line)
                enddo
            end associate
        enddo
        call fclose(funit)
    end subroutine write

    subroutine kill( self )
        class(emulation_report), intent(inout) :: self
        if( allocated(self%command) ) deallocate(self%command)
        if( allocated(self%source)  ) deallocate(self%source)
        if( allocated(self%steps)   ) deallocate(self%steps)
        self%nselected  = 0
        self%nbase      = 0
        self%naddon     = 0
        self%nchunks    = 0
        self%nstates    = 0
        self%nsteps     = 0
        self%l_rollback = .true.
    end subroutine kill

end module simple_solve3D_stream_emulation
