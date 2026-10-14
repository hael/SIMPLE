!@descr: partitioning, step projects and the report of the solve3D_chunks workflow (solve3D, then solve3D_addon)
! The split cuts the active particles (ptcl2D state > 0) in row order into a base and consecutive add-on runs
! at stack boundaries, with balanced sizes. A step project is the prepared source with later runs deselected.
! Contract: doc/policies/3D/solve3D_chunks_policy.md.
module simple_solve3D_chunks_plan
use simple_core_module_api
use simple_sp_project, only: sp_project
implicit none

public :: partition_selected_balanced, derive_chunk_run_count, chunk_nsample
public :: prepare_chunks_source, select_chunks_step
public :: chunks_step, chunks_report, CHUNKS_REPORT_FNAME
public :: NSAMPLE_TOTAL_CHUNK, NPTCLS_CHUNK_MIN, NPTCLS_CHUNK_MAX
private
#include "simple_local_flags.inc"

character(len=*), parameter :: CHUNKS_REPORT_FNAME    = 'solve3D_chunks_report.txt'
integer,          parameter :: MAX_STATE_ARTIFACTS    = 20  !< states whose map and FSC entries the source preparation drops
integer,          parameter :: VERDICT_LEN            = 12
integer,          parameter :: NSAMPLE_TOTAL_CHUNK    = 5000 !< particles solve3D_chunks samples per iteration, all states together
integer,          parameter :: NPTCLS_CHUNK_MIN       = 2 * NSAMPLE_TOTAL_CHUNK
integer,          parameter :: NPTCLS_CHUNK_MAX       = 6 * NSAMPLE_TOTAL_CHUNK

!> one executed step: the base run (kind 'base') or an add-on (kind 'addon')
type :: chunks_step
    character(len=5)      :: kind       = 'base'
    integer               :: ipart      = 0     !< chunk index: 0 for the base
    integer               :: nfrozen    = 0     !< selected particles of the adopted results the step started from (0 for the base)
    integer               :: nadded     = 0     !< selected particles the step searches
    integer               :: last_stage = 0     !< the stage the step ran last (from its manifest)
    integer               :: box_crop   = 0
    real                  :: lp         = 0.    !< that stage's low-pass limit (A)
    real                  :: seconds    = 0.    !< wall time of the step
    logical               :: l_regressed = .false. !< the add-on report holds a REGRESSED state
    character(len=STDLEN) :: projfile   = ''
    real,              allocatable :: res0143(:), res05(:)   !< per state, A (0: not available)
    real,              allocatable :: corr(:), res_cohort(:) !< per state, add-on steps (0: not compared)
    integer,           allocatable :: dshell(:)
    character(len=VERDICT_LEN), allocatable :: verdict(:)
end type chunks_step

!> the steps run so far, written as text after every step
type :: chunks_report
    private
    character(len=:), allocatable :: command, source
    integer,          allocatable :: part_sizes(:)
    integer :: nselected = 0, nstates = 0, nsample = 0, ninside = 0
    integer :: nsteps = 0
    type(chunks_step), allocatable :: steps(:)
  contains
    procedure :: new_chunks
    procedure :: add_step
    procedure :: get_nsteps
    procedure :: write
    procedure :: kill
end type chunks_report

contains

    !> The nsample that makes solve3D, which samples nsample*nstates particles per iteration, sample
    !! NSAMPLE_TOTAL_CHUNK particles whatever the number of states (rounded down; nstates >= 1).
    pure integer function chunk_nsample( nstates ) result( nsample )
        integer, intent(in) :: nstates
        nsample = NSAMPLE_TOTAL_CHUNK / nstates
    end function chunk_nsample

    !> The fewest runs that keep every run at or below the chunk cap. With the cap at three times the
    !! minimum, every part of a balanced split into two or more runs holds at least half the cap.
    subroutine derive_chunk_run_count( nactive, nruns, status, msg )
        integer,          intent(in)  :: nactive
        integer,          intent(out) :: nruns, status
        character(len=*), intent(out) :: msg
        nruns  = 0
        status = 1
        msg    = ''
        if( nactive < NPTCLS_CHUNK_MIN )then
            write(msg,'(A,I0,A,I0)') 'active particles (', nactive, ') are below the minimum per run ', NPTCLS_CHUNK_MIN
            return
        endif
        nruns  = (nactive - 1) / NPTCLS_CHUNK_MAX + 1
        status = 0
    end subroutine derive_chunk_run_count

    !> Cut the selected rows, in row order, into the fewest capped runs of consecutive whole stacks, as the
    !! stream forms its runs. Cut k goes to the stack boundary nearest its balanced target; when no boundary
    !! keeps the runs within the limits it goes inside a stack (ninside counts those). Unselected rows get -1.
    subroutine partition_selected_balanced( selected, stkind, part, nchunks, ninside, status, msg )
        logical,              intent(in)  :: selected(:)
        integer,              intent(in)  :: stkind(:)
        integer, allocatable, intent(out) :: part(:)
        integer,              intent(out) :: nchunks, ninside, status
        character(len=*),     intent(out) :: msg
        integer, allocatable :: sel_stk(:), cuts(:)
        integer :: nsel, nruns, q, r, k, i, rank, irun, target, lo, hi, a, b, cut
        allocate(part(size(selected)), source=-1)
        nchunks = 0
        ninside = 0
        if( size(stkind) /= size(selected) )then
            status = 1
            msg    = 'the stack indices and the selection differ in size'
            return
        endif
        nsel = count(selected)
        call derive_chunk_run_count(nsel, nruns, status, msg)
        if( status /= 0 ) return
        sel_stk = pack(stkind, selected) ! the stack of each selected particle, by rank
        allocate(cuts(0:nruns), source=0)
        cuts(nruns) = nsel               ! cuts(k): selected particles in runs 0..k-1
        q      = nsel / nruns
        r      = mod(nsel, nruns)
        target = 0
        do k = 1, nruns - 1
            target = target + q
            if( k <= r ) target = target + 1
            ! the cuts that keep run k-1 and the remaining runs within the limits
            a = max(cuts(k-1) + NPTCLS_CHUNK_MIN, nsel - (nruns - k) * NPTCLS_CHUNK_MAX)
            b = min(cuts(k-1) + NPTCLS_CHUNK_MAX, nsel - (nruns - k) * NPTCLS_CHUNK_MIN)
            ! the nearest stack boundaries at or below and above the target (-1: none in [a,b])
            lo = min(target, b)
            do while( lo >= a )
                if( sel_stk(lo) /= sel_stk(lo+1) ) exit
                lo = lo - 1
            enddo
            if( lo < a ) lo = -1
            hi = max(target, a)
            do while( hi <= b )
                if( sel_stk(hi) /= sel_stk(hi+1) ) exit
                hi = hi + 1
            enddo
            if( hi > b ) hi = -1
            if( lo > 0 .and. hi > 0 )then
                cut = lo
                if( hi - target < target - lo ) cut = hi
            else if( lo > 0 )then
                cut = lo
            else if( hi > 0 )then
                cut = hi
            else
                cut     = min(max(target, a), b)
                ninside = ninside + 1
            endif
            cuts(k) = cut
        enddo
        rank = 0
        irun = 0
        do i = 1, size(selected)
            if( .not. selected(i) ) cycle
            rank = rank + 1
            do while( rank > cuts(irun+1) )
                irun = irun + 1
            enddo
            part(i) = irun
        enddo
        nchunks = nruns - 1
    end subroutine partition_selected_balanced

    !> The source copy every step project is cut from: no 3D solution
    !! (alignment, update counts, class 3D, state maps and FSCs, state weights)
    !! and no run registration, so a step project is a plain project with the
    !! chosen particles selected. The 2D solution stays. ptcl3D states follow
    !! ptcl2D's selection (1 selected, 0 not).
    subroutine prepare_chunks_source( spproj )
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
    end subroutine prepare_chunks_source

    !> The project of the step that has chunks 0..kmax: the rows of later
    !! chunks, and of unselected particles, are deselected in ptcl2D and ptcl3D.
    subroutine select_chunks_step( spproj, part, kmax )
        class(sp_project), intent(inout) :: spproj
        integer,           intent(in)    :: part(:), kmax
        integer :: i
        if( size(part) /= spproj%os_ptcl2D%get_noris() .or. size(part) /= spproj%os_ptcl3D%get_noris() )then
            THROW_HARD('the partition does not cover the rows of the project; select_chunks_step')
        endif
        do i = 1, size(part)
            if( part(i) >= 0 .and. part(i) <= kmax ) cycle
            call spproj%os_ptcl2D%set_state(i, 0)
            call spproj%os_ptcl3D%set_state(i, 0)
        enddo
    end subroutine select_chunks_step

    ! ---- report ----------------------------------------------------------------

    !> Initialize the report for the automatic active-balanced workflow; ninside counts the run cuts
    !! that fell inside a stack.
    subroutine new_chunks( self, command, source, part_sizes, nstates, nsample, ninside )
        class(chunks_report), intent(inout) :: self
        character(len=*),        intent(in)    :: command, source
        integer,                 intent(in)    :: part_sizes(:), nstates, nsample, ninside
        call self%kill
        self%command    = trim(command)
        self%source     = trim(source)
        self%part_sizes = part_sizes
        self%nselected  = sum(part_sizes)
        self%nstates    = nstates
        self%nsample    = nsample
        self%ninside    = ninside
        allocate(self%steps(0))
    end subroutine new_chunks

    subroutine add_step( self, step )
        class(chunks_report), intent(inout) :: self
        type(chunks_step),    intent(in)    :: step
        type(chunks_step), allocatable :: tmp(:)
        if( .not. allocated(self%steps) ) allocate(self%steps(0))
        allocate(tmp(self%nsteps + 1))
        if( self%nsteps > 0 ) tmp(1:self%nsteps) = self%steps(1:self%nsteps)
        tmp(self%nsteps + 1) = step
        call move_alloc(tmp, self%steps)
        self%nsteps = self%nsteps + 1
    end subroutine add_step

    integer function get_nsteps( self )
        class(chunks_report), intent(in) :: self
        get_nsteps = self%nsteps
    end function get_nsteps

    !> The text report, rewritten whole: what was run, a table of the steps, then
    !! the per-state resolutions and verdicts of each step.
    subroutine write( self, fname )
        class(chunks_report), intent(in) :: self
        class(string),           intent(in) :: fname
        character(len=STDLEN) :: line
        character(len=16)     :: sec1000, kindres
        integer :: funit, io_stat, i, s
        call fopen(funit, file=fname, status='replace', action='write', iostat=io_stat)
        if( io_stat /= 0 )then
            THROW_WARN('could not write the solve3D_chunks report: '//fname%to_char())
            return
        endif
        write(funit,'(A)') 'SOLVE3D_CHUNKS REPORT'
        write(funit,'(A)') 'source project: '//trim(self%source)
        write(funit,'(A)') 'command: '//trim(self%command)
        write(funit,'(A,I0,A,I0,A,I0)') 'selected particles: ', self%nselected, '  runs: ', size(self%part_sizes), &
            &'  states: ', self%nstates
        write(funit,'(A,I0,A,I0)') 'nsample: ', self%nsample, '  sampled per iteration (nsample x states): ', &
            &self%nsample * self%nstates
        write(funit,'(A,I0)') 'runs: consecutive whole stacks in arrival order; cuts inside a stack: ', &
            &self%ninside
        write(funit,'(A)',advance='no') 'run particles:'
        do i = 1, size(self%part_sizes)
            write(funit,'(1X,I0)',advance='no') self%part_sizes(i)
        enddo
        write(funit,'(A)') ''
        write(funit,'(A)') ''
        write(funit,'(A)') 'step kind   chunk   frozen   added   union  stage  lp(A)  box   wall(s)  s/1000added  result'
        do i = 1, self%nsteps
            associate( st => self%steps(i) )
                sec1000 = '-'
                if( st%kind == 'addon' .and. st%nadded > 0 ) write(sec1000,'(F12.1)') st%seconds / (real(st%nadded) / 1000.)
                kindres = 'adopted'   ! solve3D_chunks adopts every add-on
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
        class(chunks_report), intent(inout) :: self
        if( allocated(self%command)    ) deallocate(self%command)
        if( allocated(self%source)     ) deallocate(self%source)
        if( allocated(self%part_sizes) ) deallocate(self%part_sizes)
        if( allocated(self%steps)      ) deallocate(self%steps)
        self%nselected  = 0
        self%nstates    = 0
        self%nsample    = 0
        self%ninside    = 0
        self%nsteps     = 0
    end subroutine kill

end module simple_solve3D_chunks_plan
