!@descr: unit tests of chunk planning, source preparation, step projects and the report (simple_solve3D_chunks_plan)
! The automatic plan (nsample per state count, capped runs of whole stacks) on synthetic selections; the
! 3D-free source and the step projects on a 20-row fixture; the solve3D_chunks report.
module simple_solve3D_chunks_plan_tester
use simple_defs,                    only: STDLEN
use simple_string,                  only: string
use simple_sp_project,              only: sp_project
use simple_solve3D_chunks_plan
use simple_test_utils
implicit none
private
public :: run_all_solve3D_chunks_plan_tests

integer, parameter :: NROWS = 20 !< rows of the fixture; rows 3 and 7 are not selected (18 selected particles)

contains

    subroutine run_all_solve3D_chunks_plan_tests()
        write(*,'(A)') '**** running all solve3D chunks plan tests ****'
        call test_chunk_nsample()
        call test_chunk_run_count()
        call test_balanced_partition()
        call test_source_and_step_projects()
        call test_chunks_report()
    end subroutine run_all_solve3D_chunks_plan_tests

    !> the stack index of n rows in consecutive stacks of stksz rows
    subroutine make_stacks( stkind, n, stksz )
        integer, allocatable, intent(out) :: stkind(:)
        integer,              intent(in)  :: n, stksz
        integer :: i
        allocate(stkind(n))
        do i = 1, n
            stkind(i) = (i - 1) / stksz + 1
        enddo
    end subroutine make_stacks

    !> the number of rows in chunk k
    integer function chunk_size( part, k ) result( n )
        integer, intent(in) :: part(:), k
        n = count(part == k)
    end function chunk_size

    !> solve3D samples nsample*nstates per iteration: the derived nsample keeps that total at the fixed target
    subroutine test_chunk_nsample()
        integer :: s
        logical :: l_total_kept
        write(*,'(A)') 'test_chunk_nsample'
        call assert_int(5000,  NSAMPLE_TOTAL_CHUNK, 'the total sampled per iteration is fixed')
        call assert_int(10000, NPTCLS_CHUNK_MIN,    'the particle minimum is two sampling totals')
        call assert_int(30000, NPTCLS_CHUNK_MAX,    'the particle cap is six sampling totals')
        call assert_int(5000, chunk_nsample(1), 'one state samples the whole total')
        call assert_int(2500, chunk_nsample(2), 'two states split the total')
        call assert_int(1666, chunk_nsample(3), 'three states round nsample down')
        call assert_int(1250, chunk_nsample(4), 'four states split the total')
        l_total_kept = .true.
        do s = 1, 20
            if( chunk_nsample(s) * s > NSAMPLE_TOTAL_CHUNK      ) l_total_kept = .false.
            if( chunk_nsample(s) * s <= NSAMPLE_TOTAL_CHUNK - s ) l_total_kept = .false.
        enddo
        call assert_true(l_total_kept, 'nsample*nstates stays within nstates-1 of the total for 1..20 states')
    end subroutine test_chunk_nsample

    subroutine test_chunk_run_count()
        integer :: nruns, status, n
        logical :: l_bounds_kept
        character(len=STDLEN) :: msg
        write(*,'(A)') 'test_chunk_run_count'
        call derive_chunk_run_count(100000, nruns, status, msg)
        call assert_int(0, status, '100000 particles have a legal plan')
        call assert_int(4, nruns, 'the plan uses the fewest capped runs')
        call derive_chunk_run_count(10000, nruns, status, msg)
        call assert_int(0, status, 'the exact minimum is one legal run')
        call assert_int(1, nruns, 'the minimum needs no add-on')
        call derive_chunk_run_count(30000, nruns, status, msg)
        call assert_int(0, status, 'the exact maximum is one legal run')
        call assert_int(1, nruns, 'the maximum needs no add-on')
        call derive_chunk_run_count(30001, nruns, status, msg)
        call assert_int(0, status, 'one particle over the cap has a legal balanced plan')
        call assert_int(2, nruns, 'one particle over the cap needs two runs')
        call derive_chunk_run_count(9999, nruns, status, msg)
        call assert_true(status /= 0, 'fewer particles than the minimum are refused')
        call assert_int(0, nruns, 'a refused plan has no runs')
        ! every legal count gives balanced runs inside the bounds
        l_bounds_kept = .true.
        do n = NPTCLS_CHUNK_MIN, 20 * NPTCLS_CHUNK_MAX, 997
            call derive_chunk_run_count(n, nruns, status, msg)
            if( status /= 0 )then
                l_bounds_kept = .false.
                cycle
            endif
            ! the smallest run meets the minimum, the largest the cap, and one run fewer would exceed the cap
            if( n / nruns < NPTCLS_CHUNK_MIN ) l_bounds_kept = .false.
            if( (n + nruns - 1) / nruns > NPTCLS_CHUNK_MAX ) l_bounds_kept = .false.
            if( nruns > 1 )then
                if( (n - 1) / (nruns - 1) < NPTCLS_CHUNK_MAX ) l_bounds_kept = .false.
            endif
        enddo
        call assert_true(l_bounds_kept, 'every legal count splits into the fewest runs within the bounds')
    end subroutine test_chunk_run_count

    !> Runs are consecutive whole stacks in row order: stacks of 50 or 37 rows, unselected rows inside,
    !! base-only bounds, a cut that must fall inside a single large stack, and the refusals.
    subroutine test_balanced_partition()
        logical, allocatable :: selected(:)
        integer, allocatable :: part(:), stkind(:)
        integer :: nchunks, ninside, status, i, sizes(0:2)
        logical :: l_whole
        character(len=STDLEN) :: msg
        write(*,'(A)') 'test_balanced_partition'
        ! Four runs whose balanced cuts fall exactly on stack boundaries.
        allocate(selected(100000), source=.true.)
        call make_stacks(stkind, 100000, 50)
        call partition_selected_balanced(selected, stkind, part, nchunks, ninside, status, msg)
        call assert_int(0, status, '100000 active particles split successfully')
        call assert_int(3, nchunks, 'four total runs mean three add-on chunks')
        call assert_int(0, ninside, 'every cut falls on a stack boundary')
        call assert_true(all([(chunk_size(part, i), i=0,3)] == 25000), 'all four runs hold 25000 particles')
        deallocate(selected)
        ! Stacks of 37 rows with unselected rows: whole stacks, sizes within one stack of the balanced share.
        allocate(selected(60005), source=.true.)
        call make_stacks(stkind, 60005, 37)
        selected([2, 20003, 40004, 60005]) = .false.
        call partition_selected_balanced(selected, stkind, part, nchunks, ninside, status, msg)
        call assert_int(0, status, '60001 active particles with holes split successfully')
        call assert_int(2, nchunks, 'the cap makes three total runs')
        call assert_int(0, ninside, 'no cut falls inside a stack')
        sizes   = [(chunk_size(part, i), i=0,2)]
        l_whole = .true.
        do i = 2, size(part)
            if( stkind(i) /= stkind(i-1) ) cycle
            if( part(i) >= 0 .and. part(i-1) >= 0 .and. part(i) /= part(i-1) ) l_whole = .false.
        enddo
        call assert_true(l_whole, 'no stack is split between runs')
        call assert_int(60001, sum(sizes), 'every selected row is in one run')
        call assert_true(maxval(abs(sizes - 20000)) <= 37, 'run sizes stay within one stack of the balanced share')
        call assert_true(all(part([2, 20003, 40004, 60005]) == -1), 'inactive holes belong to no chunk')
        deallocate(selected)
        ! A legal one-run plan labels the base only, at both bounds.
        allocate(selected(10000), source=.true.)
        call make_stacks(stkind, 10000, 50)
        call partition_selected_balanced(selected, stkind, part, nchunks, ninside, status, msg)
        call assert_int(0, status, 'the minimum is accepted')
        call assert_int(0, nchunks, 'the minimum is a base-only run')
        call assert_true(all(part == 0), 'every active row is in the base')
        deallocate(selected)
        allocate(selected(30000), source=.true.)
        call make_stacks(stkind, 30000, 50)
        call partition_selected_balanced(selected, stkind, part, nchunks, ninside, status, msg)
        call assert_int(0, status, 'the chunk maximum is accepted')
        call assert_int(0, nchunks, 'the maximum is a base-only run')
        call assert_int(30000, chunk_size(part, 0), 'the base holds the full capped selection')
        deallocate(selected)
        ! One particle over the cap: the cut moves to the nearest stack boundary.
        allocate(selected(30001), source=.true.)
        call make_stacks(stkind, 30001, 50)
        call partition_selected_balanced(selected, stkind, part, nchunks, ninside, status, msg)
        call assert_int(0, status, 'one particle over the cap splits successfully')
        call assert_int(1, nchunks, 'one particle over the cap makes one add-on')
        call assert_int(15000, chunk_size(part, 0), 'the base ends on the nearest stack boundary')
        call assert_int(15001, chunk_size(part, 1), 'the add-on takes the rest')
        ! A single stack has no boundary: the cut falls inside it, at the balanced share.
        stkind = 1
        call partition_selected_balanced(selected, stkind, part, nchunks, ninside, status, msg)
        call assert_int(1, ninside, 'a cut without a usable boundary falls inside the stack')
        call assert_int(15001, chunk_size(part, 0), 'the inside cut keeps the balanced share')
        deallocate(selected)
        ! A refused plan leaves every row unassigned.
        allocate(selected(9999), source=.true.)
        call make_stacks(stkind, 9999, 50)
        call partition_selected_balanced(selected, stkind, part, nchunks, ninside, status, msg)
        call assert_true(status /= 0, '9999 particles are refused')
        call assert_int(0, nchunks, 'a refused plan has no add-on chunks')
        call assert_true(all(part == -1), 'a refused plan assigns no row')
        deallocate(stkind)
        allocate(stkind(10), source=1)
        call partition_selected_balanced(selected, stkind, part, nchunks, ninside, status, msg)
        call assert_true(status /= 0, 'stack indices of another length are refused')
    end subroutine test_balanced_partition

    !> NROWS rows in one stack, rows 3 and 7 unselected, with a leftover 3D
    !! solution and run registration the source preparation removes
    subroutine make_project( spproj )
        type(sp_project), intent(inout) :: spproj
        integer :: i
        call spproj%kill
        call spproj%projinfo%new(1, is_ptcl=.false.)
        call spproj%projinfo%set(1, 'projname', 'chunks')
        call spproj%projinfo%set(1, 'sigma2_state',     '/data/run/sigma2_state.bin')
        call spproj%projinfo%set(1, 'solve3D_manifest', '/data/run/solve3D_manifest.txt')
        call spproj%projinfo%set(1, 'solve3D_run_id',   'solve3D_run')
        call spproj%os_stk%new(1, is_ptcl=.false.)
        call spproj%os_stk%set(1, 'stk',   '/data/stacks/stack_1.mrcs')
        call spproj%os_stk%set(1, 'fromp', 1)
        call spproj%os_stk%set(1, 'top',   NROWS)
        call spproj%os_stk%set(1, 'nptcls_stk', NROWS)
        call spproj%os_stk%set(1, 'box',   64)
        call spproj%os_stk%set(1, 'smpd',  1.3)
        call spproj%os_ptcl3D%new(NROWS, is_ptcl=.true.)
        do i = 1, NROWS
            call spproj%os_ptcl3D%set(i, 'stkind', 1)
            call spproj%os_ptcl3D%set(i, 'indstk', i)
            call spproj%os_ptcl3D%set_state(i, 1)
        enddo
        spproj%os_ptcl2D = spproj%os_ptcl3D
        call spproj%os_ptcl2D%set_state(3, 0)
        call spproj%os_ptcl2D%set_state(7, 0)
        ! a leftover 3D solution: the 3D state is not the selection, orientations, update counts, a resolution
        do i = 1, NROWS
            call spproj%os_ptcl3D%set_euler(i, [real(10*i), real(5*i), real(3*i)])
            call spproj%os_ptcl3D%set(i, 'updatecnt', 4)
            call spproj%os_ptcl3D%set(i, 'res',       7.5)
        enddo
        call spproj%os_ptcl3D%set_state(3, 1)
        call spproj%os_ptcl3D%set_state(5, 0)
    end subroutine make_project

    subroutine test_source_and_step_projects()
        type(sp_project)     :: source, step
        integer, allocatable :: part(:)
        integer :: i, n2D, n3D
        write(*,'(A)') 'test_source_and_step_projects'
        call make_project(source)
        call prepare_chunks_source(source)
        call assert_true(.not. source%projinfo%isthere(1, 'sigma2_state'),     'the source drops the sigma2 registration')
        call assert_true(.not. source%projinfo%isthere(1, 'solve3D_manifest'), 'the source drops the run manifest registration')
        call assert_true(.not. source%projinfo%isthere(1, 'solve3D_run_id'),   'the source drops the run identifier')
        call assert_int(0, source%os_ptcl3D%get_updatecnt(1), 'the source has no update counts')
        call assert_true(.not. source%os_ptcl3D%isthere(1, 'res'), 'the source has no resolution field')
        call assert_int(0, source%os_ptcl3D%get_state(3), 'ptcl3D follows ptcl2D: unselected row 3 is state 0')
        call assert_int(1, source%os_ptcl3D%get_state(5), 'ptcl3D follows ptcl2D: selected row 5 is state 1')
        call assert_int(NROWS, source%os_ptcl3D%get_noris(), 'the rows are kept')
        ! the 18 selected rows in order: a base of 6 (rows 1-2, 4-6, 8), then three runs of 4 (rows 9-12, 13-16, 17-20)
        allocate(part(NROWS), source=-1)
        part([1,2,4,5,6,8]) = 0
        part(9:12)          = 1
        part(13:16)         = 2
        part(17:20)         = 3
        ! the base: chunk 0 only
        step = source
        call select_chunks_step(step, part, 0)
        n2D = 0
        n3D = 0
        do i = 1, NROWS
            if( step%os_ptcl2D%get_state(i) > 0 ) n2D = n2D + 1
            if( step%os_ptcl3D%get_state(i) > 0 ) n3D = n3D + 1
        enddo
        call assert_int(6, n2D, 'the base project selects the base particles in ptcl2D')
        call assert_int(6, n3D, 'the base project selects the base particles in ptcl3D')
        call assert_int(1, step%os_ptcl2D%get_state(1), 'row 1 (base) stays selected')
        call assert_int(0, step%os_ptcl2D%get_state(9), 'row 9 (chunk 1) is deselected in ptcl2D')
        call assert_int(0, step%os_ptcl3D%get_state(9), 'row 9 (chunk 1) is deselected in ptcl3D')
        call step%kill
        ! add-on 2: chunks 0..2
        step = source
        call select_chunks_step(step, part, 2)
        n2D = 0
        do i = 1, NROWS
            if( step%os_ptcl2D%get_state(i) > 0 ) n2D = n2D + 1
        enddo
        call assert_int(14, n2D, 'the add-on 2 project selects the base and chunks 1 and 2')
        call assert_int(0, step%os_ptcl2D%get_state(3), 'an unselected row stays unselected')
        ! the source itself is not changed by a step
        n2D = 0
        do i = 1, NROWS
            if( source%os_ptcl2D%get_state(i) > 0 ) n2D = n2D + 1
        enddo
        call assert_int(18, n2D, 'a step project leaves the source selection alone')
        call step%kill
        call source%kill
    end subroutine test_source_and_step_projects

    subroutine test_chunks_report()
        type(chunks_report)  :: report
        type(string)            :: fname
        character(len=STDLEN)   :: line
        character(len=2*STDLEN) :: all_lines
        integer :: funit, io_stat
        write(*,'(A)') 'test_chunks_report'
        call report%new_chunks('simple_exec prg=solve3D_chunks', 'source.simple', [20001, 20000, 20000], 3, chunk_nsample(3), 0)
        fname = 'solve3D_chunks_report_test.txt'
        call report%write(fname)
        all_lines = ''
        open(newunit=funit, file=fname%to_char(), status='old', action='read', iostat=io_stat)
        call assert_int(0, io_stat, 'the chunks report file was written')
        do
            read(funit,'(A)',iostat=io_stat) line
            if( io_stat /= 0 ) exit
            if( index(line, 'SOLVE3D_CHUNKS REPORT') > 0 ) all_lines = trim(all_lines)//'|'//trim(line)
            if( index(line, 'selected particles:') > 0 ) all_lines = trim(all_lines)//'|'//trim(line)
            if( index(line, 'run particles:') > 0 ) all_lines = trim(all_lines)//'|'//trim(line)
            if( index(line, 'nsample:') > 0 ) all_lines = trim(all_lines)//'|'//trim(line)
            if( index(line, 'whole stacks') > 0 ) all_lines = trim(all_lines)//'|'//trim(line)
        enddo
        close(funit, status='delete')
        call assert_true(index(all_lines, 'SOLVE3D_CHUNKS REPORT') > 0, 'the chunks report has its own header')
        call assert_true(index(all_lines, 'selected particles: 60001  runs: 3  states: 3') > 0, &
            &'the chunks report summarizes the automatic plan')
        call assert_true(index(all_lines, 'nsample: 1666  sampled per iteration (nsample x states): 4998') > 0, &
            &'the chunks report records nsample and the total it samples')
        call assert_true(index(all_lines, 'run particles: 20001 20000 20000') > 0, 'the actual balanced run sizes are listed')
        call assert_true(index(all_lines, 'whole stacks in arrival order; cuts inside a stack: 0') > 0, &
            &'the chunks report names the cut rule')
        call report%kill
    end subroutine test_chunks_report

end module simple_solve3D_chunks_plan_tester
