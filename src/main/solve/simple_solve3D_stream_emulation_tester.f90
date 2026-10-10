!@descr: unit tests of the emulate_solve3D_stream partitioning, source preparation, step projects and report (simple_solve3D_stream_emulation)
! 20-row fixture with two unselected rows: the base/chunk split and the merged last chunk, every refusal,
! the 3D-free source and the deselection of later chunks, and the step table of the text report.
module simple_solve3D_stream_emulation_tester
use simple_defs,                    only: STDLEN
use simple_string,                  only: string
use simple_sp_project,              only: sp_project
use simple_solve3D_stream_emulation
use simple_test_utils
implicit none
private
public :: run_all_solve3D_stream_emulation_tests

integer, parameter :: NROWS = 20 !< rows of the fixture; rows 3 and 7 are not selected (18 selected particles)

contains

    subroutine run_all_solve3D_stream_emulation_tests()
        write(*,'(A)') '**** running all solve3D stream emulation tests ****'
        call test_partition()
        call test_partition_refusals()
        call test_source_and_step_projects()
        call test_report()
    end subroutine run_all_solve3D_stream_emulation_tests

    !> NROWS rows with rows 3 and 7 unselected
    subroutine make_selection( selected )
        logical, allocatable, intent(out) :: selected(:)
        allocate(selected(NROWS), source=.true.)
        selected(3) = .false.
        selected(7) = .false.
    end subroutine make_selection

    !> the number of rows in chunk k
    integer function chunk_size( part, k ) result( n )
        integer, intent(in) :: part(:), k
        n = count(part == k)
    end function chunk_size

    subroutine test_partition()
        logical, allocatable :: selected(:)
        integer, allocatable :: part(:)
        integer :: nchunks, status, i, rank
        character(len=STDLEN) :: msg
        write(*,'(A)') 'test_partition'
        call make_selection(selected)
        ! 18 selected: base 6, 12 left in chunks of 4 -> 3 equal chunks
        call partition_selected(selected, 6, 4, 2, part, nchunks, status, msg)
        call assert_int(0, status, 'a valid split is accepted')
        call assert_int(3, nchunks, 'three chunks of 4')
        call assert_int(6, chunk_size(part, 0), 'the base holds nptcls_base particles')
        call assert_int(4, chunk_size(part, 1), 'chunk 1 holds nptcls_addon particles')
        call assert_int(4, chunk_size(part, 2), 'chunk 2 holds nptcls_addon particles')
        call assert_int(4, chunk_size(part, 3), 'chunk 3 holds nptcls_addon particles')
        call assert_int(2, count(part == -1), 'unselected rows belong to no chunk')
        call assert_int(-1, part(3), 'row 3 is unselected')
        call assert_int(-1, part(7), 'row 7 is unselected')
        ! chunks follow the order of the selected particles
        rank = -1
        do i = 1, NROWS
            if( part(i) < 0 ) cycle
            call assert_true(part(i) >= rank, 'chunk indices never decrease along the rows')
            rank = part(i)
        enddo
        call assert_int(0, part(1), 'the first selected row is in the base')
        call assert_int(0, part(8), 'row 8 is the 6th selected particle, the last of the base')
        call assert_int(1, part(9), 'row 9 starts chunk 1')
        ! remainder 12 with chunks of 5: two chunks, the last takes the remainder (5, 7)
        call partition_selected(selected, 6, 5, 2, part, nchunks, status, msg)
        call assert_int(0, status, 'a split with a remainder is accepted')
        call assert_int(2, nchunks, 'the last short chunk is merged with the previous')
        call assert_int(5, chunk_size(part, 1), 'chunk 1 holds nptcls_addon particles')
        call assert_int(7, chunk_size(part, 2), 'the last chunk takes the remainder')
        ! fewer left than one chunk: one add-on on the remainder, the base is not enlarged
        call partition_selected(selected, 15, 5, 2, part, nchunks, status, msg)
        call assert_int(0, status, 'a remainder below one chunk is accepted')
        call assert_int(1, nchunks, 'a single add-on')
        call assert_int(15, chunk_size(part, 0), 'the base keeps nptcls_base particles')
        call assert_int(3, chunk_size(part, 1), 'the add-on takes the remainder')
        ! exactly one chunk's worth
        call partition_selected(selected, 8, 10, 2, part, nchunks, status, msg)
        call assert_int(1, nchunks, 'a remainder of exactly one chunk is one chunk')
        call assert_int(10, chunk_size(part, 1), 'that chunk holds the remainder')
        ! every row is accounted for
        call assert_int(NROWS, count(part >= 0) + count(part == -1), 'every row is in a chunk or unselected')
    end subroutine test_partition

    subroutine test_partition_refusals()
        logical, allocatable :: selected(:)
        integer, allocatable :: part(:)
        integer :: nchunks, status
        character(len=STDLEN) :: msg
        write(*,'(A)') 'test_partition_refusals'
        call make_selection(selected)
        call partition_selected(selected, 1, 4, 2, part, nchunks, status, msg)
        call assert_true(status /= 0, 'a base below the minimum is refused')
        call assert_true(index(msg, 'nptcls_base') > 0, 'the refusal names nptcls_base')
        call assert_int(0, nchunks, 'a refused split has no chunks')
        call assert_int(NROWS, count(part == -1), 'a refused split assigns no row')
        call partition_selected(selected, 6, 1, 2, part, nchunks, status, msg)
        call assert_true(status /= 0, 'a chunk below the minimum is refused')
        call assert_true(index(msg, 'nptcls_addon') > 0, 'the refusal names nptcls_addon')
        call partition_selected(selected, 18, 4, 2, part, nchunks, status, msg)
        call assert_true(status /= 0, 'a base that takes every selected particle is refused')
        call assert_true(index(msg, 'nothing to add') > 0, 'the refusal says nothing is left to add')
        call partition_selected(selected, 19, 4, 2, part, nchunks, status, msg)
        call assert_true(status /= 0, 'a base beyond the selected particles is refused')
        call partition_selected(selected, 17, 4, 2, part, nchunks, status, msg)
        call assert_true(status /= 0, 'a remainder below the minimum is refused')
        call assert_true(index(msg, 'remain') > 0, 'the refusal names the remainder')
    end subroutine test_partition_refusals

    !> NROWS rows in one stack, rows 3 and 7 unselected, with a leftover 3D
    !! solution and run registration the source preparation removes
    subroutine make_project( spproj )
        type(sp_project), intent(inout) :: spproj
        integer :: i
        call spproj%kill
        call spproj%projinfo%new(1, is_ptcl=.false.)
        call spproj%projinfo%set(1, 'projname', 'emulation')
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
        logical, allocatable :: selected(:)
        integer, allocatable :: part(:)
        integer :: nchunks, status, i, n2D, n3D
        character(len=STDLEN) :: msg
        write(*,'(A)') 'test_source_and_step_projects'
        call make_project(source)
        call prepare_emulation_source(source)
        call assert_true(.not. source%projinfo%isthere(1, 'sigma2_state'),     'the source drops the sigma2 registration')
        call assert_true(.not. source%projinfo%isthere(1, 'solve3D_manifest'), 'the source drops the run manifest registration')
        call assert_true(.not. source%projinfo%isthere(1, 'solve3D_run_id'),   'the source drops the run identifier')
        call assert_int(0, source%os_ptcl3D%get_updatecnt(1), 'the source has no update counts')
        call assert_true(.not. source%os_ptcl3D%isthere(1, 'res'), 'the source has no resolution field')
        call assert_int(0, source%os_ptcl3D%get_state(3), 'ptcl3D follows ptcl2D: unselected row 3 is state 0')
        call assert_int(1, source%os_ptcl3D%get_state(5), 'ptcl3D follows ptcl2D: selected row 5 is state 1')
        call assert_int(NROWS, source%os_ptcl3D%get_noris(), 'the rows are kept')
        selected = [(source%os_ptcl2D%get_state(i) > 0, i = 1, NROWS)]
        call partition_selected(selected, 6, 4, 2, part, nchunks, status, msg)
        call assert_int(0, status, 'the fixture splits')
        ! the base: chunk 0 only
        step = source
        call select_emulation_step(step, part, 0)
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
        call select_emulation_step(step, part, 2)
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

    subroutine test_report()
        type(emulation_report) :: report
        type(emulation_step)   :: step
        type(string)           :: fname
        character(len=STDLEN)  :: line
        character(len=2*STDLEN) :: all_lines
        integer :: funit, io_stat
        write(*,'(A)') 'test_report'
        call report%new('simple_exec prg=emulate_solve3D_stream nptcls_base=6', 'source.simple', 18, 6, 4, 3, 2, .true.)
        step%kind       = 'base'
        step%nadded     = 6
        step%last_stage = 5
        step%box_crop   = 96
        step%lp         = 6.
        step%seconds    = 100.
        allocate(step%res0143(2), source=[5.1, 5.3])
        allocate(step%res05(2),   source=[6.2, 6.4])
        call report%add_step(step)
        deallocate(step%res0143, step%res05)
        step%kind       = 'addon'
        step%ipart      = 1
        step%nfrozen    = 6
        step%nadded     = 4
        step%seconds    = 20.
        step%l_adopted  = .false.
        allocate(step%res0143(2), source=[5.0, 5.2])
        allocate(step%res05(2),   source=[6.0, 6.2])
        allocate(step%corr(2), source=[0.98, 0.97])
        allocate(step%res_cohort(2), source=[0., 0.])
        allocate(step%dshell(2), source=[0, 2])
        allocate(step%verdict(2))
        step%verdict(1) = 'UNCHANGED'
        step%verdict(2) = 'REGRESSED'
        call report%add_step(step)
        call assert_int(2, report%get_nsteps(), 'two steps are in the report')
        fname = 'emulate_solve3D_stream_report_test.txt'
        call report%write(fname)
        all_lines = ''
        open(newunit=funit, file=fname%to_char(), status='old', action='read', iostat=io_stat)
        call assert_int(0, io_stat, 'the report file was written')
        do
            read(funit,'(A)',iostat=io_stat) line
            if( io_stat /= 0 ) exit
            if( index(line, ' addon ') > 0 .or. index(line, ' base ') > 0 ) all_lines = trim(all_lines)//'|'//trim(line)
            if( index(line, 'REGRESSED') > 0 ) all_lines = trim(all_lines)//'|REGRESSED-LINE'
        enddo
        close(funit, status='delete')
        call assert_true(index(all_lines, ' base ') > 0, 'the base step is tabulated')
        call assert_true(index(all_lines, ' addon ') > 0, 'the add-on step is tabulated')
        call assert_true(index(all_lines, 'ROLLED BACK') > 0, 'a step that was not adopted is marked')
        call assert_true(index(all_lines, 'REGRESSED-LINE') > 0, 'the per-state verdict is listed')
        call report%kill
    end subroutine test_report

end module simple_solve3D_stream_emulation_tester
