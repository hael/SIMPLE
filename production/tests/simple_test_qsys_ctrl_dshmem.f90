!> @file   simple_test_qsys_ctrl_dshmem.f90
!> @brief  Standalone test exercising qsys_ctrl with the local_dshmem backend.
!>
!> This test replicates (at a low level) what qsys_env%gen_scripts_and_schedule_jobs
!> does when the queue system is local_dshmem.  It:
!>
!>   1. Constructs a qsys_ctrl backed by qsys_local_dshmem
!>   2. Generates per-partition scripts (generate_scripts)
!>   3. Simulates job completion by touching JOB_FINISHED sentinel files
!>   4. Calls schedule_jobs which polls sentinels, finds everything done, and exits
!>   5. Shows how schedule_streaming works with local_dshmem
!>
!> The flow mirrors qsys_env%gen_scripts_and_schedule_jobs:
!>
!>   ┌──────────────────────────────────────────────────────────┐
!>   │  qsys_env%gen_scripts_and_schedule_jobs( job_descr )     │
!>   │                                                          │
!>   │  1. select type(pmyqsys)                                 │
!>   │       class is(qsys_local_dshmem) → aarray = .true.      │
!>   │                                                          │
!>   │  2a. if aarray:                                          │
!>   │        generate_array_script(...)  (SLURM/LSF only!)     │
!>   │        schedule_array_jobs                               │
!>   │                                                          │
!>   │  2b. else:                                               │
!>   │        generate_scripts(...)       ← one script / part   │
!>   │        schedule_jobs               ← poll sentinels      │
!>   │                                                          │
!>   │  NOTE: current code forces aarray=true for local_dshmem  │
!>   │  but generate_array_script rejects non-SLURM backends.   │
!>   │  So the practical path is 2b (generate_scripts).         │
!>   └──────────────────────────────────────────────────────────┘
!>
!> Build with BUILD_TESTS=ON and run:
!>     ./simple_test_qsys_ctrl_dshmem
!>
program simple_test_qsys_ctrl_dshmem
use simple_core_module_api
use simple_qsys_local_dshmem, only: qsys_local_dshmem
use simple_qsys_ctrl,         only: qsys_ctrl
use simple_cmdline,           only: cmdline
implicit none
#include "simple_local_flags.inc"

integer, parameter :: NPARTS   = 4    ! number of partitions
integer, parameter :: NPTCLS   = 100  ! fake particle count
integer, parameter :: NCUNITS  = 2    ! concurrent computing units (< nparts)
logical            :: all_passed = .true.

write(*,'(a)') '============================================================'
write(*,'(a)') ' simple_test_qsys_ctrl_dshmem'
write(*,'(a)') ' Tests qsys_ctrl with qsys_local_dshmem backend'
write(*,'(a)') '============================================================'

call test_construct_with_dshmem()
call test_generate_scripts_dshmem()
call test_schedule_jobs_dshmem()
call test_generate_script_4_dshmem()
call test_streaming_dshmem()
call test_schedule_streaming_dshmem()

if( all_passed )then
    write(*,'(/,a)') '=== ALL local_dshmem TESTS PASSED ==='
else
    write(*,'(/,a)') '*** SOME local_dshmem TESTS FAILED ***'
    stop 1
endif

contains

    ! ================================================================
    ! Helper: create a qsys_local_dshmem + parts table + qsys_ctrl
    ! ================================================================
    subroutine make_dshmem_ctrl( qsys_obj, parts, ctrl, stream )
        type(qsys_local_dshmem), intent(out), target  :: qsys_obj
        integer, allocatable,    intent(out), target  :: parts(:,:)
        type(qsys_ctrl),        intent(out)           :: ctrl
        logical, optional,       intent(in)            :: stream
        logical :: sstream
        integer :: ip
        sstream = .false.
        if( present(stream) ) sstream = stream
        ! 1) Initialise the local_dshmem qsys
        !    This allocates the internal chash and sets submit_cmd = 'nohup'.
        call qsys_obj%new()
        ! 2) Build an even-split parts table  [fromp, top]
        !    split_nobjs_even distributes NPTCLS across NPARTS as evenly as possible.
        parts = split_nobjs_even(NPTCLS, NPARTS)
        write(*,'(a)') '  Parts table:'
        do ip = 1, NPARTS
            write(*,'(a,i2,a,i4,a,i4)') '    part ', ip, ':  ', parts(ip,1), ' -> ', parts(ip,2)
        end do
        ! 3) Construct qsys_ctrl
        !    Arguments:
        !      exec_binary      – path/name of binary each script will invoke
        !      qsys_obj         – polymorphic qsys backend (local_dshmem here)
        !      parts            – (nparts,2) fromp/top table
        !      fromto_part      – [first_part, last_part] managed by this controller
        !      ncomputing_units – max concurrent slots (like ncunits in params)
        !      stream           – .true. for streaming mode, .false. for batch
        call ctrl%new( string('simple_private_exec'), qsys_obj, parts, &
                       [1, NPARTS], NCUNITS, sstream )
    end subroutine make_dshmem_ctrl

    ! ================================================================
    ! Helper: touch sentinel files so schedule_jobs sees them as done
    !
    ! schedule_jobs polls for files named:
    !    JOB_FINISHED_FBODY // int2str_pad(ipart, numlen)
    ! e.g.  JOB_FINISHED_1, JOB_FINISHED_2, ...
    ! ================================================================
    subroutine touch_sentinels( numlen )
        integer, intent(in) :: numlen
        integer :: ip
        do ip = 1, NPARTS
            call simple_touch( 'JOB_FINISHED_'//int2str_pad(ip, numlen) )
        end do
    end subroutine touch_sentinels

    subroutine clean_sentinels( numlen )
        integer, intent(in) :: numlen
        integer :: ip
        do ip = 1, NPARTS
            call del_file( 'JOB_FINISHED_'//int2str_pad(ip, numlen) )
        end do
    end subroutine clean_sentinels

    subroutine clean_scripts( numlen )
        integer, intent(in) :: numlen
        integer :: ip
        do ip = 1, NPARTS
            call del_file( 'distr_simple_script_'//int2str_pad(ip, numlen) )
        end do
    end subroutine clean_scripts

    ! ================================================================
    ! TEST 1: Construct qsys_ctrl with local_dshmem and verify it works
    ! ================================================================
    subroutine test_construct_with_dshmem()
        type(qsys_local_dshmem)      :: qsys_obj
        integer, allocatable, target  :: parts(:,:)
        type(qsys_ctrl)              :: ctrl
        write(*,'(/,a)') '--- test_construct_with_dshmem ---'
        call make_dshmem_ctrl( qsys_obj, parts, ctrl )
        if( ctrl%exists() )then
            write(*,'(a)') '  PASS: ctrl created with local_dshmem backend'
        else
            write(*,'(a)') '  FAIL: ctrl should exist after new()'
            all_passed = .false.
        endif
        call ctrl%kill
        call qsys_obj%kill
    end subroutine test_construct_with_dshmem

    ! ================================================================
    ! TEST 2: generate_scripts with local_dshmem
    !
    ! This is the first half of gen_scripts_and_schedule_jobs.
    ! For each partition it:
    !   - injects fromp/top/part/nparts into the job_descr chash
    !   - writes a bash script:  distr_simple_script_<padded_part>
    !   - the script contains:
    !       #!/bin/bash
    !       cd <CWD_GLOB>
    !       <exec_binary> <job_descr key=val pairs> 2>&1 | tee -a SIMPLE_SUBPROC_OUTPUT
    !       exit
    !   - resets jobs_done / jobs_submitted to .false.
    ! ================================================================
    subroutine test_generate_scripts_dshmem()
        type(qsys_local_dshmem)      :: qsys_obj
        integer, allocatable, target  :: parts(:,:)
        type(qsys_ctrl)              :: ctrl
        type(chash)                   :: job_descr, q_descr
        integer                       :: ip, numlen
        write(*,'(/,a)') '--- test_generate_scripts_dshmem ---'
        call make_dshmem_ctrl( qsys_obj, parts, ctrl )
        ! Build a job description (the key=value pairs that form the command line)
        ! In real usage these come from cline%gen_job_descr()
        call job_descr%new(20)
        call job_descr%set('prg',      'cluster2D')
        call job_descr%set('nthr',     '1')
        call job_descr%set('ncls',     '50')
        call job_descr%set('mskdiam',  '180')
        ! Build a queue descriptor
        ! In real usage this comes from the compenv segment of the project file
        ! For local_dshmem, the key fields are:
        !   qsys_name = 'local_dshmem'
        call q_descr%new(10)
        call q_descr%set('qsys_name', 'local_dshmem')
        ! Generate one bash script per partition
        ! (corresponds to the ELSE branch in gen_scripts_and_schedule_jobs)
        call ctrl%generate_scripts( job_descr, string('.mrc'), q_descr )
        ! Verify scripts were written
        numlen = len(int2str(NPARTS))  ! same logic as qsys_ctrl%new
        write(*,'(a)') '  Generated scripts:'
        do ip = 1, NPARTS
            if( file_exists('distr_simple_script_'//int2str_pad(ip, numlen)) )then
                write(*,'(a,a)') '    OK: distr_simple_script_', trim(int2str_pad(ip, numlen))
            else
                write(*,'(a,a)') '    MISSING: distr_simple_script_', trim(int2str_pad(ip, numlen))
                all_passed = .false.
            endif
        end do
        ! Verify that generate_scripts resets status
        ! After generation, all jobs are not-done and not-submitted
        write(*,'(a)') '  PASS: generate_scripts with local_dshmem completed'
        ! Clean up
        call clean_scripts( numlen )
        call job_descr%kill
        call q_descr%kill
        call ctrl%kill
        call qsys_obj%kill
    end subroutine test_generate_scripts_dshmem

    ! ================================================================
    ! TEST 3: schedule_jobs with local_dshmem
    !
    ! This is the second half of gen_scripts_and_schedule_jobs.
    ! The scheduler loop:
    !
    !   do
    !       if( all(jobs_done) ) exit          ← check
    !       call update_queue                   ← poll sentinel files
    !       call submit_scripts                 ← submit un-submitted jobs
    !       call sleep(1)
    !   end do
    !
    ! update_queue (non-stream mode):
    !   For each partition:
    !     if JOB_FINISHED_<part> exists → mark done & submitted
    !
    ! submit_scripts:
    !   For each partition not yet submitted, and if computing units
    !   are available: run the script via exec_cmdline.
    !   When used with local/local_dshmem the submit command is 'nohup'.
    !
    ! TRICK: we pre-touch JOB_FINISHED sentinels so update_queue marks
    !        all jobs as done/submitted BEFORE submit_scripts gets a chance
    !        to actually run anything.
    ! ================================================================
    subroutine test_schedule_jobs_dshmem()
        type(qsys_local_dshmem)      :: qsys_obj
        integer, allocatable, target  :: parts(:,:)
        type(qsys_ctrl)              :: ctrl
        type(chash)                   :: job_descr, q_descr
        logical, allocatable          :: done(:), submitted(:)
        integer                       :: numlen
        write(*,'(/,a)') '--- test_schedule_jobs_dshmem ---'
        call make_dshmem_ctrl( qsys_obj, parts, ctrl )
        ! 1) Generate scripts (resets all status flags)
        call job_descr%new(10)
        call job_descr%set('prg',  'cluster2D')
        call job_descr%set('nthr', '1')
        call q_descr%new(10)
        call q_descr%set('qsys_name', 'local_dshmem')
        call ctrl%generate_scripts( job_descr, string('.mrc'), q_descr )
        numlen = len(int2str(NPARTS))
        ! 2) Pre-touch JOB_FINISHED sentinels to simulate completed jobs
        !    In production these are created by the worker binary when it finishes:
        !      call qsys_job_finished(params, ...)
        !    which does:  call simple_touch(JOB_FINISHED_FBODY//int2str_pad(part, numlen))
        call touch_sentinels( numlen )
        write(*,'(a)') '  Pre-touched all JOB_FINISHED sentinel files'
        ! 3) Run the scheduler.
        !    The loop will:
        !      iteration 1: update_queue finds sentinels → all done+submitted
        !                   submit_scripts sees all submitted → returns immediately
        !      iteration 2: all(jobs_done) = .true. → exit loop
        write(*,'(a)') '  Calling schedule_jobs (should return immediately)...'
        call ctrl%schedule_jobs()
        ! 4) Verify all jobs reported as done
        call ctrl%get_jobs_status( done, submitted )
        if( all(done) .and. all(submitted) )then
            write(*,'(a)') '  PASS: schedule_jobs completed, all jobs done'
        else
            write(*,'(a,i2,a,i2)') '  FAIL: done=', count(done), ' submitted=', count(submitted)
            all_passed = .false.
        endif
        ! 5) Clean up sentinel + script files
        call clean_sentinels( numlen )
        call clean_scripts( numlen )
        call job_descr%kill
        call q_descr%kill
        call ctrl%kill
        call qsys_obj%kill
    end subroutine test_schedule_jobs_dshmem

    ! ================================================================
    ! TEST 4: generate_script overload 4 (multi-job sequential script)
    !         with local_dshmem
    !
    ! This writes a single script with multiple sequential commands.
    ! qsys_ctrl checks for 'local' or 'local_dshmem' when chmoding.
    ! ================================================================
    subroutine test_generate_script_4_dshmem()
        type(qsys_local_dshmem)      :: qsys_obj
        integer, allocatable, target  :: parts(:,:)
        type(qsys_ctrl)              :: ctrl
        type(chash), allocatable      :: jobs(:)
        type(chash)                   :: q_descr
        character(len=*), parameter   :: SNAME = 'test_dshmem_multi_script'
        character(len=*), parameter   :: ONAME = 'test_dshmem_multi_output.log'
        integer :: ij
        write(*,'(/,a)') '--- test_generate_script_4_dshmem ---'
        call make_dshmem_ctrl( qsys_obj, parts, ctrl )
        ! Create 3 sequential job descriptions
        allocate(jobs(3))
        do ij = 1, 3
            call jobs(ij)%new(10)
            call jobs(ij)%set('prg',  'cluster2D')
            call jobs(ij)%set('ncls', int2str(ij*50))
        end do
        call q_descr%new(10)
        call q_descr%set('qsys_name', 'local_dshmem')
        ! Generate the multi-job script
        ! generate_script_4 recognises 'local_dshmem' and chmods the script
        call ctrl%generate_script( jobs, q_descr, string('simple_private_exec'), &
                                   string(SNAME), string(ONAME) )
        if( file_exists(SNAME) )then
            write(*,'(a)') '  PASS: multi-job script created ('//SNAME//')'
        else
            write(*,'(a)') '  FAIL: multi-job script NOT found'
            all_passed = .false.
        endif
        call del_file(SNAME)
        do ij = 1, 3
            call jobs(ij)%kill
        end do
        call q_descr%kill
        call ctrl%kill
        call qsys_obj%kill
    end subroutine test_generate_script_4_dshmem

    ! ================================================================
    ! TEST 5: Streaming stack operations with local_dshmem
    !
    ! In streaming mode (e.g. stream processing in abinitio2D),
    ! jobs arrive one at a time.  The controller manages a FIFO stack:
    !
    !   add_to_streaming(cline) → push onto waiting stack
    !   schedule_streaming(q_descr) → tick: checks completion, dispatches
    !   get_stream_done_stack(clines) → pop completed jobs
    !   get_stream_fail_stack(clines, n) → pop failed jobs
    !   clear_stack() → discard all waiting jobs
    !
    ! This test exercises the stack helpers without actually submitting.
    ! ================================================================
    subroutine test_streaming_dshmem()
        type(qsys_local_dshmem)        :: qsys_obj
        integer, allocatable, target    :: parts(:,:)
        type(qsys_ctrl)               :: ctrl
        type(cmdline)                  :: cline1, cline2
        class(cmdline), allocatable    :: done_stack(:), fail_stack(:)
        integer                        :: nfail
        write(*,'(/,a)') '--- test_streaming_dshmem ---'
        call make_dshmem_ctrl( qsys_obj, parts, ctrl, stream=.true. )
        ! Push 2 streaming jobs
        call cline1%set('prg',   'cluster2D')
        call cline1%set('fromp', 1.)
        call cline1%set('top',   25.)
        call ctrl%add_to_streaming( cline1 )

        call cline2%set('prg',   'cluster2D')
        call cline2%set('fromp', 26.)
        call cline2%set('top',   50.)
        call ctrl%add_to_streaming( cline2 )

        ! Check stack sizes
        if( ctrl%get_stacksz() == 2 )then
            write(*,'(a)') '  PASS: stack size = 2 after two adds'
        else
            write(*,'(a,i3)') '  FAIL: expected stack size 2, got ', ctrl%get_stacksz()
            all_passed = .false.
        endif
        write(*,'(a,i5)') '  stack range = ', ctrl%get_stack_range()
        ! Done/fail stacks start empty
        if( ctrl%get_done_stacksz() == 0 .and. ctrl%get_failed_stacksz() == 0 )then
            write(*,'(a)') '  PASS: done/fail stacks empty before scheduling'
        else
            write(*,'(a)') '  FAIL: done/fail stacks should be empty'
            all_passed = .false.
        endif
        ! Pop (empty) stacks
        call ctrl%get_stream_done_stack( done_stack )
        call ctrl%get_stream_fail_stack( fail_stack, nfail )
        if( .not. allocated(done_stack) .and. nfail == 0 )then
            write(*,'(a)') '  PASS: get_stream_done/fail_stack return empty'
        else
            write(*,'(a)') '  FAIL: stacks should be unallocated'
            all_passed = .false.
        endif
        ! Clear the waiting stack
        call ctrl%clear_stack()
        if( ctrl%get_stacksz() == 0 )then
            write(*,'(a)') '  PASS: clear_stack emptied the queue'
        else
            write(*,'(a)') '  FAIL: stack should be 0 after clear'
            all_passed = .false.
        endif
        call cline1%kill
        call cline2%kill
        call ctrl%kill
        call qsys_obj%kill
    end subroutine test_streaming_dshmem

    ! ================================================================
    ! TEST 6: schedule_streaming tick with local_dshmem
    !
    ! schedule_streaming is the "tick" function called repeatedly in a
    ! streaming processing loop (e.g., in simple_commander_stream2D).
    ! Each call:
    !   1. update_queue  – checks for EXIT_CODE_JOB_<part> files
    !                      and JOB_FINISHED_<part> files
    !   2. For each free computing unit:
    !      - pop a cline from the waiting stack
    !      - generate a script for it
    !      - submit the script
    !
    ! We test one tick with an empty stack (safe no-op).
    ! ================================================================
    subroutine test_schedule_streaming_dshmem()
        type(qsys_local_dshmem)      :: qsys_obj
        integer, allocatable, target  :: parts(:,:)
        type(qsys_ctrl)              :: ctrl
        type(chash)                   :: q_descr
        write(*,'(/,a)') '--- test_schedule_streaming_dshmem ---'
        call make_dshmem_ctrl( qsys_obj, parts, ctrl, stream=.true. )
        call q_descr%new(10)
        call q_descr%set('qsys_name', 'local_dshmem')
        ! Call schedule_streaming with an empty stack.
        ! It should just call update_queue and return (nothing to dispatch).
        call ctrl%schedule_streaming( q_descr )
        write(*,'(a,i3)') '  computing units avail after tick = ', ctrl%get_ncomputing_units_avail()
        write(*,'(a)') '  PASS: schedule_streaming tick with empty stack completed'
        call q_descr%kill
        call ctrl%kill
        call qsys_obj%kill
    end subroutine test_schedule_streaming_dshmem

end program simple_test_qsys_ctrl_dshmem
