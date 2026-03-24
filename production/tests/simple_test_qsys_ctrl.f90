!> @file   simple_test_qsys_ctrl.f90
!> @brief  Standalone test program exercising every public routine in simple_qsys_ctrl.
!>
!> This test does NOT submit real jobs.  It creates a qsys_ctrl object backed by
!> the lightweight qsys_local driver, generates scripts on disk, and exercises
!> every getter / setter / streaming helper so you can see how the pieces fit
!> together.
!>
!> Build with BUILD_TESTS=ON and run:
!>     ./simple_test_qsys_ctrl
!> The test creates (and cleans up) temporary script and sentinel files in the
!> current working directory.
program simple_test_qsys_ctrl
use simple_core_module_api
use simple_qsys_local, only: qsys_local
use simple_qsys_ctrl,  only: qsys_ctrl
use simple_cmdline,    only: cmdline
implicit none
#include "simple_local_flags.inc"

integer, parameter :: NPARTS         = 4    ! number of partitions
integer, parameter :: NPTCLS         = 100  ! fake particle count
integer, parameter :: NCUNITS        = 2    ! computing units (< nparts to see scheduling effect)
logical            :: all_passed     = .true.

! ---------- run all sub-tests ----------
call test_constructor_and_exists()
call test_get_set_jobs_status()
call test_free_all_cunits()
call test_generate_scripts()
call test_generate_script_2()
call test_generate_script_3()
call test_generate_script_4()
call test_streaming_workflow()
call test_kill()

! ---------- final verdict ----------
if( all_passed )then
    write(*,'(a)') '=== ALL qsys_ctrl TESTS PASSED ==='
else
    write(*,'(a)') '*** SOME qsys_ctrl TESTS FAILED ***'
    stop 1
endif

contains

    ! ----------------------------------------------------------------
    ! Utility: create a local qsys, an even parts array, and a qsys_ctrl
    ! ----------------------------------------------------------------
    subroutine make_ctrl( qsys_obj, parts, ctrl, stream )
        type(qsys_local),  intent(out), target :: qsys_obj
        integer, allocatable, intent(out), target :: parts(:,:)
        type(qsys_ctrl),   intent(out) :: ctrl
        logical, optional,  intent(in)  :: stream
        logical :: sstream
        integer :: ip
        sstream = .false.
        if( present(stream) ) sstream = stream
        ! initialise the local qsys (it only needs new)
        call qsys_obj%new()
        ! build a simple even-split parts table  [fromp, top]
        parts = split_nobjs_even(NPTCLS, NPARTS)
        write(*,'(a)') '  Parts table (fromp , top):'
        do ip = 1, NPARTS
            write(*,'(a,i2,a,i4,a,i4)') '    part ', ip, ':  ', parts(ip,1), ' -> ', parts(ip,2)
        end do
        ! construct qsys_ctrl
        !   exec_binary  – the path of the binary that each script would invoke
        !   qsys_obj     – polymorphic qsys backend (local here)
        !   parts        – the (nparts,2) fromp/top table
        !   fromto_part  – range of partition indices this controller manages [1,NPARTS]
        !   ncunits      – maximum concurrent computing units
        !   stream       – streaming mode flag
        call ctrl%new( string('simple_private_exec'), qsys_obj, parts, &
                       [1, NPARTS], NCUNITS, sstream )
    end subroutine make_ctrl

    ! ----------------------------------------------------------------
    ! TEST: constructor + exists
    ! ----------------------------------------------------------------
    subroutine test_constructor_and_exists()
        type(qsys_local)        :: qsys_obj
        integer, allocatable, target :: parts(:,:)
        type(qsys_ctrl)         :: ctrl
        write(*,'(/,a)') '--- test_constructor_and_exists ---'
        call make_ctrl( qsys_obj, parts, ctrl, stream=.false. )
        if( ctrl%exists() )then
            write(*,'(a)') '  PASS: ctrl%exists() returns .true. after new()'
        else
            write(*,'(a)') '  FAIL: ctrl%exists() should be .true.'
            all_passed = .false.
        endif
        call ctrl%kill
        call qsys_obj%kill
    end subroutine test_constructor_and_exists

    ! ----------------------------------------------------------------
    ! TEST: get_jobs_status / set_jobs_status / print_jobs_status
    !
    !  get_jobs_status(jobs_done, jobs_submitted)
    !      Returns allocatable copies of the internal done / submitted
    !      boolean arrays (one element per partition).
    !
    !  set_jobs_status(jobs_done, jobs_submitted)
    !      Overwrites the internal arrays with the caller's values.
    !      Useful for resuming from a checkpoint.
    !
    !  print_jobs_status()
    !      Writes a human-readable table to logfhandle.
    ! ----------------------------------------------------------------
    subroutine test_get_set_jobs_status()
        type(qsys_local)            :: qsys_obj
        integer, allocatable, target :: parts(:,:)
        type(qsys_ctrl)             :: ctrl
        logical, allocatable         :: done(:), submitted(:)
        write(*,'(/,a)') '--- test_get_set_jobs_status ---'
        call make_ctrl( qsys_obj, parts, ctrl )
        ! After construction every job is not-done and not-submitted
        call ctrl%get_jobs_status( done, submitted )
        write(*,'(a,i2)') '  #done      = ', count(done)
        write(*,'(a,i2)') '  #submitted = ', count(submitted)
        if( count(done) /= 0 .and. count(submitted) /= 0 )then
            write(*,'(a)') '  FAIL: fresh ctrl should have 0 done, 0 submitted'
            all_passed = .false.
        else
            write(*,'(a)') '  PASS: initial status is all-false'
        endif
        ! Now pretend parts 1 and 3 finished
        done(1) = .true.;  submitted(1) = .true.
        done(3) = .true.;  submitted(3) = .true.
        call ctrl%set_jobs_status( done, submitted )
        ! Re-read and verify
        deallocate(done, submitted)
        call ctrl%get_jobs_status( done, submitted )
        if( count(done) == 2 .and. count(submitted) == 2 )then
            write(*,'(a)') '  PASS: set/get round-trip works'
        else
            write(*,'(a)') '  FAIL: set/get mismatch'
            all_passed = .false.
        endif
        ! print_jobs_status writes to logfhandle (stdout by default)
        write(*,'(a)') '  Calling print_jobs_status:'
        call ctrl%print_jobs_status()
        call ctrl%kill
        call qsys_obj%kill
    end subroutine test_get_set_jobs_status

    ! ----------------------------------------------------------------
    ! TEST: free_all_cunits
    !
    !  Resets ncomputing_units_avail back to ncomputing_units.
    !  Useful after forcibly reclaiming all slots (e.g. error recovery).
    ! ----------------------------------------------------------------
    subroutine test_free_all_cunits()
        type(qsys_local)            :: qsys_obj
        integer, allocatable, target :: parts(:,:)
        type(qsys_ctrl)             :: ctrl
        write(*,'(/,a)') '--- test_free_all_cunits ---'
        call make_ctrl( qsys_obj, parts, ctrl )
        ! After new(), available units == NCUNITS.
        ! free_all_cunits explicitly resets available units to NCUNITS.
        call ctrl%free_all_cunits()
        write(*,'(a)') '  PASS: free_all_cunits completed without error'
        call ctrl%kill
        call qsys_obj%kill
    end subroutine test_free_all_cunits

    ! ----------------------------------------------------------------
    ! TEST: generate_scripts  (the main bulk script generator)
    !
    !  For each partition in [fromto_part(1) .. fromto_part(2)] it:
    !    1.  Injects  fromp, top, part, nparts  into job_descr
    !    2.  Optionally sets outfile = outfile_body // padded_part // '.mrc'
    !    3.  Optionally merges per-part parameters from part_params(:)
    !    4.  Calls generate_script_1 which writes the bash script
    !    5.  Resets ncomputing_units_avail (for non-stream mode)
    !
    !  The generated scripts live in the current directory as
    !     distr_simple_script_<padded_part>
    ! ----------------------------------------------------------------
    subroutine test_generate_scripts()
        type(qsys_local)            :: qsys_obj
        integer, allocatable, target :: parts(:,:)
        type(qsys_ctrl)             :: ctrl
        type(chash)                  :: job_descr, q_descr
        integer                      :: ip
        write(*,'(/,a)') '--- test_generate_scripts ---'
        call make_ctrl( qsys_obj, parts, ctrl )
        ! Build a minimal job description  (key-value pairs that become
        ! the command line of the launched binary)
        call job_descr%new(20)
        call job_descr%set('prg',     'cluster2D')
        call job_descr%set('nthr',    '1')
        call job_descr%set('ncls',    '50')
        ! Build a minimal queue description
        call q_descr%new(10)
        call q_descr%set('qsys_name', 'local')
        ! Generate one script per partition
        call ctrl%generate_scripts( job_descr, string('.mrc'), q_descr )
        ! Verify the script files were created on disk
        write(*,'(a)') '  Generated scripts:'
        do ip = 1, NPARTS
            if( file_exists('distr_simple_script_'//int2str(ip)) )then
                write(*,'(a,a)') '    found: distr_simple_script_', trim(int2str(ip))
            else
                write(*,'(a,a)') '    MISSING: distr_simple_script_', trim(int2str(ip))
                all_passed = .false.
            endif
        end do
        write(*,'(a)') '  PASS: generate_scripts created all script files'
        ! clean up files
        do ip = 1, NPARTS
            call del_file('distr_simple_script_'//int2str(ip))
        end do
        call job_descr%kill
        call q_descr%kill
        call ctrl%kill
        call qsys_obj%kill
    end subroutine test_generate_scripts

    ! ----------------------------------------------------------------
    ! TEST: generate_script_2  (single-job, explicit params overload)
    !
    !  Writes a single bash script that invokes exec_bin with the
    !  key-value pairs in job_descr.  Output is directed to outfile
    !  (if given) or to the global SIMPLE_SUBPROC_OUT.
    !
    !  Signature:
    !    generate_script(job_descr, q_descr, exec_bin, script_name
    !                    [, outfile] [, exit_code_fname])
    ! ----------------------------------------------------------------
    subroutine test_generate_script_2()
        type(qsys_local)            :: qsys_obj
        integer, allocatable, target :: parts(:,:)
        type(qsys_ctrl)             :: ctrl
        type(chash)                  :: job_descr, q_descr
        character(len=*), parameter  :: SNAME = 'test_single_script_2'
        character(len=*), parameter  :: ONAME = 'test_single_output_2.log'
        character(len=*), parameter  :: ENAME = 'test_exit_code_2'
        write(*,'(/,a)') '--- test_generate_script_2 ---'
        call make_ctrl( qsys_obj, parts, ctrl )
        call job_descr%new(10)
        call job_descr%set('prg',  'cluster2D')
        call job_descr%set('nthr', '2')
        call q_descr%new(10)
        call q_descr%set('qsys_name', 'local')
        ! generate_script overload 2: explicit exec binary + script name
        call ctrl%generate_script( job_descr, q_descr, string('simple_private_exec'), &
                                   string(SNAME), outfile=string(ONAME), &
                                   exit_code_fname=string(ENAME) )
        if( file_exists(SNAME) )then
            write(*,'(a)') '  PASS: script file created ('//SNAME//')'
        else
            write(*,'(a)') '  FAIL: script file NOT found'
            all_passed = .false.
        endif
        call del_file(SNAME)
        call job_descr%kill
        call q_descr%kill
        call ctrl%kill
        call qsys_obj%kill
    end subroutine test_generate_script_2

    ! ----------------------------------------------------------------
    ! TEST: generate_script_3  (single-job from cmdline object)
    !
    !  Converts a cmdline object into a job description (chash) and
    !  writes a bash script.  Useful when you already have a cmdline
    !  rather than a chash.
    !
    !  Signature:
    !    generate_script(cline, q_descr, script_name, prgoutput)
    ! ----------------------------------------------------------------
    subroutine test_generate_script_3()
        type(qsys_local)            :: qsys_obj
        integer, allocatable, target :: parts(:,:)
        type(qsys_ctrl)             :: ctrl
        type(cmdline)                :: cline
        type(chash)                  :: q_descr
        character(len=*), parameter  :: SNAME = 'test_single_script_3'
        character(len=*), parameter  :: POUT  = 'test_prg_output_3.log'
        write(*,'(/,a)') '--- test_generate_script_3 ---'
        call make_ctrl( qsys_obj, parts, ctrl )
        ! Build a cmdline (in-memory command line object)
        call cline%set('prg',  'cluster2D')
        call cline%set('nthr', 4.)
        call cline%set('ncls', 50.)
        call q_descr%new(10)
        call q_descr%set('qsys_name', 'local')
        ! generate_script overload 3: from cmdline
        call ctrl%generate_script( cline, q_descr, string(SNAME), string(POUT) )
        if( file_exists(SNAME) )then
            write(*,'(a)') '  PASS: script file created ('//SNAME//')'
        else
            write(*,'(a)') '  FAIL: script file NOT found'
            all_passed = .false.
        endif
        call del_file(SNAME)
        call q_descr%kill
        call cline%kill
        call ctrl%kill
        call qsys_obj%kill
    end subroutine test_generate_script_3

    ! ----------------------------------------------------------------
    ! TEST: generate_script_4  (multi-job sequential script)
    !
    !  Packs multiple job descriptions into a SINGLE bash script that
    !  runs each job sequentially.  Useful for chaining dependencies.
    !
    !  Signature:
    !    generate_script(jobs_descr(:), q_descr, exec_bin,
    !                    script_name, outfile [, exec_bins])
    ! ----------------------------------------------------------------
    subroutine test_generate_script_4()
        type(qsys_local)            :: qsys_obj
        integer, allocatable, target :: parts(:,:)
        type(qsys_ctrl)             :: ctrl
        type(chash), allocatable     :: jobs(:)
        type(chash)                  :: q_descr
        character(len=*), parameter  :: SNAME = 'test_multi_script_4'
        character(len=*), parameter  :: ONAME = 'test_multi_output_4.log'
        integer :: ij
        write(*,'(/,a)') '--- test_generate_script_4 ---'
        call make_ctrl( qsys_obj, parts, ctrl )
        ! Create 3 sequential job descriptions
        allocate(jobs(3))
        do ij = 1, 3
            call jobs(ij)%new(10)
            call jobs(ij)%set('prg',  'cluster2D')
            call jobs(ij)%set('ncls', int2str(ij*50))
        end do
        call q_descr%new(10)
        call q_descr%set('qsys_name', 'local')
        ! generate_script overload 4: multiple jobs → 1 script
        call ctrl%generate_script( jobs, q_descr, string('simple_private_exec'), &
                                   string(SNAME), string(ONAME) )
        if( file_exists(SNAME) )then
            write(*,'(a)') '  PASS: multi-job script file created ('//SNAME//')'
        else
            write(*,'(a)') '  FAIL: multi-job script file NOT found'
            all_passed = .false.
        endif
        call del_file(SNAME)
        do ij = 1, 3
            call jobs(ij)%kill
        end do
        call q_descr%kill
        call ctrl%kill
        call qsys_obj%kill
    end subroutine test_generate_script_4

    ! ----------------------------------------------------------------
    ! TEST: streaming workflow
    !
    !  In streaming mode the controller acts as a job pool:
    !    1. add_to_streaming(cline)   – push a cmdline onto the stack
    !    2. get_stacksz()             – how many jobs are waiting
    !    3. get_stack_range()          – total particle range in stack
    !    4. get_done_stacksz()        – how many jobs finished successfully
    !    5. get_failed_stacksz()      – how many jobs failed
    !    6. get_ncomputing_units_avail() – available compute slots
    !    7. get_stream_done_stack()   – pop all completed cmdlines
    !    8. get_stream_fail_stack()   – pop all failed cmdlines
    !    9. clear_stack()             – discard all waiting jobs
    !
    !  schedule_streaming(q_descr [,path])  is the tick function:
    !    - checks which running jobs have finished (via sentinel files)
    !    - hands waiting jobs to free computing units
    !    - generates & submits scripts automatically
    !
    !  Here we only exercise the stack helpers without submitting.
    ! ----------------------------------------------------------------
    subroutine test_streaming_workflow()
        type(qsys_local)              :: qsys_obj
        integer, allocatable, target   :: parts(:,:)
        type(qsys_ctrl)               :: ctrl
        type(cmdline)                  :: cline1, cline2, cline3
        class(cmdline), allocatable    :: done_stack(:), fail_stack(:)
        integer                        :: stacksz, done_sz, fail_sz, nfail
        write(*,'(/,a)') '--- test_streaming_workflow ---'
        ! Create in STREAMING mode (stream = .true.)
        call make_ctrl( qsys_obj, parts, ctrl, stream=.true. )
        ! ----- add_to_streaming -----
        ! Push 3 fake jobs onto the waiting stack
        call cline1%set('prg',  'cluster2D')
        call cline1%set('fromp', 1.)
        call cline1%set('top',   25.)
        call ctrl%add_to_streaming( cline1 )

        call cline2%set('prg',  'cluster2D')
        call cline2%set('fromp', 26.)
        call cline2%set('top',   50.)
        call ctrl%add_to_streaming( cline2 )

        call cline3%set('prg',  'cluster2D')
        call cline3%set('fromp', 51.)
        call cline3%set('top',   75.)
        call ctrl%add_to_streaming( cline3 )

        ! ----- get_stacksz -----
        stacksz = ctrl%get_stacksz()
        write(*,'(a,i3)') '  stack size after 3 adds   = ', stacksz
        if( stacksz == 3 )then
            write(*,'(a)') '  PASS: get_stacksz'
        else
            write(*,'(a)') '  FAIL: expected stack size 3'
            all_passed = .false.
        endif

        ! ----- get_stack_range -----
        ! range = sum(top) - sum(fromp) + stacksz = (25+50+75)-(1+26+51)+3 = 75
        write(*,'(a,i5)') '  stack range               = ', ctrl%get_stack_range()

        ! ----- get_done_stacksz / get_failed_stacksz -----
        done_sz = ctrl%get_done_stacksz()
        fail_sz = ctrl%get_failed_stacksz()
        write(*,'(a,i3)') '  done  stack size          = ', done_sz
        write(*,'(a,i3)') '  fail  stack size          = ', fail_sz
        if( done_sz == 0 .and. fail_sz == 0 )then
            write(*,'(a)') '  PASS: done/fail stacks empty before any scheduling'
        else
            write(*,'(a)') '  FAIL: unexpected done/fail stack sizes'
            all_passed = .false.
        endif

        ! ----- get_ncomputing_units_avail -----
        write(*,'(a,i3)') '  computing units avail     = ', ctrl%get_ncomputing_units_avail()

        ! ----- get_stream_done_stack / get_stream_fail_stack -----
        ! (both are empty so the arrays come back unallocated)
        call ctrl%get_stream_done_stack( done_stack )
        call ctrl%get_stream_fail_stack( fail_stack, nfail )
        if( .not. allocated(done_stack) .and. nfail == 0 )then
            write(*,'(a)') '  PASS: done/fail stacks correctly empty'
        else
            write(*,'(a)') '  FAIL: done/fail stacks should be unallocated'
            all_passed = .false.
        endif

        ! ----- clear_stack -----
        call ctrl%clear_stack()
        stacksz = ctrl%get_stacksz()
        write(*,'(a,i3)') '  stack size after clear     = ', stacksz
        if( stacksz == 0 )then
            write(*,'(a)') '  PASS: clear_stack'
        else
            write(*,'(a)') '  FAIL: stack should be 0 after clear'
            all_passed = .false.
        endif

        call cline1%kill
        call cline2%kill
        call cline3%kill
        call ctrl%kill
        call qsys_obj%kill
    end subroutine test_streaming_workflow

    ! ----------------------------------------------------------------
    ! TEST: kill  (destructor)
    !
    !  Deallocates everything and sets existence = .false.
    !  Calling kill on a dead object is safe (no-op).
    ! ----------------------------------------------------------------
    subroutine test_kill()
        type(qsys_local)            :: qsys_obj
        integer, allocatable, target :: parts(:,:)
        type(qsys_ctrl)             :: ctrl
        write(*,'(/,a)') '--- test_kill ---'
        call make_ctrl( qsys_obj, parts, ctrl )
        call ctrl%kill
        if( .not. ctrl%exists() )then
            write(*,'(a)') '  PASS: ctrl%exists() is .false. after kill'
        else
            write(*,'(a)') '  FAIL: existence should be .false.'
            all_passed = .false.
        endif
        ! second kill is safe
        call ctrl%kill
        write(*,'(a)') '  PASS: double-kill is safe'
        call qsys_obj%kill
    end subroutine test_kill

end program simple_test_qsys_ctrl
