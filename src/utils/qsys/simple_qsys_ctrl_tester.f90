!@descr: unit tests for the job controller of the queue systems (simple_qsys_ctrl) on the local backend
! Nothing is submitted: the controller is built over the local backend with four partitions of
! 100 particles and two computing units, and the tests read its state and the scripts it writes
! (partition scripts carrying fromp/top/part/nparts, single-job scripts from a chash or a
! cmdline, a multi-job script), the job-status round trip, the computing-unit bookkeeping and
! the streaming stack (size, particle range, done/fail stacks, clear). Scripts are written in
! the suite's directory and deleted.
module simple_qsys_ctrl_tester
use simple_core_module_api
use simple_qsys_local, only: qsys_local
use simple_qsys_ctrl,  only: qsys_ctrl
use simple_cmdline,    only: cmdline
use simple_test_utils
implicit none
private
public :: run_all_qsys_ctrl_tests

integer, parameter :: NPARTS  = 4    ! partitions
integer, parameter :: NPTCLS  = 100  ! particles, split evenly: 1-25, 26-50, 51-75, 76-100
integer, parameter :: NCUNITS = 2    ! computing units (< NPARTS)

contains

    subroutine run_all_qsys_ctrl_tests()
        write(*,'(A)') '**** running all qsys control tests ****'
        call test_constructor_and_kill()
        call test_jobs_status_round_trip()
        call test_computing_units()
        call test_prep_part_jobs()
        call test_single_job_scripts()
        call test_multi_job_script()
        call test_streaming_stack()
    end subroutine run_all_qsys_ctrl_tests

    !> a local backend, the even parts table and a controller over all partitions
    subroutine make_ctrl( qsys_obj, parts, ctrl, stream )
        type(qsys_local),     intent(out), target :: qsys_obj
        integer, allocatable, intent(out), target :: parts(:,:)
        type(qsys_ctrl),      intent(out) :: ctrl
        logical, optional,    intent(in)  :: stream
        logical :: sstream
        sstream = .false.
        if( present(stream) ) sstream = stream
        call qsys_obj%new()
        parts = split_nobjs_even(NPTCLS, NPARTS)
        call ctrl%new(string('simple_private_exec'), qsys_obj, parts, [1, NPARTS], NCUNITS, sstream)
    end subroutine make_ctrl

    !> the whole text of a small file, lines joined by blanks ('' when it cannot be read)
    function file_text( fname ) result( txt )
        character(len=*), intent(in) :: fname
        character(len=:), allocatable :: txt
        character(len=4096) :: line
        integer :: funit, ios
        txt = ''
        open(newunit=funit, file=fname, status='old', action='read', iostat=ios)
        if( ios /= 0 ) return
        do
            read(funit, '(a)', iostat=ios) line
            if( ios /= 0 ) exit
            txt = txt//' '//trim(line)
        end do
        close(funit)
    end function file_text

    subroutine test_constructor_and_kill()
        type(qsys_local)             :: qsys_obj
        integer, allocatable, target :: parts(:,:)
        type(qsys_ctrl)              :: ctrl
        write(*,'(A)') 'test_constructor_and_kill'
        call make_ctrl(qsys_obj, parts, ctrl)
        call assert_true(ctrl%exists(), 'the controller exists after new')
        call ctrl%kill
        call assert_false(ctrl%exists(), 'the controller does not exist after kill')
        call ctrl%kill   ! a second kill is a no-op
        call assert_false(ctrl%exists(), 'a second kill is harmless')
        call qsys_obj%kill
    end subroutine test_constructor_and_kill

    !> a fresh controller has nothing done or submitted; set/get round-trips per partition
    subroutine test_jobs_status_round_trip()
        type(qsys_local)             :: qsys_obj
        integer, allocatable, target :: parts(:,:)
        type(qsys_ctrl)              :: ctrl
        logical, allocatable         :: done(:), submitted(:)
        write(*,'(A)') 'test_jobs_status_round_trip'
        call make_ctrl(qsys_obj, parts, ctrl)
        call ctrl%get_jobs_status(done, submitted)
        call assert_int(NPARTS, size(done), 'one job status per partition')
        call assert_true(count(done) == 0 .and. count(submitted) == 0, 'a fresh controller has nothing done or submitted')
        done(1) = .true.; submitted(1) = .true.
        done(3) = .true.; submitted(3) = .true.
        call ctrl%set_jobs_status(done, submitted)
        deallocate(done, submitted)
        call ctrl%get_jobs_status(done, submitted)
        call assert_true(all(done .eqv. [.true.,.false.,.true.,.false.]), 'done flags round-trip per partition')
        call assert_true(all(submitted .eqv. [.true.,.false.,.true.,.false.]), 'submitted flags round-trip per partition')
        call ctrl%print_jobs_status()
        call ctrl%kill
        call qsys_obj%kill
    end subroutine test_jobs_status_round_trip

    subroutine test_computing_units()
        type(qsys_local)             :: qsys_obj
        integer, allocatable, target :: parts(:,:)
        type(qsys_ctrl)              :: ctrl
        write(*,'(A)') 'test_computing_units'
        call make_ctrl(qsys_obj, parts, ctrl)
        call assert_int(NCUNITS, ctrl%get_ncomputing_units_avail(), 'every computing unit is free after new')
        call ctrl%free_all_cunits()
        call assert_int(NCUNITS, ctrl%get_ncomputing_units_avail(), 'free_all_cunits leaves every unit free')
        call ctrl%kill
        call qsys_obj%kill
    end subroutine test_computing_units

    !> one script per partition carrying its particle range; the job description is restored
    subroutine test_prep_part_jobs()
        type(qsys_local)             :: qsys_obj
        integer, allocatable, target :: parts(:,:)
        type(qsys_ctrl)              :: ctrl
        type(chash)                  :: job_descr, q_descr
        character(len=:), allocatable :: txt
        integer :: ip, nmissing
        write(*,'(A)') 'test_prep_part_jobs'
        call make_ctrl(qsys_obj, parts, ctrl)
        call job_descr%new(20)
        call job_descr%set('prg',  'cluster2D')
        call job_descr%set('nthr', '1')
        call job_descr%set('ncls', '50')
        call q_descr%new(10)
        call q_descr%set('qsys_name', 'local')
        call ctrl%prep_part_jobs(job_descr, string('.mrc'), q_descr)
        nmissing = 0
        do ip = 1, NPARTS
            if( .not. file_exists('distr_simple_script_'//int2str(ip)) ) nmissing = nmissing + 1
        end do
        call assert_int(0, nmissing, 'prep_part_jobs writes one script per partition')
        txt = file_text('distr_simple_script_2')
        call assert_true(index(txt, 'simple_private_exec') > 0 .and. index(txt, 'prg=cluster2D') > 0, &
            &'a partition script runs the executable with the job description')
        call assert_true(index(txt, 'fromp=26') > 0 .and. index(txt, 'top=50') > 0 .and. index(txt, 'part=2') > 0 &
            &.and. index(txt, 'nparts=4') > 0, 'the second partition script carries particles 26-50 of four parts')
        call assert_false(job_descr%isthere('part'), 'the job description is restored after the partitions')
        call assert_int(NCUNITS, ctrl%get_ncomputing_units_avail(), 'preparing partition jobs frees every computing unit')
        do ip = 1, NPARTS
            call del_file('distr_simple_script_'//int2str(ip))
        end do
        call job_descr%kill
        call q_descr%kill
        call ctrl%kill
        call qsys_obj%kill
    end subroutine test_prep_part_jobs

    !> a single-job script from a chash (explicit executable) and from a cmdline
    subroutine test_single_job_scripts()
        type(qsys_local)             :: qsys_obj
        integer, allocatable, target :: parts(:,:)
        type(qsys_ctrl)              :: ctrl
        type(chash)                  :: job_descr, q_descr
        type(cmdline)                :: cline
        character(len=*), parameter  :: SNAME2 = 'test_single_script_2', SNAME3 = 'test_single_script_3'
        write(*,'(A)') 'test_single_job_scripts'
        call make_ctrl(qsys_obj, parts, ctrl)
        call job_descr%new(10)
        call job_descr%set('prg',  'cluster2D')
        call job_descr%set('nthr', '2')
        call q_descr%new(10)
        call q_descr%set('qsys_name', 'local')
        call ctrl%generate_script(job_descr, q_descr, string('simple_private_exec'), string(SNAME2), &
            &outfile=string('test_single_output_2.log'), exit_code_fname=string('test_exit_code_2'))
        call assert_true(file_exists(SNAME2), 'a single-job script is written from a job description')
        call assert_true(index(file_text(SNAME2), 'prg=cluster2D') > 0, 'the single-job script carries the job description')
        call cline%set('prg',  'cluster2D')
        call cline%set('nthr', 4.)
        call cline%set('ncls', 50.)
        call ctrl%generate_script(cline, q_descr, string(SNAME3), string('test_prg_output_3.log'))
        call assert_true(file_exists(SNAME3), 'a single-job script is written from a command line')
        call assert_true(index(file_text(SNAME3), 'prg=cluster2D') > 0, 'the command-line script carries the program')
        call del_file(SNAME2)
        call del_file(SNAME3)
        call job_descr%kill
        call q_descr%kill
        call cline%kill
        call ctrl%kill
        call qsys_obj%kill
    end subroutine test_single_job_scripts

    !> three job descriptions packed into one sequential script
    subroutine test_multi_job_script()
        type(qsys_local)             :: qsys_obj
        integer, allocatable, target :: parts(:,:)
        type(qsys_ctrl)              :: ctrl
        type(chash), allocatable     :: jobs(:)
        type(chash)                  :: q_descr
        character(len=*), parameter  :: SNAME = 'test_multi_script_4'
        character(len=:), allocatable :: txt
        integer :: ij
        write(*,'(A)') 'test_multi_job_script'
        call make_ctrl(qsys_obj, parts, ctrl)
        allocate(jobs(3))
        do ij = 1, 3
            call jobs(ij)%new(10)
            call jobs(ij)%set('prg',  'cluster2D')
            call jobs(ij)%set('ncls', int2str(ij*50))
        end do
        call q_descr%new(10)
        call q_descr%set('qsys_name', 'local')
        call ctrl%generate_script(jobs, q_descr, string('simple_private_exec'), string(SNAME), string('test_multi_output_4.log'))
        call assert_true(file_exists(SNAME), 'a multi-job script is written')
        txt = file_text(SNAME)
        call assert_true(index(txt, 'ncls=50') > 0 .and. index(txt, 'ncls=100') > 0 .and. index(txt, 'ncls=150') > 0, &
            &'the multi-job script runs all three jobs')
        call del_file(SNAME)
        do ij = 1, 3
            call jobs(ij)%kill
        end do
        call q_descr%kill
        call ctrl%kill
        call qsys_obj%kill
    end subroutine test_multi_job_script

    !> the streaming job pool, exercised without scheduling
    subroutine test_streaming_stack()
        type(qsys_local)             :: qsys_obj
        integer, allocatable, target :: parts(:,:)
        type(qsys_ctrl)              :: ctrl
        type(cmdline)                :: cline1, cline2, cline3
        class(cmdline), allocatable  :: done_stack(:), fail_stack(:)
        integer :: nfail
        write(*,'(A)') 'test_streaming_stack'
        call make_ctrl(qsys_obj, parts, ctrl, stream=.true.)
        call cline1%set('prg', 'cluster2D'); call cline1%set('fromp',  1.); call cline1%set('top', 25.)
        call cline2%set('prg', 'cluster2D'); call cline2%set('fromp', 26.); call cline2%set('top', 50.)
        call cline3%set('prg', 'cluster2D'); call cline3%set('fromp', 51.); call cline3%set('top', 75.)
        call ctrl%add_to_streaming(cline1)
        call ctrl%add_to_streaming(cline2)
        call ctrl%add_to_streaming(cline3)
        call assert_int(3, ctrl%get_stacksz(), 'three jobs wait on the streaming stack')
        ! (25+50+75) - (1+26+51) + 3
        call assert_int(75, ctrl%get_stack_range(), 'the stack spans 75 particles')
        call assert_int(0, ctrl%get_done_stacksz(),   'no job is done before scheduling')
        call assert_int(0, ctrl%get_failed_stacksz(), 'no job has failed before scheduling')
        call assert_int(NCUNITS, ctrl%get_ncomputing_units_avail(), 'every computing unit is free before scheduling')
        call ctrl%get_stream_done_stack(done_stack)
        call ctrl%get_stream_fail_stack(fail_stack, nfail)
        call assert_true(.not. allocated(done_stack) .and. nfail == 0, 'the done and fail stacks come back empty')
        call ctrl%clear_stack()
        call assert_int(0, ctrl%get_stacksz(), 'clear_stack empties the waiting stack')
        call cline1%kill
        call cline2%kill
        call cline3%kill
        call ctrl%kill
        call qsys_obj%kill
    end subroutine test_streaming_stack

end module simple_qsys_ctrl_tester
