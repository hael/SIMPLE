!@descr: unit tests for the queue-system environment's installation-path policy and particle partitions (simple_qsys_env, simple_map_reduce, sp_project%update_compenv)
! The installation path is runtime-local: a distributed run takes the environment's SIMPLE_PATH
! (set by CTest) for its queue description and executable, never a simple_path a project carries; an empty project gets the minimum job time. Partitions given an active-particle mask
! stay contiguous and balance the active particles (split_nobjs_active).
module simple_qsys_env_tester
use simple_core_module_api
use simple_cmdline,    only: cmdline
use simple_parameters, only: parameters
use simple_qsys_env,   only: qsys_env
use simple_sp_project, only: sp_project
use simple_test_utils
implicit none
private
public :: run_all_qsys_env_tests

character(len=*), parameter :: FOREIGN_SIMPLE_PATH = '/remote/install/that/is/not/local'

contains

    subroutine run_all_qsys_env_tests()
        write(*,'(A)') '**** running all qsys environment tests ****'
        call test_installation_path_policy()
        call test_active_balanced_split()
        call test_active_balanced_qsys_parts()
    end subroutine run_all_qsys_env_tests

    subroutine test_installation_path_policy()
        type(cmdline)    :: cline
        type(parameters) :: params
        type(qsys_env)   :: qenv
        type(sp_project) :: project
        type(string)     :: exec_bin, expected_exec_bin, job_time, local_simple_path, projfile
        integer          :: iostat
        write(*,'(A)') 'test_installation_path_policy'
        projfile = 'test_qsys_env_path_policy.simple'
        call del_file(projfile)
        ! a project carrying a foreign installation path: distributed execution must ignore it
        call project%projinfo%new(1, is_ptcl=.false.)
        call project%projinfo%set(1, 'projname', 'test_qsys_env_path_policy')
        call project%compenv%new(1, is_ptcl=.false.)
        call cline%set('qsys_name', 'local')
        call project%update_compenv(cline)
        call project%compenv%set(1, 'simple_path', FOREIGN_SIMPLE_PATH)
        call project%write(projfile)
        local_simple_path = simple_getenv('SIMPLE_PATH', iostat)
        call assert_int(0, iostat, 'SIMPLE_PATH is set in the test environment')
        if( iostat == 0 )then
            params%prg             = 'qsys_env_path_policy'
            params%projfile        = projfile
            params%qsys_name       = 'local'
            params%nparts          = 1
            params%ncunits         = 1
            params%nptcls          = 0
            params%nthr            = 1
            params%worker_priority = 'normal'
            params%worker_server   = ''
            call qenv%new(params, 1, exec_bin=string('simple_exec'))
            exec_bin          = qenv%get_exec_bin()
            expected_exec_bin = filepath(local_simple_path, string('bin'), string('simple_exec'))
            call assert_string_eq(expected_exec_bin%to_char(), exec_bin, &
                &'the executable resolves from the local SIMPLE_PATH, not the project')
            call assert_true(qenv%qdescr%isthere('simple_path'), 'the queue description carries a simple_path')
            if( qenv%qdescr%isthere('simple_path') )then
                call assert_string_eq(local_simple_path%to_char(), qenv%qdescr%get('simple_path'), &
                    &'the queue description carries the local SIMPLE_PATH, not the project one')
            endif
            call assert_true(qenv%qdescr%isthere('job_time'), 'the queue description carries a derived job time')
            if( qenv%qdescr%isthere('job_time') )then
                job_time = qenv%qdescr%get('job_time')
                call assert_string_eq('0-0:1:40', job_time, 'an empty project gets the minimum job time')
            endif
            call qenv%kill
        endif
        call project%kill
        call cline%kill
        call del_file(projfile)
    end subroutine test_installation_path_policy

    ! split_nobjs_active on the layouts the distributed 3D workflows meet
    subroutine test_active_balanced_split()
        logical, allocatable :: l_active(:)
        integer, allocatable :: parts(:,:), even(:,:)
        integer              :: i, szmax
        write(*,'(A)') 'test_active_balanced_split'
        ! a stream add-on: 49524 masked frozen rows, then 16509 new ones; the shares of
        ! split_nobjs_even(16509,4) are 4128, 4127, 4127 and 4127
        allocate(l_active(66033), source=.false.)
        l_active(49525:) = .true.
        parts = split_nobjs_active(l_active, 4, szmax)
        call assert_true(all(parts(:,1) == [1, 53653, 57780, 61907]),     'stream layout: part starts')
        call assert_true(all(parts(:,2) == [53652, 57779, 61906, 66033]), 'stream layout: part ends')
        call assert_int(4128, szmax, 'stream layout: largest active share')
        call assert_true(all(active_counts(l_active, parts) == [4128, 4127, 4127, 4127]), &
            &'stream layout: every part holds its share of the new particles')
        deallocate(l_active)
        ! interleaved (a permissive selection over a harsh one): rows 3, 6, ..., 39 active;
        ! shares 4, 3, 3, 3 end on the 4th, 7th and 10th active rows (12, 21, 30)
        allocate(l_active(40))
        l_active = [(mod(i,3) == 0, i=1,40)]
        parts = split_nobjs_active(l_active, 4, szmax)
        call assert_true(all(parts(:,1) == [1, 13, 22, 31]), 'interleaved: part starts')
        call assert_true(all(parts(:,2) == [12, 21, 30, 40]), 'interleaved: part ends')
        call assert_true(all(active_counts(l_active, parts) == [4, 3, 3, 3]), 'interleaved: active shares')
        call assert_int(4, szmax, 'interleaved: largest active share')
        deallocate(l_active)
        ! every row active: the even split
        allocate(l_active(100), source=.true.)
        parts = split_nobjs_active(l_active, 7)
        even  = split_nobjs_even(100, 7)
        call assert_true(all(parts == even), 'all active: identical to split_nobjs_even')
        deallocate(l_active)
        ! fewer active rows than parts (rows 9 and 10 of 10, four parts): every part keeps a row,
        ! the first ending where the three later parts still get one each
        allocate(l_active(10), source=.false.)
        l_active(9:10) = .true.
        parts = split_nobjs_active(l_active, 4, szmax)
        call assert_true(all(parts(:,1) == [1, 8, 9, 10]), 'few active: part starts')
        call assert_true(all(parts(:,2) == [7, 8, 9, 10]), 'few active: part ends')
        call assert_int(1, szmax, 'few active: largest active share')
        ! no active rows: the even split, nothing to do anywhere
        l_active = .false.
        parts = split_nobjs_active(l_active, 4, szmax)
        even  = split_nobjs_even(10, 4)
        call assert_true(all(parts == even), 'none active: identical to split_nobjs_even')
        call assert_int(0, szmax, 'none active: no active share')
    end subroutine test_active_balanced_split

    ! qsys_env%new hands the mask to the partitioning, so the job ranges are the balanced ones
    subroutine test_active_balanced_qsys_parts()
        type(cmdline)        :: cline
        type(parameters)     :: params
        type(qsys_env)       :: qenv
        type(sp_project)     :: project
        type(string)         :: projfile, local_simple_path
        logical, allocatable :: l_active(:)
        integer, allocatable :: expected(:,:)
        integer              :: iostat
        write(*,'(A)') 'test_active_balanced_qsys_parts'
        local_simple_path = simple_getenv('SIMPLE_PATH', iostat)
        call assert_int(0, iostat, 'SIMPLE_PATH is set in the test environment')
        if( iostat /= 0 ) return
        projfile = 'test_qsys_env_parts.simple'
        call del_file(projfile)
        call project%projinfo%new(1, is_ptcl=.false.)
        call project%projinfo%set(1, 'projname', 'test_qsys_env_parts')
        call project%compenv%new(1, is_ptcl=.false.)
        call cline%set('qsys_name', 'local')
        call project%update_compenv(cline)
        call project%write(projfile)
        params%prg             = 'qsys_env_parts'
        params%projfile        = projfile
        params%qsys_name       = 'local'
        params%nparts          = 3
        params%ncunits         = 3
        params%nptcls          = 30
        params%nthr            = 1
        params%worker_priority = 'normal'
        params%worker_server   = ''
        allocate(l_active(30), source=.false.)
        l_active(21:30) = .true.
        expected = split_nobjs_active(l_active, 3)
        call qenv%new(params, 3, exec_bin=string('simple_exec'), l_active=l_active)
        call assert_true(allocated(qenv%parts), 'the environment holds a partition table')
        if( allocated(qenv%parts) )then
            call assert_true(all(qenv%parts == expected), 'the job ranges are the active-balanced partitions')
            call assert_true(all(active_counts(l_active, qenv%parts) == [4, 3, 3]), &
                &'each job gets a share of the active particles')
        endif
        call qenv%kill
        call project%kill
        call cline%kill
        call del_file(projfile)
    end subroutine test_active_balanced_qsys_parts

    function active_counts( l_active, parts ) result( counts )
        logical, intent(in)  :: l_active(:)
        integer, intent(in)  :: parts(:,:)
        integer, allocatable :: counts(:)
        integer              :: ipart
        allocate(counts(size(parts,1)))
        do ipart = 1, size(parts,1)
            counts(ipart) = count(l_active(parts(ipart,1):parts(ipart,2)))
        end do
    end function active_counts

end module simple_qsys_env_tester
