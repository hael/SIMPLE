!@descr: for all class tests
module simple_commanders_test_class
use simple_commanders_api
use simple_stream_api
implicit none
#include "simple_local_flags.inc"

type, extends(commander_base) :: commander_test_ui_hash_test
  contains
    procedure :: execute      => exec_test_ui_hash_test
end type commander_test_ui_hash_test

type, extends(commander_base) :: commander_test_units
  contains
    procedure :: execute      => exec_test_units
end type commander_test_units

contains

subroutine exec_test_ui_hash_test( self, cline )
    use simple_ui_hash, only: test_ui_hash
    class(commander_test_ui_hash_test),  intent(inout) :: self
    class(cmdline),                      intent(inout) :: cline
    call test_ui_hash
    call simple_end('**** SIMPLE_TEST_UI_HASH_TEST NORMAL STOP ****')
end subroutine exec_test_ui_hash_test

subroutine exec_test_units( self, cline )
    use simple_imghead, only: test_imghead
    use simple_oris,    only: test_oris
    use simple_atoms,   only: test_atoms
    ! core library tester modules generated with help from chatgpt
    use simple_test_utils
    use simple_string_tester
    use simple_syslib_tester
    use simple_fileio_tester
    use simple_chash_tester
    use simple_vrefhash_tester
    use simple_hash_tester
    use simple_linked_list_tester
    use simple_binary_tree_tester
    use simple_multi_dendro_tester
    use simple_cmdline_tester
    use simple_ori_tester
    use simple_oris_tester
    use simple_rec_list_tester
    use simple_persistent_worker_message_tester, only: run_all_persistent_worker_message_tests
    use simple_persistent_worker_server_tester,  only: run_all_persistent_worker_server_tests
    ! hand-written unit tests
    use simple_ipc_tcp_socket_tester, only: run_all_ipc_tcp_socket_tests
    use simple_forked_process_tester, only: run_all_forked_process_tests
    use simple_gui_metadata_tester,   only: run_all_gui_metadata_tests
    use simple_gui_assembler_tester,  only: run_all_gui_assembler_tests
    use simple_http_post_tester,      only: run_all_http_post_tests
    use simple_aff_prop,              only: test_aff_prop
    use simple_hclust,                only: test_hclust
    use simple_ftexp_shsrch,          only: test_ftexp_shsrch
    use simple_ftiter,                only: test_ftiter
    use simple_image,                 only: test_image
    use simple_online_var,            only: test_online_var
    use simple_ui,                    only: validate_ui_json
    use simple_starfile_tester,       only: run_all_starfile_tests
    use simple_project_merge_tester,  only: run_all_project_merge_tests
    use simple_motion_gain_tester,    only: run_all_motion_gain_tests
    use simple_bspline_smoother,      only: test_bspline_smoother, test_bspline_smoother_3d
    abstract interface
        subroutine no_arg_test()
        end subroutine no_arg_test
    end interface
    class(commander_test_units),  intent(inout) :: self
    class(cmdline),               intent(inout) :: cline
    character(8)          :: datestr
    character(len=STDLEN) :: folder
    logical               :: test_failed
    type(string)          :: original_cwd, report_file
    call seed_rnd
    call date_and_time(date=datestr)
    folder = 'SIMPLE_TEST_UNITS_'//datestr
    call simple_getcwd(original_cwd)
    report_file = original_cwd//'/'//trim(folder)//'/simple_test_units_report.txt'
    call simple_mkdir(folder)
    call simple_chdir(folder)
    call reset_test_report(report_file%to_char())
    ! core library tests generated with help from chatgpt
    call run_suite('string',                    run_all_string_tests)
    call run_suite('syslib',                    run_all_syslib_tests)
    call run_suite('fileio',                    run_all_fileio_tests)
    call run_suite('character hash',            run_all_chash_tests)
    call run_suite('hash',                      run_all_hash_tests)
    call run_suite('value-reference hash',      run_all_vrefhash_tests)
    call run_suite('linked list',               run_all_list_tests)
    call run_suite('binary tree',               run_all_tree_tests)
    call run_suite('multi dendrogram',           run_all_multi_dendro_tests)
    call run_suite('command line',               run_all_cmdline_tests)
    call run_suite('orientation',                run_all_ori_tests)
    call run_suite('orientation collection',     run_all_oris_tests)
    call run_suite('record list',                run_all_rec_list_tests)
    call run_suite('IPC TCP socket',             run_all_ipc_tcp_socket_tests)
    call run_suite('forked process',             run_all_forked_process_tests)
    call run_suite('GUI metadata',               run_all_gui_metadata_tests)
    call run_suite('GUI assembler',              run_all_gui_assembler_tests)
    call run_suite('HTTP POST',                  run_all_http_post_tests)
    call run_suite('persistent worker message',  run_all_persistent_worker_message_tests)
    call run_suite('persistent worker server',   run_all_persistent_worker_server_tests)
    call run_suite('STAR file',                  run_all_starfile_tests)
    call run_suite('project merge',              run_all_project_merge_tests)
    call run_suite('motion gain',                run_all_motion_gain_tests)
    ! hand-written unit tests
    call begin_test_suite('UI JSON')
    write(*,*)'VALIDATING UI JSON FILE:'
    call validate_ui_json
    write(*,*)'PASSED UI JSON FILE TEST'
    call end_test_suite
    call run_suite('online variance',            test_online_var)
    call run_suite('image header',               test_imghead)
    call begin_test_suite('orientation data')
    call test_oris(.false.)
    call end_test_suite
    call begin_test_suite('image')
    call test_image(.false.)
    call end_test_suite
    call run_suite('Fourier shift search',        test_ftexp_shsrch)
    call run_suite('Fourier iterator',            test_ftiter)
    call begin_test_suite('B-spline smoother 2D')
    call test_bspline_smoother([64,64,1], 1.0, 0.2)
    call end_test_suite
    call begin_test_suite('B-spline smoother 3D')
    call test_bspline_smoother_3d([64,64,64], 1.0, 0.2)
    call end_test_suite
    ! local test functions
    call run_suite('multinomial random draw',     test_multinomal)
    call run_suite('Euler shift',                 test_euler_shift)
    call run_suite('straight-line fit',           simple_test_fit_line)
    call run_suite('affinity propagation',        test_aff_prop)
    call run_suite('hierarchical clustering',     test_hclust)
    call run_suite('atoms',                       test_atoms)
    call report_summary(failed=test_failed)
    call simple_chdir( "../" )
    if( test_failed ) error stop 1
    call simple_end('**** SIMPLE_TEST_UNITS NORMAL STOP ****')

    contains

        subroutine run_suite( name, test_proc )
            character(len=*), intent(in) :: name
            procedure(no_arg_test)       :: test_proc
            call begin_test_suite(name)
            call test_proc
            call end_test_suite
        end subroutine run_suite

        subroutine test_multinomal
            integer :: i, irnd
            real :: pvec(10), prob
            call seed_rnd
            pvec(1) = 0.8
            do i=2,10
                pvec(i) = 0.2/9.
            end do
            write(logfhandle,*) 'this should be one:', sum(pvec)
            prob=0.
            do i=1,1000
                if( multinomal(pvec) == 1 ) prob = prob+1.
            end do
            prob = prob/1000.
            write(logfhandle,*) 'this should be 0.8:', prob
            pvec = 0.1
            write(logfhandle,*) 'this should be one:', sum(pvec)
            prob=0.
            do i=1,1000
                irnd = multinomal(pvec)
                if( irnd == 1 ) prob = prob+1.
            end do
            prob = prob/1000.
            write(logfhandle,*) 'this should be 0.1:', prob
            write(logfhandle,'(a)') 'SIMPLE_RND: MULTINOMAL TEST COMPLETED WITHOUT TERMINAL BUGS ;-)'
        end subroutine

        subroutine test_euler_shift
            type(ori) :: o
            integer   :: i
            real      :: euls(3), euls_shifted(3)
            logical   :: doshift
            call o%new(is_ptcl=.false.)
            do i=1,100000
                euls(1) = ran3()*800.-400.
                euls(2) = ran3()*500-250.
                euls(3) = ran3()*800.-400.
                call o%set_euler(euls)
                euls_shifted = o%get_euler()
                doshift = .false.
                if( euls_shifted(1) < 0. .or. euls_shifted(1) > 360. ) doshift = .true.
                if( euls_shifted(2) < 0. .or. euls_shifted(2) > 180. ) doshift = .true.
                if( euls_shifted(3) < 0. .or. euls_shifted(3) > 360. ) doshift = .true.
                if( doshift ) THROW_HARD('euler shifting does not work!')
            end do
        end subroutine

        subroutine simple_test_fit_line
            real    :: slope, intercept, datavec(100,2), corr, x
            integer :: i, j
            do i=1,10000
                ! generate the line
                slope = 5.*ran3()
                if( ran3() < 0.5 ) slope = -slope
                intercept = 10.*ran3()
                if( ran3() < 0.5 ) intercept = -intercept
    !            write(logfhandle,*) '***********************************'
    !            write(logfhandle,*) 'Slope/Intercept:', slope, intercept
                ! generate the data
                x = -1.
                do j=1,100
                    datavec(j,1) = x
                    datavec(j,2) = slope*datavec(j,1)+intercept
                    x = x+0.02
                end do
                ! fit the data
                call fit_straight_line(100, datavec, slope, intercept, corr)
    !            write(logfhandle,*) 'Fitted Slope/Intercept:', slope, intercept
    !            write(logfhandle,*) 'Corr:', corr
                if( corr < 0.9999 )then
                    THROW_HARD('fit_straight_line failed!')
                endif
            end do
            write(logfhandle,'(a)') 'FIT_STRAIGHT_LINE UNIT TEST COMPLETED ;-)'
        end subroutine

end subroutine exec_test_units

end module simple_commanders_test_class
