!@descr: the unit-test suites: the build's fast gate (test=unit_<area>), its umbrella (test=units) and the platform-tier forked-process suite
module simple_commanders_test_class
use simple_commanders_api
use simple_test_utils,                       only: begin_test_suite, end_test_suite, reset_test_report, report_summary
! core library tester modules
use simple_string_tester,                    only: run_all_string_tests
use simple_syslib_tester,                    only: run_all_syslib_tests
use simple_fileio_tester,                    only: run_all_fileio_tests
use simple_stack_io_tester,                  only: run_all_stack_io_tests
use simple_class_sample_io_tester,           only: run_all_class_sample_io_tests
use simple_chash_tester,                     only: run_all_chash_tests
use simple_vrefhash_tester,                  only: run_all_vrefhash_tests
use simple_hash_tester,                      only: run_all_hash_tests
use simple_linked_list_tester,               only: run_all_list_tests
use simple_cmdline_tester,                   only: run_all_cmdline_tests
use simple_rec_list_tester,                  only: run_all_rec_list_tests
use simple_ori_tester,                       only: run_all_ori_tests
use simple_oris_tester,                      only: run_all_oris_tests
use simple_sym_tester,                       only: run_all_sym_tests
use simple_stat_tester,                      only: run_all_stat_tests
use simple_linalg_tester,                    only: run_all_linalg_tests
use simple_kbinterpol_tester,                only: run_all_kbinterpol_tests
use simple_srch_sort_loc_tester,             only: run_all_srch_sort_loc_tests
use simple_decay_funs_tester,                only: run_all_decay_funs_tests
use simple_pca_tester,                       only: run_all_pca_tests
use simple_opt_tester,                       only: run_all_opt_tests
use simple_lpstages_tester,                  only: run_all_lpstages_tests
use simple_image_msk_tester,                 only: run_all_mask_tests, run_all_image_bin_tests
use simple_segmentation_tester,              only: run_all_segmentation_tests
use simple_accum_blend_tester,               only: run_all_accum_blend_tests
use simple_ctf_tester,                       only: run_all_ctf_tests
use simple_starfile_tester,                  only: run_all_starfile_tests
use simple_starproject_tester,               only: run_all_starproject_tests
use simple_binoris_tester,                   only: run_all_binoris_tests
use simple_sp_project_tester,                only: run_all_sp_project_tests
use simple_project_merge_tester,             only: run_all_project_merge_tests
use simple_class_compatibility_tester,       only: run_all_class_compatibility_tests
use simple_ptcl_sieve_tester,                only: run_all_ptcl_sieve_tests
use simple_motion_gain_tester,               only: run_all_motion_gain_tests
use simple_gui_metadata_tester,              only: run_all_gui_metadata_tests
use simple_gui_assembler_tester,             only: run_all_gui_assembler_tests
use simple_ui_hash_tester,                   only: run_all_ui_hash_tests
use simple_gauran_tester,                    only: run_all_gauran_tests
use simple_rec3D_strategy_tester,            only: run_all_rec3D_strategy_tests
use simple_pcg_halfset_tester,               only: run_all_pcg_halfset_tests
use simple_ipc_tcp_socket_tester,            only: run_all_ipc_tcp_socket_tests
use simple_http_post_tester,                 only: run_all_http_post_tests
use simple_persistent_worker_message_tester, only: run_all_persistent_worker_message_tests
use simple_persistent_worker_server_tester,  only: run_all_persistent_worker_server_tests
use simple_forked_process_tester,            only: run_all_forked_process_tests
! test procedures of core types
use simple_imghead,                          only: test_imghead
use simple_oris,                             only: test_oris
use simple_image,                            only: test_image
use simple_ftiter,                           only: test_ftiter
use simple_ftexp_shsrch,                     only: test_ftexp_shsrch, test_ftexp_shsrch2
use simple_bspline_smoother,                 only: test_bspline_smoother, test_bspline_smoother_3d
use simple_online_var,                       only: test_online_var
use simple_aff_prop,                         only: test_aff_prop
use simple_hclust,                           only: test_hclust
use simple_atoms,                            only: test_atoms
use simple_srchspace_map2D_io,               only: test_srchspace_map2D_io
use simple_ui,                               only: validate_ui_json
implicit none
#include "simple_local_flags.inc"

! The fast gate is eight area suites, each one CTest entry under the label
! `fast` (doc/refactoring_notes/uniform_test_environment_refactoring.md,
! section 5.1). Every sub-suite in them makes assertions through
! simple_test_utils, needs no network beyond localhost, no download and no
! user-supplied data, and runs on one OpenMP thread.
!
!   test=unit_<area>                 one area suite in one process
!   test=unit_<area> suite=<name>    one sub-suite of it (name as in test=list,
!                                    lowercase, spaces as underscores)
!   test=units                       every area suite in sequence: a developer
!                                    convenience, not the gate CTest runs
!   test=forked_process              real child processes, clock polling:
!                                    excluded from the build, label `platform`
!   test=lib_<area>                  a library suite of the nightly extensive
!                                    tier (section 5.2.1): same shape, no
!                                    30 s budget; lib_reconstruction is the first
!
! SIMPLE_UNIT_ORDER=reverse runs a suite's table backwards; a result that
! differs from the forward run is a state leak between sub-suites.


type, extends(commander_base) :: commander_test_units
  contains
    procedure :: execute      => exec_test_units
end type commander_test_units

type, extends(commander_base) :: commander_test_unit_core
  contains
    procedure :: execute      => exec_test_unit_core
end type commander_test_unit_core

type, extends(commander_base) :: commander_test_unit_ori
  contains
    procedure :: execute      => exec_test_unit_ori
end type commander_test_unit_ori

type, extends(commander_base) :: commander_test_unit_image
  contains
    procedure :: execute      => exec_test_unit_image
end type commander_test_unit_image

type, extends(commander_base) :: commander_test_unit_numerics
  contains
    procedure :: execute      => exec_test_unit_numerics
end type commander_test_unit_numerics

type, extends(commander_base) :: commander_test_unit_project
  contains
    procedure :: execute      => exec_test_unit_project
end type commander_test_unit_project

type, extends(commander_base) :: commander_test_unit_ui
  contains
    procedure :: execute      => exec_test_unit_ui
end type commander_test_unit_ui

type, extends(commander_base) :: commander_test_unit_ipc
  contains
    procedure :: execute      => exec_test_unit_ipc
end type commander_test_unit_ipc

type, extends(commander_base) :: commander_test_unit_reconstruction
  contains
    procedure :: execute      => exec_test_unit_reconstruction
end type commander_test_unit_reconstruction

type, extends(commander_base) :: commander_test_lib_reconstruction
  contains
    procedure :: execute      => exec_test_lib_reconstruction
end type commander_test_lib_reconstruction

type, extends(commander_base) :: commander_test_forked_process
  contains
    procedure :: execute      => exec_test_forked_process
end type commander_test_forked_process

! a sub-suite: a name and a no-argument procedure that asserts through simple_test_utils
abstract interface
    subroutine no_arg_test()
    end subroutine no_arg_test
end interface

type :: unit_suite
    character(len=32) :: name = ''
    procedure(no_arg_test), pointer, nopass :: run => null()
end type unit_suite

integer, parameter :: MAX_SUITES = 64

contains

    ! ---- area tables ---------------------------------------------------------
    ! One function per area returns its sub-suites in table order. The umbrella
    ! (test=units) is the concatenation of all seven.

    subroutine suites_core( s, n )
        type(unit_suite), intent(inout) :: s(:)
        integer,          intent(inout) :: n
        call add_suite(s, n, 'string',               run_all_string_tests)
        call add_suite(s, n, 'syslib',               run_all_syslib_tests)
        call add_suite(s, n, 'fileio',               run_all_fileio_tests)
        call add_suite(s, n, 'stack I/O',            run_all_stack_io_tests)
        call add_suite(s, n, 'class sample I/O',     run_all_class_sample_io_tests)
        call add_suite(s, n, 'character hash',       run_all_chash_tests)
        call add_suite(s, n, 'hash',                 run_all_hash_tests)
        call add_suite(s, n, 'value-reference hash', run_all_vrefhash_tests)
        call add_suite(s, n, 'linked list',          run_all_list_tests)
        call add_suite(s, n, 'record list',          run_all_rec_list_tests)
        call add_suite(s, n, 'command line',         run_all_cmdline_tests)
    end subroutine suites_core

    subroutine suites_ori( s, n )
        type(unit_suite), intent(inout) :: s(:)
        integer,          intent(inout) :: n
        call add_suite(s, n, 'orientation',            run_all_ori_tests)
        call add_suite(s, n, 'orientation collection', run_all_oris_tests)
        call add_suite(s, n, 'symmetry',               run_all_sym_tests)
        call add_suite(s, n, 'orientation data',       suite_orientation_data)
        call add_suite(s, n, 'Euler shift',            test_euler_shift)
    end subroutine suites_ori

    subroutine suites_image( s, n )
        type(unit_suite), intent(inout) :: s(:)
        integer,          intent(inout) :: n
        call add_suite(s, n, 'image',                suite_image)
        call add_suite(s, n, 'image header',         test_imghead)
        call add_suite(s, n, 'Fourier iterator',     test_ftiter)
        call add_suite(s, n, 'B-spline smoother 2D', suite_bspline_2d)
        call add_suite(s, n, 'B-spline smoother 3D', suite_bspline_3d)
        call add_suite(s, n, 'masks',                run_all_mask_tests)
        call add_suite(s, n, 'binary image',         run_all_image_bin_tests)
        call add_suite(s, n, 'segmentation',         run_all_segmentation_tests)
        call add_suite(s, n, 'trailing-reconstruction blend', run_all_accum_blend_tests)
        call add_suite(s, n, 'CTF',                  run_all_ctf_tests)
    end subroutine suites_image

    subroutine suites_numerics( s, n )
        type(unit_suite), intent(inout) :: s(:)
        integer,          intent(inout) :: n
        call add_suite(s, n, 'online variance',         test_online_var)
        call add_suite(s, n, 'multinomial random draw', test_multinomal)
        call add_suite(s, n, 'straight-line fit',       test_fit_line)
        call add_suite(s, n, 'affinity propagation',    test_aff_prop)
        call add_suite(s, n, 'hierarchical clustering', test_hclust)
        call add_suite(s, n, 'statistics',              run_all_stat_tests)
        call add_suite(s, n, 'linear algebra',          run_all_linalg_tests)
        call add_suite(s, n, 'Kaiser-Bessel kernel',    run_all_kbinterpol_tests)
        call add_suite(s, n, 'search, sort, locate',    run_all_srch_sort_loc_tests)
        call add_suite(s, n, 'decay schedules',         run_all_decay_funs_tests)
        call add_suite(s, n, 'PCA',                     run_all_pca_tests)
        call add_suite(s, n, 'optimisers',              run_all_opt_tests)
        call add_suite(s, n, 'low-pass stages',         run_all_lpstages_tests)
        ! motion-correction shift search on expanded Fourier transforms (an optimiser, not an image test)
        call add_suite(s, n, 'shift search, correlator', test_ftexp_shsrch)
        call add_suite(s, n, 'shift search, optimiser',  test_ftexp_shsrch2)
    end subroutine suites_numerics

    subroutine suites_project( s, n )
        type(unit_suite), intent(inout) :: s(:)
        integer,          intent(inout) :: n
        call add_suite(s, n, 'STAR file',               run_all_starfile_tests)
        call add_suite(s, n, 'STAR project',            run_all_starproject_tests)
        call add_suite(s, n, 'binoris',                 run_all_binoris_tests)
        call add_suite(s, n, 'project records',         run_all_sp_project_tests)
        call add_suite(s, n, 'project merge',           run_all_project_merge_tests)
        call add_suite(s, n, 'class compatibility',     run_all_class_compatibility_tests)
        call add_suite(s, n, 'particle sieve',          run_all_ptcl_sieve_tests)
        call add_suite(s, n, '2D search-space map I/O', test_srchspace_map2D_io)
        call add_suite(s, n, 'motion gain',             run_all_motion_gain_tests)
        call add_suite(s, n, 'atoms',                   test_atoms)
    end subroutine suites_project

    subroutine suites_ui( s, n )
        type(unit_suite), intent(inout) :: s(:)
        integer,          intent(inout) :: n
        call add_suite(s, n, 'UI JSON',       suite_ui_json)
        call add_suite(s, n, 'GUI metadata',  run_all_gui_metadata_tests)
        call add_suite(s, n, 'GUI assembler', run_all_gui_assembler_tests)
        call add_suite(s, n, 'UI hash',       run_all_ui_hash_tests)
    end subroutine suites_ui

    subroutine suites_ipc( s, n )
        type(unit_suite), intent(inout) :: s(:)
        integer,          intent(inout) :: n
        ! localhost only, bounded; `forked process` is deliberately not here
        call add_suite(s, n, 'IPC TCP socket',            run_all_ipc_tcp_socket_tests)
        call add_suite(s, n, 'HTTP POST',                 run_all_http_post_tests)
        call add_suite(s, n, 'persistent worker message', run_all_persistent_worker_message_tests)
        call add_suite(s, n, 'persistent worker server',  run_all_persistent_worker_server_tests)
    end subroutine suites_ipc

    subroutine suites_reconstruction( s, n )
        type(unit_suite), intent(inout) :: s(:)
        integer,          intent(inout) :: n
        call add_suite(s, n, 'rec3D backend',     run_all_rec3D_strategy_tests)
        call add_suite(s, n, 'observation noise', run_all_gauran_tests)
    end subroutine suites_reconstruction

    !> nightly library suite: minutes, full boxes allowed, same assertions and runner
    subroutine suites_lib_reconstruction( s, n )
        type(unit_suite), intent(inout) :: s(:)
        integer,          intent(inout) :: n
        call add_suite(s, n, 'PCG half-set', run_all_pcg_halfset_tests)
    end subroutine suites_lib_reconstruction

    subroutine add_suite( s, n, name, proc )
        type(unit_suite),      intent(inout) :: s(:)
        integer,               intent(inout) :: n
        character(len=*),      intent(in)    :: name
        procedure(no_arg_test)               :: proc
        if( n >= size(s) ) THROW_HARD('too many unit sub-suites; raise MAX_SUITES')
        n = n + 1
        s(n)%name = name
        s(n)%run  => proc
    end subroutine add_suite

    ! ---- commanders ----------------------------------------------------------

    subroutine exec_test_units( self, cline )
        class(commander_test_units), intent(inout) :: self
        class(cmdline),              intent(inout) :: cline
        type(unit_suite) :: s(MAX_SUITES)
        integer :: n
        n = 0
        call suites_core(s, n)
        call suites_ori(s, n)
        call suites_image(s, n)
        call suites_numerics(s, n)
        call suites_project(s, n)
        call suites_ui(s, n)
        call suites_ipc(s, n)
        call suites_reconstruction(s, n)
        call run_unit_suites('units', cline, s(1:n))
    end subroutine exec_test_units

    subroutine exec_test_unit_core( self, cline )
        class(commander_test_unit_core), intent(inout) :: self
        class(cmdline),                  intent(inout) :: cline
        type(unit_suite) :: s(MAX_SUITES)
        integer :: n
        n = 0
        call suites_core(s, n)
        call run_unit_suites('unit_core', cline, s(1:n))
    end subroutine exec_test_unit_core

    subroutine exec_test_unit_ori( self, cline )
        class(commander_test_unit_ori), intent(inout) :: self
        class(cmdline),                 intent(inout) :: cline
        type(unit_suite) :: s(MAX_SUITES)
        integer :: n
        n = 0
        call suites_ori(s, n)
        call run_unit_suites('unit_ori', cline, s(1:n))
    end subroutine exec_test_unit_ori

    subroutine exec_test_unit_image( self, cline )
        class(commander_test_unit_image), intent(inout) :: self
        class(cmdline),                   intent(inout) :: cline
        type(unit_suite) :: s(MAX_SUITES)
        integer :: n
        n = 0
        call suites_image(s, n)
        call run_unit_suites('unit_image', cline, s(1:n))
    end subroutine exec_test_unit_image

    subroutine exec_test_unit_numerics( self, cline )
        class(commander_test_unit_numerics), intent(inout) :: self
        class(cmdline),                      intent(inout) :: cline
        type(unit_suite) :: s(MAX_SUITES)
        integer :: n
        n = 0
        call suites_numerics(s, n)
        call run_unit_suites('unit_numerics', cline, s(1:n))
    end subroutine exec_test_unit_numerics

    subroutine exec_test_unit_project( self, cline )
        class(commander_test_unit_project), intent(inout) :: self
        class(cmdline),                     intent(inout) :: cline
        type(unit_suite) :: s(MAX_SUITES)
        integer :: n
        n = 0
        call suites_project(s, n)
        call run_unit_suites('unit_project', cline, s(1:n))
    end subroutine exec_test_unit_project

    subroutine exec_test_unit_ui( self, cline )
        class(commander_test_unit_ui), intent(inout) :: self
        class(cmdline),                intent(inout) :: cline
        type(unit_suite) :: s(MAX_SUITES)
        integer :: n
        n = 0
        call suites_ui(s, n)
        call run_unit_suites('unit_ui', cline, s(1:n))
    end subroutine exec_test_unit_ui

    subroutine exec_test_unit_ipc( self, cline )
        class(commander_test_unit_ipc), intent(inout) :: self
        class(cmdline),                 intent(inout) :: cline
        type(unit_suite) :: s(MAX_SUITES)
        integer :: n
        n = 0
        call suites_ipc(s, n)
        call run_unit_suites('unit_ipc', cline, s(1:n))
    end subroutine exec_test_unit_ipc

    subroutine exec_test_unit_reconstruction( self, cline )
        class(commander_test_unit_reconstruction), intent(inout) :: self
        class(cmdline),                            intent(inout) :: cline
        type(unit_suite) :: s(MAX_SUITES)
        integer :: n
        n = 0
        call suites_reconstruction(s, n)
        call run_unit_suites('unit_reconstruction', cline, s(1:n))
    end subroutine exec_test_unit_reconstruction

    subroutine exec_test_lib_reconstruction( self, cline )
        class(commander_test_lib_reconstruction), intent(inout) :: self
        class(cmdline),                           intent(inout) :: cline
        type(unit_suite) :: s(MAX_SUITES)
        integer :: n
        n = 0
        call suites_lib_reconstruction(s, n)
        call run_unit_suites('lib_reconstruction', cline, s(1:n))
    end subroutine exec_test_lib_reconstruction

    subroutine exec_test_forked_process( self, cline )
        class(commander_test_forked_process), intent(inout) :: self
        class(cmdline),                       intent(inout) :: cline
        type(unit_suite) :: s(1)
        integer :: n
        n = 0
        call add_suite(s, n, 'forked process', run_all_forked_process_tests)
        call run_unit_suites('forked_process', cline, s(1:n))
    end subroutine exec_test_forked_process

    ! ---- the runner ----------------------------------------------------------

    !> Runs the sub-suites of one area in this process, in its own dated
    !! directory, accumulating failures through simple_test_utils; exits
    !! non-zero if any check failed. `suite=<name>` (from the command line)
    !! runs one sub-suite; SIMPLE_UNIT_ORDER=reverse walks the table backwards.
    subroutine run_unit_suites( label, cline, suites )
        character(len=*), intent(in)    :: label
        class(cmdline),   intent(inout) :: cline
        type(unit_suite), intent(in)    :: suites(:)
        character(8)          :: datestr
        character(len=STDLEN) :: folder
        type(string)          :: original_cwd, report_file, only_suite, order
        logical               :: test_failed, l_reverse
        integer               :: i, isuite, nrun, iostat
        call seed_rnd
        call date_and_time(date=datestr)
        folder = 'SIMPLE_TEST_'//trim(label)//'_'//datestr
        call simple_getcwd(original_cwd)
        report_file = original_cwd//'/'//trim(folder)//'/simple_test_'//trim(label)//'_report.txt'
        call simple_mkdir(folder)
        call simple_chdir(folder)
        call reset_test_report(report_file%to_char())
        only_suite = ''
        if( cline%defined('suite') )then
            only_suite = cline%get_carg('suite')
            only_suite = suite_id(only_suite%to_char())
        endif
        order     = simple_getenv('SIMPLE_UNIT_ORDER', iostat)
        l_reverse = iostat == 0 .and. order == 'reverse'
        nrun = 0
        do i = 1, size(suites)
            isuite = i
            if( l_reverse ) isuite = size(suites) + 1 - i
            if( only_suite%strlen_trim() > 0 )then
                if( .not. (only_suite == suite_id(suites(isuite)%name)) ) cycle
            endif
            call begin_test_suite(trim(suites(isuite)%name))
            call suites(isuite)%run()
            call end_test_suite
            nrun = nrun + 1
        end do
        call report_summary(failed=test_failed)
        call simple_chdir(original_cwd%to_char())
        if( nrun == 0 ) THROW_HARD('no sub-suite '//only_suite%to_char()//' in '//trim(label)//'; the names are listed by test=list')
        if( test_failed ) error stop 1
        call simple_end('**** SIMPLE_TEST_'//trim(label)//' NORMAL STOP ****')
    end subroutine run_unit_suites

    !> command-line spelling of a sub-suite name: lowercase, spaces and hyphens as underscores
    function suite_id( name ) result( id )
        character(len=*), intent(in) :: name
        character(len=len(name))     :: id
        integer :: i
        id = lowercase(name)
        do i = 1, len(id)
            if( id(i:i) == ' ' .or. id(i:i) == '-' ) id(i:i) = '_'
        end do
    end function suite_id

    ! ---- wrappers for test procedures that take arguments -----------------------

    subroutine suite_orientation_data
        call test_oris(.false.)
    end subroutine suite_orientation_data

    subroutine suite_image
        call test_image(.false.)
    end subroutine suite_image

    subroutine suite_bspline_2d
        call test_bspline_smoother([64,64,1], 1.0, 0.2)
    end subroutine suite_bspline_2d

    subroutine suite_bspline_3d
        call test_bspline_smoother_3d([64,64,64], 1.0, 0.2)
    end subroutine suite_bspline_3d

    subroutine suite_ui_json
        write(logfhandle,'(a)') 'VALIDATING UI JSON FILE:'
        call validate_ui_json
        write(logfhandle,'(a)') 'PASSED UI JSON FILE TEST'
    end subroutine suite_ui_json

    ! ---- local sub-suites (formerly contained in exec_test_units) -----------------

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
    end subroutine test_multinomal

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
    end subroutine test_euler_shift

    subroutine test_fit_line
        real    :: slope, intercept, datavec(100,2), corr, x
        integer :: i, j
        do i=1,10000
            ! generate the line
            slope = 5.*ran3()
            if( ran3() < 0.5 ) slope = -slope
            intercept = 10.*ran3()
            if( ran3() < 0.5 ) intercept = -intercept
            ! generate the data
            x = -1.
            do j=1,100
                datavec(j,1) = x
                datavec(j,2) = slope*datavec(j,1)+intercept
                x = x+0.02
            end do
            ! fit the data
            call fit_straight_line(100, datavec, slope, intercept, corr)
            if( corr < 0.9999 )then
                THROW_HARD('fit_straight_line failed!')
            endif
        end do
        write(logfhandle,'(a)') 'FIT_STRAIGHT_LINE UNIT TEST COMPLETED ;-)'
    end subroutine test_fit_line

end module simple_commanders_test_class
