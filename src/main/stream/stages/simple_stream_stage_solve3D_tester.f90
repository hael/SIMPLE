!@descr: unit tests for the steps of the stream multistate 3D stage (simple_stream_stage_solve3D)
! Each test assembles a stage in a fresh fixture directory from init_params and init_gui: no queue
! environment, no waits, and a settle time of -1. The command line names no registered program
! and carries qsys_name=local for the stage project's computing environment. The upstream is a pool
! 2D directory whose completed folder holds exports. The class-average scoring and the jobs (which
! need a queue) are left to the high-level stream tests; the import is tested on sets made in
! memory, the class averages a publication brings on fixture stacks, the reference volumes on a
! fixture solve3D_cavgs result, an assignment on a fixture matching, the local refinement's seeds
! and result on fixture rows and a fixture refine3D project, the rules on fixed numbers, and the
! volume messages on a fixture project with a volume and an FSC; a 3D snapshot is written from a
! fixture result project, on a request sent through the stage's update pipe.
module simple_stream_stage_solve3D_tester
use, intrinsic :: iso_c_binding, only: c_int
use unix,                        only: c_pipe, c_close, c_fcntl, F_GETFL, F_SETFL, O_NONBLOCK
use json_module,                 only: json_core, json_value
use simple_test_utils
use simple_defs_fname,                                only: TERM_STREAM, METADATA_EXT
use simple_refine3D_fnames,                           only: refine3D_reprojs_fname
use simple_defs_stream,                               only: DIR_STREAM_COMPLETED
use simple_string,                                    only: string
use simple_string_utils,                              only: int2str, int2str_pad
use simple_fileio,                                    only: arr2file, del_file, file_exists, simple_getcwd, simple_touch, add2fbody
use simple_syslib,                                    only: simple_mkdir, dir_exists
use simple_math_ft,                                   only: get_resarr
use simple_cmdline,                                   only: cmdline
use simple_image,                                     only: image
use simple_ori,                                       only: ori
use simple_strategy3D_cont,                           only: cont_seed_is_valid
use simple_sp_project,                                only: sp_project
use simple_rec_list,                                  only: rec_iterator, chunk_rec
use simple_gui_metadata_utils,                        only: max_metadata_size
use simple_gui_metadata_types,                        only: GUI_METADATA_STREAM_SOLVE3D_MULTISTATE_TYPE, GUI_METADATA_VOL3D_TYPE,&
                                                           &GUI_METADATA_STREAM_SOLVE3D_MULTISTATE_REPROJ_TYPE,&
                                                           &GUI_METADATA_STREAM_SOLVE3D_SNAPSHOT_TYPE, GUI_METADATA_STREAM_UPDATE_TYPE
use simple_gui_metadata_stream_solve3D_multistate, only: gui_metadata_stream_solve3D_multistate
use simple_gui_metadata_stream_snapshot,              only: gui_metadata_stream_snapshot
use simple_gui_metadata_stream_update,                only: gui_metadata_stream_update
use simple_gui_metadata_vol3D,                        only: gui_metadata_vol3D
use simple_stream_pipe,                               only: stream_pipe
use simple_stream_stage_solve3D,                   only: stream_stage_solve3D, PHASE_IMPORTING, PHASE_CAVGS3D,&
                                                           &PHASE_IDLE, PHASE_ASSIGN, JOB_NONE, JOB_ROUTE, JOB_FALLBACK,&
                                                           &JOB_ASSIGN, ROUTE_TRIES
implicit none
private
public :: run_all_stream_stage_solve3D_tests, run_all_stream_stage_solve3D_lib_tests

character(len=*), parameter :: TEST_PROJFILE = 'test_solve3D.simple'
character(len=*), parameter :: UPSTREAM      = 'pool2D' ! the pool 2D stage directory (dir_target)
integer,          parameter :: NSTATES       = 3
integer,          parameter :: VOL_BOX       = 16
integer,          parameter :: DOCK_BOX      = 64   ! the docking's smallest box (dock_vols)
real,             parameter :: VOL_SMPD      = 2.0

contains

    subroutine run_all_stream_stage_solve3D_tests()
        write(*,'(A)') '**** running all stream solve 3D stage tests ****'
        call test_init_params()
        call test_restart_removes_term_stream()
        call test_watch_order_and_mskdiam()
        call test_merge_publications()
        call test_take_cavgs()
        call test_first_publication()
        call test_mskdiam_from_each_publication()
        call test_rows_problem()
        call test_rules()
        call test_fit_rules()
        call test_reference_volumes()
        call test_failed_attempts()
        call test_assignment()
        call test_local_refinement()
        call test_retention()
        call test_send_status()
        call test_send_volumes()
        call test_send_refs()
        call test_snapshot3D()
        call test_iterate_waits()
        call test_finished()
    end subroutine run_all_stream_stage_solve3D_tests

    !> the library tests: what takes the docking's time
    subroutine run_all_stream_stage_solve3D_lib_tests()
        write(*,'(A)') '**** running all stream solve 3D stage library tests ****'
        call test_reference_docking()
    end subroutine run_all_stream_stage_solve3D_lib_tests

    !> the states of an attempt compared after docking, in both hands (at the docking's smallest
    !! box): two different structures are distinct; a structure and a rotated, mirrored copy of it,
    !! as an independently refined state in another frame and hand would be, are near-duplicates
    subroutine test_reference_docking()
        class(stream_stage_solve3D), allocatable :: stage
        type(cmdline)    :: cline
        type(sp_project) :: res
        type(string)     :: cwd_saved, root
        real, parameter  :: TRIPOD(3,3)     = reshape([12., 0., 0.,  -12., 4., 0.,  0., -6., 10.], [3,3])
        real, parameter  :: TRIPOD_AMPS(3)  = [1.0, 0.6, 0.8], TRIPOD_SIGS(3) = [4.0, 2.5, 3.0]
        real, parameter  :: CHAIN(3,4)      = reshape([-18., 0., 0.,  -6., 0., 0.,  6., 0., 0.,  18., 0., 0.], [3,4])
        real, parameter  :: CHAIN_AMPS(4)   = 1., CHAIN_SIGS(4) = 3.
        real    :: tripod_copy(3,3)
        integer :: nfail0, i
        allocate(stage)
        write(*,'(A)') 'test_reference_docking'
        nfail0 = tests_failed
        call enter_fixture('a3_stage_reference_docking', cwd_saved, root)
        call set_test_cline(cline)
        call make_test_stage(stage, cline)
        stage%l_dock_refs = .true.
        ! two different structures: a tripod and a chain of four blobs
        call make_docking_result(res, TRIPOD, TRIPOD_AMPS, TRIPOD_SIGS, CHAIN, CHAIN_AMPS, CHAIN_SIGS)
        call assert_true(stage%register_attempt(res, 1), 'two different structures are distinct states')
        ! the tripod and its copy rotated and mirrored, (x, y, z) -> (-y, -x, z): one structure twice
        do i = 1,3
            tripod_copy(:,i) = [-TRIPOD(2,i), -TRIPOD(1,i), TRIPOD(3,i)]
        enddo
        call make_docking_result(res, TRIPOD, TRIPOD_AMPS, TRIPOD_SIGS, tripod_copy, TRIPOD_AMPS, TRIPOD_SIGS)
        call assert_false(stage%register_attempt(res, 2), 'a rotated, mirrored copy is a near-duplicate')
        call res%kill
        call stage%kill
        call cline%kill
        call leave_fixture(cwd_saved, root, nfail0)

    contains

        ! solve3D_cavgs' result in memory: four class averages in its states 1, 1, 2 and 2, a
        ! volume of Gaussian blobs per state (@p c1, @p a1, @p s1 and @p c2, @p a2, @p s2: centres in
        ! voxels from the centre, amplitudes, widths in voxels) at the docking's smallest box, and
        ! halves that agree
        subroutine make_docking_result( proj, c1, a1, s1, c2, a2, s2 )
            type(sp_project), intent(inout) :: proj
            real,             intent(in)    :: c1(:,:), a1(:), s1(:), c2(:,:), a2(:), s2(:)
            integer :: j
            call proj%kill
            call proj%os_cls3D%new(4, is_ptcl=.false.)
            do j = 1,4
                call proj%os_cls3D%set(j, 'e1', 30.)
                call proj%os_cls3D%set_state(j, merge(1, 2, j <= 2))
            enddo
            call simple_mkdir('dock')
            call write_blob_volume(string('dock/state01.mrc'), c1, a1, s1)
            call write_blob_volume(string('dock/state02.mrc'), c2, a2, s2)
            call proj%add_vol2os_out(string('dock/state01.mrc'), VOL_SMPD, 1, 'vol_cavg')
            call proj%add_vol2os_out(string('dock/state02.mrc'), VOL_SMPD, 2, 'vol_cavg')
            do j = 1,proj%os_out%get_noris()
                call proj%os_out%set(j, 'eo_angdist',       12.)
                call proj%os_out%set(j, 'eo_state_overlap', 85.)
            enddo
        end subroutine make_docking_result

        ! a DOCK_BOX volume of Gaussian blobs at @p centres (voxels from the centre), amplitudes
        ! @p amps and widths @p sigs (voxels)
        subroutine write_blob_volume( fname, centres, amps, sigs )
            class(string), intent(in) :: fname
            real,          intent(in) :: centres(:,:), amps(:), sigs(:)
            type(image)       :: vol
            real, allocatable :: rmat(:,:,:)
            real    :: r(3)
            integer :: x, y, z, k
            allocate(rmat(DOCK_BOX,DOCK_BOX,DOCK_BOX), source=0.)
            do z = 1,DOCK_BOX
                do y = 1,DOCK_BOX
                    do x = 1,DOCK_BOX
                        r = real([x, y, z] - (DOCK_BOX / 2 + 1))
                        do k = 1,size(amps)
                            rmat(x,y,z) = rmat(x,y,z) + amps(k) * exp(-sum((r - centres(:,k))**2) / (2. * sigs(k)**2))
                        enddo
                    enddo
                enddo
            enddo
            call vol%new([DOCK_BOX,DOCK_BOX,DOCK_BOX], VOL_SMPD)
            call vol%set_rmat(rmat, .false.)
            call vol%write(fname)
            call vol%kill
        end subroutine write_blob_volume

    end subroutine test_reference_docking

    subroutine test_init_params()
        class(stream_stage_solve3D), allocatable :: stage
        type(cmdline)                 :: cline
        type(string)                  :: cwd_saved, root
        integer                       :: nfail0
        allocate(stage)
        write(*,'(A)') 'test_init_params'
        nfail0 = tests_failed
        call enter_fixture('a3_stage_init_params', cwd_saved, root)
        call set_test_cline(cline)
        call make_test_stage(stage, cline)
        call assert_true(allocated(stage%params), 'initialization allocates the owned parameters')
        call assert_int(0, stage%spproj%os_mic%get_noris(), 'the project starts without micrographs')
        call assert_false(stage%l_restart,                  'a fresh run is no restart')
        call assert_int(PHASE_IMPORTING, stage%phase,       'the stage starts importing')
        call assert_int(NSTATES, size(stage%state_res),     'one resolution per state')
        call assert_int(NSTATES, size(stage%state_pop),     'and one population')
        call assert_int(0, size(stage%stk_names),           'no stack yet')
        call assert_int(0, stage%n_refs(),                  'no reference volumes yet')
        call stage%kill
        call assert_false(allocated(stage%params), 'cleanup releases the owned parameters')
        call cline%kill
        call leave_fixture(cwd_saved, root, nfail0)
    end subroutine test_init_params

    !> a restart removes a leftover termination file, which would end the restarted stage at once
    subroutine test_restart_removes_term_stream()
        class(stream_stage_solve3D), allocatable :: stage
        type(cmdline)                 :: cline
        type(string)                  :: cwd_saved, root
        integer                       :: nfail0
        allocate(stage)
        write(*,'(A)') 'test_restart_removes_term_stream'
        nfail0 = tests_failed
        call enter_fixture('a3_stage_restart', cwd_saved, root)
        call simple_mkdir('previous')
        call simple_touch(TERM_STREAM)
        call set_test_cline(cline)
        call cline%set('outdir', 'previous')
        call make_test_stage(stage, cline)
        call assert_true(stage%l_restart,           'an existing output folder is a restart')
        call assert_false(file_exists(TERM_STREAM), 'the leftover termination file is removed')
        call assert_false(stage%finished(),         'so the restarted stage runs')
        call stage%kill
        call cline%kill
        call leave_fixture(cwd_saved, root, nfail0)
    end subroutine test_restart_removes_term_stream

    !> exports are taken in export order whatever order they are listed in, each once
    subroutine test_watch_order_and_mskdiam()
        class(stream_stage_solve3D), allocatable :: stage
        type(cmdline)                 :: cline
        type(rec_iterator)            :: it
        type(chunk_rec)               :: crec
        type(string)                  :: cwd_saved, root
        integer                       :: nfail0
        allocate(stage)
        write(*,'(A)') 'test_watch_order_and_mskdiam'
        nfail0 = tests_failed
        call enter_fixture('a3_stage_watch', cwd_saved, root)
        call set_test_cline(cline)
        call make_test_stage(stage, cline)
        call stage%attach_upstream()
        call assert_false(stage%l_attached,      'no pool 2D folder: not attached')
        call assert_true(stage%l_waiting_logged, 'the wait is logged')
        call simple_mkdir(UPSTREAM)
        call simple_mkdir(UPSTREAM//'/'//DIR_STREAM_COMPLETED)
        call write_export(2, 170.)
        call write_export(1, 150.)
        call stage%attach_upstream()
        call assert_true(stage%l_attached, 'attached once the folder exists')
        call stage%watch_sets()
        call assert_int(2, stage%setslist%size(), 'one record per export')
        if( stage%setslist%size() == 2 )then
            it = stage%setslist%begin()
            call it%get(crec)
            call assert_true(crec%projfile%has_substr('00001'//METADATA_EXT), 'export 1 comes first')
            call it%next()
            call it%get(crec)
            call assert_true(crec%projfile%has_substr('00002'//METADATA_EXT), 'then export 2')
        endif
        call stage%watch_sets()
        call assert_int(2, stage%setslist%size(), 'an export is taken once')
        call stage%kill
        call cline%kill
        call leave_fixture(cwd_saved, root, nfail0)
    end subroutine test_watch_order_and_mskdiam

    !> publications are merged into rows that only grow: a new stack is appended with both particle
    !! segments; a stack the stage holds is matched by name in any order, its particles take the
    !! publication's class and selection and keep their 3D parameters; a stack a publication lacks
    !! keeps its rows, deselected; the classes are the newest publication's
    subroutine test_merge_publications()
        class(stream_stage_solve3D), allocatable :: stage
        type(cmdline)                 :: cline
        type(sp_project)              :: set
        type(string)                  :: cwd_saved, root
        integer                       :: nfail0, i
        allocate(stage)
        write(*,'(A)') 'test_merge_publications'
        nfail0 = tests_failed
        call enter_fixture('a3_stage_merge', cwd_saved, root)
        call set_test_cline(cline)
        call make_test_stage(stage, cline)
        ! publication 1: stacks A (3) and B (2), the first particle deselected, class 1
        call make_set(set, ['A', 'B'], [3, 2], 2, nrejected=1, icls=1)
        call stage%merge_publication(set, 1)
        call set%kill
        call assert_int(2, stage%spproj%os_mic%get_noris(),    'publication 1: its micrographs')
        call assert_int(5, stage%spproj%os_ptcl3D%get_noris(), 'its 3D particles')
        call assert_int(5, stage%spproj%os_ptcl2D%get_noris(), 'and 2D particles')
        call assert_int(4, stage%nptcls_selected,              'the selected particles are counted')
        ! a 3D result for the first particle, which later publications must keep, and a run's
        ! multistate label on the second
        call stage%spproj%os_ptcl3D%set(1, 'e1', 33.)
        call stage%spproj%os_ptcl3D%set_state(2, 2)
        ! publication 2, as a restarted pool lists it: B, A, then the new stack C (4), class 2
        call make_set(set, ['B', 'A', 'C'], [2, 3, 4], 5, icls=2)
        call stage%merge_publication(set, 2)
        call set%kill
        call assert_int(3, stage%spproj%os_mic%get_noris(),           'publication 2: one stack appended')
        call assert_int(9, stage%spproj%os_ptcl3D%get_noris(),        'with its particles; no row renumbered')
        call assert_int(6, stage%spproj%os_stk%get_fromp(3),          'the new stack follows the rows held')
        call assert_int(3, stage%spproj%os_ptcl3D%get_int(7, 'stkind'), 'its particles point at it')
        call assert_real(3., stage%spproj%os_ptcl3D%get(7, 'x'), 1.e-4, 'and carry its parameters')
        call assert_real(1., stage%spproj%os_ptcl3D%get(1, 'x'), 1.e-4, 'row 1 is still stack A''s first particle')
        call assert_int(2, stage%spproj%os_ptcl2D%get_class(1),       'which takes the new class')
        call assert_int(1, stage%spproj%os_ptcl3D%get_state(1),       'and the new selection')
        call assert_int(2, stage%spproj%os_ptcl3D%get_state(2),       'a selected particle keeps its 3D state label')
        call assert_real(33., stage%spproj%os_ptcl3D%get(1, 'e1'), 1.e-4, 'and keeps its 3D parameters')
        call assert_int(1, stage%spproj%os_ptcl2D%get_int(1, 'stkind'), 'and its stack')
        call assert_int(5, stage%spproj%os_cls2D%get_noris(),         'the classes are the newest publication''s')
        call assert_int(9, stage%nptcls_selected,                     'every particle selected')
        ! publication 3 lacks stack B (rows 4 and 5)
        call make_set(set, ['A', 'C'], [3, 4], 5, icls=2)
        call stage%merge_publication(set, 3)
        call set%kill
        call assert_int(9, stage%spproj%os_ptcl3D%get_noris(), 'publication 3: the rows are kept')
        call assert_int(0, stage%spproj%os_ptcl3D%get_state(4), 'a stack it lacks is deselected')
        call assert_int(0, stage%spproj%os_ptcl2D%get_state(5), 'in both particle segments')
        call assert_int(7, stage%nptcls_selected,               'and no longer counted')
        ! publication 4 lists stack A's images in reverse order: each updates its own row
        call make_set(set, ['A'], [3], 13, icls=1)
        do i = 1,3
            call set%os_ptcl2D%set(i, 'indstk', 4 - i)
            call set%os_ptcl2D%set_class(i, 14 - i)
        enddo
        call stage%merge_publication(set, 4)
        call set%kill
        call assert_int(11, stage%spproj%os_ptcl2D%get_class(1), 'a particle is matched by its image index, not its row')
        call assert_int(13, stage%spproj%os_ptcl2D%get_class(3), 'whatever order the publication lists them in')
        ! a publication with the optics map's table brings it
        call make_set(set, ['A'], [3], 13, icls=1)
        call set%os_optics%new(2, is_ptcl=.false.)
        call set%os_optics%set(1, 'ogid', 4)
        call set%os_optics%set(2, 'ogid', 9)
        call stage%merge_publication(set, 5)
        call set%kill
        call assert_int(2, stage%spproj%os_optics%get_noris(), 'the stage takes the publication''s optics table')
        call stage%kill
        call cline%kill
        call leave_fixture(cwd_saved, root, nfail0)
    end subroutine test_merge_publications

    !> a publication's class averages and FRCs are copied into its quality folder and replace the
    !! stage's earlier ones, whose state volume stays; a publication without FRCs cannot be used
    subroutine test_take_cavgs()
        class(stream_stage_solve3D), allocatable :: stage
        type(cmdline)      :: cline
        type(sp_project)   :: set
        type(image)        :: img
        type(string)       :: cwd_saved, root, stk, frcs, problem
        real               :: smpd, mskdiam
        integer            :: nfail0, ncls
        integer, parameter :: NCLS_PUB = 4, NCLS_OLD = 2, BOX = 8
        allocate(stage)
        write(*,'(A)') 'test_take_cavgs'
        nfail0 = tests_failed
        call enter_fixture('a3_stage_take_cavgs', cwd_saved, root)
        call set_test_cline(cline)
        call make_test_stage(stage, cline)
        ! the stage's earlier class averages and FRCs, and a state volume
        call write_stack(string('old_cavgs.mrcs'), NCLS_OLD)
        call simple_touch('old_frcs.bin')
        call stage%spproj%add_cavgs2os_out(string('old_cavgs.mrcs'), VOL_SMPD, 'cavg')
        call stage%spproj%add_frcs2os_out(string('old_frcs.bin'), 'frc2D')
        call img%new([BOX, BOX, BOX], VOL_SMPD)
        call img%write(string('recvol_state01.mrc'))
        call img%kill
        call stage%spproj%add_vol2os_out(string('recvol_state01.mrc'), VOL_SMPD, 1, 'vol')
        ! a publication with more classes
        call simple_mkdir('pub')
        call write_stack(string('pub/00003_cavgs.mrcs'), NCLS_PUB)
        call write_stack(string('pub/00003_cavgs_even.mrcs'), NCLS_PUB)
        call write_stack(string('pub/00003_cavgs_odd.mrcs'),  NCLS_PUB)
        call set%add_cavgs2os_out(string('pub/00003_cavgs.mrcs'), VOL_SMPD, 'cavg', mskdiam=150.)
        problem = stage%publication_problem(set)
        call assert_char('its FRCs are missing', problem%to_char(), 'a publication without FRCs cannot be used')
        call simple_touch('pub/00003_frcs.bin')
        call set%add_frcs2os_out(string('pub/00003_frcs.bin'), 'frc2D')
        problem = stage%publication_problem(set)
        call assert_char('', problem%to_char(), 'one with them can')
        ! its classes, as merge_publication takes them, then its class averages and FRCs
        stage%spproj%os_cls2D = set%os_cls2D
        call stage%take_cavgs(set, string('00003'))
        call stage%spproj%get_cavgs_stk(stk, ncls, smpd)
        call assert_int(NCLS_PUB, ncls, 'the class averages are the publication''s')
        call assert_true(stk%has_substr('quality_selection/00003/00003_cavgs.mrcs'), 'copied into its quality folder')
        call assert_true(file_exists(stk), 'where the copy is')
        call stage%spproj%get_frcs(frcs, 'frc2D')
        call assert_true(frcs%has_substr('quality_selection/00003/00003_frcs.bin'), 'and so are its FRCs')
        call assert_true(file_exists(frcs), 'copied too')
        call stage%spproj%get_mskdiam('cavg', mskdiam)
        call assert_real(150., mskdiam, 1.e-4, 'with its mask diameter')
        call assert_int(NCLS_PUB, stage%spproj%os_cls2D%get_noris(), 'the classes stay the publication''s')
        call assert_true(stage%spproj%isthere_in_osout('vol', 1), 'the state volume stays')
        call assert_true(file_exists(add2fbody(stk, '.mrcs', '_even')), 'the even half is copied beside them')
        call assert_true(file_exists(add2fbody(stk, '.mrcs', '_odd')),  'and the odd half')
        ! pool 2D removes its older publications
        call del_file('pub/00003_cavgs.mrcs')
        call del_file('pub/00003_frcs.bin')
        call assert_true(file_exists(stk),  'the class averages outlive the publication')
        call assert_true(file_exists(frcs), 'and the FRCs')
        call set%kill
        call stage%kill
        call cline%kill
        call leave_fixture(cwd_saved, root, nfail0)

    contains

        subroutine write_stack( fname, n )
            class(string), intent(in) :: fname
            integer,       intent(in) :: n
            integer :: i
            call img%new([BOX, BOX, 1], VOL_SMPD)
            do i = 1,n
                call img%write(fname, i)
            enddo
            call img%kill
        end subroutine write_stack

    end subroutine test_take_cavgs

    !> the first publication is taken as the pool selects it: its selected rows are the first set,
    !! whose selection later publications do not change (also when they lack its stack); the
    !! other rows follow them
    subroutine test_first_publication()
        class(stream_stage_solve3D), allocatable :: stage
        type(cmdline)                 :: cline
        type(sp_project)              :: set
        type(string)                  :: cwd_saved, root
        integer                       :: nfail0
        allocate(stage)
        write(*,'(A)') 'test_first_publication'
        nfail0 = tests_failed
        call enter_fixture('a3_stage_first_publication', cwd_saved, root)
        call set_test_cline(cline)
        call make_test_stage(stage, cline)
        ! publication 1: stacks A (3) and B (2), the first particle never updated (deselected)
        call make_set(set, ['A', 'B'], [3, 2], 2, nrejected=1, icls=1)
        call stage%take_first_publication(set, 1)
        call set%kill
        call assert_int(5, stage%spproj%os_ptcl3D%get_noris(), 'the first publication''s rows')
        call assert_int(4, stage%nptcls_selected,              'every selected particle, without a cap')
        call assert_true(allocated(stage%first_set),           'the selected rows are the first set')
        if( allocated(stage%first_set) )then
            call assert_true(all(stage%first_set .eqv. [.false., .true., .true., .true., .true.]), 'rows 2 to 5')
        endif
        call assert_int(PHASE_IMPORTING, stage%phase,          'no job yet')
        ! publication 2 selects every particle, in class 2
        call make_set(set, ['A', 'B'], [3, 2], 2, icls=2)
        call stage%merge_publication(set, 2)
        call set%kill
        call assert_int(1, stage%spproj%os_ptcl3D%get_state(1), 'a row the first never selected follows a later publication')
        call assert_int(2, stage%spproj%os_ptcl2D%get_class(3), 'a first-set row takes the 2D parameters')
        call assert_int(5, stage%nptcls_selected,               'the selected rows are counted')
        ! publication 3 deselects the first three particles
        call make_set(set, ['A', 'B'], [3, 2], 2, nrejected=3, icls=2)
        call stage%merge_publication(set, 3)
        call set%kill
        call assert_int(1, stage%spproj%os_ptcl3D%get_state(2), 'a first-set row keeps its selection')
        call assert_int(1, stage%spproj%os_ptcl2D%get_state(3), 'in both segments')
        call assert_int(0, stage%spproj%os_ptcl3D%get_state(1), 'the other rows follow the publication')
        call assert_int(4, stage%nptcls_selected,               'so four rows stay selected')
        ! publication 4 lacks stack A (rows 1 to 3)
        call make_set(set, ['B'], [2], 2, icls=2)
        call stage%merge_publication(set, 4)
        call set%kill
        call assert_int(1, stage%spproj%os_ptcl3D%get_state(3), 'a first-set row keeps its selection when its stack is missing')
        call assert_int(4, stage%nptcls_selected,               'and the count is unchanged')
        call stage%kill
        call cline%kill
        call leave_fixture(cwd_saved, root, nfail0)
    end subroutine test_first_publication

    !> the mask diameter comes from every publication taken; a change is taken for the next run, and
    !! a publication without one leaves it
    subroutine test_mskdiam_from_each_publication()
        class(stream_stage_solve3D), allocatable :: stage
        type(cmdline)    :: cline
        type(sp_project) :: set
        type(string)     :: cwd_saved, root
        integer          :: nfail0
        allocate(stage)
        write(*,'(A)') 'test_mskdiam_from_each_publication'
        nfail0 = tests_failed
        call enter_fixture('a3_stage_mskdiam', cwd_saved, root)
        call set_test_cline(cline)
        call make_test_stage(stage, cline)
        call add_cavgs_entry(set, 150.)
        call stage%take_mskdiam(set)
        call assert_real(150., stage%mskdiam, 1.e-4, 'the first publication''s mask diameter')
        call set%kill
        call add_cavgs_entry(set, 170.)
        call stage%take_mskdiam(set)
        call assert_real(170., stage%mskdiam, 1.e-4, 'a later publication''s new one')
        call set%kill
        call add_cavgs_entry(set, 0.)
        call stage%take_mskdiam(set)
        call assert_real(170., stage%mskdiam, 1.e-4, 'one without keeps the latest')
        call set%kill
        call stage%kill
        call cline%kill
        call leave_fixture(cwd_saved, root, nfail0)
    end subroutine test_mskdiam_from_each_publication

    !> a publication whose stacks disagree with the rows is found before anything is merged
    subroutine test_rows_problem()
        class(stream_stage_solve3D), allocatable :: stage
        type(cmdline)    :: cline
        type(sp_project) :: set
        type(string)     :: cwd_saved, root, problem
        integer          :: nfail0
        allocate(stage)
        write(*,'(A)') 'test_rows_problem'
        nfail0 = tests_failed
        call enter_fixture('a3_stage_rows_problem', cwd_saved, root)
        call set_test_cline(cline)
        call make_test_stage(stage, cline)
        call make_set(set, ['A'], [3], 2, icls=1)
        call stage%merge_publication(set, 1)
        call set%kill
        call make_set(set, ['A', 'B'], [3, 2], 2, icls=1)
        problem = stage%rows_problem(set)
        call assert_int(0, problem%strlen(), 'a publication that agrees with the rows')
        call set%kill
        call make_set(set, ['A'], [4], 2, icls=1)
        problem = stage%rows_problem(set)
        call assert_true(problem%has_substr('changed size'), 'a stack that changed size is found')
        call assert_int(3, stage%spproj%os_ptcl3D%get_noris(), 'and the rows are untouched')
        call set%kill
        call make_set(set, ['A'], [3], 2, icls=1)
        call set%os_ptcl2D%set(1, 'indstk', 7)
        problem = stage%rows_problem(set)
        call assert_true(problem%has_substr('outside its stack'), 'an image index outside its stack is found')
        call set%kill
        call stage%kill
        call cline%kill
        call leave_fixture(cwd_saved, root, nfail0)
    end subroutine test_rows_problem

    !> which job starts when, the check of the reference volumes, the best attempt, and the
    !! renumbering, sub-class budget and mask diameter the route uses
    subroutine test_rules()
        class(stream_stage_solve3D), allocatable :: stage
        integer, allocatable :: states(:), newind(:)
        integer :: nstates_held
        allocate(stage)
        write(*,'(A)') 'test_rules'
        ! next_job(phase, references, route due, fallback due, assignment due, selected, states)
        call assert_int(JOB_NONE,     stage%next_job(PHASE_IMPORTING, 0, .false., .false., .false., 100, 3), 'nothing due: no job')
        call assert_int(JOB_NONE,     stage%next_job(PHASE_IMPORTING, 0, .true.,  .false., .false., 14,  3), 'fewer than 5 per state: no route')
        call assert_int(JOB_ROUTE,    stage%next_job(PHASE_IMPORTING, 0, .true.,  .false., .false., 15,  3), 'the route once 5 per state are selected')
        call assert_int(JOB_FALLBACK, stage%next_job(PHASE_IMPORTING, 0, .false., .true.,  .false., 15,  3), 'the fallback when it is due')
        call assert_int(JOB_NONE,     stage%next_job(PHASE_CAVGS3D,   0, .true.,  .true.,  .false., 100, 3), 'nothing starts while a job runs')
        call assert_int(JOB_NONE,     stage%next_job(PHASE_IDLE,      2, .false., .false., .false., 100, 3), 'references without a publication: no job')
        call assert_int(JOB_ASSIGN,   stage%next_job(PHASE_IDLE,      2, .false., .false., .true.,  10,  3), 'a publication''s assignment, 5 per reference')
        call assert_int(JOB_NONE,     stage%next_job(PHASE_IDLE,      2, .false., .false., .true.,  9,   3), 'fewer wait')
        call assert_int(JOB_ASSIGN,   stage%next_job(PHASE_IDLE,      2, .true.,  .true.,  .true.,  100, 3), 'with references neither route nor fallback')
        call assert_int(JOB_NONE,     stage%next_job(PHASE_ASSIGN,    2, .false., .false., .true.,  100, 3), 'nothing while the matching runs')
        ! the check: even/odd angular distance, state agreement, distinct states, two states of 5%
        call assert_true(stage%refs_usable(12., 85., [40, 30, 30], 0.5),   'halves that agree, three distinct states')
        call assert_true(stage%refs_usable(50., 70., [50, 50], 0.97),      'at the thresholds')
        call assert_false(stage%refs_usable(55., 85., [40, 30, 30], 0.5),  'halves too far apart')
        call assert_false(stage%refs_usable(12., 60., [40, 30, 30], 0.5),  'states that disagree between the halves')
        call assert_false(stage%refs_usable(-1., 85., [40, 30, 30], 0.5),  'an angular distance not measured')
        call assert_false(stage%refs_usable(12., -1., [40, 30, 30], 0.5),  'an agreement not measured')
        call assert_false(stage%refs_usable(12., 85., [40, 30, 30], 0.99), 'two near-duplicate states')
        call assert_false(stage%refs_usable(12., 85., [96, 4], 0.5),       'a state under 5%: one state counts')
        call assert_true(stage%refs_usable(12., 85., [95, 5], 0.5),        'two states of 5%')
        call assert_false(stage%refs_usable(12., 85., [100], -1.),         'one state')
        ! the best attempt: the smallest angular distance among those with volumes
        call assert_int(3, stage%best_attempt([-1., 30., 25., huge(1.)]), 'the smallest distance')
        call assert_int(0, stage%best_attempt([-1., -1.]),                'none without volumes')
        call assert_int(1, stage%best_attempt([huge(1.), -1.]),           'one not measured, when alone')
        ! the states renumbered from 1
        states = [0, 3, 1, 3, 0]
        call stage%compact_states(states, newind, nstates_held)
        call assert_int(2, nstates_held, 'two states are held')
        call assert_true(all(states == [0, 2, 1, 2, 0]), 'renumbered from 1 in order')
        call assert_true(all(newind == [1, 0, 2]),        'state 2, held by none, gets no number')
        states = [0, 0]
        call stage%compact_states(states, newind, nstates_held)
        call assert_int(0, nstates_held, 'none when no state is held')
        ! cls_expansion's total of sub-classes, 500 and never fewer than the classes
        call assert_int(500, stage%expansion_budget(150), '150 classes: 500 sub-classes in all')
        call assert_int(500, stage%expansion_budget(0),   'no class: 500 (cls_expansion then finds nothing to split)')
        call assert_int(600, stage%expansion_budget(600), '600 classes: one each at least')
        ! the classes cls_expansion splits: distinct labels of the selected particles
        call assert_int(2, stage%count_active_classes([1, 1, 0, 1, 1], [3, 7, 9, 3, 0]), 'classes 3 and 7; a deselected or unlabelled row is not counted')
        call assert_int(0, stage%count_active_classes([0, 0], [1, 2]),                  'none without a selected particle')
        ! the mask diameter that fits the class averages
        call assert_real(100., stage%fit_mskdiam(100.,  64, 2.0), 1.e-4, 'a mask diameter that fits is kept')
        call assert_real(116., stage%fit_mskdiam(400.,  64, 2.0), 1.e-4, 'one too large for the box gets the box default')
        call assert_real(116., stage%fit_mskdiam(0.,    64, 2.0), 1.e-4, 'none gets the box default')
        call assert_real(116., stage%fit_mskdiam(-7.8,  64, 2.0), 1.e-4, 'a negative one too')
        call check_final_publication()
    end subroutine test_rules

    !> the rejection by fit (below the median less 3 MADs, the worst first, at most 30% of the
    !! candidates, never a class that is no candidate) and the ambiguity flag
    subroutine test_fit_rules()
        class(stream_stage_solve3D), allocatable :: stage
        logical, allocatable :: l_rej(:), l_cand(:)
        real    :: med, mad, thres
        integer :: nmax
        real,    parameter   :: CORRS_TWO(10)  = [0.50, 0.52, 0.48, 0.51, 0.49, 0.50, 0.53, 0.47, 0.10, 0.20]
        real,    parameter   :: CORRS_FOUR(10) = [0.50, 0.51, 0.49, 0.50, 0.52, 0.48, 0.10, 0.12, 0.11, 0.13]
        allocate(stage)
        write(*,'(A)') 'test_fit_rules'
        allocate(l_cand(10), source=.true.)
        ! two outliers, under the cap; the statistics for the log
        call stage%reject_by_fit(CORRS_TWO, l_cand, l_rej, med, mad, thres, nmax)
        call assert_int(2, count(l_rej),           'two outliers are rejected')
        call assert_true(l_rej(9) .and. l_rej(10), 'the two that fit worst')
        call assert_real(0.495, med,   1.e-5, 'the candidates'' median')
        call assert_real(0.02,  mad,   1.e-5, 'their median absolute deviation')
        call assert_real(0.435, thres, 1.e-5, 'the threshold, 3 MADs below the median')
        call assert_int(3, nmax,                   'the cap, 30% of 10')
        ! four outliers: the cap of 3 of 10, the worst first
        call stage%reject_by_fit(CORRS_FOUR, l_cand, l_rej)
        call assert_int(3, count(l_rej),                         'at most 30% of the candidates')
        call assert_true(l_rej(7) .and. l_rej(9) .and. l_rej(8), 'the worst three')
        call assert_false(l_rej(10),                             'the fourth is kept')
        ! a class that is no candidate is never rejected, and the statistics are the candidates'
        l_cand(7) = .false.
        call stage%reject_by_fit(CORRS_FOUR, l_cand, l_rej)
        call assert_false(l_rej(7),     'a class that is no candidate is not rejected')
        call assert_int(2, count(l_rej), 'two of the nine candidates, the cap')
        call assert_true(l_rej(9) .and. l_rej(8), 'the worst two candidates')
        ! none when every class fits
        call stage%reject_by_fit([0.50, 0.51, 0.49, 0.50], [.true., .true., .true., .true.], l_rej)
        call assert_int(0, count(l_rej), 'none when every class fits')
        ! no candidate: nothing, and zero statistics
        call stage%reject_by_fit([0.5, 0.1], [.false., .false.], l_rej, med, mad, thres, nmax)
        call assert_int(0, count(l_rej), 'none without candidates')
        call assert_int(0, nmax,         'and no cap')
        ! the margin flags the kept classes below 0.02
        call assert_true(all(stage%flag_ambiguous([0.01, 0.05, 0.0, 0.019], [.true., .true., .false., .true.]) .eqv.&
            &[.true., .false., .false., .true.]), 'small margins of kept classes are ambiguous')
    end subroutine test_fit_rules

    !> solve3D_cavgs' volumes as attempts: the states that hold class averages renumbered from 1
    !! and copied with their reprojections into state_refs/try_<n>/; a failing attempt is kept
    !! aside and the route runs again; a passing one gives the references, which a restart takes
    !! (the comparison of the states after docking: test_reference_docking, a library test)
    subroutine test_reference_volumes()
        class(stream_stage_solve3D), allocatable :: stage
        type(cmdline)    :: cline
        type(sp_project) :: res
        type(string)     :: cwd_saved, root
        integer          :: nfail0
        allocate(stage)
        write(*,'(A)') 'test_reference_volumes'
        nfail0 = tests_failed
        call enter_fixture('a3_stage_reference_volumes', cwd_saved, root)
        call set_test_cline(cline)
        call make_test_stage(stage, cline)
        ! attempt 1: halves 60 degrees apart
        call make_cavgs3D_result(res, 60., 85.)
        call assert_false(stage%register_attempt(res, 1), 'halves 60 degrees apart fail the check')
        call assert_int(1, size(stage%try_angdist), 'the attempt is counted')
        if( size(stage%try_angdist) == 1 ) call assert_real(60., stage%try_angdist(1), 1.e-4, 'with its angular distance')
        call assert_true(file_exists(string('state_refs/try_1/vol_state01.mrc')), 'its volumes are kept aside')
        stage%nroute_tries = 1
        call stage%next_attempt()
        call assert_true(stage%l_route_due, 'the route runs again')
        call assert_int(0, stage%n_refs(),  'without references')
        ! attempt 2 passes: its mean over every class is 60 degrees, but the median over the classes
        ! whose halves are in one state, the distance checked, is 30; solve3D_cavgs' states 1 and 3
        ! hold class averages, its state 2 none
        call make_cavgs3D_result(res, 60., 85., angdist_median=30.)
        call assert_true(stage%register_attempt(res, 2), 'halves in one state that agree pass')
        if( size(stage%try_angdist) == 2 ) call assert_real(30., stage%try_angdist(2), 1.e-4,&
            &'the median over the halves in one state is the distance checked')
        call assert_true(file_exists(string('state_refs/try_2/vol_state01.mrc')),  'its state 1 is state 1')
        call assert_true(file_exists(string('state_refs/try_2/vol_state02.mrc')),  'its state 3 is state 2')
        call assert_false(file_exists(string('state_refs/try_2/vol_state03.mrc')), 'its state 2, empty, is left out')
        call assert_true(file_exists(string('state_refs/try_2/')//refine3D_reprojs_fname(2)), 'its state 3''s reprojections are state 2''s')
        call stage%accept_refs(2, .false.)
        call assert_int(2, stage%n_refs(),          'two references')
        call assert_int(2, stage%params%nstates,    'the stage''s states are theirs')
        call assert_int(PHASE_IDLE, stage%phase,    'the stage waits for a publication')
        call assert_int(VOL_BOX, stage%ref_box,     'their box')
        call assert_real(VOL_SMPD, stage%ref_smpd, 1.e-4, 'and pixel size')
        if( stage%n_refs() == 2 ) call assert_true(stage%ref_vols(2)%has_substr('state_refs/vol_state02.mrc'), 'in state_refs')
        call assert_true(file_exists(string('state_refs/state_refs.txt')), 'with their record')
        ! a restart takes them
        call stage%kill
        call make_test_stage(stage, cline)
        call assert_int(2, stage%n_refs(),       'a restart takes the references')
        call assert_int(2, stage%params%nstates, 'and their states')
        call assert_int(PHASE_IDLE, stage%phase, 'and skips the route')
        call res%kill
        call stage%kill
        call cline%kill
        call leave_fixture(cwd_saved, root, nfail0)
    end subroutine test_reference_volumes

    !> after the route's last attempt the fallback waits for a publication; after the fallback,
    !! or when the pool publishes no more, the attempt with the smallest angular distance gives the
    !! references; without any volumes there are none
    subroutine test_failed_attempts()
        class(stream_stage_solve3D), allocatable :: stage
        type(cmdline)    :: cline
        type(sp_project) :: res
        type(string)     :: cwd_saved, root
        integer          :: nfail0
        allocate(stage)
        write(*,'(A)') 'test_failed_attempts'
        nfail0 = tests_failed
        call enter_fixture('a3_stage_failed_attempts', cwd_saved, root)
        call set_test_cline(cline)
        call make_test_stage(stage, cline)
        ! no attempt with volumes, after the fallback: no references
        stage%try_angdist     = [-1.]
        stage%l_fallback_done = .true.
        call stage%next_attempt()
        call assert_int(0, stage%n_refs(), 'no volumes, no references')
        ! the route's last attempt failed: the fallback waits for a publication
        stage%l_fallback_done = .false.
        stage%nroute_tries    = ROUTE_TRIES
        call make_cavgs3D_result(res, 60., 85.)
        call assert_false(stage%register_attempt(res, 2), 'an attempt that fails the check')
        call stage%next_attempt()
        call assert_false(stage%l_route_due,       'the route does not run again')
        call assert_true(stage%l_fallback_pending, 'the fallback waits for the next publication')
        ! after the fallback, which failed too: the best attempt, with a warning
        stage%l_fallback_pending = .false.
        stage%l_fallback_done    = .true.
        stage%try_angdist        = [stage%try_angdist, -1.]
        call stage%next_attempt()
        call assert_int(2, stage%n_refs(),       'the attempt with volumes gives the references')
        call assert_int(PHASE_IDLE, stage%phase, 'and the stage waits for a publication')
        ! the pool's final publication is in: it is assigned at once
        call assert_false(stage%l_assign_due, 'no assignment without the final publication')
        ! a stage without references (their record gone), whose pool has published its last
        call stage%kill
        call del_file('state_refs/state_refs.txt')
        call make_test_stage(stage, cline)
        call assert_int(0, stage%n_refs(), 'no references to restore')
        stage%l_final_seen = .true.
        stage%nroute_tries = ROUTE_TRIES
        stage%try_angdist  = [-1., 60.]
        call stage%next_attempt()
        call assert_false(stage%l_fallback_pending, 'no fallback when the pool publishes no more')
        call assert_int(2, stage%n_refs(),          'the best attempt gives the references')
        call assert_true(stage%l_assign_due,        'and the final publication is assigned')
        call res%kill
        call stage%kill
        call cline%kill
        call leave_fixture(cwd_saved, root, nfail0)
    end subroutine test_failed_attempts

    !> a matching mapped onto the rows: each particle takes its sub-class's state and pose, so one
    !! pool class may split across states; a sub-class that fits badly is rejected and its
    !! particles deselected until the next publication (a first-set row is selected again by its
    !! merge and judged again); a small margin is flagged; a later matching moves particles with
    !! their sub-classes; a matching of other rows is not taken
    subroutine test_assignment()
        class(stream_stage_solve3D), allocatable :: stage
        type(cmdline)    :: cline
        type(sp_project) :: set, res
        type(string)     :: cwd_saved, root
        integer          :: nfail0
        allocate(stage)
        write(*,'(A)') 'test_assignment'
        nfail0 = tests_failed
        call enter_fixture('a3_stage_assignment', cwd_saved, root)
        call set_test_cline(cline)
        call make_test_stage(stage, cline)
        ! the first set: rows 1 to 5, pool classes 1 (rows 1 to 3) and 2 (rows 4, 5)
        call make_set(set, ['A', 'B'], [3, 2], 2, icls=1)
        call set%os_ptcl2D%set_class(4, 2)
        call set%os_ptcl2D%set_class(5, 2)
        call stage%take_first_publication(set, 1)
        call set%kill
        allocate(stage%ref_vols(2))
        ! sub-classes 1 (rows 1, 2) and 2 (row 3) of pool class 1, 3 (row 4) and 4 (row 5) of class 2;
        ! sub-class 4 fits badly, sub-class 2 has a small margin
        call make_matching(res, [1, 1, 2, 3, 4], [1, 2, 2, 1], [0.50, 0.51, 0.49, 0.10], [0.1, 0.01, 0.1, 0.1])
        call assert_true(stage%take_assignment(res), 'the matching is taken')
        call assert_int(1, stage%spproj%os_ptcl3D%get_state(1), 'sub-class 1''s particles in state 1')
        call assert_int(1, stage%spproj%os_ptcl3D%get_state(2), 'both')
        call assert_int(2, stage%spproj%os_ptcl3D%get_state(3), 'row 3, of the same pool class, in state 2')
        call assert_int(2, stage%spproj%os_ptcl3D%get_state(4), 'sub-class 3''s particle in state 2')
        call assert_int(2, stage%spproj%os_ptcl2D%get_state(3), 'in both particle segments')
        call assert_int(0, stage%spproj%os_ptcl3D%get_state(5), 'the badly fitting sub-class''s particle is deselected')
        call assert_real(45., stage%spproj%os_ptcl3D%get(1, 'e2'), 1.e-3, 'a particle takes its sub-class''s direction')
        call assert_int(0, res%os_cls3D%get_state(4),             'the rejected sub-class has no state')
        call assert_int(1, res%os_cls2D%get_int(2, 'ambiguous'),  'the small margin is flagged')
        call assert_int(0, res%os_cls2D%get_int(1, 'ambiguous'),  'a large one is not')
        call assert_int(4, stage%nptcls_at_last_run,              'four particles have a state')
        call assert_int(4, stage%nptcls_selected,                 'and are selected')
        call res%kill
        call assert_true(stage%fit_deselected(5),  'the deselection is recorded')
        call assert_false(stage%fit_deselected(1), 'an assigned row is not')
        ! a later publication: the first-set row the fit deselected is selected again
        call make_set(set, ['A', 'B'], [3, 2], 2, icls=1)
        call stage%merge_publication(set, 2)
        call set%kill
        call assert_int(1, stage%spproj%os_ptcl3D%get_state(5), 'a first-set row deselected by the fit is selected again')
        call assert_int(1, stage%spproj%os_ptcl2D%get_state(5), 'in both particle segments')
        call assert_false(stage%fit_deselected(5),              'and no longer recorded')
        ! a later matching with the sub-classes' states swapped moves their particles
        call make_matching(res, [1, 1, 2, 3, 4], [2, 1, 1, 2], [0.50, 0.50, 0.50, 0.50], [0.1, 0.1, 0.1, 0.1])
        call assert_true(stage%take_assignment(res), 'the later matching is taken')
        call assert_int(2, stage%spproj%os_ptcl3D%get_state(1), 'a particle follows its sub-class to state 2')
        call assert_int(1, stage%spproj%os_ptcl3D%get_state(3), 'and another to state 1')
        call assert_int(2, stage%spproj%os_ptcl3D%get_state(5), 'the row the fit deselected is judged again, here kept')
        call res%kill
        ! a matching of other rows changes nothing
        call make_matching(res, [1, 1, 2, 3], [1, 1, 1, 1], [0.50, 0.50, 0.50, 0.50], [0.1, 0.1, 0.1, 0.1])
        call assert_false(stage%take_assignment(res), 'a matching of other rows is not taken')
        call assert_int(2, stage%spproj%os_ptcl3D%get_state(1), 'and the rows are unchanged')
        call res%kill
        call stage%kill
        call cline%kill
        call leave_fixture(cwd_saved, root, nfail0)

    contains

        ! a matched copy of the stage's rows: particle i in sub-class @p icls(i) (its 2D alignment
        ! set, as map2ptcls needs), and per sub-class its state, best correlation and margin, its
        ! population, and a direction with e2 = 45
        subroutine make_matching( proj, icls, states, corrs, margins )
            type(sp_project), intent(inout) :: proj
            integer,          intent(in)    :: icls(:), states(:)
            real,             intent(in)    :: corrs(:), margins(:)
            integer :: j, ncls
            ncls = size(states)
            call proj%os_ptcl2D%new(size(icls), is_ptcl=.true.)
            call proj%os_ptcl3D%new(size(icls), is_ptcl=.true.)
            do j = 1,size(icls)
                call proj%os_ptcl2D%set_class(j, icls(j))
                call proj%os_ptcl2D%set(j, 'corr', 0.5)
                call proj%os_ptcl2D%set_state(j, stage%spproj%os_ptcl3D%get_state(j))
                call proj%os_ptcl3D%set_state(j, stage%spproj%os_ptcl3D%get_state(j))
            enddo
            call proj%os_cls2D%new(ncls, is_ptcl=.false.)
            call proj%os_cls3D%new(ncls, is_ptcl=.false.)
            do j = 1,ncls
                call proj%os_cls2D%set(j, 'pop', count(icls == j))
                call proj%os_cls2D%set_state(j, 1)
                call proj%os_cls3D%set_euler(j, [10. * real(j), 45., 0.])
                call proj%os_cls3D%set(j, 'corr',   corrs(j))
                call proj%os_cls3D%set(j, 'margin', margins(j))
                call proj%os_cls3D%set_state(j, states(j))
            enddo
        end subroutine make_matching

    end subroutine test_assignment


    !> the local refinement: the assigned rows are valid refine=cont seeds (state, half-set by row
    !! parity, projection direction); its result comes back under the stage's states, those the
    !! job held renumbered from 1 (an empty state left out), with the refined poses; a result of
    !! other rows is not taken
    subroutine test_local_refinement()
        class(stream_stage_solve3D), allocatable :: stage
        type(cmdline)    :: cline
        type(sp_project) :: set, res
        type(ori)        :: o
        type(image)      :: vol
        type(string)     :: cwd_saved, root, volfile, fscfile
        real, allocatable :: fsc(:)
        integer, parameter :: STATES(5) = [1, 1, 3, 3, 0], JOB_STATES(5) = [1, 1, 2, 2, 0]
        integer          :: nfail0, i, k
        logical          :: l_valid
        allocate(stage)
        write(*,'(A)') 'test_local_refinement'
        nfail0 = tests_failed
        call enter_fixture('a3_stage_local_refinement', cwd_saved, root)
        call set_test_cline(cline)
        call make_test_stage(stage, cline)
        ! assigned rows: states 1, 1, 3, 3 (none in state 2) and one deselected, with mapped poses
        call make_set(set, ['A', 'B'], [3, 2], 2, icls=1)
        call stage%take_first_publication(set, 1)
        call set%kill
        allocate(stage%ref_vols(3))
        do i = 1,5
            call stage%spproj%os_ptcl3D%set_state(i, STATES(i))
            call stage%spproj%os_ptcl3D%set_euler(i, [10. * real(i), 45., 30.])
        enddo
        call stage%prepare_refine_seeds()
        l_valid = .true.
        do i = 1,4
            call stage%spproj%os_ptcl3D%get_ori(i, o)
            l_valid = l_valid .and. cont_seed_is_valid(o)
        enddo
        call o%kill
        call assert_true(l_valid, 'every assigned row is a valid refine=cont seed')
        call assert_int(0, stage%spproj%os_ptcl3D%get_eo(1), 'half-sets by row parity: row 1 even')
        call assert_int(1, stage%spproj%os_ptcl3D%get_eo(2), 'row 2 odd')
        call assert_int(3, stage%spproj%os_ptcl3D%get_state(3), 'and keeps its state')
        ! the job's project: the rows in the job's states 1 and 2 (the stage's 1 and 3), refined
        res%os_ptcl3D = stage%spproj%os_ptcl3D
        do i = 1,5
            call res%os_ptcl3D%set_state(i, JOB_STATES(i))
            call res%os_ptcl3D%set_euler(i, [10. * real(i) + 1., 46., 31.])
        enddo
        allocate(fsc(size(get_resarr(VOL_BOX, VOL_SMPD))), source=0.5)
        do k = 1,2
            volfile = 'refine3D_vol_state'//int2str_pad(k, 2)//'.mrc'
            fscfile = 'refine3D_fsc_state'//int2str_pad(k, 2)//'.bin'
            call vol%new([VOL_BOX, VOL_BOX, VOL_BOX], VOL_SMPD)
            call vol%write(volfile)
            call vol%kill
            call arr2file(fsc, fscfile)
            call res%add_vol2os_out(volfile, VOL_SMPD, k, 'vol')
            call res%add_fsc2os_out(fscfile, k, VOL_BOX)
        enddo
        call assert_int(2, stage%take_refine_result(res), 'two state volumes are taken')
        call assert_true(stage%spproj%isthere_in_osout('vol', 1),  'the job''s state 1 is the stage''s 1')
        call assert_true(stage%spproj%isthere_in_osout('vol', 3),  'the job''s state 2 is the stage''s 3')
        call assert_false(stage%spproj%isthere_in_osout('vol', 2), 'the empty state has no volume')
        call assert_true(stage%spproj%isthere_in_osout('fsc', 3),  'with their FSCs')
        call assert_real(31., stage%spproj%os_ptcl3D%get(3, 'e1'), 1.e-3, 'a row takes its refined pose')
        call assert_real(46., stage%spproj%os_ptcl3D%get(4, 'e2'), 1.e-3, 'every refined row')
        call assert_int(3, stage%spproj%os_ptcl3D%get_state(3), 'its state stays the stage''s')
        call assert_int(0, stage%spproj%os_ptcl3D%get_state(5), 'a deselected row stays deselected')
        ! a result of other rows
        call res%os_ptcl3D%new(4, is_ptcl=.true.)
        call assert_int(0, stage%take_refine_result(res), 'a result of other rows is not taken')
        call assert_true(stage%spproj%isthere_in_osout('vol', 3), 'and the volumes stay')
        call res%kill
        call stage%kill
        call cline%kill
        call leave_fixture(cwd_saved, root, nfail0)
    end subroutine test_local_refinement

    !> the newest quality folders and those the route started from are kept, and so is the latest
    !! assignment; earlier assignments and folders set aside for unfinished jobs go
    subroutine test_retention()
        class(stream_stage_solve3D), allocatable :: stage
        type(string) :: cwd_saved, root
        integer      :: nfail0, i, funit
        allocate(stage)
        write(*,'(A)') 'test_retention'
        nfail0 = tests_failed
        call enter_fixture('a3_stage_retention', cwd_saved, root)
        call simple_mkdir('quality_selection')
        do i = 1,6
            call simple_mkdir('quality_selection/'//int2str_pad(i, 5))
        enddo
        open(newunit=funit, file='quality_selection/runs.txt', action='write')
        write(funit,'(A)') '00002'
        close(funit)
        call stage%prune_quality_dirs()
        call assert_false(dir_exists('quality_selection/00001'), 'an old quality folder goes')
        call assert_true(dir_exists('quality_selection/00002'),  'one the route started from stays')
        call assert_false(dir_exists('quality_selection/00003'), 'like every other old one')
        call assert_true(dir_exists('quality_selection/00004'),  'the newest three stay')
        call assert_true(dir_exists('quality_selection/00006'),  'up to the newest')
        call simple_mkdir('assign')
        do i = 1,3
            call simple_mkdir('assign/it_'//int2str_pad(i, 1))
        enddo
        call simple_mkdir('assign/it_2_unfinished1')
        call simple_mkdir('solve3D_cavgs_unfinished1')
        stage%nassign = 3
        call stage%prune_run_dirs()
        call assert_false(dir_exists('assign/it_1'), 'an earlier assignment goes')
        call assert_false(dir_exists('assign/it_2'), 'every one')
        call assert_true(dir_exists('assign/it_3'),  'the latest stays')
        call assert_false(dir_exists('assign/it_2_unfinished1'),   'a folder set aside for an unfinished job goes')
        call assert_false(dir_exists('solve3D_cavgs_unfinished1'), 'for solve3D_cavgs too')
        call leave_fixture(cwd_saved, root, nfail0)
    end subroutine test_retention

    !> a publication flagged pool_final=yes in its out segment is the pool's final one
    subroutine check_final_publication()
        type(sp_project) :: set
        class(stream_stage_solve3D), allocatable :: stage
        allocate(stage)
        call add_cavgs_entry(set, 150.)
        call assert_false(stage%is_final_publication(set), 'a publication without the flag is not final')
        call set%os_out%set(1, 'pool_final', 'yes')
        call assert_true(stage%is_final_publication(set), 'one with pool_final=yes is')
        call set%kill
    end subroutine check_final_publication

    !> one status message per call, with the stage name, the phase and the particle count
    subroutine test_send_status()
        class(stream_stage_solve3D), allocatable                   :: stage
        type(cmdline)                                   :: cline
        type(stream_pipe)                               :: reader
        type(gui_metadata_stream_solve3D_multistate) :: status
        character(len=:), allocatable                   :: buffer
        type(string)                                    :: cwd_saved, root, stage_name
        integer(c_int)                                  :: fds(2)
        integer :: nfail0, meta_type, solve3D_stage, refine_it, nstates_got, nimported, nlast, tlast
        logical :: l_assigned
        real    :: res
        allocate(stage)
        write(*,'(A)') 'test_send_status'
        nfail0 = tests_failed
        call enter_fixture('a3_stage_status', cwd_saved, root)
        call set_test_cline(cline)
        call make_test_stage(stage, cline)
        call open_loopback(fds)
        call stage%init_gui(-1, int(fds(2)))
        call reader%new(int(fds(1)), -1, max_metadata_size(), 'test reader')
        call stage%send_status()
        call assert_true(reader%receive(buffer), 'a status message is sent')
        if( allocated(buffer) )then
            meta_type = transfer(buffer, meta_type)
            call assert_int(GUI_METADATA_STREAM_SOLVE3D_MULTISTATE_TYPE, meta_type, 'it is a multistate 3D status')
            if( meta_type == GUI_METADATA_STREAM_SOLVE3D_MULTISTATE_TYPE )then
                status     = transfer(buffer, status)
                l_assigned = status%get(stage_name, solve3D_stage, refine_it, nstates_got, nimported, nlast, tlast, res)
                call assert_char('waiting on pool 2D', stage_name%to_char(), 'before pool 2D''s folder exists')
                call assert_int(0,       solve3D_stage, 'no 3D yet')
                call assert_int(NSTATES, nstates_got,      'the number of states')
            endif
        endif
        call assert_false(reader%receive(buffer), 'one message per call')
        call reader%kill
        call stage%kill
        call close_loopback(fds)
        call cline%kill
        call leave_fixture(cwd_saved, root, nfail0)
    end subroutine test_send_status

    !> after a reconstruction, each state with a volume is sent with its reprojection tiles, and its
    !! resolution is kept; a state without a volume is not sent, and a state without an FSC
    !! curve is sent without one (not with the previous state's)
    subroutine test_send_volumes()
        class(stream_stage_solve3D), allocatable :: stage
        type(cmdline)                 :: cline
        type(stream_pipe)             :: reader
        type(image)                   :: vol
        type(gui_metadata_vol3D)      :: meta_vol
        type(gui_metadata_stream_solve3D_multistate) :: status
        type(json_core)               :: json
        type(json_value), pointer     :: json_ptr
        character(len=:), allocatable :: buffer
        real,             allocatable :: res(:), fsc(:), invres_got(:), fsc_got(:)
        type(string)                  :: cwd_saved, root, volfile, fscfile
        integer(c_int)                :: fds(2)
        integer                       :: nfail0, meta_type, nvols, ntiles, iptcl, k, pop_got
        logical                       :: l_fsc_state(2), l_found
        allocate(stage)
        write(*,'(A)') 'test_send_volumes'
        nfail0 = tests_failed
        call enter_fixture('a3_stage_volumes', cwd_saved, root)
        call set_test_cline(cline)
        call make_test_stage(stage, cline)
        call open_loopback(fds)
        call stage%init_gui(-1, int(fds(2)))
        call reader%new(int(fds(1)), -1, max_metadata_size(), 'test reader')
        ! state 1: a volume, its FSC and its reprojections; four particles in it
        volfile = 'recvol_state01.mrc'
        call vol%new([VOL_BOX, VOL_BOX, VOL_BOX], VOL_SMPD)
        call vol%write(volfile)
        call vol%kill
        res = get_resarr(VOL_BOX, VOL_SMPD)
        allocate(fsc(size(res)))
        do k = 1,size(fsc)
            fsc(k) = max(0., 1. - real(k - 1) / real(size(fsc) - 2))
        enddo
        fscfile = 'fsc_state01.bin'
        call arr2file(fsc, fscfile)
        call simple_touch(refine3D_reprojs_fname(1))
        call stage%spproj%add_vol2os_out(volfile, VOL_SMPD, 1, 'vol')
        call stage%spproj%add_fsc2os_out(fscfile, 1, VOL_BOX)
        ! state 2: a volume only; two particles in it
        volfile = 'recvol_state02.mrc'
        call vol%new([VOL_BOX, VOL_BOX, VOL_BOX], VOL_SMPD)
        call vol%write(volfile)
        call vol%kill
        call stage%spproj%add_vol2os_out(volfile, VOL_SMPD, 2, 'vol')
        call stage%spproj%os_ptcl3D%new(6, is_ptcl=.true.)
        do iptcl = 1,6
            call stage%spproj%os_ptcl3D%set_state(iptcl, merge(1, 2, iptcl <= 4))
        enddo
        call stage%send_volumes()
        nvols       = 0
        ntiles      = 0
        l_fsc_state = .false.
        do while( reader%receive(buffer) )
            meta_type = transfer(buffer, meta_type)
            if( meta_type == GUI_METADATA_VOL3D_TYPE )then
                nvols    = nvols + 1
                meta_vol = transfer(buffer, meta_vol)
                if( meta_vol%get_state() >= 1 .and. meta_vol%get_state() <= 2 )then
                    l_fsc_state(meta_vol%get_state()) = meta_vol%get_fsc(invres_got, fsc_got)
                endif
            endif
            if( meta_type == GUI_METADATA_STREAM_SOLVE3D_MULTISTATE_REPROJ_TYPE ) ntiles = ntiles + 1
        enddo
        call assert_int(2, nvols,  'one volume message per state with a volume')
        call assert_int(3, ntiles, 'and the three reprojection tiles of state 1')
        call assert_true(l_fsc_state(1),                    'state 1 is sent with its FSC curve')
        call assert_false(l_fsc_state(2),                   'state 2, without one, is sent without a curve')
        call assert_true(stage%state_res(1) > 0.,          'the resolution of state 1 is kept for the status')
        call assert_real(0., stage%state_res(2), 1.e-6,     'a state without an FSC curve has none')
        call assert_real(0., stage%state_res(3), 1.e-6,     'nor does a state without a volume')
        call assert_int(4, stage%state_pop(1), 'the population of state 1 is kept for the status')
        call assert_int(2, stage%state_pop(2), 'and that of state 2')
        call assert_int(0, stage%state_pop(3), 'a state without a volume has none')
        ! a later merge labels a particle state 1 before any assignment: the status reports the
        ! reconstruction's populations, not the rows' labels
        call stage%spproj%os_ptcl3D%set_state(5, 1)
        allocate(stage%ref_vols(2))
        call stage%send_status()
        l_found = .false.
        pop_got = -1
        do while( reader%receive(buffer) )
            meta_type = transfer(buffer, meta_type)
            if( meta_type /= GUI_METADATA_STREAM_SOLVE3D_MULTISTATE_TYPE ) cycle
            status   = transfer(buffer, status)
            json_ptr => status%jsonise()
            call json%get(json_ptr, 'states(1).population', pop_got, l_found)
            call json%destroy(json_ptr)
        enddo
        call assert_true(l_found,    'the status carries the per-state stats')
        call assert_int(4, pop_got,  'with state 1''s population in the result')
        call reader%kill
        call stage%kill
        call close_loopback(fds)
        call cline%kill
        call leave_fixture(cwd_saved, root, nfail0)
    end subroutine test_send_volumes

    !> the reference volumes go to the GUI under their states, without FSCs, populations or
    !! orientations, with the reprojections copied beside them
    subroutine test_send_refs()
        class(stream_stage_solve3D), allocatable :: stage
        type(cmdline)                 :: cline
        type(sp_project)              :: res
        type(stream_pipe)             :: reader
        type(gui_metadata_vol3D)      :: meta_vol
        character(len=:), allocatable :: buffer
        real,             allocatable :: invres_got(:), fsc_got(:)
        type(string)                  :: cwd_saved, root
        integer(c_int)                :: fds(2)
        integer                       :: nfail0, meta_type, nvols, ntiles
        logical                       :: l_sent(3), l_fsc
        allocate(stage)
        write(*,'(A)') 'test_send_refs'
        nfail0 = tests_failed
        call enter_fixture('a3_stage_send_refs', cwd_saved, root)
        call set_test_cline(cline)
        call make_test_stage(stage, cline)
        call open_loopback(fds)
        call stage%init_gui(-1, int(fds(2)))
        call reader%new(int(fds(1)), -1, max_metadata_size(), 'test reader')
        call make_cavgs3D_result(res, 12., 85.)
        call assert_true(stage%register_attempt(res, 1), 'the attempt passes')
        call stage%accept_refs(1, .false.)
        nvols  = 0
        ntiles = 0
        l_sent = .false.
        l_fsc  = .false.
        do while( reader%receive(buffer) )
            meta_type = transfer(buffer, meta_type)
            if( meta_type == GUI_METADATA_VOL3D_TYPE )then
                nvols    = nvols + 1
                meta_vol = transfer(buffer, meta_vol)
                if( meta_vol%get_state() >= 1 .and. meta_vol%get_state() <= 3 ) l_sent(meta_vol%get_state()) = .true.
                l_fsc = l_fsc .or. meta_vol%get_fsc(invres_got, fsc_got)
            endif
            if( meta_type == GUI_METADATA_STREAM_SOLVE3D_MULTISTATE_REPROJ_TYPE ) ntiles = ntiles + 1
        enddo
        call assert_int(2, nvols,                   'one volume per reference')
        call assert_true(l_sent(1) .and. l_sent(2), 'under states 1 and 2')
        call assert_false(l_sent(3),                'none under solve3D_cavgs'' own state 3')
        call assert_false(l_fsc,                    'without an FSC curve')
        call assert_int(3, ntiles,                  'the reprojection tiles of state 2')
        call assert_int(0, stage%state_pop(1),      'without a population')
        call assert_real(0., stage%state_res(1), 1.e-6, 'or a resolution')
        call assert_false(stage%spproj%isthere_in_osout('vol', 1), 'the stage''s project does not register them')
        call res%kill
        call reader%kill
        call stage%kill
        call close_loopback(fds)
        call cline%kill
        call leave_fixture(cwd_saved, root, nfail0)
    end subroutine test_send_refs

    !> a 3D snapshot request: before a result, answered as not written; from a result, the
    !! particles of the selected states merged into state 1, the state volumes out of the project
    !! and the selected ones copied beside it, answered with the file and the particle count; each
    !! request is written once
    subroutine test_snapshot3D()
        class(stream_stage_solve3D), allocatable :: stage
        type(cmdline)                      :: cline
        type(stream_pipe)                  :: writer, reader
        type(gui_metadata_stream_update)   :: update
        type(gui_metadata_stream_snapshot) :: report
        type(sp_project)                   :: proj
        type(image)                        :: vol
        character(len=:), allocatable      :: buffer
        integer,          allocatable      :: states(:)
        type(string)                       :: cwd_saved, root, fname, volfile
        integer(c_int)                     :: upd(2), gui(2)
        integer                            :: nfail0, id, nptcls, tstamp, istate, iptcl
        integer, parameter                 :: PTCL_STATES(6) = [1, 1, 2, 2, 3, 0]
        allocate(stage)
        write(*,'(A)') 'test_snapshot3D'
        nfail0 = tests_failed
        call enter_fixture('a3_stage_snapshot3D', cwd_saved, root)
        call set_test_cline(cline)
        call make_test_stage(stage, cline)
        call open_loopback(upd)
        call open_loopback(gui)
        call stage%init_gui(int(upd(1)), int(gui(2)))
        call writer%new(-1, int(upd(2)), max_metadata_size(), 'test writer')
        call reader%new(int(gui(1)), -1, max_metadata_size(), 'test reader')
        ! before a result: answered as not written
        call send_request(1, [1])
        call stage%apply_gui_updates()
        call assert_true(receive_report(), 'a request before a result is answered')
        call assert_int(1, id,                 'with its id')
        call assert_int(0, nptcls,             'and no particles')
        call assert_char('', fname%to_char(),  'and no file')
        ! a result: six particles in states 1, 1, 2, 2, 3 and none, and a volume per state
        call proj%os_ptcl2D%new(size(PTCL_STATES), is_ptcl=.true.)
        call proj%os_ptcl3D%new(size(PTCL_STATES), is_ptcl=.true.)
        do iptcl = 1,size(PTCL_STATES)
            call proj%os_ptcl2D%set_state(iptcl, PTCL_STATES(iptcl))
            call proj%os_ptcl3D%set_state(iptcl, PTCL_STATES(iptcl))
        enddo
        do istate = 1,NSTATES
            volfile = 'recvol_state'//int2str_pad(istate, 2)//'.mrc'
            call vol%new([VOL_BOX, VOL_BOX, VOL_BOX], VOL_SMPD)
            call vol%write(volfile)
            call vol%kill
            call proj%add_vol2os_out(volfile, VOL_SMPD, istate, 'vol')
        enddo
        call proj%write(string('result.simple'))
        call proj%kill
        stage%result_projfile = 'result.simple'
        ! states 1 and 3
        call send_request(2, [1, 3])
        call stage%apply_gui_updates()
        call assert_true(receive_report(), 'the request is answered')
        call assert_int(2, id,     'with its id')
        call assert_int(3, nptcls, 'and the particles of states 1 and 3')
        states = report%get_states()
        call assert_int(2, size(states), 'and its states')
        call assert_true(file_exists(fname), 'the snapshot project is written')
        if( file_exists(fname) )then
            call proj%read(fname)
            call assert_int(3, proj%os_ptcl3D%count_state_gt_zero(),           'with the selected particles')
            call assert_int(1, maxval(proj%os_ptcl3D%get_all_asint('state')), 'merged into state 1')
            call assert_int(3, proj%os_ptcl2D%count_state_gt_zero(),           'in both particle segments')
            call assert_false(proj%isthere_in_osout('vol', 1),                 'without the state volumes')
            call proj%kill
        endif
        call assert_true(file_exists(string('snapshots/snapshot_2/vol_state01.mrc')),  'the volume of state 1 is copied')
        call assert_true(file_exists(string('snapshots/snapshot_2/vol_state03.mrc')),  'and that of state 3')
        call assert_false(file_exists(string('snapshots/snapshot_2/vol_state02.mrc')), 'but not that of state 2')
        ! the same request again
        call send_request(2, [1, 3])
        call stage%apply_gui_updates()
        call assert_false(reader%receive(buffer), 'a request is written once')
        call update%kill
        call writer%kill
        call reader%kill
        call stage%kill
        call close_loopback(upd)
        call close_loopback(gui)
        call cline%kill
        call leave_fixture(cwd_saved, root, nfail0)

    contains

        subroutine send_request( snapshot_id, selection )
            integer, intent(in) :: snapshot_id, selection(:)
            call update%kill
            call update%new(GUI_METADATA_STREAM_UPDATE_TYPE)
            call update%set_snapshot3D_update(snapshot_id, selection, string('snapshot_'//int2str(snapshot_id)//METADATA_EXT))
            call writer%send_meta(update)
        end subroutine send_request

        ! the latest 3D snapshot report the stage sent, into report, id, fname and nptcls
        logical function receive_report()
            integer :: meta_type
            receive_report = .false.
            do while( reader%receive(buffer) )
                meta_type = transfer(buffer, meta_type)
                if( meta_type /= GUI_METADATA_STREAM_SOLVE3D_SNAPSHOT_TYPE ) cycle
                report         = transfer(buffer, report)
                receive_report = report%get(id, fname, nptcls, tstamp)
            enddo
        end function receive_report

    end subroutine test_snapshot3D

    !> the public loop waits for pool 2D's folder, then attaches; no job without particles
    subroutine test_iterate_waits()
        class(stream_stage_solve3D), allocatable :: stage
        type(cmdline)                 :: cline
        type(string)                  :: cwd_saved, root
        integer                       :: nfail0
        allocate(stage)
        write(*,'(A)') 'test_iterate_waits'
        nfail0 = tests_failed
        call enter_fixture('a3_stage_iterate', cwd_saved, root)
        call set_test_cline(cline)
        call make_test_stage(stage, cline)
        call stage%iterate()
        call assert_false(stage%l_attached, 'pass 1: no pool 2D folder, not attached')
        call simple_mkdir(UPSTREAM)
        call simple_mkdir(UPSTREAM//'/'//DIR_STREAM_COMPLETED)
        call stage%iterate()
        call assert_true(stage%l_attached,            'pass 2: attached')
        call assert_int(PHASE_IMPORTING, stage%phase, 'pass 2: no particles, no job')
        call stage%kill
        call stage%kill ! idempotence
        call cline%kill
        call leave_fixture(cwd_saved, root, nfail0)
    end subroutine test_iterate_waits

    subroutine test_finished()
        class(stream_stage_solve3D), allocatable :: stage
        type(cmdline)                 :: cline
        type(string)                  :: cwd_saved, root
        integer                       :: nfail0
        allocate(stage)
        write(*,'(A)') 'test_finished'
        nfail0 = tests_failed
        call enter_fixture('a3_stage_finished', cwd_saved, root)
        call set_test_cline(cline)
        call make_test_stage(stage, cline)
        call assert_false(stage%finished(), 'a fresh stage is not finished')
        call simple_touch(TERM_STREAM)
        call assert_true(stage%finished(), 'the termination file finishes the stage')
        call del_file(TERM_STREAM)
        call stage%kill
        call cline%kill
        call leave_fixture(cwd_saved, root, nfail0)
    end subroutine test_finished

    ! ---- fixtures ------------------------------------------------------------

    subroutine set_test_cline( cline )
        type(cmdline), intent(inout) :: cline
        call cline%set('prg',        'stream_solve3D_stage_tester')
        call cline%set('mkdir',      'no')
        call cline%set('projfile',   TEST_PROJFILE)
        call cline%set('outdir',     '')
        call cline%set('nthr',       1)
        call cline%set('nparts',     1)
        call cline%set('nstates',    NSTATES)
        call cline%set('dir_target', UPSTREAM)
        call cline%set('qsys_name',  'local')
    end subroutine set_test_cline

    subroutine make_test_stage( stage, cline )
        class(stream_stage_solve3D), intent(inout) :: stage
        type(cmdline),                 intent(inout) :: cline
        call stage%init_params(cline)
        call stage%init_gui(-1, -1)
        stage%settle_s    = -1
        stage%wait_s      = 0
        stage%l_dock_refs = .false. ! the docking takes too long for a unit test (test_reference_docking)
        stage%l_exists    = .true.
    end subroutine make_test_stage

    subroutine write_export( id, mskdiam )
        integer, intent(in) :: id
        real,    intent(in) :: mskdiam
        type(sp_project) :: proj
        type(string)     :: cwd, fname
        call simple_getcwd(cwd)
        fname = cwd//'/'//UPSTREAM//'/'//DIR_STREAM_COMPLETED//int2str_pad(id, 5)//METADATA_EXT
        call proj%os_out%new(1, is_ptcl=.false.)
        call proj%os_out%set(1, 'imgkind', 'cavg')
        call proj%os_out%set(1, 'stk',     'cavgs.mrc')
        call proj%os_out%set(1, 'nptcls',  10)
        call proj%os_out%set(1, 'mskdiam', mskdiam)
        call proj%update_projinfo(fname)
        call proj%write(fname)
        call proj%kill
    end subroutine write_export

    subroutine make_set( set, stks, nptcls_stk, ncls, nrejected, icls )
        type(sp_project),  intent(inout) :: set
        character(len=1),  intent(in)    :: stks(:)
        integer,           intent(in)    :: nptcls_stk(:), ncls
        integer, optional, intent(in)    :: nrejected, icls
        integer :: istk, iptcl, fromp, nptcls
        nptcls = sum(nptcls_stk)
        call set%os_mic%new(size(stks), is_ptcl=.false.)
        call set%os_stk%new(size(stks), is_ptcl=.false.)
        call set%os_ptcl2D%new(nptcls, is_ptcl=.true.)
        call set%os_ptcl3D%new(nptcls, is_ptcl=.true.)
        call set%os_cls2D%new(ncls, is_ptcl=.false.)
        fromp = 1
        do istk = 1,size(stks)
            call set%os_mic%set_state(istk, 1)
            call set%os_stk%set(istk, 'stk',    'stack_'//stks(istk)//'.mrcs')
            call set%os_stk%set(istk, 'nptcls', nptcls_stk(istk))
            call set%os_stk%set(istk, 'fromp',  fromp)
            call set%os_stk%set(istk, 'top',    fromp + nptcls_stk(istk) - 1)
            do iptcl = fromp,fromp + nptcls_stk(istk) - 1
                call set%os_ptcl2D%set_state(iptcl, 1)
                call set%os_ptcl3D%set_state(iptcl, 1)
                call set%os_ptcl3D%set(iptcl, 'x', real(iachar(stks(istk)) - iachar('A') + 1))
                if( present(icls) ) call set%os_ptcl2D%set_class(iptcl, icls)
            enddo
            fromp = fromp + nptcls_stk(istk)
        enddo
        if( present(nrejected) )then
            do iptcl = 1,nrejected
                call set%os_ptcl2D%set_state(iptcl, 0)
                call set%os_ptcl3D%set_state(iptcl, 0)
            enddo
        endif
    end subroutine make_set

    subroutine add_cavgs_entry( set, mskdiam )
        type(sp_project), intent(inout) :: set
        real,             intent(in)    :: mskdiam
        call set%os_out%new(1, is_ptcl=.false.)
        call set%os_out%set(1, 'imgkind', 'cavg')
        call set%os_out%set(1, 'stk',     'cavgs.mrc')
        call set%os_out%set(1, 'nptcls',  10)
        if( mskdiam > 0. ) call set%os_out%set(1, 'mskdiam', mskdiam)
    end subroutine add_cavgs_entry

    ! solve3D_cavgs' result in memory: four class averages in its states 3, none, 1 and 3 (its
    ! state 2 holds none), a volume for each of its states 1 to 3 with the even/odd values
    ! @p angdist (the mean over every class) and @p overlap, with @p angdist_median (the median
    ! over the classes whose halves are in one state) when given, and the reprojections of its
    ! state 3. Each state's volume is a wave along its own axis (the fast tests compare no states:
    ! no docking).
    subroutine make_cavgs3D_result( res, angdist, overlap, angdist_median )
        type(sp_project),  intent(inout) :: res
        real,              intent(in)    :: angdist, overlap
        real, optional,    intent(in)    :: angdist_median
        integer, parameter :: CLS_STATES(4) = [3, 0, 1, 3]
        real,    parameter :: TWOPI = 6.2831853
        type(image)  :: vol
        type(string) :: volfile
        real, allocatable :: rmat(:,:,:)
        integer      :: s, i, x, y, z, r(3)
        allocate(rmat(VOL_BOX,VOL_BOX,VOL_BOX))
        call res%kill
        call res%os_cls3D%new(size(CLS_STATES), is_ptcl=.false.)
        do i = 1,size(CLS_STATES)
            call res%os_cls3D%set(i, 'e1', 30.)
            call res%os_cls3D%set_state(i, CLS_STATES(i))
        enddo
        call simple_mkdir('cavgs3D')
        do s = 1,3
            do z = 1,VOL_BOX
                do y = 1,VOL_BOX
                    do x = 1,VOL_BOX
                        r = [x, y, z] - (VOL_BOX / 2 + 1)
                        rmat(x,y,z) = cos(TWOPI * 2. * real(r(s)) / real(VOL_BOX))
                    enddo
                enddo
            enddo
            volfile = 'cavgs3D/rec_final_state'//int2str_pad(s, 2)//'.mrc'
            call vol%new([VOL_BOX, VOL_BOX, VOL_BOX], VOL_SMPD)
            call vol%set_rmat(rmat, .false.)
            call vol%write(volfile)
            call vol%kill
            call res%add_vol2os_out(volfile, VOL_SMPD, s, 'vol_cavg')
        enddo
        call simple_touch(string('cavgs3D/')//refine3D_reprojs_fname(3))
        do i = 1,res%os_out%get_noris()
            call res%os_out%set(i, 'eo_angdist',       angdist)
            call res%os_out%set(i, 'eo_state_overlap', overlap)
            if( present(angdist_median) ) call res%os_out%set(i, 'eo_angdist_median', angdist_median)
        enddo
    end subroutine make_cavgs3D_result

    subroutine open_loopback( fds )
        integer(c_int), intent(inout) :: fds(2)
        integer(c_int) :: rc, flags
        fds   = -1
        rc    = c_pipe(fds)
        call assert_int(0, int(rc), 'pipe() succeeds')
        flags = c_fcntl(fds(1), F_GETFL, 0_c_int)
        rc    = c_fcntl(fds(1), F_SETFL, ior(flags, O_NONBLOCK))
        call assert_int(0, int(rc), 'the read end is made non-blocking')
    end subroutine open_loopback

    subroutine close_loopback( fds )
        integer(c_int), intent(inout) :: fds(2)
        integer(c_int) :: rc
        rc = c_close(fds(1))
        call assert_int(0, int(rc), 'the read end closes')
        rc = c_close(fds(2))
        call assert_int(0, int(rc), 'the write end closes')
        fds = -1
    end subroutine close_loopback

end module simple_stream_stage_solve3D_tester
