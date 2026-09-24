!@descr: library tests of the in-process stream stages: optics assignment, picking references, pick and extract
! Moved from Ruben's stream test commanders (assign_optics, gen_pickrefs, pick_extract; 2026-08)
! by the stream review (plan, section 9.7): the checks are his, now assertions, so one failed
! check no longer ends the suite. Each test builds its fixture in a fresh directory under the
! suite's working directory, runs the production stage or commander with the arguments the
! stream gives it, and removes the directory when every check passed (it is kept, and named in
! the log, when one failed). The stream's preprocessing stage submits worker jobs and stays a
! workflow entry (stream_preproc); doc/refactoring_notes/stream_area_tests_handover.md says what
! these tests should pin beyond counts and files.
module simple_stream_tester
use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
use simple_core_module_api
use simple_cmdline,    only: cmdline
use simple_sp_project, only: sp_project
use simple_image,      only: image
use simple_test_utils
implicit none
private
public :: run_all_stream_optics_tests, run_all_stream_pickrefs_tests, run_all_stream_pick_extract_tests

contains

    subroutine run_all_stream_optics_tests()
        write(*,'(A)') '**** running all stream optics assignment tests ****'
        call test_assign_optics_two_shift_clusters()
    end subroutine run_all_stream_optics_tests

    subroutine run_all_stream_pickrefs_tests()
        write(*,'(A)') '**** running all stream picking reference tests ****'
        call test_make_pickrefs_expansion()
    end subroutine run_all_stream_pickrefs_tests

    subroutine run_all_stream_pick_extract_tests()
        write(*,'(A)') '**** running all stream pick and extract tests ****'
        call test_pick_extract_three_particles()
    end subroutine run_all_stream_pick_extract_tests

    !> the optics-assignment stage (stream_p02_assign_optics) on a completed preprocessing batch
    !! of five micrographs whose beam-image shifts form two clusters, near (0,0) twice and near
    !! (5,5) three times: two optics groups with those populations and centroids, every
    !! micrograph in its cluster's group, and the STAR files and the optics map written. The
    !! stage imports a project only once it is LONGTIME (60 s) old, so this takes a minute.
    subroutine test_assign_optics_two_shift_clusters()
        use simple_stream_p02_assign_optics_new, only: stream_p02_assign_optics
        character(len=*), parameter :: TEST_PROJFILE = 'test_assign_optics.simple'
        character(len=*), parameter :: TEST_OUTDIR   = 'stream_assign_optics'
        real,             parameter :: SMPD  = 1.3
        real,             parameter :: CS    = 2.7
        real,             parameter :: KV    = 300.0
        real,             parameter :: FRACA = 0.1
        real,             parameter :: SHIFT_X(STREAM_NMOVS_SET) = [0.0, 0.1, 5.0, 5.1, 4.9]
        real,             parameter :: SHIFT_Y(STREAM_NMOVS_SET) = [0.0,-0.1, 5.0, 4.9, 5.1]
        integer,          parameter :: TILT_GROUP(STREAM_NMOVS_SET) = [1, 1, 2, 2, 2]
        type(stream_p02_assign_optics) :: xassign_optics
        type(cmdline)                  :: cline_assign
        type(sp_project)               :: upstream, result
        type(string)                   :: cwd_saved, fixture_root, preproc_root, completed_dir
        type(string)                   :: batch_path, project_path, output_dir, fname
        integer                        :: i, nfail0, group_a, group_b
        logical                        :: l_groups
        write(*,'(A)') 'test_assign_optics_two_shift_clusters'
        nfail0 = tests_failed
        call enter_fixture('test_assign_optics', cwd_saved, fixture_root)
        ! the directories and the completed five-micrograph project that preprocessing emits
        preproc_root = filepath(fixture_root, 'preproc')
        call simple_mkdir(preproc_root)
        call simple_mkdir(filepath(preproc_root, 'spprojs'))
        completed_dir = filepath(preproc_root, 'spprojs_completed')
        call simple_mkdir(completed_dir)
        batch_path = filepath(completed_dir, '00001.simple')
        call upstream%os_mic%new(STREAM_NMOVS_SET, is_ptcl=.false.)
        do i = 1, STREAM_NMOVS_SET
            fname = 'synthetic_'//int2str_pad(i, 3)
            call upstream%os_mic%set(i, 'movie',       filepath(preproc_root, fname//'.mrc'))
            call upstream%os_mic%set(i, 'intg',        filepath(preproc_root, fname//'_intg.mrc'))
            call upstream%os_mic%set(i, 'mc_starfile', filepath(preproc_root, fname//'.star'))
            call upstream%os_mic%set(i, 'ctfjpg',      filepath(preproc_root, fname//'_ctf.jpg'))
            call upstream%os_mic%set(i, 'boxfile',     filepath(preproc_root, fname//'.box'))
            call upstream%os_mic%set(i, 'imgkind',     'mic')
            call upstream%os_mic%set(i, 'state',       1.0)
            call upstream%os_mic%set(i, 'importind',   real(i))
            call upstream%os_mic%set(i, 'xdim',        512.0)
            call upstream%os_mic%set(i, 'ydim',        512.0)
            call upstream%os_mic%set(i, 'nframes',     8.0)
            call upstream%os_mic%set(i, 'smpd',        SMPD)
            call upstream%os_mic%set(i, 'cs',          CS)
            call upstream%os_mic%set(i, 'kv',          KV)
            call upstream%os_mic%set(i, 'fraca',       FRACA)
            call upstream%os_mic%set(i, 'dfx',         1.5 + 0.1 * real(i))
            call upstream%os_mic%set(i, 'dfy',         1.6 + 0.1 * real(i))
            call upstream%os_mic%set(i, 'angast',      10.0 * real(i))
            call upstream%os_mic%set(i, 'phshift',     0.0)
            call upstream%os_mic%set(i, 'ctfres',      4.0)
            call upstream%os_mic%set(i, 'icefrac',     0.1)
            call upstream%os_mic%set(i, 'astig',       0.1)
            call upstream%os_mic%set(i, 'shiftx',      SHIFT_X(i))
            call upstream%os_mic%set(i, 'shifty',      SHIFT_Y(i))
            call upstream%os_mic%set(i, 'tiltgrp',     real(TILT_GROUP(i)))
        enddo
        call upstream%update_projinfo(batch_path)
        call upstream%write(batch_path)
        call upstream%kill
        ! the stage, as the stream master starts it, stopping after the five micrographs
        project_path = filepath(fixture_root, TEST_PROJFILE)
        call cline_assign%set('prg',        'assign_optics')
        call cline_assign%set('projfile',   project_path)
        call cline_assign%set('dir_target', preproc_root)
        call cline_assign%set('outdir',     TEST_OUTDIR)
        call cline_assign%set('nthr',       1)
        call cline_assign%set('nmics',      STREAM_NMOVS_SET)
        call cline_assign%set('beamtilt',   'yes')
        call cline_assign%set('tilt_thres', 0.5)
        call xassign_optics%execute(cline_assign)
        ! the stage works in its own directory
        project_path = cline_assign%get_carg('projfile')
        call simple_getcwd(output_dir)
        call simple_chdir(fixture_root)
        call assert_true(file_exists(project_path), 'the output project exists')
        if( file_exists(project_path) )then
            call result%read(project_path)
            call assert_int(STREAM_NMOVS_SET, result%os_mic%get_noris(),    'the project holds the five micrographs')
            call assert_int(2,                result%os_optics%get_noris(), 'two optics groups')
            if( result%os_mic%get_noris() == STREAM_NMOVS_SET .and. result%os_optics%get_noris() == 2 )then
                group_a  = result%os_mic%get_int(1, 'ogid')
                group_b  = result%os_mic%get_int(3, 'ogid')
                l_groups = group_a /= group_b .and. all([group_a, group_b] >= 1) .and. all([group_a, group_b] <= 2)
                call assert_true(l_groups, 'the two shift clusters are in different optics groups')
                if( l_groups )then
                    call assert_int(group_a, result%os_mic%get_int(2, 'ogid'), 'micrograph 2 is in the (0,0) group')
                    call assert_int(group_b, result%os_mic%get_int(4, 'ogid'), 'micrograph 4 is in the (5,5) group')
                    call assert_int(group_b, result%os_mic%get_int(5, 'ogid'), 'micrograph 5 is in the (5,5) group')
                    call assert_int(2, nint(result%os_optics%get(group_a, 'pop')), 'the (0,0) group holds two micrographs')
                    call assert_int(3, nint(result%os_optics%get(group_b, 'pop')), 'the (5,5) group holds three micrographs')
                    call assert_real( 0.05, result%os_optics%get(group_a, 'opcx'), 0.01, 'the (0,0) group centroid, x')
                    call assert_real(-0.05, result%os_optics%get(group_a, 'opcy'), 0.01, 'the (0,0) group centroid, y')
                    call assert_real( 5.00, result%os_optics%get(group_b, 'opcx'), 0.01, 'the (5,5) group centroid, x')
                    call assert_real( 5.00, result%os_optics%get(group_b, 'opcy'), 0.01, 'the (5,5) group centroid, y')
                endif
            endif
            call result%kill
        endif
        call assert_true(file_exists(filepath(output_dir, 'optics.star')),      'optics.star is written')
        call assert_true(file_exists(filepath(output_dir, 'micrographs.star')), 'micrographs.star is written')
        call assert_true(file_exists(filepath(output_dir, OPTICS_MAP_PREFIX//'1'//TXT_EXT)), &
            &'the optics map table is written')
        call assert_true(file_exists(filepath(output_dir, OPTICS_MAP_PREFIX//'1'//METADATA_EXT)), &
            &'the optics map project is written')
        call cline_assign%kill
        call leave_fixture(cwd_saved, fixture_root, nfail0)
    end subroutine test_assign_optics_two_shift_clusters

    !> make_pickrefs, as the stream's reference-generation stage runs it, on three asymmetric
    !! class-average-like images of 64^2: 3 x 4 rotations x 2 mirrors references, square, no larger
    !! than the source box, at the source sampling; the diameter metadata the picking stage reads
    !! (diam_max, mskdiam, box_for_pick = the reference box, box_for_extract >= box_for_pick); the
    !! source preview and the completion marker
    subroutine test_make_pickrefs_expansion()
        use simple_commanders_pick, only: commander_make_pickrefs
        character(len=*), parameter :: SOURCE_STACK   = 'source_references.mrcs'
        character(len=*), parameter :: TEST_OUTDIR    = 'stream_gen_pickrefs'
        character(len=*), parameter :: OUTPUT_STACK   = PICKREFS_FBODY//MRC_EXT
        character(len=*), parameter :: SOURCE_PREVIEW = 'pickrefs_source.jpeg'
        character(len=*), parameter :: FINISHED_FILE  = 'MAKE_PICKREFS_FINISHED'
        real,             parameter :: SMPD           = 2.0
        integer,          parameter :: SOURCE_BOX     = 64
        integer,          parameter :: NSOURCE_REFS   = 3
        integer,          parameter :: NROTS          = 4
        type(commander_make_pickrefs) :: xmake_pickrefs
        type(cmdline)                 :: cline_make_pickrefs
        type(image)                   :: reference, feature
        type(oris)                    :: moldiam_meta
        type(string)                  :: cwd_saved, fixture_root, source_path, output_dir, output_path, moldiam_path
        integer                       :: i, nfail0, ldim(3), nrefs, box_for_pick, box_for_extract
        real                          :: diam_max, mskdiam
        logical                       :: l_fields
        write(*,'(A)') 'test_make_pickrefs_expansion'
        nfail0 = tests_failed
        call enter_fixture('test_gen_pickrefs', cwd_saved, fixture_root)
        ! three compact asymmetric images: a centred square and an off-centre feature
        do i = 1, NSOURCE_REFS
            call reference%new([SOURCE_BOX, SOURCE_BOX, 1], SMPD, wthreads=.false.)
            call reference%square(5 + i)
            call feature%new([SOURCE_BOX, SOURCE_BOX, 1], SMPD, wthreads=.false.)
            call feature%square(2 + i)
            call feature%shift([real(7 + 2 * i), real(-5 + i), 0.])
            call reference%add(feature)
            call reference%write(string(SOURCE_STACK), i, del_if_exists=(i == 1))
            call feature%kill
            call reference%kill
        enddo
        source_path = simple_abspath(string(SOURCE_STACK))
        call cline_make_pickrefs%set('prg',          'make_pickrefs')
        call cline_make_pickrefs%set('pickrefs',     source_path)
        call cline_make_pickrefs%set('smpd',         SMPD)
        call cline_make_pickrefs%set('nrots',        NROTS)
        call cline_make_pickrefs%set('mirr',         'yes')
        call cline_make_pickrefs%set('trust_header', 'yes')
        call cline_make_pickrefs%set('ncls',         0)
        call cline_make_pickrefs%set('nthr',         1)
        call cline_make_pickrefs%set('mkdir',        'yes')
        call cline_make_pickrefs%set('dir_exec',     TEST_OUTDIR)
        call xmake_pickrefs%execute(cline_make_pickrefs)
        ! the commander works in its own directory
        call simple_getcwd(output_dir)
        output_path  = filepath(output_dir, OUTPUT_STACK)
        moldiam_path = filepath(output_dir, STREAM_MOLDIAM)
        ldim = 0
        call assert_true(file_exists(output_path), 'the picking-reference stack exists')
        if( file_exists(output_path) )then
            call find_ldim_nptcls(output_path, ldim, nrefs)
            call assert_int(NSOURCE_REFS * NROTS * 2, nrefs, 'every source image in every rotation, mirrored and not')
            call assert_true(ldim(1) > 0 .and. ldim(1) == ldim(2) .and. ldim(1) <= SOURCE_BOX, &
                &'the references are square and no larger than the source box')
            call assert_real(SMPD, find_img_smpd(output_path), 0.01, 'the references keep the source sampling')
        endif
        call assert_true(file_exists(moldiam_path), 'the diameter metadata exists')
        if( file_exists(moldiam_path) )then
            call moldiam_meta%new(1, .false.)
            call moldiam_meta%read(moldiam_path)
            l_fields = moldiam_meta%isthere(1, 'diam_max')     .and. moldiam_meta%isthere(1, 'mskdiam') .and. &
                      &moldiam_meta%isthere(1, 'box_for_pick') .and. moldiam_meta%isthere(1, 'box_for_extract')
            call assert_true(l_fields, 'the diameter metadata holds diam_max, mskdiam, box_for_pick, box_for_extract')
            if( l_fields )then
                diam_max        = moldiam_meta%get(1, 'diam_max')
                mskdiam         = moldiam_meta%get(1, 'mskdiam')
                box_for_pick    = moldiam_meta%get_int(1, 'box_for_pick')
                box_for_extract = moldiam_meta%get_int(1, 'box_for_extract')
                call assert_true(diam_max > 0., 'diam_max is positive')
                call assert_true(mskdiam  > 0., 'mskdiam is positive')
                call assert_int(ldim(1), box_for_pick, 'box_for_pick is the reference box')
                call assert_true(box_for_extract >= box_for_pick, 'box_for_extract is at least box_for_pick')
            endif
            call moldiam_meta%kill
        endif
        call assert_true(file_exists(filepath(output_dir, SOURCE_PREVIEW)), 'the source preview exists')
        call assert_true(file_exists(filepath(output_dir, FINISHED_FILE)),  'the completion marker exists')
        call cline_make_pickrefs%kill
        call leave_fixture(cwd_saved, fixture_root, nfail0)
    end subroutine test_make_pickrefs_expansion

    !> pick_extract, as the stream's reference-picking stage runs it on one micrograph: three
    !! copies of an asymmetric reference on weak noise in a 256^2 micrograph are picked and
    !! extracted: the micrograph kept, three particles in the micrograph record, the 2D and 3D
    !! particle segments and the extracted stack (64^2, source sampling, finite non-zero data),
    !! the box file and the thumbnail. nboxes_max=3 caps the picks at the number asserted, so the
    !! count cannot catch over-picking; the handover asks for the positions.
    subroutine test_pick_extract_three_particles()
        use simple_commanders_pick, only: commander_pick_extract
        character(len=*), parameter :: MICROGRAPH_FILE = 'synthetic_micrograph.mrc'
        character(len=*), parameter :: PICKREFS_FILE   = 'synthetic_pickrefs.mrcs'
        character(len=*), parameter :: TEST_PROJFILE   = 'test_pick_extract.simple'
        real,             parameter :: SMPD            = 2.0
        integer,          parameter :: MICROGRAPH_BOX  = 256
        integer,          parameter :: PICKREF_BOX     = 64
        integer,          parameter :: EXTRACT_BOX     = 64
        integer,          parameter :: NPARTICLES      = 3
        integer,          parameter :: PARTICLE_COORDS(2, NPARTICLES) = reshape([32, 32, 160, 40, 96, 152], [2, NPARTICLES])
        type(commander_pick_extract) :: xpick_extract
        type(cmdline)                :: cline_pick_extract
        type(image)                  :: micrograph, reference, feature, extracted
        type(sp_project)             :: fixture, result
        type(string)                 :: cwd_saved, fixture_root, micrograph_path, pickrefs_path
        type(string)                 :: project_path, boxfile_path, thumbnail_path, stack_path
        integer                      :: i, nfail0, npicked, nimages, ldim(3)
        real                         :: extracted_variance
        write(*,'(A)') 'test_pick_extract_three_particles'
        nfail0 = tests_failed
        call enter_fixture('test_pick_extract', cwd_saved, fixture_root)
        ! an asymmetric reference, three exact copies on a weak-noise micrograph
        call reference%new([PICKREF_BOX, PICKREF_BOX, 1], SMPD, wthreads=.false.)
        call reference%square(7)
        call feature%new([PICKREF_BOX, PICKREF_BOX, 1], SMPD, wthreads=.false.)
        call feature%square(3)
        call feature%shift([11., -7., 0.])
        call reference%add(feature)
        call reference%mul(5.)
        call reference%write(string(PICKREFS_FILE), 1, del_if_exists=.true.)
        call micrograph%new([MICROGRAPH_BOX, MICROGRAPH_BOX, 1], SMPD, wthreads=.false.)
        call micrograph%gauran(0., 0.02)
        do i = 1, NPARTICLES
            call micrograph%add_window(reference, PARTICLE_COORDS(:, i))
        enddo
        call micrograph%write(string(MICROGRAPH_FILE), 1, del_if_exists=.true.)
        call feature%kill
        call reference%kill
        call micrograph%kill
        micrograph_path = simple_abspath(string(MICROGRAPH_FILE))
        pickrefs_path   = simple_abspath(string(PICKREFS_FILE))
        ! the single-micrograph project a stream pick_extract job receives
        project_path = filepath(fixture_root, TEST_PROJFILE)
        call fixture%os_mic%new(1, is_ptcl=.false.)
        call fixture%os_mic%set(1, 'movie',     micrograph_path)
        call fixture%os_mic%set(1, 'intg',      micrograph_path)
        call fixture%os_mic%set(1, 'imgkind',   'mic')
        call fixture%os_mic%set(1, 'state',     1.0)
        call fixture%os_mic%set(1, 'importind', 1.0)
        call fixture%os_mic%set(1, 'xdim',      real(MICROGRAPH_BOX))
        call fixture%os_mic%set(1, 'ydim',      real(MICROGRAPH_BOX))
        call fixture%os_mic%set(1, 'nframes',   1.0)
        call fixture%os_mic%set(1, 'smpd',      SMPD)
        call fixture%update_projinfo(project_path)
        call fixture%write(project_path)
        call fixture%kill
        call cline_pick_extract%set('prg',          'pick_extract')
        call cline_pick_extract%set('projfile',     project_path)
        call cline_pick_extract%set('pickrefs',     pickrefs_path)
        call cline_pick_extract%set('dir',          fixture_root)
        call cline_pick_extract%set('stream',       'yes')
        call cline_pick_extract%set('extract',      'yes')
        call cline_pick_extract%set('picker',       'new')
        call cline_pick_extract%set('pcontrast',    'white')
        call cline_pick_extract%set('pick_roi',     'no')
        call cline_pick_extract%set('backgr_subtr', 'no')
        call cline_pick_extract%set('nboxes_max',   NPARTICLES)
        call cline_pick_extract%set('box_extract',  EXTRACT_BOX)
        call cline_pick_extract%set('fromp',        1)
        call cline_pick_extract%set('top',          1)
        call cline_pick_extract%set('nthr',         1)
        call cline_pick_extract%set('mkdir',        'no')
        call xpick_extract%execute(cline_pick_extract)
        call simple_chdir(fixture_root)
        call assert_true(file_exists(project_path), 'the output project exists')
        if( file_exists(project_path) )then
            call result%read(project_path)
            call assert_int(1, result%os_mic%get_noris(), 'the project holds the micrograph')
            if( result%os_mic%get_noris() == 1 )then
                call assert_int(1, result%os_mic%get_state(1), 'the micrograph is kept')
                npicked = result%os_mic%get_int(1, 'nptcls')
                call assert_int(NPARTICLES, npicked, 'three particles are picked')
                boxfile_path   = result%os_mic%get_str(1, 'boxfile')
                thumbnail_path = result%os_mic%get_str(1, 'thumb_den')
                call assert_true(file_exists(boxfile_path),   'the box file exists')
                call assert_true(file_exists(thumbnail_path), 'the density thumbnail exists')
                call assert_int(npicked, result%os_ptcl2D%get_noris(), 'one 2D particle record per pick')
                call assert_int(npicked, result%os_ptcl3D%get_noris(), 'one 3D particle record per pick')
                call assert_int(1, result%os_stk%get_noris(), 'one extracted stack record')
                if( result%os_stk%get_noris() == 1 )then
                    stack_path = result%os_stk%get_str(1, 'stk')
                    call assert_true(file_exists(stack_path), 'the extracted stack exists')
                    if( file_exists(stack_path) )then
                        call find_ldim_nptcls(stack_path, ldim, nimages)
                        call assert_int(npicked, nimages, 'one extracted image per pick')
                        call assert_true(all(ldim(1:2) == [EXTRACT_BOX, EXTRACT_BOX]), 'the extracted particles are 64^2')
                        call assert_real(SMPD, find_img_smpd(stack_path), 0.01, 'extraction keeps the sampling')
                        if( nimages >= 1 .and. all(ldim(1:2) == [EXTRACT_BOX, EXTRACT_BOX]) )then
                            call extracted%new([EXTRACT_BOX, EXTRACT_BOX, 1], SMPD, wthreads=.false.)
                            call extracted%read(stack_path, 1)
                            extracted_variance = extracted%variance()
                            call assert_true(ieee_is_finite(extracted_variance) .and. extracted_variance > TINY, &
                                &'the first extracted particle holds finite, non-constant data')
                            call extracted%kill
                        endif
                    endif
                endif
            endif
            call result%kill
        endif
        call cline_pick_extract%kill
        call leave_fixture(cwd_saved, fixture_root, nfail0)
    end subroutine test_pick_extract_three_particles

    ! ---- fixture directories ---------------------------------------------------

    !> makes and enters a fresh directory for one test under the working directory
    subroutine enter_fixture( tag, cwd_saved, fixture_root )
        character(len=*), intent(in)    :: tag
        type(string),     intent(inout) :: cwd_saved, fixture_root
        call simple_getcwd(cwd_saved)
        fixture_root = filepath(cwd_saved, tag//'_'//int2str(get_process_id()))
        if( dir_exists(fixture_root) ) call simple_rmdir(fixture_root)
        call simple_mkdir(fixture_root)
        call simple_chdir(fixture_root)
    end subroutine enter_fixture

    !> returns to the working directory; the fixture goes when no check failed since nfail_before
    subroutine leave_fixture( cwd_saved, fixture_root, nfail_before )
        type(string), intent(in) :: cwd_saved, fixture_root
        integer,      intent(in) :: nfail_before
        call simple_chdir(cwd_saved)
        if( tests_failed == nfail_before )then
            call simple_rmdir(fixture_root)
        else
            write(logfhandle,'(A)') '>>> a check failed; the fixture is kept: '//fixture_root%to_char()
        endif
    end subroutine leave_fixture

end module simple_stream_tester
