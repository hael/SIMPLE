!@descr: chained tests of the stream stages without the master: particle sieving to multistate 3D from simulated particles, and preprocessing to multistate 3D from simulated movies
!==============================================================================
! MODULE: simple_stream_chain_tester
!
! PURPOSE:
!   Decision 19 of the 5 October 2026 stream fix plan. Each test drives the
!   production stage types in one process, without the master: every stage is
!   made with its own new(), on its commander's command-line defaults (the
!   commanders' set_*_cline routines), in a folder of its own in the fixture.
!   The stages run under program names outside every UI table, so params%new
!   makes no folder of its own, and the driver enters a stage's folder before
!   each of its passes. Outside the master the stages' pipes are closed, so
!   they report to no one. A loop gives every stage a pass in turn until the
!   last stage has a result or a time limit passes; then each stage is
!   finalised, which cancels its jobs.
!
!   The stages submit their jobs (2D chunks, the pool's iterations, the 3D
!   runs, preprocessing, picking) to the local queue, so the SIMPLE
!   executables must be installed and SIMPLE_PATH set, as for the
!   stream_preproc workflow test. Each test takes tens of minutes.
!
!   - sieve to 3D (lib_stream, "sieve to 3D"): simulated particle stacks in
!     completed reference-picking sets, with the picking references' mask
!     diameter and reference picking's finished marker; p05, p06, p07.
!   - movies to 3D (lib_stream, "movies to 3D"): simulated movies; p01, p02,
!     p04 with picking references reprojected from both truths, p05, p06, p07.
!     Preprocessing is finalised once every movie is processed, so its
!     finished marker ends the intake downstream.
!
!   NOT YET RUN: written on 5 October 2026 without a build; the first run may
!   need its fixture or time limits adjusted.
!==============================================================================
module simple_stream_chain_tester
use simple_core_module_api
use simple_cmdline,                                  only: cmdline
use simple_sp_project,                               only: sp_project
use simple_image,                                    only: image
use simple_test_utils
use simple_ui,                                       only: make_ui
use simple_commanders_sim,                           only: commander_simulate_particles, commander_simulate_movie
use simple_commanders_reproject,                     only: commander_reproject
use simple_stream_stage_preprocess,                  only: stream_stage_preprocess
use simple_stream_stage_optics,                      only: stream_stage_optics
use simple_stream_stage_refpick,                     only: stream_stage_refpick
use simple_stream_stage_sieve,                       only: stream_stage_sieve
use simple_stream_stage_pool2D,                      only: stream_stage_pool2D
use simple_stream_stage_solve3D,                     only: stream_stage_solve3D
use simple_commanders_stream_p01_preprocess,         only: set_preprocess_stream_cline
use simple_commanders_stream_p04_refpick_extract,    only: set_refpick_cline
use simple_commanders_stream_p05_sieve_cavgs,        only: set_sieve_cline
use simple_commanders_stream_p06_pool2D,             only: set_pool2D_cline
use simple_commanders_stream_p07_solve3D_multistate, only: set_solve3D_cline
implicit none
private
public :: run_all_stream_chain_sieve_tests, run_all_stream_chain_movie_tests

! the particles
integer, parameter :: BOX        = 64    ! sieve working and movie projection box (px)
integer, parameter :: BOX_DIRECT = 108   ! direct fixture's native particle box (px)
real,    parameter :: SMPD      = 2.0   ! A
real,    parameter :: MSKDIAM   = 100.  ! A
real,    parameter :: KV        = 300.
real,    parameter :: CS        = 2.7
real,    parameter :: FRACA     = 0.1
real,    parameter :: DEFOCUS   = 1.5   ! microns
! the driver
integer, parameter :: PASS_WAIT_S   = 5         ! between passes
integer, parameter :: TIME_LIMIT_S  = 3 * 3600  ! a test that has not finished by then fails
! the stages' work, small enough for a test
integer, parameter :: NCLS2D        = 10
integer, parameter :: NPTCLS_COARSE = 300
integer, parameter :: NPTCLS_FINE   = 600
integer, parameter :: NSTATES3D     = 2
integer, parameter :: NSTAGES3D     = 2
integer, parameter :: NTHR_JOBS     = 4
integer, parameter :: MIN_PTCLS_3D  = 5 * NSTATES3D
real,    parameter :: SIM_SNR       = 0.5

contains

    subroutine run_all_stream_chain_sieve_tests()
        write(*,'(A)') '**** running all stream chain tests: sieve to 3D ****'
        call test_sieve_to_solve3D()
    end subroutine run_all_stream_chain_sieve_tests

    subroutine run_all_stream_chain_movie_tests()
        write(*,'(A)') '**** running all stream chain tests: movies to 3D ****'
        call test_movies_to_solve3D()
    end subroutine run_all_stream_chain_movie_tests

    !> simulated particle stacks, in completed reference-picking sets, through particle sieving,
    !! pool 2D and multistate 3D: the sieve accepts particles and hands on its final set, the pool
    !! runs to its final iteration and publishes, and 3D completes its first run
    subroutine test_sieve_to_solve3D()
        integer, parameter :: NMICS = 40, NPTCLS_MIC = 30, NMICS_SET = 5, NPOSES_STATE = 12
        type(stream_stage_sieve),   allocatable :: sieve
        type(stream_stage_pool2D),  allocatable :: pool
        type(stream_stage_solve3D), allocatable :: s3D
        type(cmdline) :: cline
        type(string)  :: cwd_saved, root, dir_refpick, dir_sieve, dir_pool, dir_3D
        type(string)  :: truth1, truth2, oritab, all_stk, stk1, stk2
        integer       :: nfail0, t0
        logical       :: l_done, l_cannot_start_3D
        write(*,'(A)') 'test_sieve_to_solve3D'
        nfail0 = tests_failed
        call enter_fixture('test_stream_chain_sieve', cwd_saved, root)
        call make_ui
        call set_fixed_seed(20261005)
        dir_refpick = stage_dir(root, 'reference_based_picking')
        dir_sieve   = stage_dir(root, 'particle_sieving')
        dir_pool    = stage_dir(root, 'classification_2D')
        dir_3D      = stage_dir(root, 'solve3D_multistate')
        call simple_mkdir(dir_refpick//'/'//DIR_STREAM_COMPLETED)
        ! the two-state particles, and reference picking's completed output
        call enter(root)
        truth1 = 'truth_state1.mrc'
        truth2 = 'truth_state2.mrc'
        call write_truth_volume(truth1, 1, BOX_DIRECT)
        call write_truth_volume(truth2, 2, BOX_DIRECT)
        oritab = 'simulated_particle_orientations.txt'
        call write_repeated_orientations(oritab, NMICS * NPTCLS_MIC / 2, NPOSES_STATE)
        stk1 = simulate_particle_stack(truth1, oritab, string('simulated_particles_state1.mrcs'), NMICS * NPTCLS_MIC / 2)
        stk2 = simulate_particle_stack(truth2, oritab, string('simulated_particles_state2.mrcs'), NMICS * NPTCLS_MIC / 2)
        all_stk = interleave_stacks(stk1, stk2, NMICS * NPTCLS_MIC / 2, NPOSES_STATE, BOX_DIRECT,&
            &string('simulated_particles.mrcs'))
        call write_picking_sets(all_stk, dir_refpick, NMICS, NPTCLS_MIC, NMICS_SET, BOX_DIRECT)
        call write_moldiam(dir_refpick, BOX_DIRECT)
        call simple_touch(dir_refpick//'/'//STREAM_FINISHED_MARKER)
        ! the stages
        allocate(sieve, pool, s3D)
        call enter(dir_sieve)
        call sieve_cline(cline, dir_refpick)
        call sieve%new(cline)
        sieve%wait_s = 0
        call cline%kill
        call enter(dir_pool)
        call pool_cline(cline, dir_sieve)
        call pool%new(cline)
        pool%wait_s = 0
        call cline%kill
        call enter(dir_3D)
        call solve3D_cline(cline, dir_pool)
        call s3D%new(cline)
        s3D%wait_s = 0
        call cline%kill
        ! the passes
        t0                  = simple_gettime()
        l_done              = .false.
        l_cannot_start_3D   = .false.
        do
            call enter(dir_sieve)
            call sieve%iterate()
            call enter(dir_pool)
            call pool%iterate()
            call enter(dir_3D)
            call s3D%iterate()
            l_done = s3D%frozen_projfile%strlen() > 0
            l_cannot_start_3D = .not. l_done .and. s3D%l_final_pending .and. s3D%nptcls_selected < MIN_PTCLS_3D
            if( l_done .or. l_cannot_start_3D ) exit
            if( simple_gettime() - t0 > TIME_LIMIT_S ) exit
            call sleep(PASS_WAIT_S)
        enddo
        ! what each stage made
        call assert_true(allocated(sieve%sieve), 'chain: the sieve is made')
        if( allocated(sieve%sieve) ) call assert_true(sieve%sieve%get_n_accepted_ptcls() > 0, 'chain: the sieve accepts particles')
        call assert_true(sieve%l_final,          'chain: reference picking''s finished marker sets final ingestion')
        call assert_true(pool%l_sieve_final,     'chain: the pool takes the sieve''s final set')
        call assert_true(pool%last_export_id > 1, 'chain: the pool publishes for 3D')
        call assert_false(l_cannot_start_3D,      'chain: the final publication has enough selected particles for 3D')
        call assert_true(l_done,                 'chain: 3D completes its first run within the time limit')
        if( allocated(s3D%state_res) ) call assert_true(any(s3D%state_res > 0.), 'chain: a state has a resolution')
        ! stop, cancelling what still runs
        call enter(dir_3D)
        call s3D%finalize()
        call s3D%kill()
        call enter(dir_pool)
        call pool%finalize()
        call pool%kill()
        call enter(dir_sieve)
        call sieve%finalize()
        call sieve%kill()
        deallocate(sieve, pool, s3D)
        call leave_fixture(cwd_saved, root, nfail0)
    end subroutine test_sieve_to_solve3D

    !> simulated movies through preprocessing, optics assignment, reference picking (with picking
    !! references reprojected from both truths), particle sieving, pool 2D and multistate 3D
    subroutine test_movies_to_solve3D()
        integer, parameter :: NMOVIES = 40, NPTCLS_MIC = 24, MIC_BOX = 512, NFRAMES = 8, NREFS = 20
        type(stream_stage_preprocess), allocatable :: pre
        type(stream_stage_optics),     allocatable :: optics
        type(stream_stage_refpick),    allocatable :: refpick
        type(stream_stage_sieve),      allocatable :: sieve
        type(stream_stage_pool2D),     allocatable :: pool
        type(stream_stage_solve3D),    allocatable :: s3D
        type(cmdline) :: cline
        type(string)  :: cwd_saved, root, dir_movies, dir_pre, dir_optics, dir_refpick, dir_sieve, dir_pool, dir_3D
        type(string)  :: truth1, truth2, pickrefs, pickrefs1, pickrefs2, ptcl_stk, ptcl_stk1, ptcl_stk2
        integer       :: nfail0, t0, nprocessed
        logical       :: l_done, l_pre_active, l_cannot_start_3D
        write(*,'(A)') 'test_movies_to_solve3D'
        nfail0 = tests_failed
        call enter_fixture('test_stream_chain_movies', cwd_saved, root)
        call make_ui
        call set_fixed_seed(20261006)
        dir_movies  = stage_dir(root, 'movies')
        dir_pre     = stage_dir(root, 'preprocessing')
        dir_optics  = stage_dir(root, 'optics_assignment')
        dir_refpick = stage_dir(root, 'reference_based_picking')
        dir_sieve   = stage_dir(root, 'particle_sieving')
        dir_pool    = stage_dir(root, 'classification_2D')
        dir_3D      = stage_dir(root, 'solve3D_multistate')
        ! the two truths, their reprojections for the movies and picking references, the movies
        call enter(root)
        truth1 = 'truth_state1.mrc'
        truth2 = 'truth_state2.mrc'
        call write_truth_volume(truth1, 1, BOX)
        call write_truth_volume(truth2, 2, BOX)
        ptcl_stk1 = reproject_truth(truth1, string('movie_particles_state1.mrcs'), NPTCLS_MIC / 2)
        ptcl_stk2 = reproject_truth(truth2, string('movie_particles_state2.mrcs'), NPTCLS_MIC / 2)
        ptcl_stk  = interleave_stacks(ptcl_stk1, ptcl_stk2, NPTCLS_MIC / 2, NPTCLS_MIC / 2, BOX,&
            &string('movie_particles.mrcs'))
        pickrefs1 = reproject_truth(truth1, string('pickrefs_state1.mrcs'), NREFS / 2)
        pickrefs2 = reproject_truth(truth2, string('pickrefs_state2.mrcs'), NREFS / 2)
        pickrefs  = interleave_stacks(pickrefs1, pickrefs2, NREFS / 2, NREFS / 2, BOX, string('pickrefs.mrcs'))
        call simulate_movies(ptcl_stk, dir_movies, NMOVIES, MIC_BOX, NFRAMES)
        ! the stages (each in its folder)
        allocate(pre, optics, refpick, sieve, pool, s3D)
        call enter(dir_pre)
        call preprocess_cline(cline, dir_movies)
        call pre%new(cline)
        pre%wait_s = 0
        call cline%kill
        call enter(dir_optics)
        call optics_cline(cline, dir_pre)
        call optics%new(cline)
        optics%wait_s = 0
        call cline%kill
        call enter(dir_refpick)
        call refpick_cline(cline, dir_pre, dir_optics, pickrefs)
        call refpick%new(cline)
        refpick%wait_s = 0
        call cline%kill
        call enter(dir_sieve)
        call sieve_cline(cline, dir_refpick, dir_optics)
        call sieve%new(cline)
        sieve%wait_s = 0
        call cline%kill
        call enter(dir_pool)
        call pool_cline(cline, dir_sieve, dir_optics)
        call pool%new(cline)
        pool%wait_s = 0
        call cline%kill
        call enter(dir_3D)
        call solve3D_cline(cline, dir_pool)
        call s3D%new(cline)
        s3D%wait_s = 0
        call cline%kill
        ! the passes; preprocessing stops once every movie is processed, and its finished marker
        ! ends the intake downstream
        t0                  = simple_gettime()
        l_done              = .false.
        l_pre_active        = .true.
        l_cannot_start_3D   = .false.
        do
            if( l_pre_active )then
                call enter(dir_pre)
                call pre%iterate()
                nprocessed = pre%spproj_glob%os_mic%get_noris() + pre%n_failed_jobs
                if( nprocessed >= NMOVIES )then
                    call pre%finalize()
                    l_pre_active = .false.
                endif
            endif
            call enter(dir_optics)
            call optics%iterate()
            call enter(dir_refpick)
            call refpick%iterate()
            call enter(dir_sieve)
            call sieve%iterate()
            call enter(dir_pool)
            call pool%iterate()
            call enter(dir_3D)
            call s3D%iterate()
            l_done = s3D%frozen_projfile%strlen() > 0
            l_cannot_start_3D = .not. l_done .and. s3D%l_final_pending .and. s3D%nptcls_selected < MIN_PTCLS_3D
            if( l_done .or. l_cannot_start_3D ) exit
            if( simple_gettime() - t0 > TIME_LIMIT_S ) exit
            call sleep(PASS_WAIT_S)
        enddo
        ! what each stage made
        call assert_false(l_pre_active,               'chain: preprocessing processes every movie')
        call assert_true(file_exists(dir_pre//'/'//STREAM_FINISHED_MARKER), 'chain: and leaves its finished marker')
        call assert_true(optics%map_id > 0,           'chain: optics assignment publishes a map')
        call assert_true(refpick%nptcls_glob > 0,     'chain: reference picking extracts particles')
        call assert_true(sieve%l_final,               'chain: the sieve''s intake ends')
        if( allocated(sieve%sieve) ) call assert_true(sieve%sieve%get_n_accepted_ptcls() > 0, 'chain: the sieve accepts particles')
        call assert_true(pool%last_export_id > 1,     'chain: the pool publishes for 3D')
        call assert_false(l_cannot_start_3D,          'chain: the final publication has enough selected particles for 3D')
        call assert_true(l_done,                      'chain: 3D completes its first run within the time limit')
        ! stop, cancelling what still runs
        call enter(dir_3D)
        call s3D%finalize()
        call s3D%kill()
        call enter(dir_pool)
        call pool%finalize()
        call pool%kill()
        call enter(dir_sieve)
        call sieve%finalize()
        call sieve%kill()
        call enter(dir_refpick)
        call refpick%finalize()
        call refpick%kill()
        call enter(dir_optics)
        call optics%finalize()
        call optics%kill()
        call enter(dir_pre)
        if( l_pre_active ) call pre%finalize()
        call pre%kill()
        deallocate(pre, optics, refpick, sieve, pool, s3D)
        call leave_fixture(cwd_saved, root, nfail0)
    end subroutine test_movies_to_solve3D

    ! ---- the stages' command lines: the commanders' defaults, the test's sizes ------------------

    ! common to every stage: a program name outside every UI table (params%new makes no folder),
    ! no output folder, the local queue
    subroutine base_cline( cline, prg, projfile )
        type(cmdline),    intent(inout) :: cline
        character(len=*), intent(in)    :: prg, projfile
        call cline%set('prg',       prg)
        call cline%set('projfile',  projfile)
        call cline%set('qsys_name', 'local')
        call cline%set('walltime',  3600)
    end subroutine base_cline

    ! after a commander's defaults, which ask for an output folder of the stage's own
    subroutine no_folder( cline )
        type(cmdline), intent(inout) :: cline
        call cline%set('mkdir',  'no')
        call cline%set('outdir', '')
    end subroutine no_folder

    subroutine preprocess_cline( cline, dir_movies )
        type(cmdline), intent(inout) :: cline
        class(string), intent(in)    :: dir_movies
        call base_cline(cline, 'stream_chain_preprocess', 'preprocessing.simple')
        call cline%set('dir_movies', dir_movies)
        call cline%set('smpd',       SMPD)
        call cline%set('kv',         KV)
        call cline%set('cs',         CS)
        call cline%set('fraca',      FRACA)
        call cline%set('nthr',       NTHR_JOBS)
        call cline%set('nparts',     2)
        call cline%set('pspecsz',    256)
        call set_preprocess_stream_cline(cline)
        call no_folder(cline)
    end subroutine preprocess_cline

    subroutine optics_cline( cline, dir_pre )
        type(cmdline), intent(inout) :: cline
        class(string), intent(in)    :: dir_pre
        call base_cline(cline, 'stream_chain_optics', 'optics_assignment.simple')
        call cline%set('dir_target', dir_pre)
        call cline%set('nthr',       1)
        call no_folder(cline)
    end subroutine optics_cline

    subroutine refpick_cline( cline, dir_pre, dir_optics, pickrefs )
        type(cmdline), intent(inout) :: cline
        class(string), intent(in)    :: dir_pre, dir_optics, pickrefs
        call base_cline(cline, 'stream_chain_refpick', 'reference_based_picking.simple')
        call cline%set('dir_target', dir_pre)
        call cline%set('optics_dir', dir_optics)
        call cline%set('pickrefs',   pickrefs)
        call cline%set('nthr',       NTHR_JOBS)
        call cline%set('nparts',     2)
        call set_refpick_cline(cline)
        call no_folder(cline)
    end subroutine refpick_cline

    subroutine sieve_cline( cline, dir_refpick, dir_optics )
        type(cmdline),           intent(inout) :: cline
        class(string),           intent(in)    :: dir_refpick
        class(string), optional, intent(in)    :: dir_optics
        call base_cline(cline, 'stream_chain_sieve', 'particle_sieving.simple')
        call cline%set('dir_target',     dir_refpick)
        if( present(dir_optics) ) call cline%set('optics_dir', dir_optics)
        call cline%set('nthr',           NTHR_JOBS)
        call cline%set('nparts',         1)
        call cline%set('nchunks',        2)
        call cline%set('nptcls_coarse',  NPTCLS_COARSE)
        call cline%set('nptcls_fine',    NPTCLS_FINE)
        call cline%set('nsample_coarse', NPTCLS_COARSE)
        call cline%set('nsample_fine',   NPTCLS_FINE)
        call cline%set('box_coarse',     BOX)
        call cline%set('box_fine',       BOX)
        call cline%set('ncls_coarse',    NCLS2D)
        call cline%set('ncls_fine',      NCLS2D)
        ! the hard gates only: the learned model was not trained on synthetic class averages
        call cline%set('use_model',      'no')
        call set_sieve_cline(cline)
        call no_folder(cline)
    end subroutine sieve_cline

    subroutine pool_cline( cline, dir_sieve, dir_optics )
        type(cmdline),           intent(inout) :: cline
        class(string),           intent(in)    :: dir_sieve
        class(string), optional, intent(in)    :: dir_optics
        call base_cline(cline, 'stream_chain_pool2D', 'classification_2D.simple')
        call cline%set('dir_target', dir_sieve)
        if( present(dir_optics) ) call cline%set('optics_dir', dir_optics)
        call cline%set('ncls',       NCLS2D)
        call cline%set('nthr',       NTHR_JOBS)
        call cline%set('nparts',     1)
        call set_pool2D_cline(cline)
        call no_folder(cline)
    end subroutine pool_cline

    subroutine solve3D_cline( cline, dir_pool )
        type(cmdline), intent(inout) :: cline
        class(string), intent(in)    :: dir_pool
        call base_cline(cline, 'stream_chain_solve3D', 'solve3D_multistate.simple')
        call cline%set('dir_target', dir_pool)
        call cline%set('nthr',       NTHR_JOBS)
        call cline%set('nparts',     1)
        call cline%set('nstates',    NSTATES3D)
        call cline%set('nstages',    NSTAGES3D)
        call cline%set('nthr3D',     NTHR_JOBS)
        call cline%set('nparts3D',   1)
        call set_solve3D_cline(cline)
        call no_folder(cline)
    end subroutine solve3D_cline

    ! ---- fixtures -----------------------------------------------------------------------------

    ! @p root/@p name, made; its absolute path
    function stage_dir( root, name ) result( dir )
        class(string),    intent(in) :: root
        character(len=*), intent(in) :: name
        type(string) :: dir
        dir = root//'/'//name
        call simple_mkdir(dir)
    end function stage_dir

    ! enters @p dir, where the stage it belongs to works, for the queue's scripts too
    subroutine enter( dir )
        class(string), intent(in) :: dir
        call simple_chdir(dir)
        CWD_GLOB = dir%to_char()
    end subroutine enter

    ! one of two distinct asymmetric multi-lobed objects in a @p box^3 volume
    subroutine write_truth_volume( fname, state, box )
        class(string), intent(in) :: fname
        integer,       intent(in) :: state, box
        integer, parameter :: NBLOBS = 9
        type(image)       :: vol
        real, allocatable :: rmat(:,:,:)
        integer :: i, j, k, iblob
        real    :: x(3), c, cen(3,NBLOBS), sig(NBLOBS), weight(NBLOBS)
        if( state == 1 )then
            cen = reshape([0.,0.,0.,  7.,0.,0.,  -7.,0.,0.,  0.,7.,0.,  0.,-7.,0.,  0.,0.,7.,  0.,0.,-7.,&
                &12.,8.,-5.,  -11.,9.,6.], [3,NBLOBS])
            sig    = [3.5, 2.5, 2.5, 2.5, 2.5, 2.5, 2.5, 2., 2.]
            weight = [1., 0.8, 0.7, 0.75, 0.65, 0.8, 0.7, 0.6, 0.55]
        else
            cen = reshape([-14.,-5.,-4.,  -9.,-1.,1.,  -4.,4.,5.,  2.,8.,1.,  7.,4.,-5.,  12.,-1.,-1.,&
                &15.,-7.,5.,  0.,-10.,8.,  -8.,10.,-7.], [3,NBLOBS])
            sig    = [2.5, 2.8, 2.5, 2.5, 2.5, 2.8, 2.2, 2.2, 2.2]
            weight = [0.65, 0.8, 0.7, 0.9, 0.75, 1., 0.6, 0.55, 0.5]
        endif
        allocate(rmat(box,box,box), source=0.)
        c = real(box / 2 + 1)
        do k = 1,box
            do j = 1,box
                do i = 1,box
                    x = [real(i), real(j), real(k)] - c
                    do iblob = 1,NBLOBS
                        rmat(i,j,k) = rmat(i,j,k) + WEIGHT(iblob) * exp(-sum((x - CEN(:,iblob))**2) / (2. * SIG(iblob)**2))
                    enddo
                enddo
            enddo
        enddo
        call vol%new([box,box,box], SMPD, wthreads=.false.)
        call vol%set_rmat(rmat, .false.)
        call vol%write(fname)
        call vol%kill
    end subroutine write_truth_volume

    ! @p nptcls simulated particle images of @p vol, with CTF (one defocus) and noise; the stack's
    ! absolute path
    function simulate_particle_stack( vol, oritab, outstk, nptcls ) result( stk )
        class(string), intent(in) :: vol, oritab, outstk
        integer,       intent(in) :: nptcls
        type(string) :: stk
        type(commander_simulate_particles) :: xsim
        type(cmdline) :: cline
        call cline%set('prg',      'simulate_particles')
        call cline%set('mkdir',    'no')
        call cline%set('vol1',     vol)
        call cline%set('oritab',   oritab)
        call cline%set('outstk',   outstk)
        call cline%set('outfile',  outstk//'.txt')
        call cline%set('nptcls',   nptcls)
        call cline%set('nthr',     NTHR_JOBS)
        call cline%set('smpd',     SMPD)
        call cline%set('mskdiam',  MSKDIAM)
        call cline%set('pgrp',     'c1')
        call cline%set('ctf',      'yes')
        call cline%set('snr',      SIM_SNR)
        call cline%set('kv',       KV)
        call cline%set('cs',       CS)
        call cline%set('fraca',    FRACA)
        call cline%set('defocus',  DEFOCUS)
        call cline%set('dferr',    0.)
        call cline%set('astigerr', 0.)
        call cline%set('bfac',     0.)
        call cline%set('bfacerr',  0.)
        call cline%set('sherr',    2.)
        call xsim%execute(cline)
        call cline%kill
        stk = simple_abspath(outstk)
    end function simulate_particle_stack

    ! repeat a finite spiral so every pose has enough particles to make a coherent class average
    subroutine write_repeated_orientations( fname, nptcls, nposes )
        class(string), intent(in) :: fname
        integer,       intent(in) :: nptcls, nposes
        type(oris) :: poses, particles
        type(ori)  :: pose
        integer    :: i
        call poses%new(nposes, is_ptcl=.false.)
        call poses%spiral()
        call particles%new(nptcls, is_ptcl=.true.)
        do i = 1,nptcls
            call poses%get_ori(mod(i - 1, nposes) + 1, pose)
            call particles%set_ori(i, pose)
        enddo
        call particles%write(fname, [1,nptcls])
        call pose%kill
        call particles%kill
        call poses%kill
    end subroutine write_repeated_orientations

    ! interleave two equal-size stacks, flipping their order every @p swap_stride pairs
    function interleave_stacks( stk1, stk2, n_each, swap_stride, box, outstk ) result( stk )
        class(string), intent(in) :: stk1, stk2, outstk
        integer,       intent(in) :: n_each, swap_stride, box
        type(string) :: stk
        type(image)  :: img
        integer      :: i
        call img%new([box,box,1], SMPD, wthreads=.false.)
        do i = 1,n_each
            if( mod((i - 1) / swap_stride, 2) == 0 )then
                call img%read(stk1, i)
                call img%write(outstk, 2 * i - 1)
                call img%read(stk2, i)
                call img%write(outstk, 2 * i)
            else
                call img%read(stk2, i)
                call img%write(outstk, 2 * i - 1)
                call img%read(stk1, i)
                call img%write(outstk, 2 * i)
            endif
        enddo
        call img%kill
        stk = simple_abspath(outstk)
    end function interleave_stacks

    ! reference picking's completed sets in @p dir_refpick: the particles of @p all_stk split into
    ! @p nmics stacks of @p nptcls_mic, one per micrograph, @p nmics_set micrographs per set
    subroutine write_picking_sets( all_stk, dir_refpick, nmics, nptcls_mic, nmics_set, box )
        class(string), intent(in) :: all_stk, dir_refpick
        integer,       intent(in) :: nmics, nptcls_mic, nmics_set, box
        type(image)      :: img
        type(sp_project) :: set
        type(ctfparams)  :: ctfvars
        type(string)     :: stk, fname
        integer :: imic, iset, j, jmic, nsets
        ctfvars%ctfflag = CTFFLAG_YES
        ctfvars%smpd    = SMPD
        ctfvars%kv      = KV
        ctfvars%cs      = CS
        ctfvars%fraca   = FRACA
        ctfvars%dfx     = DEFOCUS
        ctfvars%dfy     = DEFOCUS
        call img%new([box,box,1], SMPD, wthreads=.false.)
        nsets = nmics / nmics_set
        imic  = 0
        do iset = 1,nsets
            call set%os_mic%new(nmics_set, is_ptcl=.false.)
            do jmic = 1,nmics_set
                imic = imic + 1
                stk  = dir_refpick//'/mic_'//int2str_pad(imic,3)//'.mrcs'
                do j = 1,nptcls_mic
                    call img%read(all_stk, (imic - 1) * nptcls_mic + j)
                    call img%write(stk, j)
                enddo
                call set%add_stk(stk, ctfvars)
                call set%os_mic%set(jmic, 'intg',      dir_refpick//'/mic_'//int2str_pad(imic,3)//'_intg.mrc')
                call set%os_mic%set(jmic, 'imgkind',   'mic')
                call set%os_mic%set(jmic, 'state',     1)
                call set%os_mic%set(jmic, 'nptcls',    nptcls_mic)
                call set%os_mic%set(jmic, 'importind', imic)
                call set%os_mic%set(jmic, 'smpd',      SMPD)
                call set%os_mic%set(jmic, 'kv',        KV)
                call set%os_mic%set(jmic, 'cs',        CS)
                call set%os_mic%set(jmic, 'fraca',     FRACA)
                call set%os_mic%set(jmic, 'dfx',       DEFOCUS)
                call set%os_mic%set(jmic, 'dfy',       DEFOCUS)
                call set%os_mic%set(jmic, 'angast',    0.)
                call set%os_mic%set(jmic, 'ctf',       'yes')
                call set%os_mic%set(jmic, 'xdim',      1024)
                call set%os_mic%set(jmic, 'ydim',      1024)
            enddo
            fname = dir_refpick//'/'//DIR_STREAM_COMPLETED//int2str_pad(iset, 5)//METADATA_EXT
            call set%update_projinfo(fname)
            call set%write(fname)
            call set%kill
        enddo
        call img%kill
    end subroutine write_picking_sets

    ! the moldiam.txt make_pickrefs writes beside the picking references
    subroutine write_moldiam( dir_refpick, box )
        class(string), intent(in) :: dir_refpick
        integer,       intent(in) :: box
        type(oris) :: moldiam
        call moldiam%new(1, is_ptcl=.false.)
        call moldiam%set(1, 'mskdiam',         MSKDIAM)
        call moldiam%set(1, 'moldiam',         MSKDIAM / 1.2)
        call moldiam%set(1, 'box_for_extract', box)
        call moldiam%write(dir_refpick//'/'//STREAM_MOLDIAM)
        call moldiam%kill
    end subroutine write_moldiam

    ! @p n reprojections of @p vol in @p stk (no CTF); the stack's absolute path
    function reproject_truth( vol, stk, n ) result( stk_abs )
        class(string), intent(in) :: vol, stk
        integer,       intent(in) :: n
        type(string) :: stk_abs
        type(commander_reproject) :: xreproject
        type(cmdline) :: cline
        call cline%set('prg',     'reproject')
        call cline%set('mkdir',   'no')
        call cline%set('vol1',    vol)
        call cline%set('outstk',  stk)
        call cline%set('smpd',    SMPD)
        call cline%set('pgrp',    'c1')
        call cline%set('mskdiam', MSKDIAM)
        call cline%set('nspace',  n)
        call cline%set('nthr',    1)
        call xreproject%execute(cline)
        call cline%kill
        stk_abs = simple_abspath(stk)
    end function reproject_truth

    ! @p nmovies movies of @p nframes frames of @p mic_box pixels in @p dir_movies, each with the
    ! particles of @p ptcl_stk under a CTF (defocus varying over the movies) and noise
    subroutine simulate_movies( ptcl_stk, dir_movies, nmovies, mic_box, nframes )
        class(string), intent(in) :: ptcl_stk, dir_movies
        integer,       intent(in) :: nmovies, mic_box, nframes
        character(len=*), parameter :: MOVIE_FILE = 'simulate_movie.mrc'
        type(commander_simulate_movie) :: xsim
        type(cmdline) :: cline
        type(string)  :: simdir
        integer       :: i
        simdir = stage_dir(string(CWD_GLOB), 'simulate_movies')
        call enter(simdir)
        do i = 1,nmovies
            call cline%set('prg',     'simulate_movie')
            call cline%set('mkdir',   'no')
            call cline%set('stk',     ptcl_stk)
            call cline%set('xdim',    mic_box)
            call cline%set('ydim',    mic_box)
            call cline%set('nframes', nframes)
            call cline%set('smpd',    SMPD)
            call cline%set('snr',     SIM_SNR)
            call cline%set('kv',      KV)
            call cline%set('cs',      CS)
            call cline%set('fraca',   FRACA)
            call cline%set('defocus', DEFOCUS + 0.05 * real(mod(i - 1, 10)))
            call cline%set('trs',     1.0)
            call cline%set('nthr',    1)
            call xsim%execute(cline)
            call cline%kill
            call simple_rename(MOVIE_FILE, dir_movies//'/movie_'//int2str_pad(i,3)//MRC_EXT)
        enddo
    end subroutine simulate_movies

end module simple_stream_chain_tester
